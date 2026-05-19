#!/usr/bin/env python3
"""Lightweight internal validation for the current Majorana/Pauli model.

This script does not attempt paper-level reproduction. It checks the internal
consistency of the model implementation:
- Hermiticity of representative Hamiltonians
- Unitarity of the full evolution operator
- Agreement between state-propagation and propagator-based evolution
- Convergence with respect to time-step resolution
- Full-track vs reduced effective-track fidelity consistency
"""

from __future__ import annotations

import importlib.util
import os
import sys
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm


def _load_fig2_module():
    this_dir = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(this_dir, "strict_repro_fig2.py")
    spec = importlib.util.spec_from_file_location("strict_repro_fig2_mod", mod_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {mod_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


fig2 = _load_fig2_module()


@dataclass(frozen=True)
class CheckResult:
    name: str
    ok: bool
    value: float
    threshold: float
    detail: str


def _gamma_ops(n_modes: int = 3):
    c_ops = fig2.mod.fermionic_operators(n_modes)
    gammas = []
    for c in c_ops:
        cd = c.conj().T
        gammas.append(c + cd)
        gammas.append(-1j * (c - cd))
    return gammas


def _hamiltonian(gammas, params, tau, tc, t):
    t1v, t2v, t3v, Edv = fig2.mod.time_profiles(t, tau, tc)
    return fig2.mod.build_H(gammas, params, t1v, t2v, t3v, Edv)


def total_propagator(params: Dict[str, float], tau: float, tc: float, steps_per_tau: int, repeat: int = 2):
    gammas = _gamma_ops(3)
    dim = 2 ** 3
    total_time = repeat * 3.0 * tau
    steps = max(20, int(repeat * 3.0 * tau * steps_per_tau))
    dt = total_time / steps
    U = np.eye(dim, dtype=complex)
    t = 0.0
    for _ in range(steps):
        t += dt
        Ht = _hamiltonian(gammas, params, tau=tau, tc=tc, t=t)
        U = expm(-1j * Ht * dt) @ U
    return U


def state_from_propagator(U):
    psi0 = np.zeros(U.shape[0], dtype=complex)
    psi0[0] = 1.0
    return U @ psi0


def fidelity_from_state(psi):
    psi0 = np.zeros_like(psi)
    psi0[0] = 1.0
    return float(abs(np.vdot(psi0, psi)))


def check_hermiticity(params: Dict[str, float], tau: float, tc: float, times: List[float]):
    gammas = _gamma_ops(3)
    results = []
    for t in times:
        Ht = _hamiltonian(gammas, params, tau=tau, tc=tc, t=t)
        herm = float(np.linalg.norm(Ht.conj().T - Ht, ord='fro'))
        evals = np.linalg.eigvals(Ht)
        imag_max = float(np.max(np.abs(np.imag(evals))))
        results.append((t, herm, imag_max))
    return results


def check_unitarity_and_consistency(params: Dict[str, float], tau: float, tc: float, steps_per_tau: int):
    U = total_propagator(params, tau=tau, tc=tc, steps_per_tau=steps_per_tau, repeat=2)
    dim = U.shape[0]
    I = np.eye(dim, dtype=complex)
    unitarity = float(np.linalg.norm(U.conj().T @ U - I, ord='fro'))

    fid_ref, psi_ref, _, _ = fig2.mod.run_braiding(params, tau=tau, steps_per_tau=steps_per_tau, repeat=2, tc=tc)
    psi_prop = state_from_propagator(U)
    fid_prop = fidelity_from_state(psi_prop)
    state_diff = float(np.linalg.norm(psi_ref - psi_prop))
    fid_diff = float(abs(fid_ref - fid_prop))

    return {
        'unitarity': unitarity,
        'fid_ref': float(fid_ref),
        'fid_prop': fid_prop,
        'fid_diff': fid_diff,
        'state_diff': state_diff,
    }


def check_time_step_convergence(params: Dict[str, float], tau: float, tc: float, steps_list: List[int]):
    vals = []
    for steps_per_tau in steps_list:
        fid, _, _, _ = fig2.mod.run_braiding(params, tau=tau, steps_per_tau=steps_per_tau, repeat=2, tc=tc)
        vals.append((steps_per_tau, float(fid)))
    return vals


def check_full_vs_effective(anchor_idx: int, tau_grid: np.ndarray):
    anchors = fig2.Fig2Anchors()
    e1 = float(anchors.e1[anchor_idx])
    tc = float(anchors.t1[anchor_idx] * 5.88)
    alpha_grid = np.linspace(0.2, 2.5, 28)

    full_curve = np.array([
        fig2.run_braiding_full_track(e1=e1, tc=tc, tau=max(float(t), 0.2), steps_per_tau=80, repeat=2)
        for t in tau_grid
    ])
    alpha_star = fig2.calibrate_effective_alpha(
        e1=e1,
        tc=tc,
        tau_grid=tau_grid,
        full_curve=full_curve,
        alpha_grid=alpha_grid,
    )
    eff_curve = np.array([
        fig2.run_braiding_effective_track(e1=e1, tc=tc, tau=max(float(t), 0.2), steps_per_tau=160, repeat=2, alpha=alpha_star)
        for t in tau_grid
    ])

    rmse = float(np.sqrt(np.mean((full_curve - eff_curve) ** 2)))
    mae = float(np.mean(np.abs(full_curve - eff_curve)))
    corr = float(np.corrcoef(full_curve, eff_curve)[0, 1]) if tau_grid.size > 2 else float('nan')
    return {
        'B': float(anchors.b[anchor_idx]),
        'E1': e1,
        'tc': tc,
        'alpha': float(alpha_star),
        'rmse': rmse,
        'mae': mae,
        'corr': corr,
        'full': full_curve,
        'eff': eff_curve,
    }


def summarize(results: List[CheckResult]):
    lines = []
    lines.append('Internal validation summary')
    lines.append('=' * 80)
    for r in results:
        status = 'PASS' if r.ok else 'FAIL'
        lines.append(f'[{status}] {r.name}: value={r.value:.6e}, threshold={r.threshold:.6e} | {r.detail}')
    return '\n'.join(lines)


def plot_validation_summary(checks: List[CheckResult], conv, eff_checks, output_png: str):
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], width_ratios=[1.0, 1.25])

    ax0 = fig.add_subplot(gs[0, 0])
    names = [
        'Hermiticity',
        'Unitarity',
        'State consistency',
        'Step conv 80-160',
        'Step conv 40-80',
        'RMSE B=2.25T',
        'RMSE B=2.55T',
        'RMSE B=2.75T',
    ]
    values = [
        checks[0].value,
        checks[3].value,
        checks[4].value,
        checks[5].value,
        checks[6].value,
        checks[7].value,
        checks[8].value,
        checks[9].value,
    ]
    thresholds = [
        checks[0].threshold,
        checks[3].threshold,
        checks[4].threshold,
        checks[5].threshold,
        checks[6].threshold,
        checks[7].threshold,
        checks[8].threshold,
        checks[9].threshold,
    ]
    x = np.arange(len(values))
    ax0.bar(x - 0.18, values, width=0.36, label='value', color='steelblue')
    ax0.bar(x + 0.18, thresholds, width=0.36, label='threshold', color='orange', alpha=0.8)
    ax0.set_yscale('log')
    ax0.set_xticks(x)
    ax0.set_xticklabels(names, rotation=45, ha='right', fontsize=9)
    ax0.set_title('Validation magnitudes')
    ax0.set_ylabel('log scale')
    ax0.legend(frameon=False, fontsize=9)

    ax1 = fig.add_subplot(gs[0, 1])
    steps = [s for s, _ in conv]
    fids = [f for _, f in conv]
    ax1.plot(steps, fids, '-o', color='tab:green', lw=2)
    ax1.set_title('Step-resolution convergence at τ=20')
    ax1.set_xlabel('steps_per_tau')
    ax1.set_ylabel('fidelity')
    ax1.grid(alpha=0.25)

    ax2 = fig.add_subplot(gs[1, 0])
    b_labels = [f"B={eff['B']:.2f}T" for eff in eff_checks]
    rmse_vals = [eff['rmse'] for eff in eff_checks]
    corr_vals = [eff['corr'] for eff in eff_checks]
    x2 = np.arange(len(b_labels))
    ax2.bar(x2, rmse_vals, color='tab:purple', alpha=0.85, label='RMSE')
    ax2_t = ax2.twinx()
    ax2_t.plot(x2, corr_vals, '-s', color='tab:red', lw=2, label='corr')
    ax2.set_xticks(x2)
    ax2.set_xticklabels(b_labels)
    ax2.set_ylim(0, max(rmse_vals) * 1.25)
    ax2_t.set_ylim(0.95, 1.0)
    ax2.set_title('Full vs effective-track agreement')
    ax2.set_ylabel('RMSE')
    ax2_t.set_ylabel('correlation')
    ax2.grid(axis='y', alpha=0.2)

    ax3 = fig.add_subplot(gs[1, 1])
    tau_grid = np.array([5.0, 12.0, 20.0, 30.0])
    colors = ['black', 'red', 'limegreen']
    for idx, eff in enumerate(eff_checks):
        ax3.plot(tau_grid, eff['full'], '-', color=colors[idx], lw=2.0, label=f"full {b_labels[idx]}")
        ax3.plot(tau_grid, eff['eff'], '--', color=colors[idx], lw=2.0, label=f"effective {b_labels[idx]}")
    ax3.set_title('Full vs effective fidelity curves')
    ax3.set_xlabel(r'$\tau$')
    ax3.set_ylabel('fidelity')
    ax3.set_ylim(0, 1.02)
    ax3.grid(alpha=0.25)
    ax3.legend(frameon=False, fontsize=8, ncol=2)

    fig.suptitle('Internal validation summary for the current model', y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def main():
    os.makedirs('quantity', exist_ok=True)
    report_path = os.path.join('quantity', 'internal_validation_report.txt')
    figure_path = os.path.join('quantity', 'internal_validation_summary.png')

    anchors = fig2.Fig2Anchors()
    tau_probe = 20.0
    tc_probe = float(anchors.t1[1] * 5.88)
    params_probe = {'E1': float(anchors.e1[1]), 'E2': 0.0}

    checks: List[CheckResult] = []

    for t, herm, imag_max in check_hermiticity(params_probe, tau=tau_probe, tc=tc_probe, times=[0.0, tau_probe / 2.0, tau_probe]):
        checks.append(CheckResult(
            name=f'Hermiticity t={t:.2f}',
            ok=herm < 1e-10 and imag_max < 1e-10,
            value=max(herm, imag_max),
            threshold=1e-10,
            detail='Hamiltonian should be Hermitian and eigenvalues should be real',
        ))

    cons = check_unitarity_and_consistency(params_probe, tau=tau_probe, tc=tc_probe, steps_per_tau=80)
    checks.append(CheckResult(
        name='Propagator unitarity',
        ok=cons['unitarity'] < 1e-8,
        value=cons['unitarity'],
        threshold=1e-8,
        detail='U^†U should be close to identity',
    ))
    checks.append(CheckResult(
        name='State propagation consistency',
        ok=cons['state_diff'] < 1e-8 and cons['fid_diff'] < 1e-10,
        value=max(cons['state_diff'], cons['fid_diff']),
        threshold=1e-8,
        detail='Sequential state update must match total propagator action',
    ))

    conv = check_time_step_convergence(params_probe, tau=tau_probe, tc=tc_probe, steps_list=[40, 80, 160])
    fid_map = {steps: fid for steps, fid in conv}
    diff_80_160 = float(abs(fid_map[80] - fid_map[160]))
    diff_40_80 = float(abs(fid_map[40] - fid_map[80]))
    checks.append(CheckResult(
        name='Time-step convergence 80 vs 160',
        ok=diff_80_160 < 1e-4,
        value=diff_80_160,
        threshold=1e-4,
        detail='Fidelity should stabilize under step refinement',
    ))
    checks.append(CheckResult(
        name='Time-step convergence 40 vs 80',
        ok=diff_40_80 < 1e-3,
        value=diff_40_80,
        threshold=1e-3,
        detail='Coarser step should be close but not identical',
    ))

    tau_grid = np.array([5.0, 12.0, 20.0, 30.0])
    eff_checks = []
    for idx in range(3):
        eff = check_full_vs_effective(idx, tau_grid)
        eff_checks.append(eff)
        checks.append(CheckResult(
            name=f'Full vs effective RMSE B={eff["B"]:.2f}T',
            ok=eff['rmse'] < 0.20,
            value=eff['rmse'],
            threshold=0.20,
            detail=f'corr={eff["corr"]:.4f}, alpha={eff["alpha"]:.4f}',
        ))

    summary = summarize(checks)

    lines = [summary, '', 'Step-resolution fidelities at tau=20.0:']
    for steps, fid in conv:
        lines.append(f'  steps_per_tau={steps:3d} -> fidelity={fid:.8f}')

    lines.append('')
    lines.append('Full vs effective track summary:')
    for eff in eff_checks:
        lines.append(
            f'  B={eff["B"]:.2f}T | E1={eff["E1"]:.3e} | tc={eff["tc"]:.4f} | alpha={eff["alpha"]:.4f} | rmse={eff["rmse"]:.6f} | corr={eff["corr"]:.6f}'
        )

    report = '\n'.join(lines)
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)

    plot_validation_summary(checks, conv, eff_checks, figure_path)

    print(report)
    print(f'\nSaved report: {report_path}')
    print(f'Saved figure: {figure_path}')

    failed = [r for r in checks if not r.ok]
    if failed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()