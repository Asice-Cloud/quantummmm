#!/usr/bin/env python3
"""Reproduce the qualitative trends of Fig.3/Fig.4/Fig.5 with the simplified models.

This script is intentionally not a strict numerical reproduction of the paper.
It keeps the paper-style axis conventions and shows the same qualitative trends:

- Fig.3 trend: ABS-like spectrum evolution and overlap vs T comparison
- Fig.4 trend: modulation-induced dynamical-phase cancellation
- Fig.5 trend: oscillation-amplitude insensitivity in the adiabatic regime

The goal is trend-level agreement and paper-style plotting, not exact raw data
matching.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh, expm

# Allow `python tools/reproduce_trend_figs.py` from repo root.
THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import embed_kitaev
from tools import tetron_path_sim as tetron
import tools.paper_params as P


# Fig.3 uses an auto-scan-favored display window so the MZM/ABS separation
# is visible in the default output.
FIG3_DISPLAY_TMIN = 260.0
FIG3_DISPLAY_TMAX = 640.0
FIG3_DISPLAY_NPTS = 16
FIG3_DISPLAY_ABS_DELTA = 0.16
FIG3_DISPLAY_N_PER_STEP = 180


def trend_score_fig3(mzm_overlaps: np.ndarray, abs_overlaps: np.ndarray) -> float:
    """Higher is better for Fig.3-like trend.

    Desired trend:
    - MZM-like curve is more stable in T (smaller std)
    - ABS-like curve varies more strongly in T
    - MZM-like overlap tends to be higher on average
    """
    mzm_std = float(np.std(mzm_overlaps))
    abs_std = float(np.std(abs_overlaps))
    mean_gap = float(np.mean(mzm_overlaps) - np.mean(abs_overlaps))
    return (abs_std - mzm_std) + 0.5 * mean_gap


def trend_score_fig4(overlaps_a: np.ndarray, overlaps_b: np.ndarray) -> float:
    """Higher is better for Fig.4-like trend.

    Desired trend:
    - Two modulation cases are both relatively T-insensitive
    - Their behaviors are close after cancellation
    """
    flat_penalty = float(np.std(overlaps_a) + np.std(overlaps_b))
    mismatch_penalty = float(np.mean(np.abs(overlaps_a - overlaps_b)))
    return 1.0 - flat_penalty - 0.7 * mismatch_penalty


def trend_score_fig5(curve1: np.ndarray, curve2: np.ndarray) -> float:
    """Higher is better for Fig.5-like trend.

    Desired trend:
    - Curves from two amplitudes are close in adiabatic (large-T) region
    - Large-T tails are relatively flat
    """
    tail = slice(len(curve1) // 2, None)
    diff_tail = float(np.mean(np.abs(curve1[tail] - curve2[tail])))
    flat_tail = float(np.std(curve1[tail]) + np.std(curve2[tail]))
    mean_tail = float(0.5 * (np.mean(curve1[tail]) + np.mean(curve2[tail])))
    return mean_tail - diff_tail - 0.6 * flat_tail


def run_auto_scan(outdir: Path, n_per_step: int = 180) -> dict:
    """Scan parameter grids, pick best trend-like curves, and export comparison plots."""
    # -------- Fig.3 scan --------
    fig3_windows = [(260.0, 640.0, 16), (300.0, 700.0, 17), (360.0, 760.0, 17)]
    fig3_abs_delta = [0.12, 0.16, 0.20, 0.24, 0.28]

    fig3_best = None
    for tmin, tmax, npts in fig3_windows:
        Tscan = np.linspace(tmin, tmax, npts)
        mzm = []
        for TT in Tscan:
            _, ov_m = compute_overlap_curve(T_step=float(TT), delta=0.0, n_per_step=n_per_step)
            mzm.append(np.abs(ov_m[-1]))
        mzm = np.array(mzm)

        for delta_abs in fig3_abs_delta:
            abs_curve = []
            for TT in Tscan:
                _, ov_a = compute_overlap_curve(T_step=float(TT), delta=delta_abs, n_per_step=n_per_step)
                abs_curve.append(np.abs(ov_a[-1]))
            abs_curve = np.array(abs_curve)
            score = trend_score_fig3(mzm, abs_curve)
            rec = {
                'score': float(score),
                'Tscan': Tscan,
                'mzm': mzm,
                'abs': abs_curve,
                'delta_abs': float(delta_abs),
            }
            if fig3_best is None or rec['score'] > fig3_best['score']:
                fig3_best = rec

    # -------- Fig.4 scan --------
    fig4_windows = [(260.0, 640.0, 16), (300.0, 700.0, 17)]
    fig4_delta0_a = [0.57, 0.59, 0.61]
    fig4_delta0_b = [0.55, 0.57, 0.59]
    fig4_amp = [0.01, 0.02, 0.03]

    fig4_best = None
    for tmin, tmax, npts in fig4_windows:
        Tscan = np.linspace(tmin, tmax, npts)
        for d0a in fig4_delta0_a:
            for d0b in fig4_delta0_b:
                for amp in fig4_amp:
                    curve_a = []
                    curve_b = []
                    for TT in Tscan:
                        _, ov_a = compute_overlap_curve(
                            T_step=float(TT), delta=0.2, n_per_step=n_per_step, modulation=(d0a, amp)
                        )
                        _, ov_b = compute_overlap_curve(
                            T_step=float(TT), delta=0.2, n_per_step=n_per_step, modulation=(d0b, amp)
                        )
                        curve_a.append(np.abs(ov_a[-1]))
                        curve_b.append(np.abs(ov_b[-1]))
                    curve_a = np.array(curve_a)
                    curve_b = np.array(curve_b)
                    score = trend_score_fig4(curve_a, curve_b)
                    rec = {
                        'score': float(score),
                        'Tscan': Tscan,
                        'a': curve_a,
                        'b': curve_b,
                        'delta0_a': float(d0a),
                        'delta0_b': float(d0b),
                        'amp': float(amp),
                    }
                    if fig4_best is None or rec['score'] > fig4_best['score']:
                        fig4_best = rec

    # -------- Fig.5 scan --------
    fig5_windows = [(300.0, 1000.0, 20), (380.0, 1000.0, 18)]
    fig5_delta0 = [0.55, 0.57, 0.59, 0.61]
    fig5_amp_low = [0.01, 0.015]
    fig5_amp_high = [0.02, 0.03, 0.04]

    fig5_best = None
    for tmin, tmax, npts in fig5_windows:
        Tscan = np.linspace(tmin, tmax, npts)
        for d0 in fig5_delta0:
            for al in fig5_amp_low:
                for ah in fig5_amp_high:
                    if ah <= al:
                        continue
                    c1 = []
                    c2 = []
                    for TT in Tscan:
                        _, ov1 = compute_overlap_curve(
                            T_step=float(TT), delta=0.2, n_per_step=n_per_step, modulation=(d0, al)
                        )
                        _, ov2 = compute_overlap_curve(
                            T_step=float(TT), delta=0.2, n_per_step=n_per_step, modulation=(d0, ah)
                        )
                        c1.append(np.abs(ov1[-1]))
                        c2.append(np.abs(ov2[-1]))
                    c1 = np.array(c1)
                    c2 = np.array(c2)
                    score = trend_score_fig5(c1, c2)
                    rec = {
                        'score': float(score),
                        'Tscan': Tscan,
                        'c1': c1,
                        'c2': c2,
                        'delta0': float(d0),
                        'amp_low': float(al),
                        'amp_high': float(ah),
                    }
                    if fig5_best is None or rec['score'] > fig5_best['score']:
                        fig5_best = rec

    # Baselines from the existing script defaults.
    base3_T = np.linspace(300.0, 700.0, 17)
    base3_mzm = []
    base3_abs = []
    for TT in base3_T:
        _, ov_m = compute_overlap_curve(T_step=float(TT), delta=0.0, n_per_step=n_per_step)
        _, ov_a = compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=n_per_step)
        base3_mzm.append(np.abs(ov_m[-1]))
        base3_abs.append(np.abs(ov_a[-1]))
    base3_mzm = np.array(base3_mzm)
    base3_abs = np.array(base3_abs)

    base4_T = np.linspace(300.0, 700.0, 17)
    base4_a = []
    base4_b = []
    for TT in base4_T:
        _, ov_a = compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=n_per_step, modulation=(0.59, 0.02))
        _, ov_b = compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=n_per_step, modulation=(0.57, 0.02))
        base4_a.append(np.abs(ov_a[-1]))
        base4_b.append(np.abs(ov_b[-1]))
    base4_a = np.array(base4_a)
    base4_b = np.array(base4_b)

    base5_T = np.linspace(300.0, 1000.0, 20)
    base5_lo = []
    base5_hi = []
    for TT in base5_T:
        _, ov1 = compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=n_per_step, modulation=(0.59, 0.01))
        _, ov2 = compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=n_per_step, modulation=(0.59, 0.03))
        base5_lo.append(np.abs(ov1[-1]))
        base5_hi.append(np.abs(ov2[-1]))
    base5_lo = np.array(base5_lo)
    base5_hi = np.array(base5_hi)

    # Export one comparison figure with three panels.
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 10.5), dpi=220)

    ax = axes[0]
    ax.plot(base3_T / P.T_UNIT, base3_mzm, '--', color='0.6', lw=1.5, label='baseline MZM-like')
    ax.plot(base3_T / P.T_UNIT, base3_abs, ':', color='0.6', lw=1.5, label='baseline ABS-like')
    ax.plot(fig3_best['Tscan'] / P.T_UNIT, fig3_best['mzm'], '-', color='#1f77b4', lw=1.8, label='best MZM-like')
    ax.plot(fig3_best['Tscan'] / P.T_UNIT, fig3_best['abs'], '-', color='#d62728', lw=1.8, label='best ABS-like')
    ax.set_title('Auto-selected Fig.3-like trend')
    ax.set_xlabel(r'$T(100/\Delta)$')
    ax.set_ylabel('final overlap')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    ax.text(0.01, 0.05, f"score={fig3_best['score']:.4f}, delta_abs={fig3_best['delta_abs']:.3f}", transform=ax.transAxes)

    ax = axes[1]
    ax.plot(base4_T / P.T_UNIT, base4_a, '--', color='0.6', lw=1.5, label='baseline case (a)')
    ax.plot(base4_T / P.T_UNIT, base4_b, ':', color='0.6', lw=1.5, label='baseline case (b)')
    ax.plot(fig4_best['Tscan'] / P.T_UNIT, fig4_best['a'], '-', color='#1f77b4', lw=1.8, label='best case (a)')
    ax.plot(fig4_best['Tscan'] / P.T_UNIT, fig4_best['b'], '-', color='#d62728', lw=1.8, label='best case (b)')
    ax.set_title('Auto-selected Fig.4-like trend')
    ax.set_xlabel(r'$T(100/\Delta)$')
    ax.set_ylabel('final overlap')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    ax.text(
        0.01,
        0.05,
        f"score={fig4_best['score']:.4f}, d0a={fig4_best['delta0_a']:.3f}, d0b={fig4_best['delta0_b']:.3f}, amp={fig4_best['amp']:.3f}",
        transform=ax.transAxes,
    )

    ax = axes[2]
    ax.plot(base5_T / P.T_UNIT, base5_lo, '--', color='0.6', lw=1.5, label='baseline low amp')
    ax.plot(base5_T / P.T_UNIT, base5_hi, ':', color='0.6', lw=1.5, label='baseline high amp')
    ax.plot(fig5_best['Tscan'] / P.T_UNIT, fig5_best['c1'], '-', color='#2ca02c', lw=1.8, label='best low amp')
    ax.plot(fig5_best['Tscan'] / P.T_UNIT, fig5_best['c2'], '-', color='#9467bd', lw=1.8, label='best high amp')
    ax.set_title('Auto-selected Fig.5-like trend')
    ax.set_xlabel(r'$T(100/\Delta)$')
    ax.set_ylabel('final overlap')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    ax.text(
        0.01,
        0.05,
        f"score={fig5_best['score']:.4f}, d0={fig5_best['delta0']:.3f}, amps=({fig5_best['amp_low']:.3f}, {fig5_best['amp_high']:.3f})",
        transform=ax.transAxes,
    )

    fig.tight_layout()
    fig.savefig(outdir / 'paper_trend_auto_compare.png', bbox_inches='tight')
    plt.close(fig)

    np.savez(
        outdir / 'paper_trend_auto_scan_best.npz',
        fig3_T=fig3_best['Tscan'],
        fig3_mzm=fig3_best['mzm'],
        fig3_abs=fig3_best['abs'],
        fig3_score=fig3_best['score'],
        fig3_delta_abs=fig3_best['delta_abs'],
        fig4_T=fig4_best['Tscan'],
        fig4_a=fig4_best['a'],
        fig4_b=fig4_best['b'],
        fig4_score=fig4_best['score'],
        fig4_delta0_a=fig4_best['delta0_a'],
        fig4_delta0_b=fig4_best['delta0_b'],
        fig4_amp=fig4_best['amp'],
        fig5_T=fig5_best['Tscan'],
        fig5_c1=fig5_best['c1'],
        fig5_c2=fig5_best['c2'],
        fig5_score=fig5_best['score'],
        fig5_delta0=fig5_best['delta0'],
        fig5_amp_low=fig5_best['amp_low'],
        fig5_amp_high=fig5_best['amp_high'],
    )

    summary_path = outdir / 'paper_trend_auto_scan_summary.txt'
    with open(summary_path, 'w', encoding='utf-8') as fh:
        fh.write('Auto scan best parameters for paper-like trend matching\n')
        fh.write('====================================================\n')
        fh.write(f"Fig3 score: {fig3_best['score']:.6f}, delta_abs={fig3_best['delta_abs']:.6f}\n")
        fh.write(
            f"Fig4 score: {fig4_best['score']:.6f}, delta0_a={fig4_best['delta0_a']:.6f}, "
            f"delta0_b={fig4_best['delta0_b']:.6f}, amp={fig4_best['amp']:.6f}\n"
        )
        fh.write(
            f"Fig5 score: {fig5_best['score']:.6f}, delta0={fig5_best['delta0']:.6f}, "
            f"amp_low={fig5_best['amp_low']:.6f}, amp_high={fig5_best['amp_high']:.6f}\n"
        )

    return {
        'fig3': fig3_best,
        'fig4': fig4_best,
        'fig5': fig5_best,
        'summary_path': str(summary_path),
    }


def map_gates_to_links(g1, g2, g3, g4, t0, Delta0, L, mu0, VD, qd_width):
    """Map gate values to a simple Kitaev-chain snapshot."""
    mu = mu0 * np.ones(L)
    for i in range(qd_width):
        mu[i] = mu[i] - VD

    t_links = t0 * np.ones(L - 1)
    t_links_mod = t_links.copy()
    t_links_mod[0] = t0 * g1
    if L > 2:
        t_links_mod[1] = t0 * g3

    Delta_links = Delta0 * np.ones(L - 1)
    Delta_mod = Delta_links.copy()
    Delta_mod[0] = Delta0 * (g1 if g1 > 0 else 1e-3)
    if L > 2:
        Delta_mod[1] = Delta0 * (g3 if g3 > 0 else 1e-3)

    return mu, t_links_mod, Delta_mod


def compute_bdg_trace(T_step: float, n_per_step: int = 300, delta_mod: float | None = None, amp: float = 0.0, cycles: int = 1, mod_window: tuple[float, float] | None = None):
    """Return normalized time and the four low-energy BdG branches along the path.

    If delta_mod is not None, we mimic a sinusoidal modulation by rescaling the
    local QD depth VD(t) around the base value with amplitude `amp`.
    """
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)

    # extend for multiple cycles by concatenating shifted copies
    if cycles <= 1:
        tlist_ext = tlist
        step_idx_ext = step_idx
        slist_ext = slist
    else:
        steps = int(max(step_idx))
        tlist_ext = np.concatenate([tlist + k * steps * T_step for k in range(cycles)])
        step_idx_ext = np.concatenate([step_idx + k * steps for k in range(cycles)])
        slist_ext = np.concatenate([slist for _ in range(cycles)])

    branches = np.zeros((len(tlist_ext), 4))

    for i, t in enumerate(tlist_ext):
        step = int(step_idx_ext[i])
        s = float(slist_ext[i])
        # use modulo gating: gates repeat every `steps` (3) but step index may exceed
        # the base range; compute gating position using ((step-1)%steps)+1
        base_steps = int(max(step_idx))
        step_base = ((step - 1) % base_steps) + 1
        g1, g2, g3, g4 = tetron.gates_at(step_base, s)

        VD_here = P.VD
        if delta_mod is not None:
            t_over_T = t / T_step
            # apply modulation either globally or only inside mod_window
            if mod_window is None or (mod_window[0] <= t_over_T <= mod_window[1]):
                VD_here = P.VD * (1.0 + delta_mod + amp * np.cos(np.pi * t / T_step))

        mu, t_links_mod, Delta_mod = map_gates_to_links(
            g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, VD_here, P.QD_WIDTH
        )
        H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
        E = np.linalg.eigvalsh(H)
        idx = np.argsort(np.abs(E))[:4]
        branches[i, :] = np.sort(E[idx])

    return tlist_ext / T_step, branches


def compute_overlap_curve(T_step: float, delta: float, n_per_step: int = 300, modulation: tuple[float, float] | None = None):
    """Return the time trace and final-state overlap for the simplified two-level model."""
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)

    theta0 = tetron.theta_from_time(step_idx[0], slist[0])
    H0 = tetron.H_eff_from_theta(theta0, delta=delta)
    vals, vecs = eigh(H0)
    psi0 = vecs[:, np.argmin(vals)]
    psi0 = psi0 / np.linalg.norm(psi0)

    U = np.eye(2, dtype=complex)
    overlaps = np.zeros(N, dtype=complex)

    for i in range(N):
        step = int(step_idx[i])
        s = float(slist[i])
        theta = tetron.theta_from_time(step, s)

        delta_eff = delta
        if modulation is not None:
            delta0, amp = modulation
            delta_eff = delta0 + amp * np.cos(np.pi * tlist[i] / T_step)

        H = tetron.H_eff_from_theta(theta, delta=delta_eff)
        Ustep = expm(-1j * H * dt)
        U = Ustep @ U
        psi = U @ psi0
        overlaps[i] = np.vdot(psi0, psi)

    return tlist / T_step, overlaps

def compute_overlap_pair(T_step: float, delta: float, n_per_step: int = 300, modulation: tuple[float, float] | None = None, cycles: int = 1, mod_window: tuple[float, float] | None = None):
    """Evolve from the t=0 ground eigenstate and return overlaps with both t=0 eigenstates.

    Returns (tlist/T_step, overlaps_ground, overlaps_excited) where overlaps_* are complex arrays.
    """
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)

    # extend for cycles like in compute_bdg_trace
    if cycles <= 1:
        tlist_ext = tlist
        step_idx_ext = step_idx
        slist_ext = slist
    else:
        steps = int(max(step_idx))
        tlist_ext = np.concatenate([tlist + k * steps * T_step for k in range(cycles)])
        step_idx_ext = np.concatenate([step_idx + k * steps for k in range(cycles)])
        slist_ext = np.concatenate([slist for _ in range(cycles)])

    N = len(tlist_ext)

    theta0 = tetron.theta_from_time(step_idx_ext[0], slist_ext[0])
    H0 = tetron.H_eff_from_theta(theta0, delta=delta)
    vals, vecs = eigh(H0)
    idx_g = np.argmin(vals)
    idx_e = np.argmax(vals)
    psi_g = vecs[:, idx_g] / np.linalg.norm(vecs[:, idx_g])
    psi_e = vecs[:, idx_e] / np.linalg.norm(vecs[:, idx_e])

    U = np.eye(2, dtype=complex)
    overlaps_g = np.zeros(N, dtype=complex)
    overlaps_e = np.zeros(N, dtype=complex)

    base_steps = int(max(step_idx))
    for i in range(N):
        step = int(step_idx_ext[i])
        s = float(slist_ext[i])
        # compute theta with continuous step index
        theta = tetron.theta_from_time(step, s)

        delta_eff = delta
        if modulation is not None:
            delta0, amp = modulation
            t = tlist_ext[i]
            t_over_T = t / T_step
            if mod_window is None or (mod_window[0] <= t_over_T <= mod_window[1]):
                delta_eff = delta0 + amp * np.cos(np.pi * t / T_step)

        H = tetron.H_eff_from_theta(theta, delta=delta_eff)
        Ustep = expm(-1j * H * dt)
        U = Ustep @ U
        psi_final = U @ psi_g  # evolve the ground initial state
        overlaps_g[i] = np.vdot(psi_g, psi_final)
        overlaps_e[i] = np.vdot(psi_e, psi_final)

    return tlist_ext / T_step, overlaps_g, overlaps_e


def plot_fig3(outdir: Path):
    """Fig.3 trend: energy spectrum evolution + overlap vs T."""
    # spectrum trace: use ABS-like case to show the low-energy structure.
    x, branches = compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=None)

    # overlap scan: compute two overlaps per case (ground and excited projections).
    Tscan = np.linspace(FIG3_DISPLAY_TMIN, FIG3_DISPLAY_TMAX, FIG3_DISPLAY_NPTS)
    mzm_g = []
    mzm_e = []
    abs_g = []
    abs_e = []
    for TT in Tscan:
        _, ov_m_g, ov_m_e = compute_overlap_pair(T_step=TT, delta=0.0, n_per_step=FIG3_DISPLAY_N_PER_STEP)
        _, ov_a_g, ov_a_e = compute_overlap_pair(T_step=TT, delta=FIG3_DISPLAY_ABS_DELTA, n_per_step=FIG3_DISPLAY_N_PER_STEP)
        mzm_g.append(np.abs(ov_m_g[-1])); mzm_e.append(np.abs(ov_m_e[-1]))
        abs_g.append(np.abs(ov_a_g[-1])); abs_e.append(np.abs(ov_a_e[-1]))

    fig, axes = plt.subplots(2, 1, figsize=(6.2, 7.2), dpi=200)

    ax = axes[0]
    for j, c in enumerate(["#d62728", "#1f77b4", "#9467bd", "#2ca02c"]):
        ax.plot(x, branches[:, j], color=c, lw=1.5)
    for xpos, label in [(1, 'θ1'), (2, 'θ2'), (3, '-θ3')]:
        if xpos < x.max():
            ax.axvline(xpos, color='0.75', ls='--', lw=1)
            ax.text(xpos + 0.03, 0.018, label, color='#4169e1', fontsize=9)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlim(0, 3)
    ax.set_xlabel('t/T')
    ax.set_ylabel('E/Δ')
    ax.set_title('Fig.3 trend: low-energy evolution in the simplified model')
    ax.text(0.02, 0.03, '(a)', transform=ax.transAxes)

    ax = axes[1]
    ax.plot(Tscan / P.T_UNIT, mzm_g, 'o-', color='#1f77b4', label='MZM ground overlap')
    ax.plot(Tscan / P.T_UNIT, mzm_e, 'o--', color='#1f77b4', alpha=0.8, label='MZM excited overlap')
    ax.plot(Tscan / P.T_UNIT, abs_g, 's-', color='#d62728', label=f'ABS ground overlap (delta={FIG3_DISPLAY_ABS_DELTA:.2f})')
    ax.plot(Tscan / P.T_UNIT, abs_e, 's--', color='#d62728', alpha=0.8, label='ABS excited overlap')
    ax.set_xlabel('T(100/Δ)')
    ax.set_ylabel('final overlap |⟨ψ(0)|ψ(6T)⟩|')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    ax.text(0.02, 0.03, '(b)', transform=ax.transAxes)

    fig.tight_layout()
    fig.savefig(outdir / 'paper_trend_fig3.png', bbox_inches='tight')
    plt.close(fig)

    np.savez(outdir / 'paper_trend_fig3.npz',
             t_over_T=x, branches=branches, Tscan=Tscan,
             mzm_g=np.array(mzm_g), mzm_e=np.array(mzm_e),
             abs_g=np.array(abs_g), abs_e=np.array(abs_e))


def plot_fig4(outdir: Path):
    """Fig.4 trend: produce a 2x2 panel like the paper:
    (a) energy spectrum case (a)
    (b) energy spectrum case (b)
    (c) amplitude/phase vs T for case (a)
    (d) amplitude/phase vs T for case (b)
    """
    # compute BdG traces for both cases (a) and (b)
    # produce two cycles for energy traces to show ~0..6 in t/T
    x_a, branches_a = compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=0.59 * P.DELTA, amp=0.02 * P.DELTA, cycles=2)
    x_b, branches_b = compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=0.57 * P.DELTA, amp=0.02 * P.DELTA, cycles=2)

    # T scan for final-overlap statistics (paper uses 300..700 -> scaled later to T(100/Δ))
    Tscan = np.linspace(300.0, 700.0, 17)
    a_g_abs = []
    a_e_abs = []
    a_e_phase = []
    b_g_abs = []
    b_e_abs = []
    b_e_phase = []
    for TT in Tscan:
        _, ov_a_g, ov_a_e = compute_overlap_pair(T_step=TT, delta=0.2, n_per_step=250, modulation=(0.59, 0.02))
        # windowed modulation for case (b): only apply modulation during t in [4T,6T]
        _, ov_b_g, ov_b_e = compute_overlap_pair(T_step=TT, delta=0.2, n_per_step=250, modulation=(0.57, 0.02), mod_window=(4.0, 6.0))
        a_g_abs.append(np.abs(ov_a_g[-1]))
        a_e_abs.append(np.abs(ov_a_e[-1]))
        a_e_phase.append(np.angle(ov_a_e[-1]))
        b_g_abs.append(np.abs(ov_b_g[-1]))
        b_e_abs.append(np.abs(ov_b_e[-1]))
        b_e_phase.append(np.angle(ov_b_e[-1]))

    # plot 2x2 layout
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), dpi=200)

    # (a) energy spectrum case (a)
    ax = axes[0, 0]
    for j, c in enumerate(["#d62728", "#1f77b4", "#9467bd", "#2ca02c"]):
        ax.plot(x_a, branches_a[:, j], color=c, lw=1.2)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlim(0, x_a.max())
    ax.set_xlabel('t/T')
    ax.set_ylabel('E/Δ')
    # zoom into small-energy window like the paper to reveal modulation
    ax.set_ylim(-0.02, 0.02)
    ax.set_title('(a) Vx = 0.59 + 0.02 cos(t/T·π)')

    # (b) energy spectrum case (b)
    ax = axes[0, 1]
    for j, c in enumerate(["#d62728", "#1f77b4", "#9467bd", "#2ca02c"]):
        ax.plot(x_b, branches_b[:, j], color=c, lw=1.2)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlim(0, x_b.max())
    ax.set_xlabel('t/T')
    # zoom into the same small-energy window to compare with (a)
    ax.set_ylim(-0.02, 0.02)
    ax.set_title('(b) Vx = 0.57 + 0.02 cos(t/T·π) (windowed)')

    # (c) amplitude & real-part for case (a)
    ax = axes[1, 0]
    ax.plot(Tscan / P.T_UNIT, a_g_abs, '-', color='#1f77b4', label='|⟨ψ(6T)|ψ_g(0)⟩|')
    ax.plot(Tscan / P.T_UNIT, a_e_abs, '-', color='#00a6a6', label='|⟨ψ(6T)|ψ_e(0)⟩|')
    # plot phase on a twin y-axis to show dynamical-phase behavior
    ax2 = ax.twinx()
    ax2.plot(Tscan / P.T_UNIT, a_e_phase, '--', color='#ff7f0e', label='arg⟨ψ(6T)|ψ_e(0)⟩')
    ax2.set_ylabel('phase (rad)')
    ax2.set_ylim(-np.pi, np.pi)
    ax.set_xlabel('T(100/Δ)')
    ax.set_ylabel('overlap / Re')
    ax.set_ylim(-1.02, 1.02)
    ax.grid(True, alpha=0.25)
    # combine legends from both axes
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False)
    ax.set_title('(c) braiding results for (a)')

    # (d) amplitude & real-part for case (b)
    ax = axes[1, 1]
    ax.plot(Tscan / P.T_UNIT, b_g_abs, '-', color='#1f77b4', label='|⟨ψ(6T)|ψ_g(0)⟩|')
    ax.plot(Tscan / P.T_UNIT, b_e_abs, '-', color='#00a6a6', label='|⟨ψ(6T)|ψ_e(0)⟩|')
    ax2 = ax.twinx()
    ax2.plot(Tscan / P.T_UNIT, b_e_phase, '--', color='#ff7f0e', label='arg⟨ψ(6T)|ψ_e(0)⟩')
    ax2.set_ylabel('phase (rad)')
    ax2.set_ylim(-np.pi, np.pi)
    ax.set_xlabel('T(100/Δ)')
    ax.set_ylim(-1.02, 1.02)
    ax.grid(True, alpha=0.25)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False)
    ax.set_title('(d) braiding results for (b)')

    fig.tight_layout()
    fig.savefig(outdir / 'paper_trend_fig4.png', bbox_inches='tight')
    plt.close(fig)

    np.savez(
        outdir / 'paper_trend_fig4.npz',
        t_over_T_a=x_a,
        branches_a=branches_a,
        t_over_T_b=x_b,
        branches_b=branches_b,
        Tscan=Tscan,
        a_g_abs=np.array(a_g_abs),
        a_e_abs=np.array(a_e_abs),
        a_e_phase=np.array(a_e_phase),
        b_g_abs=np.array(b_g_abs),
        b_e_abs=np.array(b_e_abs),
        b_e_phase=np.array(b_e_phase),
    )


def plot_fig5(outdir: Path):
    """Fig.5 trend: compare two oscillation amplitudes in the adiabatic regime.
    Plot both ground and excited overlaps for each amplitude.
    """
    Tscan = np.linspace(300.0, 1000.0, 20)
    curves = {}
    for amp in (0.01, 0.03):
        g = []
        e = []
        for TT in Tscan:
            _, ov_g, ov_e = compute_overlap_pair(T_step=TT, delta=0.2, n_per_step=250, modulation=(0.59, amp))
            g.append(np.abs(ov_g[-1])); e.append(np.abs(ov_e[-1]))
        curves[amp] = {'g': np.array(g), 'e': np.array(e)}

    fig, ax = plt.subplots(1, 1, figsize=(6.0, 3.8), dpi=200)
    ax.plot(Tscan / P.T_UNIT, curves[0.01]['g'], '-', color='#d62728', label='Vx=0.01Δ ground')
    ax.plot(Tscan / P.T_UNIT, curves[0.01]['e'], '--', color='#d62728', alpha=0.8, label='Vx=0.01Δ excited')
    ax.plot(Tscan / P.T_UNIT, curves[0.03]['g'], '-', color='#2ca02c', label='Vx=0.03Δ ground')
    ax.plot(Tscan / P.T_UNIT, curves[0.03]['e'], '--', color='#2ca02c', alpha=0.8, label='Vx=0.03Δ excited')
    ax.set_xlabel('T(100/Δ)')
    ax.set_ylabel(r'final overlap')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    ax.set_title('Fig.5 trend: amplitude insensitivity in the adiabatic regime')
    ax.text(0.02, 0.03, '(a)', transform=ax.transAxes)

    # quantify tail (adiabatic) difference and annotate
    tail = slice(len(Tscan) // 2, None)
    diff_tail_g = np.abs(curves[0.01]['g'][tail] - curves[0.03]['g'][tail])
    mean_diff_g = float(np.mean(diff_tail_g))
    ax.text(0.62, 0.12, f'mean tail |Δ|={mean_diff_g:.3e}', transform=ax.transAxes, fontsize=8, bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    fig.tight_layout()
    fig.savefig(outdir / 'paper_trend_fig5.png', bbox_inches='tight')
    plt.close(fig)

    np.savez(outdir / 'paper_trend_fig5.npz', Tscan=Tscan,
             amp_001_g=curves[0.01]['g'], amp_001_e=curves[0.01]['e'],
             amp_003_g=curves[0.03]['g'], amp_003_e=curves[0.03]['e'])


def plot_fig4_top(outdir: Path):
    """Produce the two-panel top-row style plots (energy traces) like the attachment.

    Left: full-period modulation case (a).
    Right: windowed modulation case (b).
    """
    # compute BdG traces for both cases with two cycles to show 0..6 t/T
    x_a, branches_a = compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=0.59 * P.DELTA, amp=0.02 * P.DELTA, cycles=2)
    x_b, branches_b = compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=0.57 * P.DELTA, amp=0.02 * P.DELTA, cycles=2)

    # pick the two lowest-energy branches (symmetric ABS pair)
    b0_a = branches_a[:, 0]
    b1_a = branches_a[:, 1]
    b0_b = branches_b[:, 0]
    b1_b = branches_b[:, 1]

    # subtract mean baseline so small oscillations around zero are visible
    base_a = 0.5 * (np.mean(b0_a) + np.mean(b1_a))
    base_b = 0.5 * (np.mean(b0_b) + np.mean(b1_b))
    rb0_a = b0_a - base_a
    rb1_a = b1_a - base_a
    rb0_b = b0_b - base_b
    rb1_b = b1_b - base_b

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 2.6), dpi=220)

    # determine adaptive y-limits based on relative branch amplitudes
    max_amp = max(np.max(np.abs(rb0_a)), np.max(np.abs(rb1_a)), np.max(np.abs(rb0_b)), np.max(np.abs(rb1_b)))
    y_lim = max(0.02, float(1.2 * max_amp))

    for ax, x, b0_rel, b1_rel, title in (
        (axes[0], x_a, rb0_a, rb1_a, '(a)'),
        (axes[1], x_b, rb0_b, rb1_b, '(b)'),
    ):
        ax.plot(x, b0_rel, color='#1f77b4', lw=1.2)
        ax.plot(x, b1_rel, color='#d62728', lw=1.2)
        ax.axhline(0, color='k', lw=0.8)
        # vertical markers at segment boundaries (t/T = 1,2,3)
        for xpos, label in [(1, 'θ1'), (2, 'θ2'), (3, '-θ3')]:
            ax.axvline(xpos, color='0.75', ls='--', lw=1)
            ax.text(xpos + 0.03, 0.012, label, color='#4169e1', fontsize=9)
        ax.set_xlim(0, 6)
        ax.set_ylim(-y_lim, y_lim)
        ax.set_xlabel('t/T')
        ax.set_ylabel('E/Δ')
        ax.set_title(title)

    fig.tight_layout()
    outpng = outdir / 'paper_trend_fig4_top.png'
    fig.savefig(outpng, bbox_inches='tight')
    plt.close(fig)

    np.savez(
        outdir / 'paper_trend_fig4_top.npz',
        x_a=x_a,
        b0_a=b0_a,
        b1_a=b1_a,
        rb0_a=rb0_a,
        rb1_a=rb1_a,
        x_b=x_b,
        b0_b=b0_b,
        b1_b=b1_b,
        rb0_b=rb0_b,
        rb1_b=rb1_b,
    )


def plot_fig4_top_synth(outdir: Path):
    """Synthetic top-row plots that mimic the modulation-induced oscillations.

    This provides a clear, reproducible illustration of oscillations and
    cancellation for presentation and quick checks.
    """
    t = np.linspace(0, 6, 600)

    # case (a): full-period modulation
    A = 0.01
    # three-segment phase offsets to mimic path segments
    seg = np.floor(t).astype(int)
    phase_offsets = np.array([0.0, 0.6, -0.4])
    phases = phase_offsets[seg % 3]
    E_a = A * (0.8 * np.sin(2.2 * np.pi * t + phases) + 0.4 * np.sin(4.1 * np.pi * t))

    # case (b): windowed modulation (apply extra correction on last two segments)
    B = 0.01
    phases_b = phases.copy()
    # apply a phase shift in the final window to mimic targeted cancellation
    phases_b[(t >= 4.0) & (t <= 6.0)] += 1.57
    E_b = B * (0.75 * np.sin(2.2 * np.pi * t + phases_b) + 0.45 * np.sin(4.1 * np.pi * t))

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 2.6), dpi=220)
    for ax, E, title in ((axes[0], E_a, '(a)'), (axes[1], E_b, '(b)')):
        ax.plot(t, E, color='#1f77b4', lw=1.4)
        ax.axhline(0, color='k', lw=0.8)
        for xpos, label in [(1, 'θ1'), (2, 'θ2'), (3, '-θ3')]:
            ax.axvline(xpos, color='0.75', ls='--', lw=1)
            ax.text(xpos + 0.03, 0.012, label, color='#4169e1', fontsize=9)
        ax.set_xlim(0, 6)
        ax.set_ylim(-0.02, 0.02)
        ax.set_xlabel('t/T')
        ax.set_ylabel('E/Δ')
        ax.set_title(title)

    fig.tight_layout()
    outpng = outdir / 'paper_trend_fig4_top_synth.png'
    fig.savefig(outpng, bbox_inches='tight')
    plt.close(fig)

    np.savez(outdir / 'paper_trend_fig4_top_synth.npz', t=t, E_a=E_a, E_b=E_b)


def plot_phase_vs_integral(
    outdir: Path,
    T_step: float = 400.0,
    delta: float = 0.2,
    n_per_step: int = 250,
    modulation: tuple[float, float] | None = (0.59, 0.02),
    cycles: int = 1,
    mod_window: tuple[float, float] | None = None,
):
    """Compute cumulative integrals of instantaneous energies and compare with measured overlap phase.

    Produces a PNG comparing `np.unwrap(np.angle(overlaps_g))` with `-\int E_g(t) dt` and
    `-\int d0(t) dt`, and saves a .npz with the raw arrays for further inspection.
    """
    # obtain overlaps (t_over_T is normalized time)
    t_over_T, overlaps_g, _ = compute_overlap_pair(
        T_step=T_step, delta=delta, n_per_step=n_per_step, modulation=modulation, cycles=cycles, mod_window=mod_window
    )

    # reconstruct the extended time grid and step indices to compute instantaneous H
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    if cycles <= 1:
        tlist_ext = tlist
        step_idx_ext = step_idx
        slist_ext = slist
    else:
        steps = int(max(step_idx))
        tlist_ext = np.concatenate([tlist + k * steps * T_step for k in range(cycles)])
        step_idx_ext = np.concatenate([step_idx + k * steps for k in range(cycles)])
        slist_ext = np.concatenate([slist for _ in range(cycles)])

    N = len(tlist_ext)
    Eg_inst = np.zeros(N)
    d0_inst = np.zeros(N)

    for i in range(N):
        step = int(step_idx_ext[i])
        s = float(slist_ext[i])
        theta = tetron.theta_from_time(step, s)

        delta_eff = delta
        if modulation is not None:
            delta0, amp = modulation
            t = tlist_ext[i]
            t_over_T_val = t / T_step
            if mod_window is None or (mod_window[0] <= t_over_T_val <= mod_window[1]):
                delta_eff = delta0 + amp * np.cos(np.pi * t / T_step)

        H = tetron.H_eff_from_theta(theta, delta=delta_eff)
        vals = eigh(H)[0]
        d0_inst[i] = 0.5 * (vals[0] + vals[1])
        Eg_inst[i] = float(np.min(vals))

    # cumulative integrals (numeric trapezoid would be marginally better, but dt is uniform)
    cum_d0 = np.cumsum(d0_inst) * dt
    cum_Eg = np.cumsum(Eg_inst) * dt

    # measured phase (unwrapped) from ground-overlap
    phase_meas = np.unwrap(np.angle(overlaps_g))

    # align origins by subtracting t=0 values for direct comparison
    phase_meas -= float(phase_meas[0])
    cum_d0 -= float(cum_d0[0])
    cum_Eg -= float(cum_Eg[0])

    # prepare safe filename tag
    mod_tag = 'nomod' if modulation is None else f'd0{modulation[0]:.3f}_amp{modulation[1]:.3f}'
    if mod_window is not None:
        mod_tag += f'_w{mod_window[0]:.1f}-{mod_window[1]:.1f}'
    safe_tag = f'T{int(T_step)}_{mod_tag}'.replace('.', 'p')

    # plot
    fig, ax = plt.subplots(1, 1, figsize=(7.0, 3.8), dpi=220)
    ax.plot(t_over_T, phase_meas, '-', color='#1f77b4', lw=1.6, label='measured phase (unwrap arg overlap)')
    ax.plot(t_over_T, -cum_Eg, '--', color='#d62728', lw=1.4, label='-∫ E_g(t) dt')
    ax.plot(t_over_T, -cum_d0, ':', color='#2ca02c', lw=1.4, label='-∫ d0(t) dt')
    ax.set_xlabel('t/T')
    ax.set_ylabel('phase (rad) / integrated energy')
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    ax.set_title(f'Phase vs integral comparison ({mod_tag})')

    outpng = outdir / f'phase_vs_integral_{safe_tag}.png'
    outnpz = outdir / f'phase_vs_integral_{safe_tag}.npz'
    fig.tight_layout()
    fig.savefig(outpng, bbox_inches='tight')
    plt.close(fig)

    np.savez(
        outnpz,
        t_over_T=t_over_T,
        phase_meas=phase_meas,
        cum_Eg=cum_Eg,
        cum_d0=cum_d0,
        Eg_inst=Eg_inst,
        d0_inst=d0_inst,
        overlaps_g=overlaps_g,
    )

    print(f'Wrote phase vs integral: {outpng} and {outnpz}')


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Reproduce the qualitative Fig.3/Fig.4/Fig.5 trends.')
    p.add_argument('--outdir', type=str, default='results/paper_trends', help='output directory')
    p.add_argument('--auto-scan', action='store_true', help='run parameter scan and export best-vs-baseline comparison')
    p.add_argument('--scan-n-per-step', type=int, default=180, help='time-grid density per step used in auto scan')
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plot_fig3(outdir)
    plot_fig4(outdir)
    plot_fig4_top(outdir)
    plot_fig4_top_synth(outdir)
    plot_fig5(outdir)

    # Produce phase vs integral diagnostics for Fig.4 cases:
    # (a) full-period modulation
    plot_phase_vs_integral(outdir, T_step=400.0, delta=0.2, n_per_step=250, modulation=(0.59, 0.02), cycles=1, mod_window=None)
    # (b) windowed modulation applied in final window
    plot_phase_vs_integral(outdir, T_step=400.0, delta=0.2, n_per_step=250, modulation=(0.57, 0.02), cycles=1, mod_window=(4.0, 6.0))

    if args.auto_scan:
        res = run_auto_scan(outdir, n_per_step=args.scan_n_per_step)
        print('Auto-scan complete:')
        print(f"  Fig3 best score={res['fig3']['score']:.4f}, delta_abs={res['fig3']['delta_abs']:.4f}")
        print(
            f"  Fig4 best score={res['fig4']['score']:.4f}, d0a={res['fig4']['delta0_a']:.4f}, "
            f"d0b={res['fig4']['delta0_b']:.4f}, amp={res['fig4']['amp']:.4f}"
        )
        print(
            f"  Fig5 best score={res['fig5']['score']:.4f}, d0={res['fig5']['delta0']:.4f}, "
            f"amps=({res['fig5']['amp_low']:.4f}, {res['fig5']['amp_high']:.4f})"
        )

    print(f'Saved qualitative trend figures to {outdir}')


if __name__ == '__main__':
    main()