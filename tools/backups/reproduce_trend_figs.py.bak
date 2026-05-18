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


def compute_bdg_trace(T_step: float, n_per_step: int = 300, delta_mod: float | None = None, amp: float = 0.0):
    """Return normalized time and the four low-energy BdG branches along the path.

    If delta_mod is not None, we mimic a sinusoidal modulation by rescaling the
    local QD depth VD(t) around the base value with amplitude `amp`.
    """
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    branches = np.zeros((len(tlist), 4))

    for i, t in enumerate(tlist):
        step = int(step_idx[i])
        s = float(slist[i])
        g1, g2, g3, g4 = tetron.gates_at(step, s)

        VD_here = P.VD
        if delta_mod is not None:
            # Paper-like sinusoidal modulation: keep the x-axis label and trend,
            # while using the simplified QD-depth proxy in the toy BdG model.
            VD_here = P.VD * (1.0 + delta_mod + amp * np.cos(np.pi * t / T_step))

        mu, t_links_mod, Delta_mod = map_gates_to_links(
            g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, VD_here, P.QD_WIDTH
        )
        H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
        E = np.linalg.eigvalsh(H)
        idx = np.argsort(np.abs(E))[:4]
        branches[i, :] = np.sort(E[idx])

    return tlist / T_step, branches


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


def plot_fig3(outdir: Path):
    """Fig.3 trend: energy spectrum evolution + overlap vs T."""
    # spectrum trace: use ABS-like case to show the low-energy structure.
    x, branches = compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=None)

    # overlap scan: use the clearer auto-selected display window by default.
    Tscan = np.linspace(FIG3_DISPLAY_TMIN, FIG3_DISPLAY_TMAX, FIG3_DISPLAY_NPTS)
    mzm_overlaps = []
    abs_overlaps = []
    for TT in Tscan:
        _, ov_m = compute_overlap_curve(T_step=TT, delta=0.0, n_per_step=FIG3_DISPLAY_N_PER_STEP)
        _, ov_a = compute_overlap_curve(T_step=TT, delta=FIG3_DISPLAY_ABS_DELTA, n_per_step=FIG3_DISPLAY_N_PER_STEP)
        mzm_overlaps.append(np.abs(ov_m[-1]))
        abs_overlaps.append(np.abs(ov_a[-1]))

    fig, axes = plt.subplots(2, 1, figsize=(6.2, 7.2), dpi=200)

    ax = axes[0]
    for j, c in enumerate(["#d62728", "#1f77b4", "#9467bd", "#2ca02c"]):
        ax.plot(x, branches[:, j], color=c, lw=1.5)
    for xpos, label in [(1, r'$\theta_1$'), (2, r'$\theta_2$'), (3, r'$-\theta_3$')]:
        if xpos < x.max():
            ax.axvline(xpos, color='0.75', ls='--', lw=1)
            ax.text(xpos + 0.03, 0.018, label, color='#4169e1', fontsize=9)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlim(0, 3)
    ax.set_xlabel(r'$t/T$')
    ax.set_ylabel(r'$E/\Delta$')
    ax.set_title('Fig.3 trend: low-energy evolution in the simplified model')
    ax.text(0.02, 0.03, '(a)', transform=ax.transAxes)

    ax = axes[1]
    ax.plot(Tscan / P.T_UNIT, mzm_overlaps, 'o-', color='#1f77b4', label='MZM-like')
    ax.plot(Tscan / P.T_UNIT, abs_overlaps, 's--', color='#d62728', label=f'ABS-like (delta={FIG3_DISPLAY_ABS_DELTA:.2f})')
    ax.set_xlabel(r'$T(100/\Delta)$')
    ax.set_ylabel(r'final overlap $|\langle\psi(0)|\psi(6T)\rangle|$')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    ax.text(0.02, 0.03, '(b)', transform=ax.transAxes)

    fig.tight_layout()
    fig.savefig(outdir / 'paper_trend_fig3.png', bbox_inches='tight')
    plt.close(fig)

    np.savez(outdir / 'paper_trend_fig3.npz', t_over_T=x, branches=branches, Tscan=Tscan, mzm_overlaps=mzm_overlaps, abs_overlaps=abs_overlaps)


def plot_fig4(outdir: Path):
    """Fig.4 trend: sinusoidal modulation and time-independence after phase cancellation."""
    x, branches_a = compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=0.59 * P.DELTA, amp=0.02 * P.DELTA)
    x2, branches_b = compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=0.57 * P.DELTA, amp=0.02 * P.DELTA)

    Tscan = np.linspace(300.0, 700.0, 17)
    overlaps_a = []
    overlaps_b = []
    for TT in Tscan:
        _, ov_a = compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=250, modulation=(0.59, 0.02))
        _, ov_b = compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=250, modulation=(0.57, 0.02))
        overlaps_a.append(np.abs(ov_a[-1]))
        overlaps_b.append(np.abs(ov_b[-1]))

    fig, axes = plt.subplots(2, 1, figsize=(6.2, 7.2), dpi=200)

    ax = axes[0]
    for j, c in enumerate(["#d62728", "#1f77b4", "#9467bd", "#2ca02c"]):
        ax.plot(x, branches_a[:, j], color=c, lw=1.4)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlim(0, 3)
    ax.set_xlabel(r'$t/T$')
    ax.set_ylabel(r'$E/\Delta$')
    ax.set_title(r'Fig.4 trend: modulation-cancelled evolution (case a)')
    ax.text(0.02, 0.03, '(a)', transform=ax.transAxes)

    ax = axes[1]
    ax.plot(Tscan / P.T_UNIT, overlaps_a, 'o-', color='#1f77b4', label=r'case (a): $0.59+0.02\cos$')
    ax.plot(Tscan / P.T_UNIT, overlaps_b, 's--', color='#d62728', label=r'case (b): $0.57+0.02\cos$')
    ax.set_xlabel(r'$T(100/\Delta)$')
    ax.set_ylabel(r'final overlap')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    ax.text(0.02, 0.03, '(b)', transform=ax.transAxes)

    fig.tight_layout()
    fig.savefig(outdir / 'paper_trend_fig4.png', bbox_inches='tight')
    plt.close(fig)

    np.savez(outdir / 'paper_trend_fig4.npz', t_over_T=x, branches=branches_a, Tscan=Tscan, overlaps_a=overlaps_a, overlaps_b=overlaps_b)


def plot_fig5(outdir: Path):
    """Fig.5 trend: compare two oscillation amplitudes in the adiabatic regime."""
    Tscan = np.linspace(300.0, 1000.0, 20)
    curves = {}
    for amp in (0.01, 0.03):
        ovs = []
        for TT in Tscan:
            _, ov = compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=250, modulation=(0.59, amp))
            ovs.append(np.abs(ov[-1]))
        curves[amp] = np.array(ovs)

    fig, ax = plt.subplots(1, 1, figsize=(6.0, 3.8), dpi=200)
    ax.plot(Tscan / P.T_UNIT, curves[0.01], '-', color='#d62728', label=r'$V_{x1}=0.01\Delta$')
    ax.plot(Tscan / P.T_UNIT, curves[0.03], '--', color='#2ca02c', label=r'$V_{x1}=0.03\Delta$')
    ax.set_xlabel(r'$T(100/\Delta)$')
    ax.set_ylabel(r'final overlap')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    ax.set_title('Fig.5 trend: amplitude insensitivity in the adiabatic regime')
    ax.text(0.02, 0.03, '(a)', transform=ax.transAxes)

    fig.tight_layout()
    fig.savefig(outdir / 'paper_trend_fig5.png', bbox_inches='tight')
    plt.close(fig)

    np.savez(outdir / 'paper_trend_fig5.npz', Tscan=Tscan, amp_001=curves[0.01], amp_003=curves[0.03])


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
    plot_fig5(outdir)

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