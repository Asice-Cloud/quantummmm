#!/usr/bin/env python3
"""Ultra-high-clarity Fig.4 rendering: larger figure, thicker lines, clearer labels.
Saves both PNG and PDF to results/paper_trends for highest quality.
"""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import tools.reproduce_trend_figs as rt
import tools.paper_params as P


def render(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    # use same data generation as before
    x_a, branches_a = rt.compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=0.59 * P.DELTA, amp=0.02 * P.DELTA, cycles=2)
    x_b, branches_b = rt.compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=0.57 * P.DELTA, amp=0.02 * P.DELTA, cycles=2)

    Tscan = np.linspace(300.0, 700.0, 17)
    a_g_abs = []
    a_e_abs = []
    a_e_phase = []
    b_g_abs = []
    b_e_abs = []
    b_e_phase = []
    for TT in Tscan:
        _, ov_a_g, ov_a_e = rt.compute_overlap_pair(T_step=TT, delta=0.2, n_per_step=250, modulation=(0.59, 0.02))
        _, ov_b_g, ov_b_e = rt.compute_overlap_pair(T_step=TT, delta=0.2, n_per_step=250, modulation=(0.57, 0.02), mod_window=(4.0, 6.0))
        a_g_abs.append(np.abs(ov_a_g[-1]))
        a_e_abs.append(np.abs(ov_a_e[-1]))
        a_e_phase.append(np.angle(ov_a_e[-1]))
        b_g_abs.append(np.abs(ov_b_g[-1]))
        b_e_abs.append(np.abs(ov_b_e[-1]))
        b_e_phase.append(np.angle(ov_b_e[-1]))

    # pick two lowest |E| branches per time
    def pick_two(branches):
        abs_br = np.abs(branches)
        idx = np.argsort(abs_br, axis=1)[:, :2]
        N = branches.shape[0]
        e1 = np.zeros(N)
        e2 = np.zeros(N)
        for i in range(N):
            e1[i] = branches[i, idx[i, 0]]
            e2[i] = branches[i, idx[i, 1]]
        return e1, e2

    e1_a, e2_a = pick_two(branches_a)
    e1_b, e2_b = pick_two(branches_b)

    # baseline subtraction
    base_a = 0.5 * (np.mean(e1_a) + np.mean(e2_a))
    base_b = 0.5 * (np.mean(e1_b) + np.mean(e2_b))
    rb1_a = e1_a - base_a
    rb2_a = e2_a - base_a
    rb1_b = e1_b - base_b
    rb2_b = e2_b - base_b

    # high-clarity plotting params
    figsize = (14, 10)
    dpi = 600
    lw = 2.6
    fontsz = 18
    ticksz = 12
    legend_font = 14

    plt.rcParams.update({'font.size': fontsz})

    fig, axes = plt.subplots(2, 2, figsize=figsize, dpi=dpi)

    # (a)
    ax = axes[0, 0]
    ax.plot(x_a, rb1_a, color='#1f77b4', lw=lw)
    ax.plot(x_a, rb2_a, color='#d62728', lw=lw)
    ax.axhline(0, color='k', lw=1.0)
    for xpos, label in [(1, 'θ1'), (2, 'θ2'), (3, '-θ3')]:
        ax.axvline(xpos, color='0.5', ls='--', lw=1.2)
        ax.text(xpos + 0.03, 0.012, label, color='#1f77b4', fontsize=14)
    ax.set_xlim(0, 6)
    ax.set_ylim(-0.02, 0.02)
    ax.set_xlabel('t/T', fontsize=fontsz)
    ax.set_ylabel('E/Δ', fontsize=fontsz)
    ax.set_title('(a) Vx = 0.59 + 0.02 cos(t/T · π)', fontsize=fontsz)
    ax.tick_params(labelsize=ticksz)

    # (b)
    ax = axes[0, 1]
    ax.plot(x_b, rb1_b, color='#1f77b4', lw=lw)
    ax.plot(x_b, rb2_b, color='#d62728', lw=lw)
    ax.axhline(0, color='k', lw=1.0)
    for xpos, label in [(1, 'θ1'), (2, 'θ2'), (3, '-θ3')]:
        ax.axvline(xpos, color='0.5', ls='--', lw=1.2)
        ax.text(xpos + 0.03, 0.012, label, color='#1f77b4', fontsize=14)
    ax.set_xlim(0, 6)
    ax.set_ylim(-0.02, 0.02)
    ax.set_xlabel('t/T', fontsize=fontsz)
    ax.set_title('(b) Vx = 0.57 + 0.02 cos(t/T · π) (windowed)', fontsize=fontsz)
    ax.tick_params(labelsize=ticksz)

    # (c)
    ax = axes[1, 0]
    ax.plot(Tscan / P.T_UNIT, a_g_abs, '-', color='#1f77b4', lw=lw, label='|⟨ψ(6T)|ψ_g(0)⟩|')
    ax.plot(Tscan / P.T_UNIT, a_e_abs, '-', color='#00a6a6', lw=lw, label='|⟨ψ(6T)|ψ_e(0)⟩|')
    ax.set_xlabel('T (100/Δ)', fontsize=fontsz)
    ax.set_ylabel('overlap magnitude', fontsize=fontsz)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.tick_params(labelsize=ticksz)

    ax2 = ax.twinx()
    ax2.plot(Tscan / P.T_UNIT, a_e_phase, '--', color='#ff7f0e', lw=lw, label='phase (excited)')
    ax2.set_ylabel('phase (rad)', fontsize=fontsz)
    ax2.set_ylim(-np.pi, np.pi)
    ax2.tick_params(labelsize=ticksz)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=legend_font)
    ax.set_title('(c) braiding results for (a)', fontsize=fontsz)

    # (d)
    ax = axes[1, 1]
    ax.plot(Tscan / P.T_UNIT, b_g_abs, '-', color='#1f77b4', lw=lw, label='|⟨ψ(6T)|ψ_g(0)⟩|')
    ax.plot(Tscan / P.T_UNIT, b_e_abs, '-', color='#00a6a6', lw=lw, label='|⟨ψ(6T)|ψ_e(0)⟩|')
    ax.set_xlabel('T (100/Δ)', fontsize=fontsz)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.tick_params(labelsize=ticksz)

    ax2 = ax.twinx()
    ax2.plot(Tscan / P.T_UNIT, b_e_phase, '--', color='#ff7f0e', lw=lw, label='phase (excited)')
    ax2.set_ylabel('phase (rad)', fontsize=fontsz)
    ax2.set_ylim(-np.pi, np.pi)
    ax2.tick_params(labelsize=ticksz)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=legend_font)
    ax.set_title('(d) braiding results for (b)', fontsize=fontsz)

    # add big panel letters
    axes[0,0].text(-0.08, 1.05, '(a)', transform=axes[0,0].transAxes, fontsize=22, weight='bold')
    axes[0,1].text(-0.08, 1.05, '(b)', transform=axes[0,1].transAxes, fontsize=22, weight='bold')
    axes[1,0].text(-0.08, 1.05, '(c)', transform=axes[1,0].transAxes, fontsize=22, weight='bold')
    axes[1,1].text(-0.08, 1.05, '(d)', transform=axes[1,1].transAxes, fontsize=22, weight='bold')

    fig.tight_layout(pad=2.0)
    outpng = outdir / 'paper_trend_fig4_ultraclear.png'
    outpdf = outdir / 'paper_trend_fig4_ultraclear.pdf'
    fig.savefig(outpng, bbox_inches='tight')
    fig.savefig(outpdf, bbox_inches='tight')
    plt.close(fig)
    print('Wrote', outpng, outpdf)


if __name__ == '__main__':
    render(Path('results/paper_trends'))
