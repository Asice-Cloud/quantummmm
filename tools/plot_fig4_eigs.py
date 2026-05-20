#!/usr/bin/env python3
"""Plot Fig.4-like instantaneous low-energy eigenvalues vs t/T for two modulation cases.

Generates a 1x2 figure: (a) full-period modulation, (b) windowed modulation in [4,6] (normalized t/T).
Saves PNG to results/paper_trends/fig4_eigs_comparison.png
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


def main():
    outdir = Path('results/paper_trends')
    outdir.mkdir(parents=True, exist_ok=True)

    T_step = 400.0
    n_per_step = 300

    # case (a): full-period modulation
    t_a, branches_a = rt.compute_bdg_trace(T_step=T_step, n_per_step=n_per_step, delta_mod=0.59 * P.DELTA, amp=0.02 * P.DELTA, cycles=1, mod_window=None)

    # case (b): windowed modulation applied in final window [4,6]
    t_b, branches_b = rt.compute_bdg_trace(T_step=T_step, n_per_step=n_per_step, delta_mod=0.57 * P.DELTA, amp=0.02 * P.DELTA, cycles=1, mod_window=(4.0, 6.0))

    # pick two lowest |E| branches (sorted per time by absolute value)
    def pick_two(branches):
        # branches is (N,4) sorted ascending by abs then sorted by value
        # we'll pick the two with smallest absolute value at each time
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

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), dpi=220)
    ax = axes[0]
    ax.plot(t_a, e1_a, '-', color='#1f77b4', lw=1.6, label='branch 1')
    ax.plot(t_a, e2_a, '-', color='#d62728', lw=1.6, label='branch 2')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlabel('t/T')
    ax.set_ylabel('E/Δ')
    ax.set_title('(a) full-period modulation')
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(t_b, e1_b, '-', color='#1f77b4', lw=1.6, label='branch 1')
    ax.plot(t_b, e2_b, '-', color='#d62728', lw=1.6, label='branch 2')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlabel('t/T')
    ax.set_title('(b) windowed modulation (t/T in [4,6])')
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)

    fig.tight_layout()
    outpng = outdir / 'fig4_eigs_comparison.png'
    fig.savefig(outpng, bbox_inches='tight')
    plt.close(fig)
    print('Wrote', outpng)


if __name__ == '__main__':
    main()
