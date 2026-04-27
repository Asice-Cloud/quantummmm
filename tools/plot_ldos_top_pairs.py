#!/usr/bin/env python3
"""Plot LDOS for top-N candidates: show mode1, mode2, and combined LDOS.

Saves files: results/ldos_top{i}_mode1.png, _mode2.png, _sum.png and combined multi-panel results/ldos_top_pairs_combined.png
"""
import os
import json
import numpy as np
import importlib.util
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location('scan_all_mixtures', os.path.join(here, 'scan_all_mixtures.py'))
scan_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan_mod)
build_bdg_from_params = scan_mod.build_bdg_from_params


def ldos_from_vec(vec):
    half = vec.shape[0]//2
    u = vec[:half]
    v = vec[half:]
    return np.abs(u)**2 + np.abs(v)**2


def compute_modes_ldos(cand, N=120):
    t = complex(cand['t'])
    Delta = complex(cand['Delta'])
    mu = complex(cand['mu'])
    H = build_bdg_from_params(N, t, Delta, mu)
    vals, vecs = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))
    v1 = vecs[:, idx[0]]
    v2 = vecs[:, idx[1]]
    l1 = ldos_from_vec(v1)
    l2 = ldos_from_vec(v2)
    return l1, l2, vals[idx[0]], vals[idx[1]]


def main(topk=6, N=120):
    os.makedirs('results', exist_ok=True)
    data = json.load(open('results/scan_all_mixtures_validated.json'))
    data_sorted = sorted(data, key=lambda d: d.get('min_abs', 1e9))
    chosen = data_sorted[:topk]
    combined_fig, axs = plt.subplots(topk, 3, figsize=(12, 2.5*topk))
    if topk==1:
        axs = np.expand_dims(axs, 0)
    for i, cand in enumerate(chosen, start=1):
        l1, l2, e1, e2 = compute_modes_ldos(cand, N=N)
        x = np.arange(len(l1))
        # individual plots
        fig1, ax1 = plt.subplots(figsize=(8,2.5))
        ax1.plot(x, l1, '-o', markersize=3)
        ax1.set_title(f'top{i} mode1 E={e1:.2e}')
        ax1.set_xlabel('site')
        ax1.set_ylabel('LDOS')
        fig1.tight_layout()
        fig1.savefig(f'results/ldos_top{i}_mode1.png')
        plt.close(fig1)

        fig2, ax2 = plt.subplots(figsize=(8,2.5))
        ax2.plot(x, l2, '-o', markersize=3)
        ax2.set_title(f'top{i} mode2 E={e2:.2e}')
        ax2.set_xlabel('site')
        ax2.set_ylabel('LDOS')
        fig2.tight_layout()
        fig2.savefig(f'results/ldos_top{i}_mode2.png')
        plt.close(fig2)

        fig3, ax3 = plt.subplots(figsize=(8,2.5))
        s = l1 + l2
        ax3.plot(x, s, '-o', markersize=3)
        ax3.set_title(f'top{i} combined (mode1+mode2)')
        ax3.set_xlabel('site')
        ax3.set_ylabel('LDOS')
        fig3.tight_layout()
        fig3.savefig(f'results/ldos_top{i}_sum.png')
        plt.close(fig3)

        # add to combined figure
        axs[i-1,0].plot(x, l1)
        axs[i-1,0].set_title(f'top{i} m1 E={e1:.2e}')
        axs[i-1,1].plot(x, l2)
        axs[i-1,1].set_title(f'top{i} m2 E={e2:.2e}')
        axs[i-1,2].plot(x, s)
        axs[i-1,2].set_title(f'top{i} sum')

    for axcol in axs.flatten():
        axcol.set_xlim(0, N-1)
    combined_fig.tight_layout()
    combined_fig.savefig('results/ldos_top_pairs_combined.png')
    print('Saved paired LDOS plots for top', topk)


if __name__ == '__main__':
    main()
