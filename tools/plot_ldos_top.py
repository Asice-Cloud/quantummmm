#!/usr/bin/env python3
"""Plot LDOS for top-N candidates from validated scan.

Saves PNG files to results/ldos_top{i}.png and a combined figure results/ldos_top_combined.png
"""
import os
import json
import numpy as np
import importlib.util
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# load build_bdg_from_params from scan_all_mixtures
here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location('scan_all_mixtures', os.path.join(here, 'scan_all_mixtures.py'))
scan_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan_mod)
build_bdg_from_params = scan_mod.build_bdg_from_params


def compute_ldos_for_candidate(cand, N=120):
    t = complex(cand['t'])
    Delta = complex(cand['Delta'])
    mu = complex(cand['mu'])
    H = build_bdg_from_params(N, t, Delta, mu)
    vals, vecs = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))
    v0 = vecs[:, idx[0]]
    half = H.shape[0]//2
    u = v0[:half]
    v = v0[half:]
    ldos = np.abs(u)**2 + np.abs(v)**2
    return ldos, vals[idx[0]]


def main(topk=6, N=120):
    os.makedirs('results', exist_ok=True)
    data = json.load(open('results/scan_all_mixtures_validated.json'))
    data_sorted = sorted(data, key=lambda d: d.get('min_abs', 1e9))
    chosen = data_sorted[:topk]
    figs = []
    for i, cand in enumerate(chosen, start=1):
        ldos, ev = compute_ldos_for_candidate(cand, N=N)
        fig, ax = plt.subplots(figsize=(8,3))
        ax.plot(np.arange(len(ldos)), ldos, '-o', markersize=3)
        ax.set_title(f"top{i}: combo={'+'.join(cand.get('combo',[]))}, E={ev:.2e}")
        ax.set_xlabel('site')
        ax.set_ylabel('LDOS')
        plt.tight_layout()
        out = f'results/ldos_top{i}.png'
        fig.savefig(out)
        plt.close(fig)
        figs.append(out)

    # combined figure
    n = len(figs)
    fig, axs = plt.subplots(n, 1, figsize=(8, 2.5*n))
    if n==1:
        axs = [axs]
    for j, cand in enumerate(chosen):
        ldos, ev = compute_ldos_for_candidate(cand, N=N)
        axs[j].plot(np.arange(len(ldos)), ldos, '-o', markersize=3)
        axs[j].set_title(f"top{j+1}: { '+'.join(cand.get('combo',[])) }, E={ev:.2e}")
        axs[j].set_ylabel('LDOS')
    axs[-1].set_xlabel('site')
    plt.tight_layout()
    fig.savefig('results/ldos_top_combined.png')
    print('Saved LDOS plots for top', topk)


if __name__ == '__main__':
    main()
