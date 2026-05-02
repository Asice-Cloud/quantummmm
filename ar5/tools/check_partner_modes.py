#!/usr/bin/env python3
"""Check partner (opposite-end) near-zero modes for top candidates.

Loads candidates, diagonalizes BdG for N=200, finds two smallest-|E| modes,
computes max_site, max_weight, maj_sim for each, and saves report+plots.
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

from verify_mzm import map_c_to_params, build_kitaev_bdg


def analyze_candidate(cxy, cyx, baseline=0.05, N=200):
    c_xx = float(baseline)
    c_yy = float(baseline)
    c_xy = complex(cxy)
    c_yx = complex(cyx)
    c_zz = 0.0
    c_z0 = 0.0
    c_0z = 0.0
    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    H = build_kitaev_bdg(N, t, Delta, mu)
    w, vecs = eigh(H)
    order = np.argsort(np.abs(w))
    modes = []
    for idx in order[:4]:
        val = float(w[idx])
        vec = vecs[:, idx]
        u = vec[:N]
        vpart = vec[N:]
        weights = (np.abs(u)**2 + np.abs(vpart)**2)
        max_site = int(np.argmax(weights))
        max_weight = float(np.max(weights))
        num = np.linalg.norm(u - np.conjugate(vpart))
        denom = np.linalg.norm(u) + np.linalg.norm(vpart) + 1e-12
        maj_norm = float(num / denom)
        maj_sim = float(1.0 - maj_norm)
        modes.append({'E': val, 'max_site': max_site, 'max_weight': max_weight, 'maj_sim': maj_sim, 'weights': weights.tolist()})
    return modes


def main(candidates_path='results/abs_scan_multi_candidates.json', top_k=10):
    os.makedirs('results', exist_ok=True)
    with open(candidates_path, 'r') as f:
        cand = json.load(f)
    cand_sorted = sorted(cand, key=lambda x: x.get('min_abs', 1e9))
    top = cand_sorted[:top_k]
    report = []
    for i, c in enumerate(top, start=1):
        cxy = float(c['cxy'])
        cyx = float(c['cyx'])
        print(f'Analyzing top {i}: cxy={cxy}, cyx={cyx}')
        modes = analyze_candidate(cxy, cyx, baseline=0.05, N=200)
        report.append({'index': i, 'cxy': cxy, 'cyx': cyx, 'modes': modes[:2]})

        # plot weights of first two modes
        fig, ax = plt.subplots(2, 1, figsize=(8,4), sharex=True)
        for j in range(2):
            weights = np.array(modes[j]['weights'])
            ax[j].plot(weights)
            ax[j].set_ylabel(f'mode{j+1} |psi|^2 E={modes[j]["E"]:.2e}')
        ax[-1].set_xlabel('site')
        fig.suptitle(f'top{i} cxy={cxy} cyx={cyx}')
        out = f'results/abs_partner_top{i}_modes.png'
        fig.tight_layout()
        fig.savefig(out)
        plt.close(fig)

    outjson = 'results/abs_partner_report.json'
    with open(outjson, 'w') as f:
        json.dump(report, f, indent=2)
    print('Partner analysis saved to', outjson)


if __name__ == '__main__':
    main()
