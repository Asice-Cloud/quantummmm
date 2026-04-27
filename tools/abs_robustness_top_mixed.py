#!/usr/bin/env python3
"""Robustness and scaling for top-k mixed+baseline candidates.

Loads results/abs_scan_mixed_with_baseline_candidates.json, picks top-k by
min_abs, computes BdG for multiple N, performs small random perturbations
around (cxy,cyx) and saves plots and a JSON report in results/.
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

from verify_mzm import map_c_to_params, build_kitaev_bdg


def load_candidates(path):
    with open(path, 'r') as f:
        return json.load(f)


def compute_min_mode(cxy, cyx, baseline=0.05, N=200):
    c_xx = float(baseline)
    c_yy = float(baseline)
    c_xy = complex(cxy)
    c_yx = complex(cyx)
    c_zz = 0.0
    c_z0 = 0.0
    c_0z = 0.0
    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    H = build_kitaev_bdg(N, t, Delta, mu)
    w, v = eigh(H)
    idx = int(np.argmin(np.abs(w)))
    min_abs = float(np.abs(w[idx]))
    vec = v[:, idx]
    u = vec[:N]
    vpart = vec[N:]
    weights = (np.abs(u)**2 + np.abs(vpart)**2)
    max_site = int(np.argmax(weights))
    max_weight = float(max(weights))
    num = np.linalg.norm(u - np.conjugate(vpart))
    denom = np.linalg.norm(u) + np.linalg.norm(vpart) + 1e-12
    maj_norm = float(num / denom)
    maj_sim = float(1.0 - maj_norm)
    return dict(min_abs=min_abs, max_site=max_site, max_weight=max_weight, maj_sim=maj_sim, weights=weights.tolist())


def main(top_k=3, candidates_path='results/abs_scan_mixed_with_baseline_candidates.json'):
    os.makedirs('results', exist_ok=True)
    cand = load_candidates(candidates_path)
    if not cand:
        print('No candidates found at', candidates_path)
        return
    cand_sorted = sorted(cand, key=lambda x: x.get('min_abs', 1e9))
    top = cand_sorted[:top_k]

    report = []
    Ns = [100, 200, 400]
    perturb_N = 200
    perturb_samples = 100
    perturb_sigma = 0.01

    for i, c in enumerate(top, start=1):
        cxy = float(c['cxy'])
        cyx = float(c['cyx'])
        print(f'Processing top {i}: cxy={cxy}, cyx={cyx}')
        entry = {'cxy': cxy, 'cyx': cyx, 'per_N': {}, 'perturb': {}}

        # N scaling
        mins = []
        for N in Ns:
            res = compute_min_mode(cxy, cyx, baseline=0.05, N=N)
            entry['per_N'][str(N)] = res
            mins.append(res['min_abs'])

        # plot min_abs vs N
        fig = plt.figure(figsize=(5,3))
        plt.plot(Ns, mins, 'o-')
        plt.yscale('log')
        plt.xlabel('N')
        plt.ylabel('min |E|')
        plt.title(f'top{i} min|E| vs N')
        out1 = f'results/abs_top{i}_minE_vs_N.png'
        fig.tight_layout()
        fig.savefig(out1)
        plt.close(fig)

        # LDOS for N = perturb_N
        res200 = entry['per_N'][str(perturb_N)]
        weights = np.array(res200['weights'])
        fig = plt.figure(figsize=(8,3))
        plt.plot(weights)
        plt.xlabel('Site')
        plt.ylabel('LDOS')
        plt.title(f'top{i} LDOS (N={perturb_N}) site {res200["max_site"]}')
        out2 = f'results/abs_top{i}_ldos_N{perturb_N}.png'
        fig.tight_layout()
        fig.savefig(out2)
        plt.close(fig)

        # Perturbation robustness: random gaussian noise on cxy/cyx
        pert_mins = []
        for s in range(perturb_samples):
            dcxy = np.random.normal(scale=perturb_sigma)
            dcyx = np.random.normal(scale=perturb_sigma)
            res = compute_min_mode(cxy + dcxy, cyx + dcyx, baseline=0.05, N=perturb_N)
            pert_mins.append(res['min_abs'])
        entry['perturb']['mins'] = pert_mins

        # plot perturbation histogram
        fig = plt.figure(figsize=(5,3))
        plt.hist(pert_mins, bins=30)
        plt.xlabel('min |E|')
        plt.ylabel('counts')
        plt.xscale('log')
        plt.title(f'top{i} perturbation min|E| histogram')
        out3 = f'results/abs_top{i}_perturb_hist.png'
        fig.tight_layout()
        fig.savefig(out3)
        plt.close(fig)

        report.append(entry)

    out_json = 'results/abs_robustness_top_mixed_report.json'
    with open(out_json, 'w') as f:
        json.dump(report, f, indent=2)

    print('Robustness analysis complete. Report:', out_json)


if __name__ == '__main__':
    main()
