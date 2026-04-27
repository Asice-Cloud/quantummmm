#!/usr/bin/env python3
"""Targeted mu_well scan for top candidates.

Loads results/abs_scan_multi_candidates.json, picks top_k, and for each scans
mu_site values placing a local mu at site 0 (mapping c_z0=c_0z=-mu/4), computes
min|E|, IPR, LDOS, maj_sim, and saves report and plots under results/.
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

from verify_mzm import map_c_to_params, build_kitaev_bdg


def analyze_with_mu(cxy, cyx, mu_site, baseline=0.05, N=200):
    c_xx = float(baseline)
    c_yy = float(baseline)
    c_xy = complex(cxy)
    c_yx = complex(cyx)
    c_zz = 0.0
    # mapping used elsewhere: c_z0 = c_0z = -mu/4
    v = -float(mu_site) / 4.0
    c_z0 = v
    c_0z = v
    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    H = build_kitaev_bdg(N, t, Delta, mu)
    w, vecs = eigh(H)
    idx = int(np.argmin(np.abs(w)))
    min_abs = float(np.abs(w[idx]))
    vec = vecs[:, idx]
    u = vec[:N]
    vpart = vec[N:]
    weights = np.abs(u)**2 + np.abs(vpart)**2
    max_site = int(np.argmax(weights))
    max_weight = float(np.max(weights))
    # IPR
    ipr = float(np.sum(weights**2))
    num = np.linalg.norm(u - np.conjugate(vpart))
    denom = np.linalg.norm(u) + np.linalg.norm(vpart) + 1e-12
    maj_norm = float(num / denom)
    maj_sim = float(1.0 - maj_norm)
    return dict(min_abs=min_abs, max_site=max_site, max_weight=max_weight, ipr=ipr, maj_sim=maj_sim, weights=weights.tolist())


def main(candidates_path='results/abs_scan_multi_candidates.json', top_k=10):
    os.makedirs('results', exist_ok=True)
    with open(candidates_path, 'r') as f:
        cand = json.load(f)
    cand_sorted = sorted(cand, key=lambda x: x.get('min_abs', 1e9))
    top = cand_sorted[:top_k]

    mu_sites = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    report = []
    for i, c in enumerate(top, start=1):
        cxy = float(c['cxy'])
        cyx = float(c['cyx'])
        print(f'Scanning top{i}: cxy={cxy}, cyx={cyx}')
        entry = {'index': i, 'cxy': cxy, 'cyx': cyx, 'mu_scan': {}}
        for mu in mu_sites:
            res = analyze_with_mu(cxy, cyx, mu, baseline=0.05, N=200)
            entry['mu_scan'][str(mu)] = res

        # Save per-candidate JSON and plots
        outj = f'results/abs_top{i}_mu_scan.json'
        with open(outj, 'w') as fj:
            json.dump(entry, fj, indent=2)

        # plot min_abs vs mu
        mus = np.array([float(k) for k in entry['mu_scan'].keys()])
        mins = np.array([entry['mu_scan'][str(m)]['min_abs'] for m in mus])
        fig = plt.figure(figsize=(5,3))
        plt.plot(mus, mins, 'o-')
        plt.yscale('log')
        plt.xlabel('mu_site')
        plt.ylabel('min |E|')
        plt.title(f'top{i} min|E| vs mu_site')
        fig.tight_layout()
        fig.savefig(f'results/abs_top{i}_minE_vs_mu.png')
        plt.close(fig)

        # plot LDOS for best mu (smallest min_abs)
        best_mu = mus[np.argmin(mins)]
        best = entry['mu_scan'][str(float(best_mu))]
        weights = np.array(best['weights'])
        fig = plt.figure(figsize=(8,3))
        plt.plot(weights)
        plt.xlabel('site')
        plt.ylabel('LDOS')
        plt.title(f'top{i} best mu={best_mu} site {best["max_site"]}')
        fig.tight_layout()
        fig.savefig(f'results/abs_top{i}_best_mu_ldos.png')
        plt.close(fig)

        report.append(entry)

    out = 'results/abs_targeted_mu_scan_report.json'
    with open(out, 'w') as f:
        json.dump(report, f, indent=2)
    print('Targeted mu_well scan complete. Report:', out)


if __name__ == '__main__':
    main()
