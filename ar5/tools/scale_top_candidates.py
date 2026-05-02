#!/usr/bin/env python3
"""Compute finite-size scaling for top candidates (min|E| vs N) and fit exponential.

Reads `results/abs_partner_report.json` for candidate list, runs Ns and saves plots and JSON report.
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

from verify_mzm import map_c_to_params, build_kitaev_bdg


def compute_min_energy(cxy, cyx, baseline=0.05, N=200):
    c_xx = float(baseline)
    c_yy = float(baseline)
    c_xy = complex(cxy)
    c_yx = complex(cyx)
    c_zz = 0.0
    c_z0 = 0.0
    c_0z = 0.0
    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    H = build_kitaev_bdg(N, t, Delta, mu)
    w, _ = eigh(H)
    return float(np.min(np.abs(w)))


def fit_exponential(Ns, vals):
    # ignore zeros or extremely small values
    Ns = np.array(Ns)
    vals = np.array(vals)
    mask = vals > 0
    if np.sum(mask) < 2:
        return None
    x = Ns[mask]
    y = np.log(vals[mask])
    A = np.polyfit(x, y, 1)
    slope, intercept = A[0], A[1]
    xi = -1.0 / slope if slope != 0 else None
    A0 = np.exp(intercept)
    return {'slope': float(slope), 'intercept': float(intercept), 'xi': xi if xi is None else float(xi), 'A0': float(A0)}


def main(partner_report='results/abs_partner_report.json', outjson='results/abs_scaling_report.json', top_k=6):
    os.makedirs('results', exist_ok=True)
    with open(partner_report, 'r') as f:
        cand = json.load(f)
    top = cand[:top_k]
    Ns = [100, 200, 400, 800]
    report = []
    for c in top:
        i = c['index']
        cxy = c['cxy']; cyx = c['cyx']
        print(f'Scaling top{i}: cxy={cxy}, cyx={cyx}')
        vals = []
        for N in Ns:
            val = compute_min_energy(cxy, cyx, baseline=0.05, N=N)
            vals.append(val)
        fit = fit_exponential(Ns, vals)
        report.append({'index': i, 'cxy': cxy, 'cyx': cyx, 'Ns': Ns, 'mins': vals, 'fit': fit})

        # plot
        fig = plt.figure(figsize=(5,3))
        plt.plot(Ns, vals, 'o-')
        if fit is not None:
            xs = np.linspace(min(Ns), max(Ns), 100)
            ys = np.exp(fit['intercept'] + fit['slope'] * xs)
            plt.plot(xs, ys, '--', label=f"xi={fit['xi']:.2f}")
            plt.legend()
        plt.yscale('log')
        plt.xlabel('N')
        plt.ylabel('min |E|')
        plt.title(f'top{i} scaling')
        fig.tight_layout()
        fig.savefig(f'results/abs_top{i}_scaling.png')
        plt.close(fig)

    with open(outjson, 'w') as f:
        json.dump(report, f, indent=2)
    print('Scaling report saved to', outjson)


if __name__ == '__main__':
    main()
