#!/usr/bin/env python3
"""Build a centered-interface BdG chain from XYZ demo parameters and save spectrum+density.

Usage: run from project root: python3 tools/run_centered_xyz_demo.py
"""
import os
import sys
import numpy as np
from numpy.linalg import eigh

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from tools.pulse_abs import build_chain_from_params


def run(N=200, a1=0.0, b1=1.0, c1=None, d1=0.5, rho=1.0, interface_width=4):
    if c1 is None:
        c1 = a1
    t_val = b1 / rho
    Delta_val = d1 / rho
    mu_val = 0.0

    # left: paired, right: normal
    t_left = t_val
    Delta_left = Delta_val
    mu_left = mu_val
    t_right = t_val
    Delta_right = 0.0
    mu_right = mu_val

    os.makedirs('results', exist_ok=True)

    H, t_hop, Delta_bond, mu_site = build_chain_from_params(N, t_left, Delta_left, mu_left, t_right, Delta_right, mu_right, interface_width=interface_width)
    vals, vecs = eigh(H)

    # save spectrum plot
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6,4))
        plt.plot(np.sort(vals), '.', markersize=3)
        plt.title('Centered interface spectrum')
        plt.tight_layout()
        plt.savefig('results/xyz_centered_spectrum.png')
        plt.close()
    except Exception:
        pass

    # lowest-by-abs eigenstate
    idx = np.argsort(np.abs(vals))
    low_idx = idx[0]
    vec0 = vecs[:, low_idx]
    Nloc = vec0.shape[0] // 2
    u = vec0[:Nloc]
    v = vec0[Nloc:]
    dens = np.abs(u)**2 + np.abs(v)**2
    np.savetxt('results/xyz_centered_density0.txt', dens)

    peak = int(np.argmax(dens))
    peak_val = float(dens[peak])
    right = peak
    while right < N and dens[right] >= peak_val/np.e:
        right += 1
    xi_est = (right - peak) if right < N else float('nan')

    with open('results/xyz_centered_report.txt', 'w') as f:
        f.write(f'params: a1=c1={a1}, b1={b1}, d1={d1}, rho={rho}\n')
        f.write(f't={t_val}, Delta_left={Delta_left}, Delta_right={Delta_right}, mu=0\n')
        f.write('\nLowest energies (abs-sorted first 6):\n')
        for i in idx[:6]:
            f.write(f'{vals[i]:.6e}\n')
        f.write(f'\npeak_idx={peak}\npeak_val={peak_val:.6e}\nxi_est={xi_est}\n')

    print('Saved results/xyz_centered_spectrum.png, results/xyz_centered_density0.txt, results/xyz_centered_report.txt')


if __name__ == '__main__':
    run()
