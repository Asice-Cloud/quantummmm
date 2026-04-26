#!/usr/bin/env python3
"""Run a simple XYZ-derived BdG demo and save spectrum + report.
Example: a1=c1=0, b1=1.0 -> t=1.0; d1=0.5 -> Delta=0.5 on left, right Delta=0.
"""
import os
import sys
import numpy as np
from numpy.linalg import eigh

# ensure repo root is on sys.path so `import tools.*` works when running from project root
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from tools.pulse_abs import build_chain_from_params, plot_spectrum


def run_demo():
    os.makedirs('results', exist_ok=True)
    rho = 1.0
    a1 = 0.0
    c1 = a1
    b1 = 1.0
    d1 = 0.5
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

    N = 200
    H, t_hop, Delta_bond, mu_site = build_chain_from_params(N, t_left, Delta_left, mu_left, t_right, Delta_right, mu_right, interface_width=4)
    vals, vecs = eigh(H)
    vals_sorted = np.sort(np.real_if_close(vals))
    plot_spectrum(vals_sorted, 'results/xyz_demo_spectrum.png', title='XYZ demo S-N spectrum')

    # find lowest-energy eigenstates by absolute value
    idx = np.argsort(np.abs(vals))
    low_idx = idx[:6]

    # compute densities for lowest state
    def density_from_vec(vec):
        N2 = vec.shape[0]
        Nloc = N2 // 2
        u = vec[:Nloc]
        v = vec[Nloc:]
        dens = np.abs(u)**2 + np.abs(v)**2
        return dens

    with open('results/xyz_demo_report.txt', 'w') as f:
        f.write(f'params: a1=c1={a1}, b1={b1}, d1={d1}, rho={rho}\n')
        f.write(f't={t_val}, Delta_left={Delta_left}, Delta_right={Delta_right}, mu=0\n')
        f.write('\nLowest energies (by value):\n')
        for i in low_idx:
            f.write(f'{vals[i]:.6e}\n')
        f.write('\nDensities for lowest eigenstate:\n')
        dens0 = density_from_vec(vecs[:, low_idx[0]])
        # save density sampled
        np.savetxt('results/xyz_demo_density0.txt', dens0)
        f.write('density0 saved to results/xyz_demo_density0.txt\n')

    print('Saved results/xyz_demo_spectrum.png and results/xyz_demo_report.txt')


if __name__ == '__main__':
    run_demo()
