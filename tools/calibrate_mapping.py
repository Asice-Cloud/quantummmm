#!/usr/bin/env python3
"""Calibrate mapping constants A0 and C0 by fitting BdG low-energy splitting.

Procedure:
- Measure lowest positive BdG eigenvalue E(t_left) for L=2 chain while varying t_left.
- Fit E = A0 * t_left (linear) to estimate A0.
- Measure E(mu_diff) for L=2 chain while varying onsite potential difference
  between site0 and site1; fit E = C0 * mu_diff to estimate C0.

Results are saved to `results/mapping_calibration.npz` and diagnostic plots saved
to `results/`.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

# make repo root importable
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools.embed_kitaev import build_bdg


def lowest_positive_energy(H, tol=1e-8):
    E = eigh(H, eigvals_only=True)
    # take positive eigenvalues
    Epos = np.array([e for e in E if e > tol])
    if len(Epos) == 0:
        return 0.0
    return float(np.min(Epos))


def measure_vs_t_left(t_list, L=2, Delta0=0.3, mu0=0.0):
    Evals = []
    for t_left in t_list:
        mu = mu0 * np.ones(L)
        t_links = np.zeros(L - 1)
        t_links[0] = t_left
        Delta_links = Delta0 * np.ones(L - 1)
        H = build_bdg(mu, t_links, Delta_links)
        E = lowest_positive_energy(H)
        Evals.append(E)
    return np.array(Evals)


def measure_vs_mu_diff(mu_diff_list, L=2, t0=0.3, Delta0=0.3):
    Evals = []
    for mu_diff in mu_diff_list:
        # set mu0 = +mu_diff/2 on site0, -mu_diff/2 on site1
        mu = np.zeros(L)
        mu[0] = mu_diff / 2.0
        mu[1] = -mu_diff / 2.0
        t_links = t0 * np.ones(L - 1)
        Delta_links = Delta0 * np.ones(L - 1)
        H = build_bdg(mu, t_links, Delta_links)
        E = lowest_positive_energy(H)
        Evals.append(E)
    return np.array(Evals)


def fit_slope(x, y):
    # linear fit y = m x (+ b). Prefer forced-through-origin fit if intercept small.
    m, b = np.polyfit(x, y, 1)
    # also slope through origin
    m0 = np.sum(x * y) / np.sum(x * x)
    return {'m': float(m), 'b': float(b), 'm0': float(m0)}


def main():
    os.makedirs('results', exist_ok=True)

    # A0 calibration
    t_list = np.linspace(0.05, 1.0, 10)
    E_t = measure_vs_t_left(t_list, L=2, Delta0=0.3, mu0=0.0)
    fit_t = fit_slope(t_list, E_t)

    plt.figure()
    plt.plot(t_list, E_t, 'o', label='data')
    m_val = fit_t['m']
    b_val = fit_t['b']
    m0_val = fit_t['m0']
    plt.plot(t_list, m_val * t_list + b_val, '-', label=f'fit m={m_val:.3f}, b={b_val:.3e}')
    plt.plot(t_list, m0_val * t_list, '--', label=f'origin fit m0={m0_val:.3f}')

    plt.xlabel('t_left')
    plt.ylabel('lowest positive BdG E')
    plt.legend()
    out1 = 'results/calibrate_A0.png'
    plt.title('Calibration of A0 (E vs t_left)')
    plt.savefig(out1)
    plt.close()

    # C0 calibration
    mu_diff_list = np.linspace(0.05, 1.0, 10)
    E_mu = measure_vs_mu_diff(mu_diff_list, L=2, t0=0.1, Delta0=0.3)
    fit_mu = fit_slope(mu_diff_list, E_mu)

    plt.figure()
    plt.plot(mu_diff_list, E_mu, 'o', label='data')
    m_val2 = fit_mu['m']
    b_val2 = fit_mu['b']
    m0_val2 = fit_mu['m0']
    plt.plot(mu_diff_list, m_val2 * mu_diff_list + b_val2, '-', label=f'fit m={m_val2:.3f}, b={b_val2:.3e}')
    plt.plot(mu_diff_list, m0_val2 * mu_diff_list, '--', label=f'origin fit m0={m0_val2:.3f}')
    plt.xlabel('mu_diff (mu1 - mu2)')
    plt.ylabel('lowest positive BdG E')
    plt.legend()
    out2 = 'results/calibrate_C0.png'
    plt.title('Calibration of C0 (E vs mu_diff)')
    plt.savefig(out2)
    plt.close()

    np.savez('results/mapping_calibration.npz', t_list=t_list, E_t=E_t, fit_t=fit_t,
             mu_diff_list=mu_diff_list, E_mu=E_mu, fit_mu=fit_mu)
    print('Saved results/mapping_calibration.npz and plots')


if __name__ == '__main__':
    main()
