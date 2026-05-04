#!/usr/bin/env python3
"""Scan chain length L and measure lowest BdG eigenvalue splitting.

Compute the smallest positive eigenvalue (or splitting between lowest
pair) as a function of L and fit to an exponential decay to estimate
localization length.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# make repo importable
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import embed_kitaev
import tools.paper_params as P


def lowest_positive_eig(H, tol=1e-12):
    E = np.linalg.eigvalsh(H)
    Epos = np.array([e for e in E if e > tol])
    if len(Epos) == 0:
        return 0.0
    return float(np.min(Epos))


def exp_fit(L, a, xi):
    return a * np.exp(-L / xi)


def main():
    os.makedirs('results', exist_ok=True)
    L_list = np.array([40, 60, 80, 100, 140, 200])
    Es = []
    for L in L_list:
        # build simple chain with QD on left
        mu = P.mu0 * np.ones(L)
        t_links = P.t0 * np.ones(L - 1)
        Delta_links = P.DELTA * np.ones(L - 1)
        # QD depth VD at left qd_width
        for i in range(min(P.QD_WIDTH, L)):
            mu[i] = mu[i] - P.VD
        H = embed_kitaev.build_bdg(mu, t_links, Delta_links)
        Epos = lowest_positive_eig(H)
        Es.append(Epos)
        print(f'L={L}, E_minpos={Epos:.6e}')

    Es = np.array(Es)
    # fit log-linear or exponential
    try:
        popt, pcov = curve_fit(exp_fit, L_list, Es, p0=(1e-1, 20.0))
        a_fit, xi_fit = popt
        print(f'Exponential fit: a={a_fit:.3e}, xi={xi_fit:.3f}')
    except Exception as e:
        print('Exponential fit failed:', e)
        a_fit, xi_fit = np.nan, np.nan

    np.savez('results/eig_splitting_vs_L.npz', L_list=L_list, Es=Es, a_fit=a_fit, xi_fit=xi_fit)

    plt.figure()
    plt.semilogy(L_list, Es, 'o', label='data')
    if not np.isnan(a_fit):
        Lfine = np.linspace(L_list.min(), L_list.max(), 200)
        plt.semilogy(Lfine, exp_fit(Lfine, a_fit, xi_fit), '-', label=f'fit xi={xi_fit:.2f}')
    plt.xlabel('L')
    plt.ylabel('lowest positive BdG E')
    plt.title('Zero-mode splitting vs chain length')
    plt.legend()
    plt.tight_layout()
    out = 'results/eig_splitting_vs_L.png'
    plt.savefig(out)
    print('Saved', out)


if __name__ == '__main__':
    main()
