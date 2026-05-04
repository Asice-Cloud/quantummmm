#!/usr/bin/env python3
"""Optimize mapping constants A0 and C0 to match BdG energies along the tetron path.

Minimize: sum_t ( sqrt( (A0*S(t))^2 + (C0*M(t))^2 ) - E_bdg(t) )^2
where S(t)=t_left(t)+t_right(t), M(t)=mu1(t)-mu2(t), and E_bdg(t) is the
lowest absolute BdG eigenenergy at time t computed from the full chain.
"""
import os
import sys
import numpy as np
from scipy.optimize import least_squares

# ensure repo root importable
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import embed_kitaev
from tools import tetron_path_sim as tetron
from tools.reproduce_figs import map_gates_to_links
import tools.paper_params as P


def compute_bdg_minE(T_step=200.0, n_per_step=200):
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)
    E_bdg = np.zeros(N)
    S = np.zeros(N)
    M = np.zeros(N)

    for i in range(N):
        step = int(step_idx[i])
        s = float(slist[i])
        g1, g2, g3, g4 = tetron.gates_at(step, s)
        mu, t_links_mod, Delta_mod = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, P.VD, P.QD_WIDTH)
        # BdG Hamiltonian
        H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
        E = np.linalg.eigvalsh(H)
        E_bdg[i] = np.min(np.abs(E))
        # summary statistics for mapping
        t_left = t_links_mod[0] if len(t_links_mod) > 0 else 0.0
        t_right = t_links_mod[1] if len(t_links_mod) > 1 else 0.0
        S[i] = t_left + t_right
        M[i] = mu[0] - mu[1] if len(mu) > 1 else 0.0

    return tlist, E_bdg, S, M


def residuals(x, S, M, E_bdg):
    A0, C0 = x
    pred = np.sqrt((A0 * S) ** 2 + (C0 * M) ** 2)
    return pred - E_bdg


def main():
    os.makedirs('results', exist_ok=True)
    # select a representative T (paper uses e.g., T=100..500); use P.FIG3_T
    T_step = P.FIG3_T
    n_per_step = 300
    print('Computing BdG min energies along path (this may take a while)...')
    tlist, E_bdg, S, M = compute_bdg_minE(T_step=T_step, n_per_step=n_per_step)
    print('Done computing BdG energies. Running least-squares fit for A0,C0...')

    # initial guess from earlier calibration
    x0 = np.array([0.6, 0.16])
    res = least_squares(residuals, x0, args=(S, M, E_bdg), bounds=(0, np.inf), xtol=1e-8, ftol=1e-8)

    A0_fit, C0_fit = res.x
    print(f'Fit results: A0={A0_fit:.6f}, C0={C0_fit:.6f}')

    np.savez('results/mapping_fit.npz', tlist=tlist, E_bdg=E_bdg, S=S, M=M, A0_fit=A0_fit, C0_fit=C0_fit, res=res.x)
    print('Saved results/mapping_fit.npz')


if __name__ == '__main__':
    main()
