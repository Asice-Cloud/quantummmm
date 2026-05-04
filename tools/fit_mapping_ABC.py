#!/usr/bin/env python3
"""Joint fit for mapping constants A0, B0, C0 using BdG min energies.

Fits: E_pred(t) = sqrt( (A0*S(t))^2 + (B0*sin(theta(t)))^2 + (C0*M(t))^2 )
to the lowest-absolute BdG eigenenergy along the tetron path.

S(t) = t_left + t_right
M(t) = mu1 - mu2
theta(t) from `tools.tetron_path_sim.theta_from_time`.
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


def compute_bdg_minE_with_theta(T_step=200.0, n_per_step=200):
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)
    E_bdg = np.zeros(N)
    S = np.zeros(N)
    M = np.zeros(N)
    Theta = np.zeros(N)

    for i in range(N):
        step = int(step_idx[i])
        s = float(slist[i])
        theta = tetron.theta_from_time(step, s)
        Theta[i] = theta
        g1, g2, g3, g4 = tetron.gates_at(step, s)
        mu, t_links_mod, Delta_mod = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, P.VD, P.QD_WIDTH)
        H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
        E = np.linalg.eigvalsh(H)
        E_bdg[i] = np.min(np.abs(E))
        t_left = t_links_mod[0] if len(t_links_mod) > 0 else 0.0
        t_right = t_links_mod[1] if len(t_links_mod) > 1 else 0.0
        S[i] = t_left + t_right
        M[i] = mu[0] - mu[1] if len(mu) > 1 else 0.0

    return tlist, E_bdg, S, M, Theta, dt


def residuals(x, S, M, Theta, E_bdg):
    A0, B0, C0 = x
    pred = np.sqrt((A0 * S) ** 2 + (B0 * np.sin(Theta)) ** 2 + (C0 * M) ** 2)
    return pred - E_bdg


def main():
    os.makedirs('results', exist_ok=True)
    T_step = P.FIG3_T
    n_per_step = max(200, P.FIG3_N_PER_STEP)
    print('Computing BdG min energies along path (this may take a while)...')
    tlist, E_bdg, S, M, Theta, dt = compute_bdg_minE_with_theta(T_step=T_step, n_per_step=n_per_step)
    print('Done computing BdG energies. Running least-squares fit for A0,B0,C0...')

    # initial guess: use previous A0,C0 and small B0
    x0 = np.array([0.00265, 0.005, 0.16])
    bounds = (0.0, np.inf)
    res = least_squares(residuals, x0, args=(S, M, Theta, E_bdg), bounds=bounds, xtol=1e-8, ftol=1e-8)

    A0_fit, B0_fit, C0_fit = res.x
    print(f'Fit results: A0={A0_fit:.6g}, B0={B0_fit:.6g}, C0={C0_fit:.6g}')

    np.savez('results/mapping_fit_ABC.npz', tlist=tlist, E_bdg=E_bdg, S=S, M=M, Theta=Theta,
             A0_fit=A0_fit, B0_fit=B0_fit, C0_fit=C0_fit, res=res.x)
    print('Saved results/mapping_fit_ABC.npz')


if __name__ == '__main__':
    main()
