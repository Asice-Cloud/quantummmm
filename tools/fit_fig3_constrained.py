#!/usr/bin/env python3
"""Constrained per-figure fit for Fig.3: fit A0,B0,C0,ts with C0 <= 0.5.

Saves results to `results/mapping_fit_fig3_constrained_C0le0.5.npz`.
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


def compute_bdg_minE_with_theta(T_step=100.0, n_per_step=300):
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
    A0, B0, C0, ts = x
    pred = np.sqrt((A0 * S) ** 2 + (B0 * np.sin(Theta * ts)) ** 2 + (C0 * M) ** 2)
    return pred - E_bdg


def main():
    os.makedirs('results', exist_ok=True)
    T_step = P.FIG3_T
    n_per_step = max(200, P.FIG3_N_PER_STEP)
    print('Computing BdG min energies along path for Fig.3...')
    tlist, E_bdg, S, M, Theta, dt = compute_bdg_minE_with_theta(T_step=T_step, n_per_step=n_per_step)
    print('Done. Running constrained least-squares fit (A0,B0,C0,ts) with C0<=0.5 ...')

    # initial guess from previous fit if available
    x0 = np.array([0.00265, 0.01, 0.16, 1.0])
    try:
        # prefer global ABC fit if present
        d = np.load('results/mapping_fit_ABC.npz')
        if 'res' in d.files:
            r = d['res']
            x0[0] = float(r[0]) if len(r) > 0 else x0[0]
            x0[1] = float(r[1]) if len(r) > 1 else x0[1]
            x0[2] = float(r[2]) if len(r) > 2 else x0[2]
    except Exception:
        pass
    try:
        d2 = np.load('results/mapping_fit_fig3.npz')
        x0[0] = float(d2['A0_fit']) if 'A0_fit' in d2.files else x0[0]
        x0[1] = float(d2['B0_fit']) if 'B0_fit' in d2.files else x0[1]
        x0[2] = float(d2['C0_fit']) if 'C0_fit' in d2.files else x0[2]
    except Exception:
        pass

    lb = np.array([0.0, 0.0, 0.0, 0.2])
    ub = np.array([np.inf, np.inf, 0.5, 5.0])  # constrain C0 <= 0.5
    # ensure initial guess lies within bounds
    x0 = np.clip(x0, lb, ub)

    res = least_squares(residuals, x0, args=(S, M, Theta, E_bdg), bounds=(lb.tolist(), ub.tolist()), xtol=1e-8, ftol=1e-8)

    A0_fit, B0_fit, C0_fit, ts_fit = res.x
    print(f'Constrained fit results: A0={A0_fit:.6g}, B0={B0_fit:.6g}, C0={C0_fit:.6g}, ts={ts_fit:.6g}')

    outpath = 'results/mapping_fit_fig3_constrained_C0le0.5.npz'
    np.savez(outpath, tlist=tlist, E_bdg=E_bdg, S=S, M=M, Theta=Theta,
             A0_fit=A0_fit, B0_fit=B0_fit, C0_fit=C0_fit, ts_fit=ts_fit, res=res.x)
    print('Saved', outpath)


if __name__ == '__main__':
    main()
