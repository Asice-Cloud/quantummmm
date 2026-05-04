#!/usr/bin/env python3
"""Batch per-figure mapping fits for Fig.2, Fig.4 and Fig.5.

For each figure this script computes BdG minimal energies along the tetron
path (with modulation when appropriate) and fits A0,B0,C0,ts to match the
Pauli-model prediction.

Outputs: results/mapping_fit_fig2.npz, results/mapping_fit_fig4.npz,
         results/mapping_fit_fig5_<amp>.npz
"""
import os
import sys
import numpy as np
from scipy.optimize import least_squares

# repo import
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import tetron_path_sim as tetron
from tools.reproduce_figs import map_gates_to_links
from tools import embed_kitaev
import tools.paper_params as P


def compute_bdg_minE(T_step=100.0, n_per_step=300, mod_fn=None):
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
        # apply modulation if given (e.g., VD_t scaling)
        if mod_fn is None:
            VD_t = P.VD
        else:
            t = tlist[i]
            VD_t = mod_fn(t)
        mu, t_links_mod, Delta_mod = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, VD_t, P.QD_WIDTH)
        H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
        E = np.linalg.eigvalsh(H)
        E_bdg[i] = np.min(np.abs(E))
        t_left = t_links_mod[0] if len(t_links_mod) > 0 else 0.0
        t_right = t_links_mod[1] if len(t_links_mod) > 1 else 0.0
        S[i] = t_left + t_right
        M[i] = mu[0] - mu[1] if len(mu) > 1 else 0.0

    return tlist, E_bdg, S, M, Theta, dt


def do_fit(tlist, E_bdg, S, M, Theta, x0=None, bounds=None):
    if x0 is None:
        x0 = np.array([0.00265, 0.01, 0.16, 1.0])
    if bounds is None:
        lb = [0.0, 0.0, 0.0, 0.2]
        ub = [np.inf, np.inf, np.inf, 5.0]
    else:
        lb, ub = bounds

    def residuals(x, S, M, Theta, E_bdg):
        A0, B0, C0, ts = x
        pred = np.sqrt((A0 * S) ** 2 + (B0 * np.sin(Theta * ts)) ** 2 + (C0 * M) ** 2)
        return pred - E_bdg

    res = least_squares(residuals, x0, args=(S, M, Theta, E_bdg), bounds=(lb, ub), xtol=1e-8, ftol=1e-8)
    return res


def fit_fig2():
    outf = 'results/mapping_fit_fig2.npz'
    T_list = P.FIG2_TS
    results = {}
    for T in T_list:
        print('Fig2: computing BdG for T=', T)
        tlist, E_bdg, S, M, Theta, dt = compute_bdg_minE(T_step=T, n_per_step=300)
        print('Fitting mapping for T=', T)
        res = do_fit(tlist, E_bdg, S, M, Theta)
        results[f'T{int(T)}'] = res.x
    np.savez(outf, **results)
    print('Saved', outf)


def fit_fig4():
    outf = 'results/mapping_fit_fig4.npz'
    # modulation used in reproduce_figs: VD_t = P.VD * (1.0 + mod) where mod = Vx0 + Vx1*cos(pi*t/T)
    Vx0 = P.FIG4_VX0
    Vx1 = P.FIG4_VX1
    T = 200.0

    def mod_fn(t):
        mod = Vx0 + Vx1 * np.cos(np.pi * t / T)
        return P.VD * (1.0 + mod)

    print('Fig4: computing BdG with modulation')
    tlist, E_bdg, S, M, Theta, dt = compute_bdg_minE(T_step=T, n_per_step=300, mod_fn=mod_fn)
    res = do_fit(tlist, E_bdg, S, M, Theta)
    np.savez(outf, res=res.x)
    print('Saved', outf)


def fit_fig5():
    outs = []
    for amp in P.FIG5_VX1_OPTIONS:
        name = f'results/mapping_fit_fig5_amp{amp:.6g}.npz'
        Vx0 = P.FIG4_VX0
        Vx1 = amp
        T = 200.0

        def mod_fn(t):
            mod = Vx0 + Vx1 * np.cos(np.pi * t / T)
            return P.VD * (1.0 + mod)

        print('Fig5: computing BdG with modulation amp=', amp)
        tlist, E_bdg, S, M, Theta, dt = compute_bdg_minE(T_step=T, n_per_step=300, mod_fn=mod_fn)
        res = do_fit(tlist, E_bdg, S, M, Theta)
        np.savez(name, res=res.x)
        outs.append(name)
        print('Saved', name)
    return outs


def main():
    os.makedirs('results', exist_ok=True)
    fit_fig2()
    fit_fig4()
    fit_fig5()


if __name__ == '__main__':
    main()
