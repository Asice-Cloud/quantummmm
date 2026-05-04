#!/usr/bin/env python3
"""Fit mapping parameters by matching overlap(T) between BdG and mapped Pauli models (Fig.3).

This script computes, for a set of T values, the final overlap |<psi0|psi(T)>|
from (A) the full-chain BdG evolution projected to the low-energy subspace,
and (B) the mapped 2x2 Pauli evolution. It fits A0,B0,C0,ts to minimize
the L2 difference across the T scan.

Warning: full-chain time evolution is compute-heavy (matrix exponentials).
Use moderate `n_per_step` and `Tscan` to control runtime.
"""
import os
import sys
import time
import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import expm, eigh

# make repo importable
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import tetron_path_sim as tetron
from tools.reproduce_figs import map_gates_to_links
from tools import embed_kitaev
import tools.paper_params as P


def compute_bdg_overlap(T_step=100.0, n_per_step=120):
    """Evolve full BdG U and return overlap |<e1|Ueff|e1>| where e1 is initial low-energy eigenvector."""
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    L = P.L
    dim = 2 * L

    # initial H and Vlow0
    g1, g2, g3, g4 = tetron.gates_at(int(step_idx[0]), float(slist[0]))
    mu0, t_links0, Delta0 = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, P.VD, P.QD_WIDTH)
    H0 = embed_kitaev.build_bdg(mu0, t_links0, Delta0)
    vals0, vecs0 = eigh(H0)
    # choose two lowest-abs eigenvectors as low-energy subspace
    idx = np.argsort(np.abs(vals0))[:2]
    Vlow0 = vecs0[:, idx]

    U = np.eye(dim, dtype=complex)
    for i in range(len(tlist)):
        step = int(step_idx[i])
        s = float(slist[i])
        g1, g2, g3, g4 = tetron.gates_at(step, s)
        mu, t_links_mod, Delta_mod = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, P.VD, P.QD_WIDTH)
        H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
        Ustep = expm(-1j * H * dt)
        U = Ustep @ U

    Ueff = Vlow0.conj().T @ U @ Vlow0
    # initial basis e1 corresponds to first column of Vlow0 -> overlap = |Ueff[0,0]|
    ov = abs(Ueff[0, 0])
    return ov


def compute_mapped_overlap(T_step=100.0, n_per_step=120, A0=0.00265, B0=0.02, C0=0.16, ts=1.0):
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    # Pauli matrices
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sy = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
    sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

    # initial mapped H and eigenvector
    step0 = int(step_idx[0])
    s0 = float(slist[0])
    g1, g2, g3, g4 = tetron.gates_at(step0, s0)
    mu0, t_links0, Delta0 = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, P.VD, P.QD_WIDTH)
    t_left0 = t_links0[0] if len(t_links0) > 0 else 0.0
    t_right0 = t_links0[1] if len(t_links0) > 1 else 0.0
    mu_diff0 = mu0[0] - mu0[1] if len(mu0) > 1 else 0.0
    theta0 = tetron.theta_from_time(step0, s0)
    dx0 = A0 * (t_left0 + t_right0)
    dy0 = B0 * np.sin(theta0 * ts)
    dz0 = C0 * mu_diff0
    H0 = dx0 * sx + dy0 * sy + dz0 * sz
    vals0, vecs0 = eigh(H0)
    psi0 = vecs0[:, np.argmin(vals0)]

    U2 = np.eye(2, dtype=complex)
    for i in range(len(tlist)):
        step = int(step_idx[i])
        s = float(slist[i])
        g1, g2, g3, g4 = tetron.gates_at(step, s)
        mu, t_links_mod, Delta_mod = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, P.VD, P.QD_WIDTH)
        t_left = t_links_mod[0] if len(t_links_mod) > 0 else 0.0
        t_right = t_links_mod[1] if len(t_links_mod) > 1 else 0.0
        mu_diff = mu[0] - mu[1] if len(mu) > 1 else 0.0
        theta = tetron.theta_from_time(step, s)
        dx = A0 * (t_left + t_right)
        dy = B0 * np.sin(theta * ts)
        dz = C0 * mu_diff
        H = dx * sx + dy * sy + dz * sz
        Ustep = expm(-1j * H * dt)
        U2 = Ustep @ U2

    psi_final = U2 @ psi0
    ov = abs(np.vdot(psi0, psi_final))
    return ov


def residuals(x, Tscan, n_per_step):
    A0, B0, C0, ts = x
    res = []
    for T in Tscan:
        print(f'Computing T={T:.1f} (A0={A0:.4g},B0={B0:.4g},C0={C0:.4g},ts={ts:.4g})')
        t0 = time.time()
        ov_bdg = compute_bdg_overlap(T_step=T, n_per_step=n_per_step)
        ov_map = compute_mapped_overlap(T_step=T, n_per_step=n_per_step, A0=A0, B0=B0, C0=C0, ts=ts)
        dt = time.time() - t0
        print(f'  BdG ov={ov_bdg:.6f}, mapped ov={ov_map:.6f} (done {dt:.1f}s)')
        res.append(ov_map - ov_bdg)
    return np.array(res)


def main():
    os.makedirs('results', exist_ok=True)
    # choose Tscan similar to reproduce_figs (use reduced points to save time)
    Tscan = np.linspace(50, 1000, 10)
    n_per_step = 120

    # initial guess from mapping_fit_ABC if available
    x0 = np.array([0.0023, 0.02, 0.16, 1.0])
    try:
        d = np.load('results/mapping_fit_ABC.npz')
        x0[0] = float(d.get('A0_fit', x0[0]))
        x0[1] = float(d.get('B0_fit', x0[1]))
        x0[2] = float(d.get('C0_fit', x0[2]))
    except Exception:
        pass

    lb = [0.0, 0.0, 0.0, 0.2]
    ub = [np.inf, np.inf, np.inf, 5.0]

    print('Starting overlap-based fit (this may take significant time)...')
    res = least_squares(lambda x: residuals(x, Tscan, n_per_step), x0, bounds=(lb, ub), xtol=1e-6, ftol=1e-6)

    A0f, B0f, C0f, tsf = res.x
    print('Fit finished:', res.x)

    # compute final diagnostics
    ov_bdg_list = []
    ov_map_list = []
    for T in Tscan:
        ov_bdg = compute_bdg_overlap(T_step=T, n_per_step=n_per_step)
        ov_map = compute_mapped_overlap(T_step=T, n_per_step=n_per_step, A0=A0f, B0=B0f, C0=C0f, ts=tsf)
        ov_bdg_list.append(ov_bdg)
        ov_map_list.append(ov_map)

    np.savez('results/mapping_fit_overlap_fig3.npz', Tscan=Tscan, ov_bdg=np.array(ov_bdg_list), ov_map=np.array(ov_map_list),
             A0_fit=A0f, B0_fit=B0f, C0_fit=C0f, ts_fit=tsf, res=res.x)
    print('Saved results/mapping_fit_overlap_fig3.npz')


if __name__ == '__main__':
    main()
