#!/usr/bin/env python3
"""Run tetron two-level simulation using mapping constants (A0, B0, C0).

This script maps gate values to effective Pauli coefficients using the
fitted mapping constants (from results/mapping_fit.npz) and evolves the
two-level system for the tetron path. Outputs saved to `results/`.
"""
import os
import sys
import numpy as np
from scipy.linalg import expm, eigh
import matplotlib.pyplot as plt

# import repo modules
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import tetron_path_sim as tetron
from tools.reproduce_figs import map_gates_to_links
import tools.paper_params as P


def H_from_mapped(dx, dy, dz):
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sy = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
    sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    return dx * sx + dy * sy + dz * sz


def run_mapped(T_step=200.0, n_per_step=300, A0=0.00265, B0=0.0, C0=0.16, save_prefix='results/mapped'):
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)
    psi_traj = np.zeros((N, 2), dtype=complex)
    bloch_traj = np.zeros((N, 3))

    # initial H using first step
    theta0 = tetron.theta_from_time(step_idx[0], slist[0])
    # initial gate mapping
    g1, g2, g3, g4 = tetron.gates_at(int(step_idx[0]), float(slist[0]))
    mu0, t_links0, Delta0 = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, P.VD, P.QD_WIDTH)
    t_left0 = t_links0[0] if len(t_links0) > 0 else 0.0
    t_right0 = t_links0[1] if len(t_links0) > 1 else 0.0
    mu_diff0 = mu0[0] - mu0[1] if len(mu0) > 1 else 0.0
    dx0 = A0 * (t_left0 + t_right0)
    dy0 = B0 * 0.0
    dz0 = C0 * mu_diff0
    H0 = H_from_mapped(dx0, dy0, dz0)
    vals, vecs = eigh(H0)
    psi0 = vecs[:, np.argmin(vals)]

    U = np.eye(2, dtype=complex)
    for i in range(N):
        step = int(step_idx[i])
        s = float(slist[i])
        g1, g2, g3, g4 = tetron.gates_at(step, s)
        mu, t_links_mod, Delta_mod = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, P.VD, P.QD_WIDTH)
        t_left = t_links_mod[0] if len(t_links_mod) > 0 else 0.0
        t_right = t_links_mod[1] if len(t_links_mod) > 1 else 0.0
        mu_diff = mu[0] - mu[1] if len(mu) > 1 else 0.0
        dx = A0 * (t_left + t_right)
        dy = B0 * 0.0
        dz = C0 * mu_diff
        H = H_from_mapped(dx, dy, dz)
        Ustep = expm(-1j * H * dt)
        U = Ustep @ U
        psi = U @ psi0
        psi_traj[i, :] = psi
        # bloch vector
        sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
        sy = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
        sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
        bx = np.vdot(psi, sx @ psi).real
        by = np.vdot(psi, sy @ psi).real
        bz = np.vdot(psi, sz @ psi).real
        bloch_traj[i, :] = np.array([bx, by, bz])

    os.makedirs(os.path.dirname(save_prefix), exist_ok=True)
    np.save(save_prefix + f"_T{int(T_step)}_mapped.npy", {"t": tlist, "bloch": bloch_traj})
    # simple plot: final bloch components
    bfinal = bloch_traj[-1]
    plt.figure()
    plt.plot(bloch_traj[:,0], label='bx')
    plt.plot(bloch_traj[:,1], label='by')
    plt.plot(bloch_traj[:,2], label='bz')
    plt.legend()
    plt.title(f'mapped final bloch (T={int(T_step)})')
    out = save_prefix + f"_T{int(T_step)}_mapped.png"
    plt.savefig(out)
    plt.close()
    print('Saved', out)


if __name__ == '__main__':
    # load fitted constants if available
    fitf = 'results/mapping_fit.npz'
    if os.path.exists(fitf):
        d = np.load(fitf)
        A0 = float(d['A0_fit'])
        C0 = float(d['C0_fit'])
    else:
        A0 = 0.00265382
        C0 = 0.16
    B0 = 0.0
    for T in [400, 450, 500]:
        run_mapped(T_step=T, n_per_step=400, A0=A0, B0=B0, C0=C0, save_prefix=f'results/mapped_tetron_T{T}')
