#!/usr/bin/env python3
"""Simulate piecewise tetron gate path on a two-level (Bloch) effective Hamiltonian.

This script builds the 3-step gate sequence described in ybe_ver.md,
constructs H(t)=d_x(t) sigma_x + d_y(t) sigma_y + d_z(t) sigma_z,
and time-evolves an initial eigenstate. It saves results and plots.
"""
import os
import numpy as np
from scipy.linalg import expm, eigh
import matplotlib.pyplot as plt


def make_time_grid(T_step=200.0, n_per_step=200):
    steps = 3
    dt = T_step / float(n_per_step)
    tlist = []
    step_idx = []
    slist = []
    for i in range(steps):
        for k in range(n_per_step):
            s = k / float(n_per_step)
            t = (i + s) * T_step
            tlist.append(t)
            step_idx.append(i + 1)
            slist.append(s)
    return np.array(tlist), np.array(step_idx), np.array(slist), dt


def gates_at(step, s):
    # piecewise gate definitions from documentation
    if step == 1:
        g1 = 1.0 - s
        g3 = s
        g2 = 1.0
        g4 = 0.0
    elif step == 2:
        g2 = 1.0 - s
        g1 = s
        g3 = 1.0
        g4 = 0.0
    elif step == 3:
        g3 = 1.0 - s
        g2 = s
        g1 = 1.0
        g4 = 0.0
    else:
        g1 = g2 = g3 = g4 = 0.0
    return g1, g2, g3, g4


def theta_from_time(step, s):
    # normalized phase over 3 steps
    frac = ((step - 1) + s) / 3.0
    return 2.0 * np.pi * frac


def H_eff_from_theta(theta, delta=0.0):
    dx = np.cos(theta)
    dy = np.sin(theta)
    dz = delta
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sy = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
    sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    return dx * sx + dy * sy + dz * sz


def initial_eigenstate(H0):
    vals, vecs = eigh(H0)
    # pick the ground (lowest) eigenvector
    idx = np.argmin(vals)
    psi0 = vecs[:, idx]
    return psi0


def run_sim(T_step=200.0, n_per_step=200, delta=0.0, save_prefix="results/tetron"):
    tlist, step_idx, slist, dt = make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)
    # prepare arrays
    probs = np.zeros((N, 2))
    overlaps = np.zeros(N, dtype=complex)
    thetas = np.zeros(N)

    # initial H and initial state
    theta0 = theta_from_time(step_idx[0], slist[0])
    H0 = H_eff_from_theta(theta0, delta=delta)
    psi = initial_eigenstate(H0)
    psi = psi / np.linalg.norm(psi)
    psi0 = psi.copy()

    # time evolution
    U = np.eye(2, dtype=complex)
    for i in range(N):
        step = int(step_idx[i])
        s = float(slist[i])
        g1, g2, g3, g4 = gates_at(step, s)
        theta = theta_from_time(step, s)
        thetas[i] = theta
        H = H_eff_from_theta(theta, delta=delta)
        # short-time propagator
        Ustep = expm(-1j * H * dt)
        U = Ustep @ U
        psi = U @ psi0
        probs[i, 0] = np.abs(psi[0]) ** 2
        probs[i, 1] = np.abs(psi[1]) ** 2
        overlaps[i] = np.vdot(psi0, psi)

    # save results and plot
    os.makedirs(os.path.dirname(save_prefix), exist_ok=True)
    np.save(save_prefix + f"_T{int(T_step)}_delta{delta}.npy", {"t": tlist, "probs": probs, "overlaps": overlaps, "thetas": thetas})

    plt.figure(figsize=(8, 4))
    plt.plot(tlist, probs[:, 0], label="p(site1)")
    plt.plot(tlist, probs[:, 1], label="p(site2)")
    plt.xlabel('time')
    plt.ylabel('probability')
    plt.legend()
    plt.tight_layout()
    pngfile = save_prefix + f"_T{int(T_step)}_delta{delta}.png"
    plt.savefig(pngfile)
    plt.close()

    print(f"Saved results to {save_prefix}_T{int(T_step)}_delta{delta}.npy and {pngfile}")


if __name__ == '__main__':
    # run a quick scan for a few T_step values and two delta cases
    T_list = [50, 100, 200, 400]
    for T in T_list:
        run_sim(T_step=T, n_per_step=300, delta=0.0, save_prefix=f"results/tetron_MZM")
        run_sim(T_step=T, n_per_step=300, delta=0.2, save_prefix=f"results/tetron_ABS")
