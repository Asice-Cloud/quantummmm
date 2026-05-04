#!/usr/bin/env python3
"""Simulate Bloch rotations induced by a time-dependent coupling vector and ABS hybridization.

This script builds a 3-step piecewise path for a coupling vector delta(t) on the Bloch sphere
and optionally adds a time-dependent splitting d(t) (ABS hybridization) along sigma_z.
It time-evolves the initial eigenstate and records the final Bloch vector as a function of T_step,
showing the time-dependence of the rotation.
"""
import os
import numpy as np
from scipy.linalg import expm, eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def make_time_grid(T_step=200.0, n_per_step=300):
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


def phi_from_time(step, s):
    frac = ((step - 1) + s) / 3.0
    return 2.0 * np.pi * frac


def delta_vector(theta, alpha=0.8, A=1.0):
    # theta = azimuthal angle on equator, alpha = polar inclination
    dx = A * np.sin(alpha) * np.cos(theta)
    dy = A * np.sin(alpha) * np.sin(theta)
    dz = A * np.cos(alpha)
    return np.array([dx, dy, dz])


def H_from(delta_vec, d_split=0.0):
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sy = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
    sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    return delta_vec[0] * sx + delta_vec[1] * sy + (delta_vec[2] + d_split) * sz


def initial_eigenstate(H0, which='plus'):
    vals, vecs = eigh(H0)
    # choose eigenstate: plus=highest energy, minus=lowest
    if which == 'plus':
        idx = np.argmax(vals)
    else:
        idx = np.argmin(vals)
    psi0 = vecs[:, idx]
    return psi0


def bloch_vector(psi):
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sy = np.array([[0.0, -1j], [1j, 0.0]], dtype=complex)
    sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    bx = np.vdot(psi, sx @ psi).real
    by = np.vdot(psi, sy @ psi).real
    bz = np.vdot(psi, sz @ psi).real
    return np.array([bx, by, bz])


def run_for_T(T_step=200.0, n_per_step=300, alpha=0.8, A=1.0, d_profile=None):
    tlist, step_idx, slist, dt = make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)
    psi_traj = np.zeros((N, 2), dtype=complex)
    bloch_traj = np.zeros((N, 3))

    # initial H
    theta0 = phi_from_time(step_idx[0], slist[0])
    delta0 = delta_vector(theta0, alpha=alpha, A=A)
    H0 = H_from(delta0, d_split=0.0 if d_profile is None else d_profile(0.0))
    psi0 = initial_eigenstate(H0, which='minus')
    psi = psi0.copy()

    U = np.eye(2, dtype=complex)
    for i in range(N):
        step = int(step_idx[i])
        s = float(slist[i])
        theta = phi_from_time(step, s)
        delta_vec = delta_vector(theta, alpha=alpha, A=A)
        # determine d_split at this time
        t = tlist[i]
        d_split = 0.0 if d_profile is None else d_profile(t)
        H = H_from(delta_vec, d_split=d_split)
        Ustep = expm(-1j * H * dt)
        U = Ustep @ U
        psi = U @ psi0
        psi_traj[i, :] = psi
        bloch_traj[i, :] = bloch_vector(psi)

    return tlist, bloch_traj, psi0, psi_traj


def scan_T_and_plot(T_list, d_profile_fn, label, outprefix):
    os.makedirs('results', exist_ok=True)
    final_bx = []
    final_by = []
    final_bz = []
    for T in T_list:
        t, traj, psi0, psi_traj = run_for_T(T_step=T, n_per_step=400, alpha=0.8, A=1.0, d_profile=d_profile_fn)
        b_final = traj[-1]
        final_bx.append(b_final[0])
        final_by.append(b_final[1])
        final_bz.append(b_final[2])
    final_bx = np.array(final_bx)
    final_by = np.array(final_by)
    final_bz = np.array(final_bz)

    plt.figure()
    plt.plot(T_list, final_bx, '-o', label='bx')
    plt.plot(T_list, final_by, '-o', label='by')
    plt.plot(T_list, final_bz, '-o', label='bz')
    plt.xlabel('T_step')
    plt.ylabel('final Bloch vector components')
    plt.title(label)
    plt.legend()
    out = f'results/{outprefix}_final_bloch_vs_T.png'
    plt.savefig(out)
    print('Saved', out)

    # also save sample trajectory for longest T
    t, traj, psi0, psi_traj = run_for_T(T_step=T_list[-1], n_per_step=400, alpha=0.8, A=1.0, d_profile=d_profile_fn)
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(traj[:,0], traj[:,1], traj[:,2], '-', lw=1)
    ax.scatter([traj[0,0]],[traj[0,1]],[traj[0,2]], color='green', s=50, label='start')
    ax.scatter([traj[-1,0]],[traj[-1,1]],[traj[-1,2]], color='red', s=50, label='end')
    ax.set_xlabel('bx'); ax.set_ylabel('by'); ax.set_zlabel('bz')
    ax.set_title(f'{label} trajectory (T={T_list[-1]})')
    ax.legend()
    out2 = f'results/{outprefix}_traj_T{T_list[-1]}.png'
    plt.savefig(out2)
    print('Saved', out2)


if __name__ == '__main__':
    # define T list
    T_list = [50, 100, 200, 400, 800]

    # Case A: MZM, d(t)=0
    d_none = lambda t: 0.0
    scan_T_and_plot(T_list, d_none, label='MZM (d=0)', outprefix='bloch_MZM')

    # Case B: ABS, constant hybridization d0 during entire protocol
    d_const = lambda t: 0.3
    scan_T_and_plot(T_list, d_const, label='ABS (d=0.3)', outprefix='bloch_ABS')

    # Case C: ABS pulsed during certain step (e.g., nonzero only in Step2)
    def d_pulse(t):
        # total time = 3*T_step => we need ratio; choose T_step=400 baseline or inspect t range
        # we'll consider normalized time fraction
        # but approximate: assume maximum t ~ 3*800 for safety; use heuristic pulse in middle third
        # find total duration by reading environment: instead we implement pulse on time window
        # Here choose t in [200,400] roughly pulse
        if 200.0 <= t <= 400.0:
            return 0.3
        return 0.0

    scan_T_and_plot(T_list, d_pulse, label='ABS pulse mid', outprefix='bloch_ABS_pulse')
