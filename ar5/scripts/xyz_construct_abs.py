#!/usr/bin/env python3
"""Construct an inhomogeneous Kitaev chain with a quantum-dot like potential
to create an ABS localized at one end, then report lowest modes and LDOS.
"""
import numpy as np
import argparse
import os
import sys

# allow importing tools when run from repo root
ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from tools.verify_mzm import build_kitaev_bdg
from tools.xxz_R_and_H import expand_on_pauli


def build_inhom_bdg(N, t, Delta, mu_bulk, mu_dot, dot_width=4):
    mu = np.ones(N, dtype=float) * mu_bulk
    # place dot at left end spanning dot_width sites
    mu[:dot_width] = mu_dot
    # construct A,B with site-dependent mu
    A = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)
    for i in range(N):
        A[i, i] = -mu[i]
        if i + 1 < N:
            A[i, i + 1] = -t
            A[i + 1, i] = -t
            B[i, i + 1] = Delta
            B[i + 1, i] = -Delta
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    Hbdg = np.concatenate((top, bottom), axis=0)
    return Hbdg


def run_example():
    N = 80
    t = 1.0
    Delta = 0.15
    mu_bulk = 2.0
    mu_dot = -0.5
    H = build_inhom_bdg(N, t, Delta, mu_bulk, mu_dot, dot_width=8)
    w, v = np.linalg.eigh(H)
    idx = np.argsort(np.abs(w))
    print('Lowest 8 eigenvalues:')
    for i in range(8):
        print(f'{w[idx[i]]:.6e}')
    psi0 = v[:, idx[0]]
    u = psi0[:N]
    ldos = np.abs(u)**2
    os.makedirs('results', exist_ok=True)
    np.savez('results/xyz_abs_example.npz', eigvals=w, psi0=psi0, ldos=ldos)
    print('Saved results/xyz_abs_example.npz')


def run_from_time_seq(time_seq_path, L=80, t_index=0, mu_dot=-0.5, dot_width=8):
    d = np.load(time_seq_path)
    key = f'H{t_index}'
    if key not in d.files:
        raise SystemExit(f'{key} not in {time_seq_path}')
    H2 = d[key]
    mapping = expand_on_pauli(H2)
    # map to Kitaev params
    c_xx = mapping.get('X_X', 0.0)
    c_yy = mapping.get('Y_Y', 0.0)
    c_xy = mapping.get('X_Y', 0.0)
    c_yx = mapping.get('Y_X', 0.0)
    c_zz = mapping.get('Z_Z', 0.0)
    c_z0 = mapping.get('Z_I', 0.0)
    c_0z = mapping.get('I_Z', 0.0)
    from tools.verify_mzm import map_c_to_params
    t_val, Delta_val, mu_val, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    # build inhom BdG using these uniform bond params, with a local dot at left
    H = build_inhom_bdg(L, t_val, Delta_val, mu_val, mu_dot, dot_width=dot_width)
    w, v = np.linalg.eigh(H)
    idx = np.argsort(np.abs(w))
    print('Lowest 8 eigenvalues:')
    for i in range(8):
        print(f'{w[idx[i]]:.6e}')
    psi0 = v[:, idx[0]]
    u = psi0[:L]
    ldos = np.abs(u)**2
    os.makedirs('results', exist_ok=True)
    np.savez('results/xyz_abs_example_from_model.npz', eigvals=w, psi0=psi0, ldos=ldos)
    print('Saved results/xyz_abs_example_from_model.npz')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--time_seq', type=str, default=None,
                        help='Path to time_H_sequence.npz to derive base params')
    parser.add_argument('--L', type=int, default=80)
    parser.add_argument('--t_index', type=int, default=0)
    parser.add_argument('--mu_dot', type=float, default=-0.5)
    parser.add_argument('--dot_width', type=int, default=8)
    args = parser.parse_args()
    if args.time_seq:
        run_from_time_seq(args.time_seq, L=args.L, t_index=args.t_index, mu_dot=args.mu_dot, dot_width=args.dot_width)
    else:
        run_example()
