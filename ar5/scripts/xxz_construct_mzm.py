#!/usr/bin/env python3
"""Build a Kitaev chain from XXZ-derived Pauli coefficients and inspect low-energy modes.

Usage: adjust c_xx,c_yy,... or call from CLI.
"""
import os
import sys
ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np
import argparse
from tools.verify_mzm import map_c_to_params, build_kitaev_bdg
from tools.xxz_R_and_H import expand_on_pauli


def run_example_from_coeffs(N=40, c_xx=1.0, c_yy=0.0, c_xy=0.0, c_yx=0.0, c_zz=0.0, c_z0=0.0, c_0z=0.0):
    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    Hbdg = build_kitaev_bdg(N, t, Delta, mu)
    w, v = np.linalg.eigh(Hbdg)
    w = np.real_if_close(w)
    idx = np.argsort(np.abs(w))
    print('Lowest 6 eigenvalues:')
    for i in range(6):
        print(f'{w[idx[i]]:.6e}')
    psi0 = v[:, idx[0]]
    os.makedirs('results', exist_ok=True)
    np.savez('results/xxz_mzm_example.npz', eigvals=w, psi0=psi0)
    print('Saved results/xxz_mzm_example.npz')


def run_from_time_seq(time_seq_path, N=40, t_index=0):
    d = np.load(time_seq_path)
    # pick Ht at t_index
    key = f'H{t_index}'
    if key not in d.files:
        raise SystemExit(f'{key} not found in {time_seq_path}')
    H2 = d[key]
    mapping = expand_on_pauli(H2)
    c_xx = mapping.get('X_X', 0.0)
    c_yy = mapping.get('Y_Y', 0.0)
    c_xy = mapping.get('X_Y', 0.0)
    c_yx = mapping.get('Y_X', 0.0)
    c_zz = mapping.get('Z_Z', 0.0)
    c_z0 = mapping.get('Z_I', 0.0)
    c_0z = mapping.get('I_Z', 0.0)
    print('Mapped Pauli coeffs from time_seq:', {k: mapping.get(k,0) for k in ['X_X','Y_Y','X_Y','Y_X','Z_Z','Z_I','I_Z']})
    run_example_from_coeffs(N=N, c_xx=c_xx, c_yy=c_yy, c_xy=c_xy, c_yx=c_yx, c_zz=c_zz, c_z0=c_z0, c_0z=c_0z)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--N', type=int, default=40)
    parser.add_argument('--c_xx', type=float, default=1.0)
    parser.add_argument('--c_yy', type=float, default=0.0)
    parser.add_argument('--time_seq', type=str, default=None,
                        help='Path to time_H_sequence.npz to derive coefficients from model')
    parser.add_argument('--t_index', type=int, default=0,
                        help='Time index within time sequence to use')
    args = parser.parse_args()
    if args.time_seq:
        run_from_time_seq(args.time_seq, N=args.N, t_index=args.t_index)
    else:
        run_example_from_coeffs(N=args.N, c_xx=args.c_xx, c_yy=args.c_yy)
