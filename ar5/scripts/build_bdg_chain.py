#!/usr/bin/env python3
"""Map a sequence of 4x4 local two-site operators H2(t) into site-resolved BdG H(t).

Reads a time-local sequence NPZ (times, H0..HN) where each Hk is 4x4 acting
on sites (i,i+1). Produces a sequence of BdG single-particle matrices sized 2N
and saves to an output NPZ.
"""
import os
import sys
ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np
import argparse

from tools.xxz_R_and_H import expand_on_pauli
from tools.verify_mzm import map_c_to_params


def build_kitaev_bdg_general(N, t_bond, Delta_bond, mu_site):
    A = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)
    for i in range(N):
        A[i, i] = -mu_site[i]
        if i + 1 < N:
            t = t_bond[i]
            Delta = Delta_bond[i]
            A[i, i + 1] = -t
            A[i + 1, i] = -t
            B[i, i + 1] = Delta
            B[i + 1, i] = -Delta
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    Hbdg = np.concatenate((top, bottom), axis=0)
    return Hbdg


def load_local_sequence(path):
    d = np.load(path)
    times = d['times']
    Hkeys = sorted([k for k in d.files if k.startswith('H')], key=lambda s: int(s[1:]))
    Hs = [d[k] for k in Hkeys]
    return times, Hs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', dest='inp', required=True)
    parser.add_argument('--L', type=int, required=True)
    parser.add_argument('--out', default='results/H_bdg_chain.npz')
    args = parser.parse_args()

    times, H2s = load_local_sequence(args.inp)
    N = args.L
    Hbdg_list = []

    for H2 in H2s:
        mapping = expand_on_pauli(H2)
        c_xx = mapping.get('X_X', 0.0)
        c_yy = mapping.get('Y_Y', 0.0)
        c_xy = mapping.get('X_Y', 0.0)
        c_yx = mapping.get('Y_X', 0.0)
        c_zz = mapping.get('Z_Z', 0.0)
        c_z0 = mapping.get('Z_I', 0.0)
        c_0z = mapping.get('I_Z', 0.0)

        t_val, Delta_val, mu_val, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)

        # distribute bond parameters into chain arrays (assume uniform along chain)
        t_bond = np.full(N - 1, t_val, dtype=complex)
        Delta_bond = np.full(N - 1, Delta_val, dtype=complex)
        # distribute mu_val equally to the two sites of the bond (heuristic)
        mu_site = np.zeros(N, dtype=complex)
        mu_site += mu_val / 2.0

        Hbdg = build_kitaev_bdg_general(N, t_bond, Delta_bond, mu_site)
        Hbdg_list.append(Hbdg)

    # save sequence
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    arrs = {'times': times}
    for i, H in enumerate(Hbdg_list):
        arrs[f'H{i}'] = H
    np.savez(args.out, **arrs)
    print('Saved BdG chain sequence to', args.out)


if __name__ == '__main__':
    main()
