#!/usr/bin/env python3
"""Embed a time sequence of 2-site local operators into an L-site spin chain.

Produces a .npz with full-chain Hamiltonians for each time step: H0..HN.

Usage example:
  python3 scripts/embed_local_to_chain.py --in results/time_H_sequence.npz --L 6 --out results/H_spin_chain_L6.npz --apply_all
"""
import numpy as np
import argparse
import os


def load_local_sequence(path):
    d = np.load(path)
    times = d['times']
    Hkeys = sorted([k for k in d.files if k.startswith('H')], key=lambda s: int(s[1:]))
    Hs = [d[k] for k in Hkeys]
    return times, Hs


def embed_on_bond(H2, L, bond):
    # H2 is 4x4 acting on sites (bond, bond+1), bond in [0, L-2]
    if bond < 0 or bond > L - 2:
        raise ValueError('bond out of range')
    left_dim = 2 ** bond
    right_dim = 2 ** (L - bond - 2)
    full = np.kron(np.kron(np.eye(left_dim, dtype=complex), H2), np.eye(right_dim, dtype=complex))
    return full


def build_chain_sequence(times, H2_list, L, apply_all=False, bond=0):
    H_chain_list = []
    for H2 in H2_list:
        if apply_all:
            Hfull = np.zeros((2**L, 2**L), dtype=complex)
            for b in range(L - 1):
                Hfull += embed_on_bond(H2, L, b)
        else:
            Hfull = embed_on_bond(H2, L, bond)
        # ensure Hermitian
        Hfull = 0.5 * (Hfull + Hfull.conj().T)
        H_chain_list.append(Hfull)
    return times, H_chain_list


def save_sequence(path, times, H_list):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    arrs = {f'H{i}': H_list[i] for i in range(len(H_list))}
    arrs['times'] = times
    np.savez(path, **arrs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', dest='inp', required=True)
    parser.add_argument('--L', type=int, default=6)
    parser.add_argument('--bond', type=int, default=0,
                        help='Bond index to embed (0-based)')
    parser.add_argument('--apply_all', action='store_true',
                        help='Apply the local operator on all bonds and sum')
    parser.add_argument('--out', type=str, default='results/H_spin_chain.npz')
    args = parser.parse_args()

    times, H2s = load_local_sequence(args.inp)
    times, H_chain_list = build_chain_sequence(times, H2s, args.L, apply_all=args.apply_all, bond=args.bond)
    save_sequence(args.out, times, H_chain_list)
    print('Saved chain H sequence to', args.out, 'L=', args.L, 'apply_all=', args.apply_all)


if __name__ == '__main__':
    main()
