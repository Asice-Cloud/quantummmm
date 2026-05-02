#!/usr/bin/env python3
"""Create a time-dependent hybridization H_hyb(t)=i d(t) gamma1 gamma2 projected
into the subspace of two chosen BdG eigenmodes.

Usage: provide a saved BdG eigenvector NPZ (psi0, psi1) or point to results from other scripts.
This script will build H_full(t) by embedding a 2x2 diag(-eps,+eps) in the subspace.
"""
import numpy as np
import argparse
import os
import sys

# ensure repo root on sys.path
ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def load_modes(npz_path, idx0=0, idx1=1):
    d = np.load(npz_path)
    # expect d to contain 'eigvecs' or individual 'psi0' entries
    if 'eigvecs' in d.files:
        vecs = d['eigvecs']
        return vecs[:, idx0], vecs[:, idx1]
    if 'psi0' in d.files and 'psi1' in d.files:
        return d['psi0'], d['psi1']
    # fallback: take two lowest eigenvectors from 'psi' like saved format
    keys = [k for k in d.files if k.startswith('psi')]
    keys = sorted(keys)
    if len(keys) >= 2:
        return d[keys[0]], d[keys[1]]
    raise SystemExit('Could not find mode vectors in ' + npz_path)


def build_hybrid_sequence(times, v1, v2, eps_func):
    # build subspace projector
    V = np.column_stack([v1, v2])  # shape (M,2)
    Hs = []
    for t in times:
        eps = eps_func(t)
        Hsub = np.diag([-eps, +eps])
        Hfull = V @ Hsub @ V.conj().T
        # ensure Hermitian
        Hfull = 0.5 * (Hfull + Hfull.conj().T)
        Hs.append(Hfull)
    return Hs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--modes', required=True, help='NPZ file with two mode vectors')
    parser.add_argument('--out', default='results/hybrid_sequence.npz')
    parser.add_argument('--T', type=float, default=1.0)
    parser.add_argument('--n_steps', type=int, default=200)
    parser.add_argument('--eps_max', type=float, default=0.1)
    args = parser.parse_args()

    v1, v2 = load_modes(args.modes)
    times = np.linspace(0.0, args.T, args.n_steps + 1)

    def eps_func(t):
        # simple sinusoidal pulse from 0 to eps_max and back
        return args.eps_max * np.sin(np.pi * t / args.T)**2

    Hs = build_hybrid_sequence(times, v1, v2, eps_func)
    # save
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    arrs = {'times': times}
    for i, H in enumerate(Hs):
        arrs[f'H{i}'] = H
    np.savez(args.out, **arrs)
    print('Saved hybrid sequence to', args.out)


if __name__ == '__main__':
    main()
