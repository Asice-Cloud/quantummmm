#!/usr/bin/env python3
"""Evolve a time-sequence of BdG single-particle Hamiltonians H(t)
by piecewise-constant propagation and extract the effective unitary
on the lowest two eigenmodes.

Input: NPZ with 'times' and H0..HN (each 2N x 2N).
Output: saves results including U_eff (2x2) and eigenphases.
"""
import numpy as np
import argparse
import os


def load_sequence(path):
    d = np.load(path)
    times = d['times']
    Hkeys = sorted([k for k in d.files if k.startswith('H')], key=lambda s: int(s[1:]))
    Hs = [d[k] for k in Hkeys]
    return times, Hs


def step_propagator(H, dt):
    # diagonalize and build exp(-i H dt)
    w, v = np.linalg.eigh(H)
    U = (v @ np.diag(np.exp(-1j * w * dt)) @ v.conj().T)
    return U


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', dest='inp', required=True)
    parser.add_argument('--out', dest='out', default='results/evolution_bdg.npz')
    parser.add_argument('--subspace_dim', type=int, default=2,
                        help='Dimension of low-energy subspace to project onto')
    args = parser.parse_args()

    times, Hs = load_sequence(args.inp)
    dt_list = np.diff(times)
    if np.any(dt_list <= 0):
        raise SystemExit('Non-increasing times in input')

    U = np.eye(Hs[0].shape[0], dtype=complex)
    for i in range(len(Hs) - 1):
        dt = times[i+1] - times[i]
        Ustep = step_propagator(Hs[i], dt)
        U = Ustep @ U

    # project onto lowest subspace at t=0
    w0, v0 = np.linalg.eigh(Hs[0])
    idx = np.argsort(np.abs(w0))[:args.subspace_dim]
    P = v0[:, idx]
    U_eff = P.conj().T @ U @ P

    # diagonalize U_eff for eigenphases
    uu, sv = np.linalg.eig(U_eff)
    phases = np.angle(uu)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    np.savez(args.out, U_eff=U_eff, phases=phases)
    print('Saved effective evolution to', args.out)
    print('Eigenphases (rad):', phases)


if __name__ == '__main__':
    main()
