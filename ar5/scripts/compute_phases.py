#!/usr/bin/env python3
"""Compute dynamic and geometric (Berry) phases from a time-dependent H(t).

Supports two ways to get the total phase:
 - Instantaneous eigenstate path (gauge-fixed continuity)
 - Direct time-ordered evolution U(T) via stepwise diagonalization

Usage:
  python3 scripts/compute_phases.py --in results/time_H_sequence.npz
"""
import numpy as np
import argparse
import os


def load_H_sequence(path):
    d = np.load(path)
    times = d['times']
    # collect keys H0..Hn in order
    Hkeys = sorted([k for k in d.files if k.startswith('H')], key=lambda s: int(s[1:]))
    Hs = [d[k] for k in Hkeys]
    return times, Hs


def gauge_fix_vectors(vecs):
    """Enforce phase continuity: choose phase of vecs[t] so <vecs[t-1]|vecs[t]> is real positive."""
    vecs = [v.copy() for v in vecs]
    for i in range(1, len(vecs)):
        overlap = np.vdot(vecs[i-1], vecs[i])
        if np.abs(overlap) < 1e-16:
            # orthogonal (level crossing); leave as is
            continue
        phase = np.angle(overlap)
        vecs[i] *= np.exp(-1j * phase)
    return vecs


def compute_instantaneous_phases(times, Hs, eig_index=0):
    n = len(Hs)
    eigvals = []
    eigvecs = []
    for H in Hs:
        w, v = np.linalg.eigh(H)
        eigvals.append(w)
        eigvecs.append(v[:, eig_index])

    eigvals = np.array(eigvals)
    eigvecs = gauge_fix_vectors(eigvecs)

    # dynamic phase: -∫ <ψ|H|ψ> dt (trapezoid)
    expect = np.array([np.real(np.vdot(psi, H @ psi)) for psi, H in zip(eigvecs, Hs)])
    dt = np.diff(times)
    # trapezoidal integration
    dyn = 0.0
    for i in range(len(dt)):
        dyn += -0.5 * (expect[i] + expect[i+1]) * dt[i]

    # total phase via overlap between initial and final (gauged) states
    total = np.angle(np.vdot(eigvecs[0], eigvecs[-1]))
    geo = total - dyn
    # normalize to (-pi, pi]
    def wrap(x):
        return (x + np.pi) % (2 * np.pi) - np.pi

    return {
        'dynamic': dyn,
        'total': wrap(total),
        'geometric': wrap(geo),
        'eigvals': eigvals,
    }


def evolve_via_propagator(times, Hs, psi0=None):
    """Build step propagators U_k = V exp(-i E dt) V^† and apply to psi0.
    If psi0 is None, use instantaneous ground state at t=0.
    Returns final state psiT and full propagator U.
    """
    n = len(Hs)
    if psi0 is None:
        w0, v0 = np.linalg.eigh(Hs[0])
        psi0 = v0[:, 0]

    U = np.eye(Hs[0].shape[0], dtype=complex)
    psi = psi0.copy()
    for i in range(n - 1):
        dt = times[i+1] - times[i]
        w, v = np.linalg.eigh(Hs[i])
        Ustep = (v @ np.diag(np.exp(-1j * w * dt)) @ v.conj().T)
        U = Ustep @ U
        psi = Ustep @ psi
    return psi, U


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', dest='inp', required=True)
    parser.add_argument('--eig_index', type=int, default=0,
                        help='Index of instantaneous eigenstate to follow (0=lowest)')
    parser.add_argument('--use_propagator', action='store_true')
    args = parser.parse_args()

    times, Hs = load_H_sequence(args.inp)
    res_inst = compute_instantaneous_phases(times, Hs, eig_index=args.eig_index)

    print('Instantaneous eigenstate results:')
    print('  dynamic phase:', res_inst['dynamic'])
    print('  total phase:', res_inst['total'])
    print('  geometric phase:', res_inst['geometric'])

    if args.use_propagator:
        psiT, U = evolve_via_propagator(times, Hs)
        # overlap
        w0, v0 = np.linalg.eigh(Hs[0])
        psi0 = v0[:, args.eig_index]
        overlap = np.vdot(psi0, psiT)
        tot = np.angle(overlap)
        # dynamic via integrated expectation
        expect = np.array([np.real(np.vdot(psi0, H @ psi0)) for H in Hs])
        # note: this is not state-following dynamic — for evolved state use energy expectation along path
        print('Propagator evolution overlap total phase:', tot)


if __name__ == '__main__':
    main()
