#!/usr/bin/env python3
"""Diagnose and reconcile U/U_eff between two result files.

Usage:
  python3 tools/diagnose_projection.py --a results/fileA.npz --b results/fileB.npz

This script:
- Loads two .npz result files (defaults chosen to existing runs).
- Chooses a canonical full-space `U` and a canonical `Vlow0`.
- Projects `U` -> `U_eff`, computes fidelity against ideal braid,
  compares to stored `Ueff` entries, and maps min-gap step to snapshots.
"""
import argparse
import numpy as np
import numpy.linalg as la


def ideal_braid_2x2():
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    B = np.cos(np.pi/4.0) * np.eye(2, dtype=complex) + 1j * np.sin(np.pi/4.0) * sy
    return B


def load_npz(path):
    return np.load(path, allow_pickle=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--a', default='results/baxterize_demo_h0.150_L160.npz')
    p.add_argument('--b', default='results/braid_sim_record_h0.150_L160.npz')
    p.add_argument('--save', action='store_true', help='Save numeric diagnostics to results/')
    args = p.parse_args()

    A = load_npz(args.a)
    B = load_npz(args.b)
    print('A keys:', A.files)
    print('B keys:', B.files)

    # choose canonical full U: prefer B.U_final else A.U
    U_full = B.get('U_final', None)
    if U_full is None:
        U_full = A.get('U', None)
    print('Using U_full from', 'B.U_final' if 'U_final' in B.files else 'A.U' if 'U' in A.files else 'none')

    # choose canonical Vlow: prefer B then A
    V = B.get('Vlow0', None)
    if V is None:
        V = A.get('Vlow0', None)
    print('Using Vlow from', 'B' if 'Vlow0' in B.files else 'A' if 'Vlow0' in A.files else 'none')

    if U_full is None or V is None:
        print('Missing U_full or Vlow0; abort')
        return

    # project and evaluate
    Ueff = V.conj().T.dot(U_full).dot(V)
    B2 = ideal_braid_2x2()
    frob = la.norm(Ueff - B2)
    fid_like = np.real(np.trace(Ueff.conj().T @ B2)) / 2.0
    phi = np.angle(np.vdot(Ueff.ravel(), B2.ravel()))
    fid_phase_invariant = np.real(np.vdot(Ueff.ravel(), B2.ravel()) * np.exp(-1j*phi)) / 2.0

    print('Projected Ueff norm diff to B2 (frob)=', frob)
    print('fid_like (canonical)=', fid_like)
    print('fid (phase-invariant inner prod /2)=', fid_phase_invariant)

    # compare to stored Ueffs if present
    if 'Ueff' in A.files:
        print('norm(Ueff - A.Ueff_stored)=', la.norm(Ueff - A['Ueff']))
    if 'Ueff_final' in B.files:
        print('norm(Ueff - B.Ueff_final)=', la.norm(Ueff - B['Ueff_final']))

    # compare Vlow overlap and condition numbers
    V1 = A.get('Vlow0', None)
    V2 = B.get('Vlow0', None)
    if V1 is not None and V2 is not None:
        print('cond V1', la.cond(V1), 'cond V2', la.cond(V2))
        S = V1.conj().T.dot(V2)
        print('overlap singular values', la.svd(S, compute_uv=False))

    # map min gap to snapshot if available
    if 'gaps' in B.files:
        g = B['gaps']
        imin = int(np.argmin(g))
        print('min gap value=', float(g[imin]), 'at step index=', imin)
        if 'times' in B.files and 'Ueff_list' in B.files:
            times = B['times']
            Ulist = B['Ueff_list']
            save_every = len(g) // len(times)
            snap = imin // save_every
            print('save_every=', save_every, 'nearest snapshot=', snap, 'snapshot time=', float(times[snap]))
            rng = range(max(0, snap-2), min(len(Ulist), snap+3))
            for s in rng:
                print('\nsnapshot', s, 'Ueff_list[s]=')
                print(Ulist[s])

    if args.save:
        out = dict(frob=frob, fid_like=fid_like, fid_phase_invariant=fid_phase_invariant)
        np.savez('results/diagnostics_projection.npz', **out)
        print('Saved diagnostics to results/diagnostics_projection.npz')


if __name__ == '__main__':
    main()
