#!/usr/bin/env python3
"""Optimize two single-qubit phases applied on the projected 2x2 space
to maximize Re Tr((D Ueff)^
                \dagger B)/2, where D=diag(e^{i a}, e^{i b}).

Usage: python3 tools/phase_compensate.py path/to/baxterize_demo_....npz
"""
import sys
import numpy as np
from pathlib import Path


def ideal_braid_2x2():
    sy = np.array([[0, -1j],[1j, 0]], dtype=complex)
    B = np.cos(np.pi/4.0) * np.eye(2, dtype=complex) + 1j * np.sin(np.pi/4.0) * sy
    return B


def fidelity_like(Ueff, B):
    return np.real(np.trace(Ueff.conj().T @ B)) / 2.0


def optimize_phases(Ueff, B, na=181):
    # grid search on [-pi,pi] for two angles
    a_vals = np.linspace(-np.pi, np.pi, na)
    b_vals = np.linspace(-np.pi, np.pi, na)
    best = {'fid': -1e9, 'a': 0.0, 'b': 0.0, 'Ueff': Ueff}
    for a in a_vals:
        ea = np.exp(1j * a)
        for b in b_vals:
            eb = np.exp(1j * b)
            D = np.diag([ea, eb])
            Uc = D @ Ueff
            fid = fidelity_like(Uc, B)
            if fid > best['fid']:
                best.update({'fid': fid, 'a': a, 'b': b, 'Ueff': Uc})
    return best


def lift_and_apply(U_full, Vlow, D):
    # Lift D (2x2) to full space via Vlow: V = Vlow @ D @ Vlow^dagger
    V = Vlow @ D @ Vlow.conj().T
    return V @ U_full


def main():
    if len(sys.argv) < 2:
        print('Usage: phase_compensate.py path/to/baxterize_demo_...npz')
        sys.exit(1)
    p = Path(sys.argv[1])
    data = np.load(p, allow_pickle=True)
    Ueff = data['Ueff']
    U = data['U'] if 'U' in data else None
    Vlow = data['Vlow0'] if 'Vlow0' in data else None

    B = ideal_braid_2x2()

    base_fid = fidelity_like(Ueff, B)
    print('Base fid_like=', base_fid)

    best = optimize_phases(Ueff, B, na=121)
    print('Best projected fid_like=', best['fid'])
    print('Best phases (a,b)=', best['a'], best['b'])

    # lifted test if full U and Vlow available
    if U is not None and Vlow is not None:
        D = np.diag([np.exp(1j*best['a']), np.exp(1j*best['b'])])
        U_new = lift_and_apply(U, Vlow, D)
        Ueff_new = Vlow.conj().T @ U_new @ Vlow
        fid_new = fidelity_like(Ueff_new, B)
        print('Lifted fid_like after applying lifted D=', fid_new)
    else:
        print('Full U or Vlow0 not available in npz; lifted test skipped')


if __name__ == '__main__':
    main()
