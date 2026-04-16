#!/usr/bin/env python3
"""Majorana / Jordan-Wigner mapping for 2-site diagonal-K model.

Computes JW fermionic operators c_j (as 4x4 matrices), Majorana operators
gamma_p, extracts the Majorana antisymmetric matrix A via commutators
and reports eigenvalues / zero modes for representative (Jx,Jy,Jz) points.

Usage: run without args — it will load representative points from
`scripts/representative_spectra.csv` if present, otherwise uses a small
default list.
"""
import csv
from pathlib import Path
import numpy as np
from math import cos, sin
from itertools import product

import sys

try:
    from sympy import Matrix, I, kronecker_product as kron
except Exception:
    print('Requires sympy; please activate the project venv.')
    raise


sigma0 = Matrix([[1, 0], [0, 1]])
sigma1 = Matrix([[0, 1], [1, 0]])
sigma2 = Matrix([[0, -I], [I, 0]])
sigma3 = Matrix([[1, 0], [0, -1]])


def build_K(Jx, Jy, Jz, c00=0):
    Jx = float(Jx)
    Jy = float(Jy)
    Jz = float(Jz)
    K = (c00 * kron(sigma0, sigma0)
         + Jx * kron(sigma1, sigma1)
         + Jy * kron(sigma2, sigma2)
         + Jz * kron(sigma3, sigma3))
    return Matrix(K)


def jw_operators():
    """Return c0,c1 (annihilation) as 4x4 matrices using JW mapping."""
    # site 0 operators
    X0 = kron(sigma1, sigma0)
    Y0 = kron(sigma2, sigma0)
    Z0 = kron(sigma3, sigma0)
    X1 = kron(sigma0, sigma1)
    Y1 = kron(sigma0, sigma2)
    Z1 = kron(sigma0, sigma3)

    c0 = (X0 - I * Y0) / 2
    # c1 = Z0 * (X1 - i Y1)/2
    c1 = (Z0 * (X1 - I * Y1)) / 2
    return c0, c1


def majorana_from_c(cmat):
    cd = cmat.conjugate().T
    gamma0 = (cmat + cd)
    gamma1 = -I * (cmat - cd)
    return gamma0, gamma1


def assemble_A(H, gammas):
    """Given Hamiltonian H (4x4 sympy Matrix) and list of 4 Majorana
    matrices, compute real antisymmetric A (4x4) defined by
    H = (i/4) sum_{p,q} A_{pq} gamma_p gamma_q + const.

    We use [H, gamma_r] = (i/2) sum_p A_{pr} gamma_p to solve for A.
    """
    import numpy.linalg as la

    G = [np.array(g.tolist(), dtype=complex) for g in gammas]
    Hm = np.array(H.tolist(), dtype=complex)
    n = len(G)
    A = np.zeros((n, n), dtype=float)
    # precompute overlaps Tr(g_p g_q)
    overlaps = np.array([[np.trace(G[p] @ G[q]) for q in range(n)] for p in range(n)])
    for r in range(n):
        comm = Hm @ G[r] - G[r] @ Hm
        # express comm as linear combination of G[p]
        coeffs = np.zeros(n, dtype=complex)
        for p in range(n):
            coeffs[p] = np.trace(G[p].conj().T @ comm) / overlaps[p, p]
        # from relation comm = (i/2) sum_p A_{pr} gamma_p
        A_col = (-2j) * coeffs
        # A should be real; take real part
        A[:, r] = np.real(A_col)
    # enforce antisymmetry
    A = 0.5 * (A - A.T)
    return A


def analyze_point(Jx, Jy, Jz):
    K = build_K(Jx, Jy, Jz)
    c0, c1 = jw_operators()
    g0, g1 = majorana_from_c(c0)
    g2, g3 = majorana_from_c(c1)
    gammas = [g0, g1, g2, g3]
    H = K  # treat K as Hamiltonian matrix
    A = assemble_A(H, gammas)
    evals = np.linalg.eigvals(A)
    # eigenvalues of A come in ± pairs; zero eigenvalues indicate Majorana zero modes
    zero_count = np.sum(np.isclose(evals, 0.0, atol=1e-8))
    return {
        'Jx': float(Jx), 'Jy': float(Jy), 'Jz': float(Jz),
        'A': A, 'A_eigs': np.sort(np.real(evals)), 'zero_modes': int(zero_count)
    }


def main():
    rep = []
    csvp = Path('scripts/representative_spectra.csv')
    if csvp.exists():
        with csvp.open() as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    rep.append((float(row['Jx']), float(row['Jy']), float(row['Jz'])))
                except Exception:
                    continue
    if not rep:
        # fallback: a few representative points
        rep = [(0.0, 0.0, 0.0), (0.78539816339, 0.78539816339, 0.78539816339), (0.0, 0.78539816339, 0.3926990817)]

    out = []
    for Jx, Jy, Jz in rep:
        res = analyze_point(Jx, Jy, Jz)
        out.append(res)
        print(f"Jx={Jx:.6f}, Jy={Jy:.6f}, Jz={Jz:.6f} -> zero Majorana modes: {res['zero_modes']}")
        print('A matrix:\n', res['A'])
        print('A eigenvalues:', res['A_eigs'])
        print('---')


if __name__ == '__main__':
    main()
