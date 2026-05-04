#!/usr/bin/env python3
"""
tools/verify_bloch_model.py

Verify the Bloch-sphere effective model described in ybe223.md.

This script:
- builds the two-spin Hamiltonian H(u) from Pauli terms
- projects to the subspace S={|01>,|10>} to get H_eff (2x2)
- extracts the Bloch vector d(u) from H_eff
- checks whether the origin lies inside the convex hull of d(u)
  (a practical numerical proxy for the "path encloses origin" criterion)

Example: eight-vertex toy model in ybe223.md with
  H(u) = cos(u) XX + sin(u) YY + delta ZZ

Usage: python3 tools/verify_bloch_model.py
"""
import numpy as np
import math
import sys

try:
    from scipy.optimize import linprog
except Exception:
    linprog = None


def paulis():
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return {'I': I, 'X': X, 'Y': Y, 'Z': Z}


def build_H_from_coeffs(hcoeffs):
    """Construct 4x4 two-spin Hamiltonian from dict hcoeffs[(a,b)] -> value
    where a,b in {'I','X','Y','Z'}.
    """
    P = paulis()
    H = np.zeros((4, 4), dtype=complex)
    for (a, b), v in hcoeffs.items():
        H += complex(v) * np.kron(P[a], P[b])
    return H


def project_to_subspace(H4):
    """Project H4 (4x4) onto subspace S = {|01>,|10>} and return 2x2 Heff
    in basis [|01>,|10>]."""
    psi01 = np.array([0, 1, 0, 0], dtype=complex)
    psi10 = np.array([0, 0, 1, 0], dtype=complex)
    states = [psi01, psi10]
    Heff = np.zeros((2, 2), dtype=complex)
    for i, psi_i in enumerate(states):
        for j, psi_j in enumerate(states):
            Heff[i, j] = np.vdot(psi_i, H4.dot(psi_j))
    return Heff


def d_from_Heff(Heff):
    """Return (d0, dx, dy, dz) for Heff = d0 I + d·σ"""
    P = paulis()
    d0 = 0.5 * np.trace(Heff)
    dx = 0.5 * np.trace(Heff.dot(P['X']))
    dy = 0.5 * np.trace(Heff.dot(P['Y']))
    dz = 0.5 * np.trace(Heff.dot(P['Z']))
    return np.real_if_close(d0), np.real_if_close(dx), np.real_if_close(dy), np.real_if_close(dz)


def origin_in_convex_hull(points, tol=1e-9):
    """Check if origin (0,0,0) belongs to convex hull of `points` using
    a linear programming feasibility test:
      find w >= 0, sum(w)=1, and sum_i w_i * points[i] = 0.

    Falls back to a simple bounding-box test if scipy.linprog not available.
    """
    pts = np.asarray(points, dtype=float)
    N = pts.shape[0]
    if N == 0:
        return False
    if linprog is None:
        # crude fallback: check if origin within axis-aligned bounding box
        mins = pts.min(axis=0)
        maxs = pts.max(axis=0)
        return np.all(mins - tol <= 0) and np.all(maxs + tol >= 0)

    # variables w (length N)
    c = np.zeros(N)
    A_eq = np.vstack([pts.T, np.ones(N)])
    b_eq = np.zeros(pts.shape[1] + 1)
    b_eq[-1] = 1.0
    bounds = [(0.0, None) for _ in range(N)]
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
    return bool(res.success)


def test_eight_vertex(deltas=(0.0, 0.3), N=400):
    results = {}
    us = np.linspace(0.0, 2.0 * math.pi, N, endpoint=False)
    for delta in deltas:
        dpoints = []
        for u in us:
            hcoeffs = {}
            # eight-vertex model toy: H = cos(u) XX + sin(u) YY + delta ZZ
            hcoeffs[('X', 'X')] = math.cos(u)
            hcoeffs[('Y', 'Y')] = math.sin(u)
            hcoeffs[('Z', 'Z')] = delta
            H4 = build_H_from_coeffs(hcoeffs)
            Heff = project_to_subspace(H4)
            _, dx, dy, dz = d_from_Heff(Heff)
            dpoints.append([float(dx), float(dy), float(dz)])
        dpoints = np.array(dpoints)
        inside = origin_in_convex_hull(dpoints)
        min_norm = np.min(np.linalg.norm(dpoints, axis=1))
        results[delta] = {'inside_convex_hull': inside, 'min|d|': float(min_norm), 'dpoints': dpoints}
    return results


def main():
    deltas = (0.0, 0.3)
    print('Verifying eight-vertex Bloch model from ybe223.md')
    res = test_eight_vertex(deltas=deltas, N=400)
    for delta, info in res.items():
        print('\ndelta =', delta)
        print('  origin in convex hull of d(u)? ->', info['inside_convex_hull'])
        print('  min |d(u)| ->', info['min|d|'])
        # print a small sample of d(u)
        print('  sample d(u) [first 6] ->')
        print(info['dpoints'][:6])


if __name__ == '__main__':
    main()
