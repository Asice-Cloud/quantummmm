#!/usr/bin/env python3
"""Classify small two-site Pauli coefficients c_{αβ} into XXX/XXZ/XYZ/free-fermion.

Usage:
  python tools/classify_R_constraints.py --demo
  python tools/classify_R_constraints.py --json my_coeffs.json

The script builds h = sum_{α,β} c_{αβ} σ^α ⊗ σ^β and runs tests:
 - SU(2) (XXX): [h,S^a_total]=0 for a=x,y,z
 - U(1)  (XXZ): [h,S^z_total]=0 and c_xx≈c_yy
 - XYZ: fallback when anisotropic
 - Free-fermion: whether h can be expressed as a Majorana bilinear i/4 sum A_{jk} γ_j γ_k

This is intended as a diagnostic tool with numeric tolerances.
"""
from __future__ import annotations
import json
import argparse
import math
import numpy as np
from typing import Dict, Tuple


PAULI = {
    '0': np.array([[1.0, 0.0], [0.0, 1.0]], dtype=complex),
    'x': np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    'y': np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    'z': np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
}


def build_h(c: Dict[str, float]) -> np.ndarray:
    """Build 4x4 Hamiltonian matrix from c dict with keys like 'xx','0z', etc."""
    H = np.zeros((4, 4), dtype=complex)
    for k, v in c.items():
        if len(k) != 2:
            continue
        a, b = k[0], k[1]
        pa = PAULI[a]
        pb = PAULI[b]
        H += v * np.kron(pa, pb)
    return H


def frob_norm(A: np.ndarray) -> float:
    return np.linalg.norm(A, ord='fro')


def comm_norm(A: np.ndarray, B: np.ndarray) -> float:
    return frob_norm(A @ B - B @ A)


def check_su2(H: np.ndarray, tol: float = 1e-8) -> Tuple[bool, Dict[str, float]]:
    sx = np.kron(PAULI['x'], PAULI['0']) + np.kron(PAULI['0'], PAULI['x'])
    sy = np.kron(PAULI['y'], PAULI['0']) + np.kron(PAULI['0'], PAULI['y'])
    sz = np.kron(PAULI['z'], PAULI['0']) + np.kron(PAULI['0'], PAULI['z'])
    n1 = comm_norm(H, sx)
    n2 = comm_norm(H, sy)
    n3 = comm_norm(H, sz)
    return (n1 < tol and n2 < tol and n3 < tol), {'comm_x': n1, 'comm_y': n2, 'comm_z': n3}


def check_u1(H: np.ndarray, tol: float = 1e-8) -> Tuple[bool, float]:
    sz = 0.5 * (np.kron(PAULI['z'], PAULI['0']) + np.kron(PAULI['0'], PAULI['z']))
    cn = comm_norm(H, sz)
    return cn < tol, cn


def check_free_fermion(H: np.ndarray, tol: float = 1e-8) -> Tuple[bool, float]:
    """Test if H can be expressed as i/4 sum_{j<k} A_{jk} γ_j γ_k + const.

    Use JW ordering: γ1=σ^x⊗I, γ2=σ^y⊗I, γ3=σ^z⊗σ^x, γ4=σ^z⊗σ^y.
    Form basis B_m = (i/4) γ_j γ_k for j<k (6 matrices). Solve linear system.
    """
    # Define Majoranas
    g1 = np.kron(PAULI['x'], PAULI['0'])
    g2 = np.kron(PAULI['y'], PAULI['0'])
    g3 = np.kron(PAULI['z'], PAULI['x'])
    g4 = np.kron(PAULI['z'], PAULI['y'])
    gammas = [g1, g2, g3, g4]
    basis = []
    labels = []
    for j in range(4):
        for k in range(j + 1, 4):
            B = 1j / 4.0 * (gammas[j] @ gammas[k])
            basis.append(B)
            labels.append((j + 1, k + 1))
    # Flatten basis to vectors
    M = np.column_stack([B.reshape(16) for B in basis])
    h_vec = H.reshape(16)
    # Least squares
    coeffs, residuals, rank, s = np.linalg.lstsq(M, h_vec, rcond=None)
    approx = M @ coeffs
    res_norm = np.linalg.norm(h_vec - approx)
    rel = res_norm / max(1e-12, np.linalg.norm(h_vec))
    return (rel < tol), float(rel)


def summarize(c: Dict[str, float]) -> None:
    H = build_h(c)
    normH = frob_norm(H)
    print("Hamiltonian Frobenius norm:", normH)
    is_su2, su2_info = check_su2(H)
    is_u1, u1_comm = check_u1(H)
    is_ff, ff_rel = check_free_fermion(H)

    print('\nDiagnostics:')
    print(' - SU(2) (XXX) check:', is_su2, su2_info)
    print(' - U(1) (XXZ) check (comm with Sz):', is_u1, 'comm_norm=', u1_comm)
    print(' - Free-fermion (Majorana bilinear) check:', is_ff, 'rel_res=', ff_rel)

    # Additional quick checks from coefficients
    c_xx = c.get('xx', 0.0)
    c_yy = c.get('yy', 0.0)
    c_zz = c.get('zz', 0.0)
    print('\nCoefficient quick view:')
    print(' c_xx, c_yy, c_zz =', c_xx, c_yy, c_zz)
    if is_su2:
        print('\nClassification: XXX (SU(2) symmetric)')
    elif is_u1 and abs(c_xx - c_yy) < 1e-6:
        print('\nClassification: XXZ (U(1) symmetric, c_xx≈c_yy)')
    elif is_ff:
        print('\nClassification: Free-fermion subfamily (expressible as Majorana bilinear)')
    else:
        print('\nClassification: XYZ or anisotropic two-site Hamiltonian (no full SU(2)/U(1)/free-fermion match)')


def demo() -> None:
    eta = 0.6
    c = {
        'xx': 1.0 / (2.0 * math.sin(eta)),
        'yy': 1.0 / (2.0 * math.sin(eta)),
        'zz': 0.5 * (1.0 / math.tan(eta)),
    }
    print('Demo: XXZ-like coefficients with eta=', eta)
    summarize(c)


def load_json(path: str) -> Dict[str, float]:
    with open(path, 'r') as f:
        data = json.load(f)
    return {k: float(v) for k, v in data.items()}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--demo', action='store_true')
    p.add_argument('--json', type=str, help='Path to JSON with c_{uv} dict')
    args = p.parse_args()
    if args.demo:
        demo()
        return
    if args.json:
        c = load_json(args.json)
        summarize(c)
        return
    p.print_help()


if __name__ == '__main__':
    main()
