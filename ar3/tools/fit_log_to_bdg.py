#!/usr/bin/env python3
"""Fit matrix-log generator of target braid into local BdG basis.

Usage: python3 tools/fit_log_to_bdg.py --h 0.15 --L 120 --t 2.1 --Delta 0.1
"""
import argparse
import numpy as np
from scipy.linalg import eigh, logm
from pathlib import Path


def build_bdg(L, t, Delta, mu_vec):
    H0 = np.zeros((L, L), dtype=complex)
    for i in range(L-1):
        H0[i, i+1] = -t
        H0[i+1, i] = -t
    for i in range(L):
        H0[i, i] += -mu_vec[i]
    D = np.zeros((L, L), dtype=complex)
    for i in range(L-1):
        D[i, i+1] = Delta
        D[i+1, i] = -Delta
    top = np.hstack([H0, D])
    bottom = np.hstack([-D.conj(), -H0.T])
    return np.vstack([top, bottom])


def make_basis_ops(L):
    ops = []
    # on-site mu variations
    for i in range(L):
        mu = np.zeros(L)
        mu[i] = 1.0
        ops.append((f'mu_{i}', build_bdg(L, t=0.0, Delta=0.0, mu_vec=mu)))
    # hopping terms (t) on bonds
    for i in range(L-1):
        # construct Hbdg with t=1 on bond i
        H0 = np.zeros((L,L), dtype=complex)
        H0[i,i+1] = -1.0
        H0[i+1,i] = -1.0
        D = np.zeros((L,L), dtype=complex)
        top = np.hstack([H0, D])
        bottom = np.hstack([-D.conj(), -H0.T])
        ops.append((f't_{i}', np.vstack([top,bottom])))
    # pairing terms Delta on bonds
    for i in range(L-1):
        H0 = np.zeros((L,L), dtype=complex)
        D = np.zeros((L,L), dtype=complex)
        D[i,i+1] = 1.0
        D[i+1,i] = -1.0
        top = np.hstack([H0, D])
        bottom = np.hstack([-D.conj(), -H0.T])
        ops.append((f'D_{i}', np.vstack([top,bottom])))
    return ops


def ideal_braid_2x2():
    # choose B = exp(i * (pi/4) * sigma_y)
    sy = np.array([[0, -1j],[1j, 0]], dtype=complex)
    B = np.cos(np.pi/4.0) * np.eye(2, dtype=complex) + 1j * np.sin(np.pi/4.0) * sy
    return B


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h', type=float, required=True)
    p.add_argument('--L', type=int, default=120)
    p.add_argument('--t', type=float, default=2.1)
    p.add_argument('--Delta', type=float, default=0.1)
    args = p.parse_args()

    L = args.L
    mid = L//2
    mu0 = 0.0
    mu_vec = np.zeros(L) + mu0
    mu_vec[mid] += 2.0 * args.h
    Hbdg = build_bdg(L, args.t, args.Delta, mu_vec)

    evals, evecs = eigh(Hbdg)
    idx = np.argsort(np.abs(evals))
    low_idx = idx[:2]
    Vlow = evecs[:, low_idx]

    # target braid in low subspace
    B2 = ideal_braid_2x2()
    # embed into full single-particle space
    Bfull = Vlow @ B2 @ Vlow.conj().T

    # compute generator G = -i log(Bfull) (Hermitian)
    Lm = logm(Bfull)
    Gfull = -1j * Lm
    # project to 2x2 basis
    Gproj = Vlow.conj().T @ Gfull @ Vlow

    # build basis operators (BdG matrices)
    basis = make_basis_ops(L)
    Mk = []
    names = []
    for name, Op in basis:
        M = Vlow.conj().T @ Op @ Vlow
        Mk.append(M)
        names.append(name)

    # Flatten complex 2x2 Hermitian matrices to real vector (4 real dims -> use real+imag)
    def mat2vec(A):
        # returns real vector of length 8: real and imag flattened
        return np.concatenate([A.real.ravel(), A.imag.ravel()])

    A = np.column_stack([mat2vec(M) for M in Mk])
    b = mat2vec(Gproj)

    # least squares solve for coefficients
    coeffs, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    approx = A @ coeffs
    err = np.linalg.norm(approx - b)
    normb = np.linalg.norm(b)

    out = Path('results')
    out.mkdir(exist_ok=True)
    with open(out / f'fit_log_h{args.h:.3f}_L{L}.txt','w') as f:
        f.write(f'Gproj norm: {normb}\n')
        f.write(f'residual norm: {err}\n')
        f.write(f'relative residual: {err/normb if normb>0 else np.nan}\n')
        f.write('\nTop coefficients (abs):\n')
        order = np.argsort(-np.abs(coeffs))
        for i in order[:20]:
            f.write(f'{names[i]}: {coeffs[i]}\n')
    print('Saved fit results to', out / f'fit_log_h{args.h:.3f}_L{L}.txt')
    print('Relative residual:', err/normb)


if __name__ == '__main__':
    main()
