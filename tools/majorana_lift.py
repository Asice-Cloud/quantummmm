#!/usr/bin/env python3
"""Lift spin Hamiltonians (from R(u)) to Majorana space.
- compute_Hs_from_R(delta,N) from verify_from_R
- build JW annihilation operators c_j, Majorana gammas
- project each spin H (4x4) onto bilinear basis B_ab = (i/2) gamma_a gamma_b
  to obtain real antisymmetric Theta (4x4)
- construct O = Texp(\int Theta du) and U_mb = Texp(-i\int H du)
- verify U_mb gamma U_mb^dagger = O gamma (within numerical tolerance)

Saves results to results/majorana_delta{delta:.3g}.npz and prints diagnostics.
"""
import os
import argparse
import numpy as np
from scipy.linalg import expm, norm

# import helper from tools
import verify_from_R as vr

# Pauli matrices
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def jw_annihilation(L, j):
    """Return matrix for annihilation c_j on 2^L Hilbert space (1-based j)."""
    op = None
    for k in range(1, L + 1):
        if k < j:
            term = Z
        elif k == j:
            term = 0.5 * (X - 1j * Y)
        else:
            term = I2
        op = term if op is None else np.kron(op, term)
    return op


def majorana_list(L):
    cs = [jw_annihilation(L, j + 1) for j in range(L)]
    gam = []
    for c in cs:
        cd = c.conj().T
        gam.append(c + cd)               # gamma_{2j-1}
        gam.append(-1j * (c - cd))      # gamma_{2j}
    return gam


def build_B_basis(gammas):
    n = len(gammas)
    d = gammas[0].shape[0]
    pairs = []
    Bmats = []
    for a in range(n):
        for b in range(a + 1, n):
            B = 0.5j * (gammas[a] @ gammas[b])  # (i/2) gamma_a gamma_b
            pairs.append((a, b))
            Bmats.append(B)
    Bmats = np.array(Bmats)
    return pairs, Bmats


def H_to_Theta(H, gammas, tol=1e-12):
    """Project H (dxd) onto basis B_ab to extract Theta (n x n antisymmetric).
    Returns Theta (real antisymmetric), coefficients x, recon_error.
    """
    pairs, Bmats = build_B_basis(gammas)
    d = H.shape[0]
    m = Bmats.shape[0]
    A = Bmats.reshape((m, d * d)).T  # cols are flattened B_k
    h = H.reshape((d * d,))
    # solve least squares A x = h
    x, *_ = np.linalg.lstsq(A, h, rcond=None)
    # reconstruct
    Hrec = (A @ x).reshape((d, d))
    recon_err = float(norm(H - Hrec, ord='fro'))
    # build Theta
    n = len(gammas)
    Theta = np.zeros((n, n), dtype=float)
    for k, (a, b) in enumerate(pairs):
        val = np.real_if_close(x[k])
        val = float(np.real(val))
        Theta[a, b] = val
        Theta[b, a] = -val
    # small cleanup
    Theta[np.abs(Theta) < tol] = 0.0
    return Theta, x, recon_err


def compute_O_from_Theta_list(Theta_list, dt):
    O = np.eye(Theta_list[0].shape[0])
    for Th in Theta_list:
        O = expm(Th * dt) @ O
    return O


def compute_Umb_from_Hs(Hs, dt):
    U = np.eye(Hs[0].shape[0], dtype=complex)
    for H in Hs:
        U = expm(-1j * H * dt) @ U
    return U


def verify_action(U_mb, gammas, O):
    n = len(gammas)
    d = gammas[0].shape[0]
    errs = []
    for a in range(n):
        LHS = U_mb @ gammas[a] @ U_mb.conj().T
        RHS = sum(float(O[a, b]) * gammas[b] for b in range(n))
        errs.append(norm(LHS - RHS, ord='fro'))
    errs = np.array(errs)
    return errs


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--delta', type=float, default=0.015)
    p.add_argument('--N', type=int, default=600)
    p.add_argument('--L', type=int, default=2)
    args = p.parse_args()

    delta = args.delta
    N = args.N
    L = args.L
    os.makedirs('results', exist_ok=True)

    print('Computing H(u) from R(u) for delta=', delta)
    us, du, H4s = vr.compute_Hs_from_R(delta, N=N)
    dt = du

    print('Building Majorana operators for L=', L)
    gammas = majorana_list(L)
    nmaj = len(gammas)

    print('Projecting each H to Theta basis (this may take a moment)')
    Theta_list = []
    recon_errs = []
    for i, H in enumerate(H4s):
        Th, x, rerr = H_to_Theta(H, gammas)
        Theta_list.append(Th)
        recon_errs.append(rerr)
    Theta_list = np.array(Theta_list)
    recon_errs = np.array(recon_errs)

    print('Constructing O from Theta list')
    O = compute_O_from_Theta_list(Theta_list, dt)

    print('Constructing many-body U from H4s')
    U_mb = compute_Umb_from_Hs(H4s, dt)

    print('Verifying conjugation action on Majoranas')
    errs = verify_action(U_mb, gammas, O)
    max_err = float(np.max(errs))
    mean_err = float(np.mean(errs))

    # orthogonality of O
    ortho_err = float(norm(O.T @ O - np.eye(nmaj), ord='fro'))
    detO = float(np.linalg.det(O))

    outnpz = f'results/majorana_delta{delta:.3g}.npz'
    np.savez(outnpz,
             us=us, du=du, Theta_list=Theta_list, O=O, U_mb=U_mb,
             recon_errs=recon_errs, verify_errs=errs)

    print('Saved results to', outnpz)
    print(f'recon_errs: min={recon_errs.min():.3e}, max={recon_errs.max():.3e}')
    print(f'verify_errs: mean={mean_err:.3e}, max={max_err:.3e}')
    print(f'O orthogonality err (Fro): {ortho_err:.3e}, det(O)={detO:.6f}')


if __name__ == '__main__':
    main()
