#!/usr/bin/env python3
"""Verify MZM existence from c_{mu,nu} coefficients (XX,YY,XY,YX types).
Maps c to t,Delta,mu, builds BdG single-particle Hamiltonian, diagonalizes,
and reports lowest positive eigenvalues (edge modes).
"""
import numpy as np

# mapping from doc: C = [c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z, ...]
# t = c_xx + c_yy + i (c_xy - c_yx)
# Delta = c_xx - c_yy - i (c_xy + c_yx)
# U = 4 c_zz
# mu = 4 c_zz - 2( c_z0 + c_0z )

def map_c_to_params(c_xx, c_yy, c_xy=0.0, c_yx=0.0, c_zz=0.0, c_z0=0.0, c_0z=0.0):
    t = c_xx + c_yy + 1j*(c_xy - c_yx)
    Delta = c_xx - c_yy - 1j*(c_xy + c_yx)
    U = 4.0 * c_zz
    mu = 4.0*c_zz - 2.0*(c_z0 + c_0z)
    return t, Delta, mu, U

def build_kitaev_bdg(N, t, Delta, mu):
    # Construct A and B following convention: H = c^dag A c + 1/2 (c^dag B c^dag + h.c.)
    A = np.zeros((N,N), dtype=complex)
    B = np.zeros((N,N), dtype=complex)
    for i in range(N):
        A[i,i] = -mu
        if i+1 < N:
            A[i,i+1] = -t
            A[i+1,i] = -t
            B[i,i+1] = Delta
            B[i+1,i] = -Delta
    # BdG matrix
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    Hbdg = np.concatenate((top, bottom), axis=0)
    return Hbdg

if __name__ == '__main__':
    # Example: choose coefficients with only XX,YY nonzero to produce t=1, Delta=1, mu=0
    c_xx = 1.0
    c_yy = 0.0
    c_xy = 0.0
    c_yx = 0.0
    c_zz = 0.0
    c_z0 = 0.0
    c_0z = 0.0

    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    print(f'Parameters from c: t={t}, Delta={Delta}, mu={mu}, U={U}')

    # choose chain length
    N = 80
    Hbdg = build_kitaev_bdg(N, t, Delta, mu)

    # eigenvalues
    eigs = np.linalg.eigvals(Hbdg)
    eigs = np.sort(np.real_if_close(eigs))
    # we expect eigenvalues in +/- pairs; get positive ones
    pos = eigs[eigs >= -1e-8]
    # print smallest magnitudes
    idx = np.argsort(np.abs(eigs))
    print('\nLowest 10 eigenvalues (sorted by magnitude):')
    for i in range(10):
        val = eigs[idx[i]]
        print(f'{val:.6e}')

    # Count near-zero modes
    nz_tol = 1e-6
    zero_modes = np.sum(np.abs(eigs) < nz_tol)
    print(f'\nNumber of (near-)zero eigenvalues (|E|<{nz_tol}): {zero_modes}')

    # For clarity, check gap (smallest non-zero positive eigenvalue)
    positive = np.sort(eigs[eigs > 1e-12])
    if positive.size>0:
        gap = positive[0]
        print(f'Smallest positive eigenvalue (gap): {gap:.6e}')
    else:
        print('No positive eigenvalues found (unexpected).')
