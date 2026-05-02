#!/usr/bin/env python3
"""Check that R(u)=exp(i u H_P) realizes the braiding B(gamma1,gamma4)
for the adjacent Majorana pair mapping.
"""
import numpy as np

# Pauli
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

def kron(A,B):
    return np.kron(A,B)

def mat_exp(M):
    # use eigh for Hermitian; M may be non-Hermitian but here should be unitary generator
    lam, V = np.linalg.eig(M)
    return (V @ np.diag(np.exp(lam)) @ np.linalg.inv(V))

def exp_iuH(H, u=1.0):
    # compute exp(i u H) using eigendecomposition for Hermitian H
    lam, V = np.linalg.eigh(H)
    return V @ np.diag(np.exp(1j * u * lam)) @ V.conj().T

if __name__ == '__main__':
    # target: B = exp(pi/4 * gamma1 gamma4)
    # mapping: gamma1 gamma4 = -i (sx\otimes sx - sy\otimes sy)
    M_gamma14 = -1j * (kron(sx,sx) - kron(sy,sy))
    B_target = mat_exp((np.pi/4)*M_gamma14)

    # construct H_P such that exp(i H_P) == B_target
    # H_P = -(pi/4) * (sx sx - sy sy)
    H_P = - (np.pi/4.0) * (kron(sx,sx) - kron(sy,sy))
    R_from_HP = exp_iuH(H_P, u=1.0)

    # compare
    diff = np.linalg.norm(R_from_HP - B_target)
    print(f'||R_from_HP - B_target|| = {diff:.3e}')

    # Print whether they are equal up to numerical tolerance
    tol = 1e-10
    if diff < tol:
        print('Matched: R(u=1) equals B_target (within tolerance).')
    else:
        print('Mismatch: difference above tolerance.')

    # Also print unitarity checks
    e1 = np.linalg.norm(R_from_HP.conj().T @ R_from_HP - np.eye(4))
    e2 = np.linalg.norm(B_target.conj().T @ B_target - np.eye(4))
    print(f'Unitarity errors: R: {e1:.3e}, B: {e2:.3e}')
