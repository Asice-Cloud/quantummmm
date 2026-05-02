#!/usr/bin/env python3
"""构造XXZ（6-vertex）R(u)，计算局域哈密顿 H_local = -i R^{-1} dR/du |_{u=0}
并将 H_local 展开到 Pauli 基上；同时数值检验参数 YBE。
"""
import numpy as np

# Pauli matrices
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

def kron(*mats):
    out = np.array([1], dtype=complex)
    for M in mats:
        out = np.kron(out, M)
    return out

def R_xxz(u, eta):
    # trigonometric 6-vertex R-matrix (up to overall factor)
    s = np.sin
    a = s(u + eta)
    b = s(u)
    c = s(eta)
    # basis ordering: 00,01,10,11
    R = np.zeros((4,4), dtype=complex)
    R[0,0] = a
    R[1,1] = b
    R[1,2] = c
    R[2,1] = c
    R[2,2] = b
    R[3,3] = a
    return R

def dRdu(u, eta):
    c = np.cos
    da = c(u + eta)
    db = c(u)
    dc = 0.0  # c does not depend on u
    dR = np.zeros((4,4), dtype=complex)
    dR[0,0] = da
    dR[1,1] = db
    dR[1,2] = dc
    dR[2,1] = dc
    dR[2,2] = db
    dR[3,3] = da
    return dR


def canonical_H_local(u, eta):
    """Return canonical local generator H_local(u) = -i R^{-1} dR/du (Hermitian symmetrized)."""
    R = R_xxz(u, eta)
    dR = dRdu(u, eta)
    if abs(np.linalg.det(R)) < 1e-12:
        raise RuntimeError('R(u) is singular or nearly singular; cannot compute canonical generator')
    Rinv = np.linalg.inv(R)
    H_local = -1j * (Rinv @ dR)
    H_local = 0.5 * (H_local + H_local.conj().T)
    return H_local

# embedding and YBE check

def swap_23():
    dim = 8
    S = np.zeros((dim,dim), dtype=complex)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                in_idx = (a<<2) | (b<<1) | c
                out_idx = (a<<2) | (c<<1) | b
                S[out_idx, in_idx] = 1.0
    return S


def embed_R_three(R):
    R12 = np.kron(R, I2)
    R23 = np.kron(I2, R)
    S23 = swap_23()
    R13 = S23 @ R12 @ S23
    return R12, R13, R23


def expand_on_pauli(H):
    # expand 4x4 H on basis {I,x,y,z} tensor {I,x,y,z}
    basis = [I2, sx, sy, sz]
    B = []
    labels = []
    for i,A in enumerate(basis):
        for j,Bb in enumerate(basis):
            M = np.kron(A, Bb)
            B.append(M.reshape(16))
            labels.append((i,j))
    Bmat = np.stack(B, axis=1)  # 16 x 16
    coeffs, *_ = np.linalg.lstsq(Bmat, H.reshape(16), rcond=None)
    # map indices to labels
    mapping = {}
    for idx, (i,j) in enumerate(labels):
        key = ['I','X','Y','Z'][i] + '_' + ['I','X','Y','Z'][j]
        mapping[key] = coeffs[idx]
    return mapping

if __name__ == '__main__':
    eta = 0.6
    R0 = R_xxz(0.0, eta)
    dR0 = dRdu(0.0, eta)
    # check invertibility
    if abs(np.linalg.det(R0)) < 1e-12:
        print('R(0) is singular or nearly singular; choose eta not multiple of pi.')
    R0inv = np.linalg.inv(R0)
    H_local = -1j * (R0inv @ dR0)
    print('H_local (4x4 matrix) =')
    np.set_printoptions(precision=4, suppress=True)
    print(H_local)

    # expand on Pauli basis
    mapping = expand_on_pauli(H_local)
    print('\nExpansion coefficients on Pauli basis (scale factors):')
    for k in sorted(mapping.keys()):
        val = mapping[k]
        if abs(val) > 1e-8:
            print(f'{k}: {val}')

    # Show relation to XXZ form: compute coefficients for XX+YY and ZZ
    coef_xx = mapping.get('X_X', 0)
    coef_yy = mapping.get('Y_Y', 0)
    coef_zz = mapping.get('Z_Z', 0)
    coef_ii = mapping.get('I_I', 0)
    print('\nSelected combinations:')
    print(f'XX coeff = {coef_xx}')
    print(f'YY coeff = {coef_yy}')
    print(f'ZZ coeff = {coef_zz}')
    print(f'Identity coeff = {coef_ii}')

    # Verify YBE numerically for this R(u)
    us = [0.1, 0.3, 0.7]
    vs = [0.2, 0.5, 1.0]
    print('\nYBE numerical check (||LHS-RHS||):')
    for u in us:
        for v in vs:
            R_uv = R_xxz(u - v, eta)
            R_u = R_xxz(u, eta)
            R_v = R_xxz(v, eta)
            R12_uv, R13_uv, R23_uv = embed_R_three(R_uv)
            R12_u, R13_u, R23_u = embed_R_three(R_u)
            R12_v, R13_v, R23_v = embed_R_three(R_v)
            LHS = R12_uv @ R13_u @ R23_v
            RHS = R23_v @ R13_u @ R12_uv
            d = np.linalg.norm(LHS - RHS)
            print(f'u={u:.3f}, v={v:.3f}, ||LHS-RHS|| = {d:.3e}')

    # Also show that H_local is local two-body operator (Hermitian up to numerical noise)
    herm_err = np.linalg.norm(H_local - H_local.conj().T)
    print(f"\nHermiticity error ||H - H^\\dagger|| = {herm_err:.3e}")

    # Construct permutation (swap) operator P on two qubits: P|a,b> = |b,a>
    P = np.zeros((4,4), dtype=complex)
    for a in range(2):
        for b in range(2):
            in_idx = (a<<1) | b
            out_idx = (b<<1) | a
            P[out_idx, in_idx] = 1.0

    # Historically the repo sometimes used the heuristic h = P * dR/du / rho
    # (evaluated at a regular point) to obtain a familiar XXZ density. That
    # construction is kept here for reference but is NOT the canonical
    # time-ordered generator. Prefer `canonical_H_local(u,eta)` for strict
    # time-evolution reconstruction:
    rho = np.sin(eta)
    h_local = P @ dR0 / rho
    herm_err_h = np.linalg.norm(h_local - h_local.conj().T)
    print(f"\n(Deprecated/heuristic) h = P dR/du|_0 / sin(eta) hermiticity error: {herm_err_h:.3e}")
    print('\n(Deprecated/heuristic) h_local (Pauli expansion):')
    mapping_h = expand_on_pauli(h_local)
    for k in sorted(mapping_h.keys()):
        val = mapping_h[k]
        if abs(val) > 1e-8:
            print(f'{k}: {val}')

    # Show canonical generator (preferred) for direct time-ordered evolution
    H_local_canonical = canonical_H_local(0.0, eta)
    print('\nCanonical H_local = -i R^{-1} dR/du (Pauli expansion):')
    mapping_c = expand_on_pauli(H_local_canonical)
    for k in sorted(mapping_c.keys()):
        val = mapping_c[k]
        if abs(val) > 1e-8:
            print(f'{k}: {val}')
