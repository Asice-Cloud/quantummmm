#!/usr/bin/env python3
"""检查谱参数 R(u)=exp(i u H_P) 是否满足三体 YBE：
R12(u-v) R13(u) R23(v) = R23(v) R13(u) R12(u-v)

只依赖 numpy（通过本征分解计算矩阵指数，适用于厄米 H_P）。
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

def mat_exp_iut(H, u):
    # H is Hermitian; use eigh for stable expm: exp(i u H) = V diag(exp(i u lam)) V^dagger
    lam, V = np.linalg.eigh(H)
    return (V @ np.diag(np.exp(1j * u * lam)) @ V.conj().T)

def make_R(H2, u):
    return mat_exp_iut(H2, u)

def swap_23():
    # swap qubit 2 and 3 in 3-qubit space (basis |q1 q2 q3>, q in {0,1})
    dim = 8
    S = np.zeros((dim,dim), dtype=complex)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                in_idx = (a<<2) | (b<<1) | c
                out_idx = (a<<2) | (c<<1) | b
                S[out_idx, in_idx] = 1.0
    return S

def embed_R_on_123(R):
    # R acts on qubit pair (1,2) by default for construction convenience
    R12 = np.kron(R, I2)
    R23 = np.kron(I2, R)
    S23 = swap_23()
    R13 = S23 @ R12 @ S23
    return R12, R13, R23

def check_ybe(H2, u, v, tol=1e-10):
    R_uv = make_R(H2, u - v)  # R(u-v)
    R_u = make_R(H2, u)
    R_v = make_R(H2, v)
    R12_a, R13_a, R23_a = embed_R_on_123(R_uv)
    # but embed_R_on_123 expects R for each; rebuild properly
    R12_uv, R13_uv, R23_uv = embed_R_on_123(R_uv)
    R12_u,  R13_u,  R23_u  = embed_R_on_123(R_u)
    R12_v,  R13_v,  R23_v  = embed_R_on_123(R_v)

    LHS = R12_uv @ R13_u @ R23_v
    RHS = R23_v @ R13_u @ R12_uv
    diff = np.linalg.norm(LHS - RHS)
    return diff

if __name__ == '__main__':
    # Example: H_P = a_x XX + a_y YY + a_z ZZ (these two-body terms commute pairwise)
    ax, ay, az = 1.0, 0.8, 0.5
    H2 = ax * kron(sx, sx) + ay * kron(sy, sy) + az * kron(sz, sz)

    us = [0.1, 0.5, 1.2]
    vs = [0.2, 0.3, 0.7]

    for u in us:
        for v in vs:
            d = check_ybe(H2, u, v)
            print(f"u={u:.3f}, v={v:.3f}, ||LHS-RHS|| = {d:.3e}")

    # If desired, try a noncommuting example to see deviations
    H2_nc = kron(sx, sy) + 0.3 * kron(sz, sx)
    print('\nNoncommuting example:')
    for u in [0.4, 1.0]:
        for v in [0.2, 0.9]:
            d = check_ybe(H2_nc, u, v)
            print(f"u={u:.3f}, v={v:.3f}, ||LHS-RHS|| = {d:.3e}")
