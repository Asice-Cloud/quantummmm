#!/usr/bin/env python3
"""Verify whether R=exp(i H_P) acts like a Majorana braid in the low-energy subspace.

Simple workflow (L=3 spins by default):
- build Pauli operators and H_P from per-bond c_{mu,nu} (uses defaults if none parsed)
- compute R = expm( i H_P )
- build JW fermion operators and Majoranas gamma_a
- construct ideal braid B = exp(pi/4 * gamma_p gamma_q)
- diagonalize H_P, form projector P onto lowest 2 states, compute U_eff = P R P
- compute operator fidelity ||P B P - U_eff|| / ||P B P||

Usage: python3 tools/abs_r_verify.py
"""
import re
import numpy as np
from scipy.linalg import expm, eigh, norm


def kron(*mats):
    res = np.array([[1.0]])
    for m in mats:
        res = np.kron(res, m)
    return res


X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def pauli_op(pauli, site, L):
    ops = []
    for i in range(L):
        if i == site:
            ops.append({'x':X,'y':Y,'z':Z,'0':I2}[pauli])
        else:
            ops.append(I2)
    return kron(*ops)


def build_spin_ops(L):
    sx = [pauli_op('x', i, L) for i in range(L)]
    sy = [pauli_op('y', i, L) for i in range(L)]
    sz = [pauli_op('z', i, L) for i in range(L)]
    return sx, sy, sz


def read_c_from_ybe(path):
    # Try to parse simple assignments like c_xx = 0.5 in ybe.md, fallback to None
    try:
        text = open(path, 'r', encoding='utf-8').read()
    except Exception:
        return None
    pattern = re.compile(r'c_([a-z0-9]+)\s*=\s*([-+0-9.eE]+)')
    vals = {m.group(1): float(m.group(2)) for m in pattern.finditer(text)}
    return vals if vals else None


def jw_c_op(j, L, sx, sy, sz):
    # c_j = (prod_{k<j} Z_k) * (X_j - i Y_j)/2
    # Use single-site Pauli matrices (X,Y,Z,I2) to build the Kronecker product
    ops = []
    for k in range(L):
        if k < j:
            ops.append(Z)
        elif k == j:
            ops.append((X - 1j * Y) / 2)
        else:
            ops.append(I2)
    return kron(*ops)


def make_majoranas(L, sx, sy, sz):
    c = [jw_c_op(j, L, sx, sy, sz) for j in range(L)]
    gam = []
    for j in range(L):
        gam.append(c[j] + c[j].conj().T)            # gamma_{2j+1}
        gam.append(-1j*(c[j] - c[j].conj().T))     # gamma_{2j+2}
    return gam


def default_coeffs():
    # Example: two bonds (1-2) and (2-3). keys as '1_2:c_xx' etc.
    # We pick values that give a nontrivial quadratic model but not large interactions.
    return {
        '1_2': {'c_xx': 1.0, 'c_yy': 1.0, 'c_xy': 0.0, 'c_yx': 0.0, 'c_zz': 0.0},
        '2_3': {'c_xx': 0.8, 'c_yy': 1.2, 'c_xy': 0.0, 'c_yx': 0.0, 'c_zz': 0.0},
    }


def build_H_P(L, coeffs, sx, sy, sz):
    H = np.zeros((2**L,2**L), dtype=complex)
    # bonds: (0,1) and (1,2)
    bonds = [(0,1,'1_2'), (1,2,'2_3')]
    for a,b,name in bonds:
        c = coeffs.get(name, {})
        # enumerate mu,nu
        # include terms c_xx sigma^x_a sigma^x_b, etc.
        if 'c_xx' in c:
            H += c['c_xx'] * (sx[a] @ sx[b])
        if 'c_yy' in c:
            H += c['c_yy'] * (sy[a] @ sy[b])
        if 'c_xy' in c:
            H += c['c_xy'] * (sx[a] @ sy[b])
        if 'c_yx' in c:
            H += c['c_yx'] * (sy[a] @ sx[b])
        if 'c_zz' in c:
            H += c['c_zz'] * (sz[a] @ sz[b])
    return H.real


def main():
    L = 3
    sx, sy, sz = build_spin_ops(L)
    # try to read c from ybe.md, else use defaults
    vals = read_c_from_ybe('../ybe.md')
    if vals is None:
        coeffs = default_coeffs()
    else:
        # naive mapping: expect keys like c_xx_1_2 etc (best-effort)
        coeffs = default_coeffs()
    H_P = build_H_P(L, coeffs, sx, sy, sz)
    R = expm(1j * H_P)

    # build Majoranas
    gam = make_majoranas(L, sx, sy, sz)
    # debug: show shapes/types of Majorana operators
    print("DEBUG: gam[0].shape=", getattr(gam[0], 'shape', None), "dtype=", type(gam[0]))

    # diagonalize H_P and form projector onto two lowest eigenstates
    evals, evecs = eigh(H_P)
    idx = np.argsort(evals)
    low_idx = idx[:2]
    P = evecs[:, low_idx] @ evecs[:, low_idx].conj().T
    U_eff = P @ R @ P

    # compute each Majorana's support in the low-energy subspace and pick top two
    weights = []
    for j, g in enumerate(gam):
        # projected matrix in low-energy basis (2x2)
        Gp = evecs[:, low_idx].conj().T @ (g @ evecs[:, low_idx])
        w = norm(Gp)
        weights.append((w, j))
    weights.sort(reverse=True)
    (w1, p), (w2, q) = weights[0], weights[1]
    print(f"Selected Majorana pair indices: p={p}, q={q} with weights {w1:.6f}, {w2:.6f}")
    B = expm((np.pi/4.0) * (gam[p] @ gam[q]))

    # debug shapes before computing fidelity
    PB = None
    try:
        PB = P @ B
        PBP = PB @ P
        num = norm(PBP - U_eff)
        den = norm(PBP)
        fid = num / (den + 1e-12)
    except Exception as e:
        print("Shape debug: P.shape=", P.shape)
        print("Shape debug: B.shape=", B.shape)
        print("Shape debug: U_eff.shape=", U_eff.shape)
        if PB is not None:
            print("Shape debug: (P@B).shape=", PB.shape)
        raise

    print("Eigenvalues (H_P):", np.round(evals,6))
    print("Operator fidelity ||PBP - U_eff||/||PBP|| =", fid)


if __name__ == '__main__':
    main()
