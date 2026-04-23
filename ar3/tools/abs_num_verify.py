#!/usr/bin/env python3
"""Construct H (spin) with local potential to create ABS, compute R=exp(iH),
project to low-energy subspace and compare U_eff with ideal Majorana braid.
"""
import numpy as np
from scipy.linalg import expm, eigh, norm
from pathlib import Path


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


def jw_c_op(j, L):
    ops = []
    for k in range(L):
        if k < j:
            ops.append(Z)
        elif k == j:
            ops.append((X - 1j * Y) / 2)
        else:
            ops.append(I2)
    res = kron(*ops)
    return res


def make_majoranas(L):
    c = [jw_c_op(j, L) for j in range(L)]
    gam = []
    for j in range(L):
        gam.append(c[j] + c[j].conj().T)
        gam.append(-1j*(c[j] - c[j].conj().T))
    return gam


def build_H_total(L, bond_coeffs, onsite_h):
    sx, sy, sz = build_spin_ops(L)
    H = np.zeros((2**L,2**L), dtype=complex)
    # bonds: j and j+1 for j=0..L-2
    for j in range(L-1):
        c = bond_coeffs.get(f'{j}_{j+1}', {})
        if 'c_xx' in c:
            H += c['c_xx'] * (sx[j] @ sx[j+1])
        if 'c_yy' in c:
            H += c['c_yy'] * (sy[j] @ sy[j+1])
        if 'c_xy' in c:
            H += c['c_xy'] * (sx[j] @ sy[j+1])
        if 'c_yx' in c:
            H += c['c_yx'] * (sy[j] @ sx[j+1])
        if 'c_zz' in c:
            H += c['c_zz'] * (sz[j] @ sz[j+1])
    # onsite h * sz
    for j,h in onsite_h.items():
        H += h * sz[j]
    return H.real


def select_majorana_pair_by_support(gam, low_evecs):
    weights = []
    for j,g in enumerate(gam):
        Gp = low_evecs.conj().T @ (g @ low_evecs)
        w = norm(Gp)
        weights.append((w,j))
    weights.sort(reverse=True)
    return weights[0][1], weights[1][1], weights


def main():
    # construct a chain L=5, with bond quadratic coefficients chosen to be in topological regime
    L = 5
    # set bond coefficients uniform (producing t>0, Delta>0)
    # use Jx=1.1, Jy=1.0 => t=Jx+Jy=2.1, Delta=Jx-Jy=0.1
    bond_coeffs = {}
    for j in range(L-1):
        bond_coeffs[f'{j}_{j+1}'] = {'c_xx':1.1, 'c_yy':1.0, 'c_xy':0.0, 'c_yx':0.0, 'c_zz':0.0}
    # add a local potential at site 2 to create ABS when desired
    onsite_h = {2: 2.0}  # strong local Zeeman -> local chemical potential effect

    H = build_H_total(L, bond_coeffs, onsite_h)
    R = expm(1j * H)

    # diagonalize H to get low-energy subspace
    evals, evecs = eigh(H)
    idx = np.argsort(evals)
    # take two lowest states
    low_idx = idx[:2]
    low_evecs = evecs[:, low_idx]
    P = low_evecs @ low_evecs.conj().T

    gam = make_majoranas(L)
    p, q, weights = select_majorana_pair_by_support(gam, low_evecs)
    print('Selected Majorana pair:', p, q)

    B = expm((np.pi/4.0) * (gam[p] @ gam[q]))
    U_eff = P @ R @ P
    PB = P @ B @ P
    fid = norm(PB - U_eff) / (norm(PB) + 1e-12)

    print('Lowest eigenvalues:', np.round(evals[idx[:4]],6))
    print('Operator fidelity ||PBP - U_eff||/||PBP|| =', fid)

    # save some diagnostics: LDOS of lowest states (particle components u)
    # compute u components from JW c operators
    # here we reuse jw representation: particle prob = |u_j|^2 from c_j^
    c_ops = [jw_c_op(j, L) for j in range(L)]
    # project c_ops into low subspace
    for k in range(2):
        vec = low_evecs[:,k]
        print(f'state{k} energy {evals[idx[k]]:.6f}')
        probs = [np.sum(np.abs(low_evecs.conj().T @ (c @ vec))**2) for c in c_ops]
        print(' particle_probs per site:', np.round(probs,4))


if __name__ == '__main__':
    main()
