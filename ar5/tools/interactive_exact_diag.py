#!/usr/bin/env python3
"""Exact diagonalization for small chains including extra Pauli terms (interactive check).

Builds full spin Hamiltonian for N sites from local 2-site h (from R_xxz mapping)
and optional extra Pauli coefficients supplied on CLI. Diagonalizes lowest eigenstates,
computes local fermion density via JW mapping `n_j = (1 - Z_j)/2`, and single-particle
matrix elements <psi_exc| c_j^† |psi_gs> to detect single-particle-like localized excitations.
"""
import argparse
import sys
import pathlib
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from math import log
# ensure project root is on sys.path so imports like `tools.xxz_R_and_H` work
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from tools.xxz_R_and_H import sx, sy, sz, I2, R_xxz, dRdu, expand_on_pauli


def kron_n(op_list):
    out = sp.csr_matrix([1.0])
    for M in op_list:
        out = sp.kron(out, sp.csr_matrix(M), format='csr')
    return out


def one_site_ops(N):
    ops = {'X':[], 'Y':[], 'Z':[], 'I':[]}
    for i in range(N):
        pre = [I2]*i
        post = [I2]*(N-i-1)
        ops['X'].append(kron_n(pre + [sx] + post))
        ops['Y'].append(kron_n(pre + [sy] + post))
        ops['Z'].append(kron_n(pre + [sz] + post))
        ops['I'].append(kron_n(pre + [I2] + post))
    return ops


def build_chain_h(N, h_bond_coeffs, extra_coeffs=None):
    # h_bond_coeffs: mapping like 'X_X'->val from expand_on_pauli for a single bond
    # extra_coeffs: dict of additional c_ab keys to add (applied uniformly)
    ops = one_site_ops(N)
    H = sp.csr_matrix((2**N, 2**N), dtype=complex)
    for j in range(N-1):
        hb = sp.csr_matrix((2**N, 2**N), dtype=complex)
        for key, val in h_bond_coeffs.items():
            a,b = key.split('_')
            A = ops[a][j]
            B = ops[b][j+1]
            hb += val * (A.dot(B))
        if extra_coeffs:
            for key, val in extra_coeffs.items():
                a,b = key.split('_')
                A = ops[a][j]
                B = ops[b][j+1]
                hb += val * (A.dot(B))
        H += hb
    return H


def jw_c_ops(N):
    # return list of c_j operators in spin basis
    ops = one_site_ops(N)
    c_ops = []
    for j in range(N):
        # build Z string for k<j (use full-dimension identity for j==0)
        dim = 2**N
        if j == 0:
            Zstr = sp.identity(dim, format='csr')
        else:
            Zstr = ops['Z'][0]
            for k in range(1, j):
                Zstr = Zstr.dot(ops['Z'][k])
        # sigma^- = (X - i Y)/2
        s_minus = (ops['X'][j] - 1j * ops['Y'][j]) * 0.5
        c_j = Zstr.dot(s_minus)
        c_ops.append(c_j)
    return c_ops


def compute_local_density(psi, N):
    # n_j = (1 - Z_j)/2
    ops = one_site_ops(N)
    dens = np.zeros(N)
    for j in range(N):
        Z = ops['Z'][j]
        exp = np.vdot(psi, (Z.dot(psi)))
        dens[j] = 0.5*(1 - exp.real)
    return dens


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--N', type=int, default=8)
    p.add_argument('--eta', type=float, default=0.6)
    p.add_argument('--extra', nargs='*', help='Extra Pauli coeffs like X_Z=0.2')
    p.add_argument('--k', type=int, default=6, help='number of eigenpairs to compute')
    args = p.parse_args()

    N = args.N
    eta = args.eta

    # get local h from xxz
    R0 = R_xxz(0.0, eta)
    dR0 = dRdu(0.0, eta)
    rho = np.sin(eta)
    # permutation P
    P = np.zeros((4,4), dtype=complex)
    for a in range(2):
        for b in range(2):
            in_idx = (a<<1) | b
            out_idx = (b<<1) | a
            P[out_idx, in_idx] = 1.0
    h_local = P @ dR0 / rho
    mapping = expand_on_pauli(h_local)
    # mapping keys 'X_X' etc

    extra = {}
    if args.extra:
        for s in args.extra:
            if '=' in s:
                k,v = s.split('=')
                k = k.replace('-', '_').replace(',', '_')
                extra[k] = float(v)

    H = build_chain_h(N, mapping, extra_coeffs=extra)
    H = (H + H.getH())*0.5

    print(f'Constructed H for N={N}, dim={2**N}')
    k = min(args.k, max(2, 2**N-2))
    eigvals, eigvecs = spla.eigsh(H, k=k, which='SA')
    idx = np.argsort(eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    print('Lowest eigenvalues:')
    for i,e in enumerate(eigvals[:6]):
        print(f'{i}: {e:.6f}')

    psi0 = eigvecs[:,0]
    psi1 = eigvecs[:,1]
    dens0 = compute_local_density(psi0, N)
    dens1 = compute_local_density(psi1, N)
    print('\nLocal densities (ground state):', np.round(dens0,4))
    print('Local densities (first excited):', np.round(dens1,4))

    # build JW c ops and compute overlaps
    c_ops = jw_c_ops(N)
    overlaps = np.zeros(N, dtype=complex)
    for j in range(N):
        val = np.vdot(psi1, c_ops[j].dot(psi0))
        overlaps[j] = val
    print('\nSingle-particle overlaps |<psi1|c_j^\u2020|psi0>|:')
    for j in range(N):
        print(f'j={j}: {abs(overlaps[j]):.4e}')

if __name__ == '__main__':
    main()
