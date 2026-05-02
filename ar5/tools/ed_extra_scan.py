#!/usr/bin/env python3
"""Scan extra Pauli couplings (X_Z, Z_X) and record max single-particle overlap.

Saves results/ed_extra_scan.json and results/ed_extra_scan_max_overlap.png
"""
import os
import sys
# ensure repo root is on sys.path so `from tools.xxz_R_and_H import ...` works when run as script
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
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
    ops = one_site_ops(N)
    c_ops = []
    dim = 2**N
    for j in range(N):
        if j == 0:
            Zstr = sp.identity(dim, format='csr')
        else:
            Zstr = ops['Z'][0]
            for k in range(1, j):
                Zstr = Zstr.dot(ops['Z'][k])
        s_minus = (ops['X'][j] - 1j * ops['Y'][j]) * 0.5
        c_j = Zstr.dot(s_minus)
        c_ops.append(c_j)
    return c_ops


def compute_local_density(psi, N):
    ops = one_site_ops(N)
    dens = np.zeros(N)
    for j in range(N):
        Z = ops['Z'][j]
        exp = np.vdot(psi, (Z.dot(psi)))
        dens[j] = 0.5*(1 - exp.real)
    return dens


def run_scan(N=8, eta=0.6, vals=np.linspace(0,0.6,7)):
    # prepare base mapping from R
    R0 = R_xxz(0.0, eta)
    dR0 = dRdu(0.0, eta)
    rho = np.sin(eta)
    P = np.zeros((4,4), dtype=complex)
    for a in range(2):
        for b in range(2):
            in_idx = (a<<1) | b
            out_idx = (b<<1) | a
            P[out_idx, in_idx] = 1.0
    h_local = P @ dR0 / rho
    mapping = expand_on_pauli(h_local)

    results = {
        'N': N,
        'eta': eta,
        'vals': vals.tolist(),
        'max_overlap': [],
        'max_pos': []
    }

    for xz in vals:
        row_overlap = []
        row_pos = []
        for zx in vals:
            extra = {'X_Z': xz, 'Z_X': zx}
            H = build_chain_h(N, mapping, extra_coeffs=extra)
            H = (H + H.getH())*0.5
            try:
                eigvals, eigvecs = spla.eigsh(H, k=6, which='SA')
            except Exception:
                eigvals, eigvecs = spla.eigsh(H, k=4, which='SA')
            idx = np.argsort(eigvals)
            eigvals = eigvals[idx]
            eigvecs = eigvecs[:, idx]
            psi0 = eigvecs[:,0]
            psi1 = eigvecs[:,1]
            c_ops = jw_c_ops(N)
            overlaps = np.array([abs(np.vdot(psi1, c_ops[j].dot(psi0))) for j in range(N)])
            row_overlap.append(float(overlaps.max()))
            row_pos.append(int(overlaps.argmax()))
        results['max_overlap'].append(row_overlap)
        results['max_pos'].append(row_pos)

    # save json
    os.makedirs('results', exist_ok=True)
    outjson = 'results/ed_extra_scan.json'
    with open(outjson, 'w') as f:
        json.dump(results, f, indent=2)

    # plot heatmap of max_overlap
    mat = np.array(results['max_overlap'])
    plt.figure(figsize=(6,5))
    plt.imshow(mat, origin='lower', extent=(vals[0], vals[-1], vals[0], vals[-1]), aspect='auto')
    plt.colorbar(label='max |<psi1|c_j^\u2020|psi0>|')
    plt.xlabel('Z_X')
    plt.ylabel('X_Z')
    plt.title(f'ED max overlap N={N}, eta={eta}')
    outpng = 'results/ed_extra_scan_max_overlap.png'
    plt.savefig(outpng, dpi=150)
    print('Wrote', outjson, 'and', outpng)


if __name__ == '__main__':
    run_scan()
