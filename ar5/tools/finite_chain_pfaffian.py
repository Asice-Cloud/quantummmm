#!/usr/bin/env python3
"""Compute finite-chain Majorana Pfaffian for top candidates.

Transforms BdG -> Majorana basis and computes Pfaffian of the
resulting real skew-symmetric matrix. Outputs JSON with Pfaffian
values for the top-N (by min_abs) candidates.
"""
import os
import json
import numpy as np
import importlib.util

# load build_bdg_from_params from scan_all_mixtures
here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location('scan_all_mixtures', os.path.join(here, 'scan_all_mixtures.py'))
scan_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan_mod)
build_bdg_from_params = scan_mod.build_bdg_from_params


def pfaffian(A):
    A = A.copy().astype(complex)
    n = A.shape[0]
    assert n % 2 == 0
    pf = 1.0 + 0j
    for k in range(0, n-1, 2):
        # pivot at (k,k+1)
        if abs(A[k, k+1]) < 1e-16:
            found = False
            for j in range(k+2, n):
                if abs(A[k, j]) > 1e-16:
                    # swap cols k+1 and j
                    A[:, [k+1, j]] = A[:, [j, k+1]]
                    A[[k+1, j], :] = A[[j, k+1], :]
                    pf *= -1
                    found = True
                    break
            if not found:
                return 0.0 + 0j
        a = A[k, k+1]
        pf *= a
        inv = 1.0 / a
        # Schur complement update
        for i in range(k+2, n):
            for j in range(i+1, n):
                A[i, j] = A[i, j] - A[k, i] * A[k+1, j] * inv + A[k, j] * A[k+1, i] * inv
                A[j, i] = -A[i, j]
    return pf


def majorana_matrix_from_bdG(H):
    # H is 2N x 2N BdG in Nambu (c, c^dagger) ordering
    n2 = H.shape[0]
    assert n2 % 2 == 0
    N = n2 // 2
    I = np.eye(N, dtype=complex)
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    # U maps Nambu -> Majorana: gamma = U^ psi, with U = 1/sqrt2 [[I, I], [-iI, iI]]
    U = np.zeros((2*N, 2*N), dtype=complex)
    U[0:N, 0:N] = I * inv_sqrt2
    U[0:N, N:2*N] = I * inv_sqrt2
    U[N:2*N, 0:N] = -1j * I * inv_sqrt2
    U[N:2*N, N:2*N] = 1j * I * inv_sqrt2
    # Transform
    M = -1j * U.conj().T @ H @ U
    # Force real skew-symmetric (numerical cleanup)
    M = 0.5 * (M - M.T)
    M_real = np.real_if_close(M, tol=1e-8)
    return M_real


def main():
    os.makedirs('results', exist_ok=True)
    data = json.load(open('results/scan_all_mixtures_validated.json'))
    # sort by min_abs ascending
    data_sorted = sorted(data, key=lambda d: d.get('min_abs', 1e9))
    topk = 10
    report = []
    for i, d in enumerate(data_sorted[:topk]):
        t = complex(d['t'])
        Delta = complex(d['Delta'])
        mu = complex(d['mu'])
        N = 120
        H = build_bdg_from_params(N, t, Delta, mu)
        M = majorana_matrix_from_bdG(H)
        pf = pfaffian(M)
        report.append({'index': i, 'candidate': d, 'pfaffian': complex(pf), 'pf_abs': float(abs(pf))})
    with open('results/finite_chain_pfaffian_top10.json', 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print('Pfaffian computed for top', topk)


if __name__ == '__main__':
    main()
