#!/usr/bin/env python3
"""Numeric grid search for full constant YBE solutions R=exp(iK) with diagonal K.

Scans Jx,Jy,Jz on a specified grid and writes results to CSV.
"""
import csv
import numpy as np

def paulis():
    I2 = np.array([[1,0],[0,1]], dtype=complex)
    X = np.array([[0,1],[1,0]], dtype=complex)
    Y = np.array([[0,-1j],[1j,0]], dtype=complex)
    Z = np.array([[1,0],[0,-1]], dtype=complex)
    return I2, X, Y, Z

def exp_iK(c00, Jx, Jy, Jz):
    I2, X, Y, Z = paulis()
    K = c00 * np.kron(I2, I2) + Jx * np.kron(X, X) + Jy * np.kron(Y, Y) + Jz * np.kron(Z, Z)
    # K is Hermitian for real Jx,Jy,Jz,c00 — use eigh for stability
    vals, vecs = np.linalg.eigh(K)
    R = (vecs @ np.diag(np.exp(1j * vals)) @ vecs.conj().T)
    return R

def compute_residual_norm(R):
    I2 = np.eye(2, dtype=complex)
    R12 = np.kron(R, I2)
    R23 = np.kron(I2, R)
    # build permutation P23 to get R13 from R12
    P23 = np.zeros((8,8), dtype=complex)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                src = a*4 + b*2 + c
                dst = a*4 + c*2 + b
                P23[dst, src] = 1
    R13 = P23 @ R12 @ P23
    M1 = R12 @ R13 @ R23
    M2 = R23 @ R13 @ R12
    res = M1 - M2
    return np.linalg.norm(res)

def main():
    c00 = 0.0
    # coarse grid: multiples of pi/4 (8 points)
    pts = np.arange(0, 2*np.pi, np.pi/4)
    out_rows = []
    best = (None, 1e9)
    for Jx in pts:
        for Jy in pts:
            for Jz in pts:
                R = exp_iK(c00, Jx, Jy, Jz)
                nrm = compute_residual_norm(R)
                out_rows.append((Jx, Jy, Jz, nrm))
                if nrm < best[1]:
                    best = ((Jx, Jy, Jz), nrm)

    # write CSV
    fn = 'scripts/ybe_grid_results_coarse.csv'
    with open(fn, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Jx','Jy','Jz','residual_norm'])
        for r in out_rows:
            w.writerow(r)

    print('Wrote', fn)
    print('Best point:', best)

    # Print solutions with tiny residual
    tol = 1e-10
    sols = [r for r in out_rows if r[3] < tol]
    print('Solutions found (res <', tol, '):', len(sols))
    for s in sols:
        print(s)

if __name__ == '__main__':
    main()
