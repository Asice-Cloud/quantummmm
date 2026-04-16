#!/usr/bin/env python3
"""Compute spectra of K and R=exp(iK) for representative J points.

Writes results to CSV for later inspection.
"""
import csv
import numpy as np

def paulis():
    I2 = np.array([[1,0],[0,1]], dtype=complex)
    X = np.array([[0,1],[1,0]], dtype=complex)
    Y = np.array([[0,-1j],[1j,0]], dtype=complex)
    Z = np.array([[1,0],[0,-1]], dtype=complex)
    return I2, X, Y, Z

def K_matrix(c00, Jx, Jy, Jz):
    I2, X, Y, Z = paulis()
    return c00 * np.kron(I2, I2) + Jx * np.kron(X, X) + Jy * np.kron(Y, Y) + Jz * np.kron(Z, Z)

def R_from_K(K):
    vals, vecs = np.linalg.eigh(K)
    R = vecs @ np.diag(np.exp(1j * vals)) @ vecs.conj().T
    return R

def analyze_point(c00, Jx, Jy, Jz):
    K = K_matrix(c00, Jx, Jy, Jz)
    evals_K = np.linalg.eigvalsh(K)
    R = R_from_K(K)
    evals_R = np.linalg.eigvals(R)
    # also compute residual norm for verification
    I2 = np.eye(2, dtype=complex)
    R12 = np.kron(R, I2)
    R23 = np.kron(I2, R)
    P23 = np.zeros((8,8), dtype=complex)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                src = a*4 + b*2 + c
                dst = a*4 + c*2 + b
                P23[dst, src] = 1
    R13 = P23 @ R12 @ P23
    res_norm = np.linalg.norm(R12 @ R13 @ R23 - R23 @ R13 @ R12)
    return evals_K, evals_R, res_norm

def main():
    c00 = 0.0
    pts = [
        (np.pi/4, np.pi/4, np.pi/4),
        (np.pi/2, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ]
    out = []
    for (Jx, Jy, Jz) in pts:
        eK, eR, res = analyze_point(c00, Jx, Jy, Jz)
        out.append((Jx, Jy, Jz, res, list(np.round(eK,12)), [complex(round(x.real,12), round(x.imag,12)) for x in eR]))

    fn = 'scripts/representative_spectra.csv'
    with open(fn, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Jx','Jy','Jz','residual_norm','eigvals_K','eigvals_R'])
        for row in out:
            w.writerow(row)

    print('Wrote', fn)
    for row in out:
        print(row)

if __name__ == '__main__':
    main()
