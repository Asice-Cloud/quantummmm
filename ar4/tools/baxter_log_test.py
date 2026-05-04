#!/usr/bin/env python3
"""Numerical Baxterization / log test.

Build K = sum_{mu,nu} c_{mu,nu} sigma^mu ⊗ sigma^nu (4x4), form R(u)=exp(i u K),
compute log(R) and compare reconstructed generator to K. Output Pauli-basis
coefficients for inspection.
"""
import numpy as np
from scipy.linalg import expm, logm, eigvalsh
from pathlib import Path


pauli = {
    '0': np.array([[1.0,0.0],[0.0,1.0]], dtype=complex),
    'x': np.array([[0.0,1.0],[1.0,0.0]], dtype=complex),
    'y': np.array([[0.0,-1j],[1j,0.0]], dtype=complex),
    'z': np.array([[1.0,0.0],[0.0,-1.0]], dtype=complex),
}

order = ['x','y','z','0']


def build_K_from_c(c):
    # c is dict with keys like 'xx','yy',... also '00','x0','0x', etc.
    K = np.zeros((4,4), dtype=complex)
    for a in ['0','x','y','z']:
        for b in ['0','x','y','z']:
            key = f'{a}{b}'
            val = c.get(key, 0.0)
            if abs(val) == 0.0:
                continue
            pa = pauli[a]
            pb = pauli[b]
            K += val * np.kron(pa, pb)
    # enforce Hermiticity
    K = 0.5 * (K + K.conj().T)
    return K


def pauli_basis_coeffs(H):
    # coefficients a_{ab} = 1/4 Tr[(sigma_a ⊗ sigma_b)^ H]
    coeffs = {}
    for a in ['x','y','z','0']:
        for b in ['x','y','z','0']:
            pa = pauli[a]
            pb = pauli[b]
            M = np.kron(pa, pb)
            coeff = 0.25 * np.trace(M.conj().T @ H)
            coeffs[f'{a}{b}'] = coeff
    return coeffs


def pretty_print_coeffs(coeffs):
    keys = sorted(coeffs.keys())
    for k in keys:
        v = coeffs[k]
        if abs(v) < 1e-12:
            v = 0.0
        print(f'{k:>3}: {v.real: .6f} {v.imag:+.6f}j')


def main():
    # example c from ybe.md
    c = {
        'xx': 1.1,
        'yy': 1.0,
        'xy': 0.0,
        'yx': 0.0,
        'zz': 0.0,
        'z0': 0.0,
        '0z': 0.0,
        '00': 0.0,
    }

    K = build_K_from_c(c)
    print('K eigenvalues:', eigvalsh(K))

    # test various u values
    us = [0.1, 0.5, 1.0, np.pi/4.0]
    outdir = Path('results')
    outdir.mkdir(exist_ok=True)

    for u in us:
        R = expm(1j * u * K)
        L = logm(R)
        # attempt to reconstruct generator: H_rec = -i * L / u
        Hrec = (-1j * L) / u
        # force Hermitian symmetrize
        Hrec = 0.5 * (Hrec + Hrec.conj().T)
        diff_norm = np.linalg.norm(Hrec - K)
        print(f'u={u:.6g}: ||Hrec - K||_F = {diff_norm:.6e}')
        coeffs = pauli_basis_coeffs(Hrec)
        print('Recovered pauli coefficients (Hrec):')
        pretty_print_coeffs(coeffs)
        print(' ')

    # save last R,Hrec
    np.savez(outdir / 'baxter_log_test.npz', K=K, Hrec=Hrec)


if __name__ == '__main__':
    main()
