#!/usr/bin/env python3
"""Construct R(u)=exp(i u K) from coefficients in ybe.md, compute log and
expand recovered generator in Pauli⊗Pauli basis to test locality.

Saves results to results/baxter_log_ybe.npz and prints a concise report.
"""
import re
from pathlib import Path
import numpy as np
from scipy.linalg import expm, logm, eigh

PAULIS = {
    '0': np.eye(2, dtype=complex),
    'x': np.array([[0, 1], [1, 0]], dtype=complex),
    'y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'z': np.array([[1, 0], [0, -1]], dtype=complex),
}

def pauli_tensor(mu, nu):
    return np.kron(PAULIS[mu], PAULIS[nu])


def parse_coeffs_from_ybe(path):
    text = Path(path).read_text()
    # regex for patterns like c_xx = 1.1 or c_xx: 1.1
    pat = re.compile(r'c_([0xyz]{1,2})\s*[:=]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)')
    found = dict()
    for m in pat.finditer(text):
        key = m.group(1)
        val = float(m.group(2))
        found[key] = val
    # provide defaults if missing (example defaults from ybe.md)
    C = {k: 0.0 for k in [a+b for a in ['0','x','y','z'] for b in ['0','x','y','z']]}
    defaults = {'xx':1.1, 'yy':1.0}
    for k,v in defaults.items():
        if k in found:
            C[k] = found[k]
        else:
            C[k] = defaults[k]
    # override with any parsed
    for k,v in found.items():
        C[k] = v
    return C


def build_K_from_C(C):
    # K = sum_{mu,nu} c_{mu
    K = np.zeros((4,4), dtype=complex)
    for mu in ['0','x','y','z']:
        for nu in ['0','x','y','z']:
            key = mu+nu
            c = C.get(key, 0.0)
            if abs(c) > 0:
                K += c * pauli_tensor(mu, nu)
    # ensure Hermitian
    K = 0.5*(K + K.conj().T)
    return K


def expand_in_pauli(H):
    coeffs = {}
    for mu in ['0','x','y','z']:
        for nu in ['0','x','y','z']:
            P = pauli_tensor(mu, nu)
            val = 0.25 * np.trace(P.conj().T @ H)
            coeffs[mu+nu] = val
    return coeffs


def locality_metrics(coeffs):
    # nonlocal components per ybe.md: single-side x/y and mixed with z
    nonlocal_keys = ['x0','y0','0x','0y','xz','yz','zx','zy']
    total_norm = np.sqrt(sum(abs(v)**2 for v in coeffs.values()))
    nonlocal_norm = np.sqrt(sum(abs(coeffs[k])**2 for k in nonlocal_keys))
    return total_norm, nonlocal_norm, nonlocal_keys


def main():
    repo = Path(__file__).resolve().parents[1]
    ybe_path = repo / 'ybe.md'
    if not ybe_path.exists():
        print('ybe.md not found at', ybe_path)
        return

    C = parse_coeffs_from_ybe(ybe_path)
    K = build_K_from_C(C)

    # report K eigenvalues
    evals, _ = eigh(K)
    print('K eigenvalues:', np.round(evals, 6))

    results = {}
    us = [0.1, 0.5, 1.0, np.pi/4.0]
    for u in us:
        R = expm(1j * u * K)
        # compute log and recover generator: Krec = -i log(R)/u
        L = logm(R)
        Hrec = -1j * L / u
        # measure reconstruction error vs original K (up to numerical)
        err = np.linalg.norm(Hrec - K)
        coeffs = expand_in_pauli(Hrec)
        total_norm, nonlocal_norm, nonlocal_keys = locality_metrics(coeffs)
        print(f'u={u:.6f}: ||Hrec - K||_F = {err:.6e}')
        print('Recovered pauli coefficients (Hrec):')
        for k in sorted(coeffs.keys()):
            v = coeffs[k]
            print(f' {k:2s}: {v.real: .6f} {v.imag:+.6f}j')
        print('nonlocal norm / total norm =', nonlocal_norm, '/', total_norm)
        print('---')
        results[f'u={u}'] = {'err': err, 'coeffs': coeffs, 'total_norm': total_norm, 'nonlocal_norm': nonlocal_norm}

    out = repo / 'results'
    out.mkdir(exist_ok=True)
    np.savez(out / 'baxter_log_ybe.npz', K=K, C=C, results=results)
    print('Saved results to', out / 'baxter_log_ybe.npz')


if __name__ == '__main__':
    main()
