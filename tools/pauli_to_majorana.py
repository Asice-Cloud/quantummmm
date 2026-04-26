#!/usr/bin/env python3
"""Convert a BdG eigenvector (u,v) into a Pauli-string expansion of a Majorana operator Gamma.

Usage:
  python3 tools/pauli_to_majorana.py --vec path/to/vector.npy --top 20
  python3 tools/pauli_to_majorana.py --test 6 --top 12

The script prints the largest Pauli-string coefficients and writes a JSON file `results/gamma_pauli.json`.
"""
import argparse
import json
import numpy as np
import os


def build_pauli_string(N, site, op):
    # op in {'X','Y'}
    s = ['I'] * N
    for k in range(site):
        s[k] = 'Z'
    s[site] = op
    return ''.join(s)


def gamma_to_pauli(u, v):
    # u,v are length-N complex arrays (BdG eigenvector split)
    N = len(u)
    terms = []
    for j in range(N):
        # coefficients derived from Gamma = gamma + gamma^dagger
        x_coeff = float(np.real(u[j] + v[j]))
        y_coeff = float(np.imag(v[j]) - np.imag(u[j]))
        if abs(x_coeff) > 0:
            terms.append((x_coeff, build_pauli_string(N, j, 'X')))
        if abs(y_coeff) > 0:
            terms.append((y_coeff, build_pauli_string(N, j, 'Y')))
    return terms


def normalize_terms(terms):
    # L2 normalize coefficients
    coeffs = np.array([t[0] for t in terms], dtype=float)
    norm = np.linalg.norm(coeffs)
    if norm == 0:
        return terms, 0.0
    scaled = [(float(t[0] / norm), t[1]) for t in terms]
    return scaled, norm


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--vec', help='npy file containing eigenvector of length 2N')
    p.add_argument('--test', type=int, help='run quick random test with N sites')
    p.add_argument('--top', type=int, default=20, help='number of top Pauli strings to show')
    args = p.parse_args()

    if not os.path.isdir('results'):
        os.makedirs('results')

    if args.vec:
        vec = np.load(args.vec)
        L = len(vec)
        if L % 2 != 0:
            raise SystemExit('vector length must be even (2N)')
        N = L // 2
        u = vec[:N]
        v = vec[N:]
    elif args.test:
        N = args.test
        # random test vector
        u = np.random.randn(N) + 1j * np.random.randn(N)
        v = np.random.randn(N) + 1j * np.random.randn(N)
        # normalize in BdG sense
        vec = np.concatenate([u, v])
        vec = vec / np.linalg.norm(vec)
        u = vec[:N]
        v = vec[N:]
    else:
        raise SystemExit('provide --vec or --test')

    terms = gamma_to_pauli(u, v)
    terms, norm = normalize_terms(terms)
    terms.sort(key=lambda t: abs(t[0]), reverse=True)

    out = {
        'N': N,
        'norm_before': float(norm),
        'terms': [{'coeff': c, 'pauli': p} for c, p in terms]
    }
    out_path = 'results/gamma_pauli.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)

    print(f'Wrote {out_path} with {len(terms)} terms (N={N}).')
    top = min(args.top, len(terms))
    for i in range(top):
        c, pauli = terms[i]
        print(f'{i+1:2d}. {c:+.6f}  {pauli}')


if __name__ == '__main__':
    main()
