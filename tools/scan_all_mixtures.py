#!/usr/bin/env python3
"""Scan single/pair/triple Pauli-term mixtures and save candidates.

This script reuses `xxz_R_and_H.expand_on_pauli` mapping format and
`verify_mzm.map_c_to_params` to produce BdG parameters and looks for
near-zero modes on finite open chains.
"""
import os
import json
import itertools
import numpy as np
from numpy.linalg import eigvals
from xxz_R_and_H import R_xxz, dRdu, expand_on_pauli
from verify_mzm import map_c_to_params


def build_bdg_from_params(N, t, Delta, mu):
    A = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)
    for i in range(N):
        A[i, i] = -mu
        if i < N - 1:
            A[i, i + 1] = -t
            A[i + 1, i] = -t
            B[i, i + 1] = Delta
            B[i + 1, i] = -Delta
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    H = np.concatenate((top, bottom), axis=0)
    return H


def pauli_map_from_mapping(mapping):
    c_xx = mapping.get('X_X', 0.0)
    c_yy = mapping.get('Y_Y', 0.0)
    c_xy = mapping.get('X_Y', 0.0)
    c_yx = mapping.get('Y_X', 0.0)
    c_zz = mapping.get('Z_Z', 0.0)
    c_z0 = mapping.get('Z_I', 0.0)
    c_0z = mapping.get('I_Z', 0.0)
    return map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)


def main():
    os.makedirs('results', exist_ok=True)
    eta = 0.6
    u0 = 0.0
    N = 120
    threshold = 1e-4

    # build base mapping from XXZ local h
    R0 = R_xxz(u0, eta)
    dR0 = dRdu(u0, eta)
    P = np.zeros((4, 4), dtype=complex)
    for a in range(2):
        for b in range(2):
            in_idx = (a << 1) | b
            out_idx = (b << 1) | a
            P[out_idx, in_idx] = 1.0
    rho = np.sin(eta)
    h_local = P @ dR0 / rho
    base_mapping = expand_on_pauli(h_local)

    # Pauli keys to try (including cross terms)
    keys = ['X_X', 'Y_Y', 'Z_Z', 'X_Y', 'Y_X', 'X_Z', 'Y_Z', 'Z_X', 'Z_Y']

    # values grid
    vals = np.linspace(0.0, 0.3, 7)

    out_csv = os.path.join('results', 'scan_all_mixtures.csv')
    candidates = []
    with open(out_csv, 'w') as outf:
        outf.write('combo,values,min_abs_E\n')
        # single, pairs, triples
        for r in [1, 2, 3]:
            for combo in itertools.combinations(keys, r):
                # iterate grid for nonzero amplitude (same amplitude for all keys in combo)
                for v in vals:
                    mapping_mod = dict(base_mapping)
                    for k in combo:
                        mapping_mod[k] = mapping_mod.get(k, 0.0) + complex(v)
                    t_mod, Delta_mod, mu_mod, U_mod = pauli_map_from_mapping(mapping_mod)
                    H = build_bdg_from_params(N, t_mod, Delta_mod, mu_mod)
                    eigs = eigvals(H)
                    min_abs = float(np.min(np.abs(np.real_if_close(eigs))))
                    outf.write(f'{"+".join(combo)},{v:.6g},{min_abs:.6e}\n')
                    if min_abs < threshold:
                        candidates.append({'combo': list(combo), 'v': float(v), 'min_abs': min_abs,
                                           't': complex(t_mod), 'Delta': complex(Delta_mod), 'mu': complex(mu_mod)})

    with open(os.path.join('results', 'scan_all_mixtures_candidates.json'), 'w') as f:
        json.dump(candidates, f, indent=2, default=str)

    print('Scan complete. CSV:', out_csv)
    print('Candidates found:', len(candidates))


if __name__ == '__main__':
    main()
