#!/usr/bin/env python3
"""Multi-dimensional ABS scan over c_XY, c_YX, c_XZ (coarse grid).

Saves CSV results to `results/abs_scan_multi.csv` and candidates to
`results/abs_scan_multi_candidates.json` and prints a short summary.
"""
import os
import json
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
    # coarse grids
    cxy_vals = np.linspace(0.0, 0.3, 7)
    cyx_vals = np.linspace(0.0, 0.3, 7)
    cxz_vals = np.linspace(0.0, 0.2, 5)
    eta = 0.6
    u0 = 0.0
    N = 100  # smaller N for multi-run
    threshold = 1e-3

    R0 = R_xxz(u0, eta)
    dR0 = dRdu(u0, eta)
    # build P
    P = np.zeros((4, 4), dtype=complex)
    for a in range(2):
        for b in range(2):
            in_idx = (a << 1) | b
            out_idx = (b << 1) | a
            P[out_idx, in_idx] = 1.0
    rho = np.sin(eta)
    h_local = P @ dR0 / rho
    base_mapping = expand_on_pauli(h_local)

    out_csv = os.path.join('results', 'abs_scan_multi.csv')
    out_cands = []
    with open(out_csv, 'w') as outf:
        outf.write('cxy,cyx,cxz,min_abs_E\n')
        for cxy in cxy_vals:
            for cyx in cyx_vals:
                for cxz in cxz_vals:
                    mapping_mod = dict(base_mapping)
                    mapping_mod['X_Y'] = mapping_mod.get('X_Y', 0.0) + complex(cxy)
                    mapping_mod['Y_X'] = mapping_mod.get('Y_X', 0.0) + complex(cyx)
                    # add X_Z as potential extra Pauli (may not affect map_c_to_params)
                    mapping_mod['X_Z'] = mapping_mod.get('X_Z', 0.0) + complex(cxz)
                    t_mod, Delta_mod, mu_mod, U_mod = pauli_map_from_mapping(mapping_mod)
                    H = build_bdg_from_params(N, t_mod, Delta_mod, mu_mod)
                    eigs = eigvals(H)
                    min_abs = float(np.min(np.abs(np.real_if_close(eigs))))
                    outf.write(f'{cxy:.6g},{cyx:.6g},{cxz:.6g},{min_abs:.6e}\n')
                    if min_abs < threshold:
                        out_cands.append({'cxy': float(cxy), 'cyx': float(cyx), 'cxz': float(cxz), 'min_abs': min_abs})

    with open(os.path.join('results', 'abs_scan_multi_candidates.json'), 'w') as f:
        json.dump(out_cands, f, indent=2)

    print('Multi-scan complete. CSV:', out_csv)
    print('Candidates found:', len(out_cands))


if __name__ == '__main__':
    main()
