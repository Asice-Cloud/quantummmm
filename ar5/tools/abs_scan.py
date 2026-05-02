#!/usr/bin/env python3
"""ABS scan: first sweep FORCE_CXY, if no candidates then sweep MU_WELL.

Saves summaries to results/abs_scan_cxy.txt and results/abs_scan_mu.txt and
writes a short report results/abs_scan_report.txt.
"""
import os
import json
import numpy as np
from numpy.linalg import eigh
from xxz_R_and_H import R_xxz, dRdu, expand_on_pauli
from verify_mzm import map_c_to_params


def build_bdg_from_arrays(N, t_hop, Delta_bond, mu_site):
    A = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)
    for i in range(N):
        A[i, i] = -mu_site[i]
        if i < N - 1:
            t = t_hop[i]
            A[i, i + 1] = -t
            A[i + 1, i] = -t
            Delta = Delta_bond[i]
            B[i, i + 1] = Delta
            B[i + 1, i] = -Delta
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    Hbdg = np.concatenate((top, bottom), axis=0)
    return Hbdg


def pauli_map_from_h(h_mapping):
    c_xx = h_mapping.get('X_X', 0.0)
    c_yy = h_mapping.get('Y_Y', 0.0)
    c_xy = h_mapping.get('X_Y', 0.0)
    c_yx = h_mapping.get('Y_X', 0.0)
    c_zz = h_mapping.get('Z_Z', 0.0)
    c_z0 = h_mapping.get('Z_I', 0.0)
    c_0z = h_mapping.get('I_Z', 0.0)
    return map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)


def scan_cxy(eta=0.6, u0=0.0, N=200, cxy_list=None, threshold=1e-3):
    os.makedirs('results', exist_ok=True)
    if cxy_list is None:
        cxy_list = np.linspace(0.01, 0.5, 10)
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
    mapping_h = expand_on_pauli(h_local)

    results = []
    candidates = []
    for cxy in cxy_list:
        mapping_mod = dict(mapping_h)
        mapping_mod['X_Y'] = mapping_mod.get('X_Y', 0.0) + cxy
        t_mod, Delta_mod, mu_mod, U_mod = pauli_map_from_h(mapping_mod)
        t_hop = np.full(N - 1, t_mod, dtype=complex)
        Delta_bond = np.full(N - 1, Delta_mod, dtype=complex)
        mu_site = np.full(N, float(mu_mod), dtype=float)
        H = build_bdg_from_arrays(N, t_hop, Delta_bond, mu_site)
        eigs = np.linalg.eigvals(H)
        min_abs = np.min(np.abs(np.real_if_close(eigs)))
        results.append((float(cxy), float(min_abs)))
        if min_abs < threshold:
            candidates.append({'cxy': float(cxy), 'min_abs': float(min_abs)})

    out_txt = 'results/abs_scan_cxy.txt'
    with open(out_txt, 'w') as f:
        f.write('cxy, min_abs_E\n')
        for cxy, m in results:
            f.write(f'{cxy:.6g}, {m:.6e}\n')
        f.write('\nCandidates (|E|<%g):\n' % threshold)
        for c in candidates:
            f.write(json.dumps(c) + '\n')

    return candidates


def scan_mu_well(eta=0.6, u0=0.0, N=200, mu_list=None, threshold=1e-3):
    os.makedirs('results', exist_ok=True)
    if mu_list is None:
        mu_list = np.linspace(0.0, 6.0, 13)
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
    mapping_h = expand_on_pauli(h_local)
    t_val, Delta_val, mu_val, U = pauli_map_from_h(mapping_h)

    results = []
    candidates = []
    for mu_w in mu_list:
        for sign in (+1, -1):
            mu_site = np.full(N, float(mu_val), dtype=float)
            mid = N // 2
            mu_site[mid] = float(mu_val) + sign * float(mu_w)
            t_hop = np.full(N - 1, t_val, dtype=complex)
            Delta_bond = np.full(N - 1, Delta_val, dtype=complex)
            H = build_bdg_from_arrays(N, t_hop, Delta_bond, mu_site)
            eigs = np.linalg.eigvals(H)
            min_abs = np.min(np.abs(np.real_if_close(eigs)))
            results.append((float(sign*mu_w), float(min_abs)))
            if min_abs < threshold:
                candidates.append({'mu_well': float(sign*mu_w), 'min_abs': float(min_abs)})

    out_txt = 'results/abs_scan_mu.txt'
    with open(out_txt, 'w') as f:
        f.write('mu_well_sign, min_abs_E\n')
        for mw, m in results:
            f.write(f'{mw:.6g}, {m:.6e}\n')
        f.write('\nCandidates (|E|<%g):\n' % threshold)
        for c in candidates:
            f.write(json.dumps(c) + '\n')

    return candidates


def main():
    # first sweep cxy
    print('Scanning c_XY...')
    cxy_vals = np.linspace(0.01, 0.5, 10)
    candidates = scan_cxy(cxy_list=cxy_vals)
    report = {'cxy_candidates': candidates}
    if not candidates:
        print('No c_XY candidates found; switching to mu_well scan')
        mu_vals = np.linspace(0.0, 6.0, 13)
        mu_cands = scan_mu_well(mu_list=mu_vals)
        report['mu_candidates'] = mu_cands
    else:
        print('Found c_XY candidates; skipping mu_well scan')

    with open('results/abs_scan_report.txt', 'w') as f:
        json.dump(report, f, indent=2)
    print('Scan complete. Reports in results/')


if __name__ == '__main__':
    main()
