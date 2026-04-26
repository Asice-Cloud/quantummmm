#!/usr/bin/env python3
"""Refined multi-dimensional ABS scan around promising c_XY,c_YX regions.

Centers default to (0.15,0.15) and (0.20,0.20). Produces CSV and JSON
with diagnostics (min |E|, localization site and weight, Majorana score).
"""
import os
import json
import numpy as np
from numpy.linalg import eigvals

from xxz_R_and_H import R_xxz, dRdu, expand_on_pauli
from verify_mzm import map_c_to_params, build_kitaev_bdg


def build_bdg_from_params(N, t, Delta, mu):
    return build_kitaev_bdg(N, t, Delta, mu)


def main():
    os.makedirs('results', exist_ok=True)
    # refined grids around strong candidates
    centers = [(0.15, 0.15), (0.20, 0.20)]
    span = 0.06
    pts = 13
    cxz_vals = np.linspace(0.0, 0.05, 5)
    N = 200
    threshold = 1e-6

    # prepare base mapping from R
    u0 = 0.0
    eta = 0.6
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

    out_rows = []
    candidates = []

    for center in centers:
        cxy_c, cyx_c = center
        cxy_vals = np.linspace(max(0.0, cxy_c - span), cxy_c + span, pts)
        cyx_vals = np.linspace(max(0.0, cyx_c - span), cyx_c + span, pts)
        for cxy in cxy_vals:
            for cyx in cyx_vals:
                for cxz in cxz_vals:
                    mapping_mod = dict(base_mapping)
                    mapping_mod['X_Y'] = mapping_mod.get('X_Y', 0.0) + complex(cxy)
                    mapping_mod['Y_X'] = mapping_mod.get('Y_X', 0.0) + complex(cyx)
                    mapping_mod['X_Z'] = mapping_mod.get('X_Z', 0.0) + complex(cxz)
                    # map pauli mapping to params
                    c_xx = mapping_mod.get('X_X', 0.0)
                    c_yy = mapping_mod.get('Y_Y', 0.0)
                    c_xy = mapping_mod.get('X_Y', 0.0)
                    c_yx = mapping_mod.get('Y_X', 0.0)
                    c_zz = mapping_mod.get('Z_Z', 0.0)
                    c_z0 = mapping_mod.get('Z_I', 0.0)
                    c_0z = mapping_mod.get('I_Z', 0.0)
                    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)

                    H = build_bdg_from_params(N, t, Delta, mu)
                    # eigen-decompose (full), get smallest magnitude eigenvalue and vector
                    from scipy.linalg import eigh
                    w, v = eigh(H)
                    idx = np.argmin(np.abs(w))
                    min_abs = float(np.abs(w[idx]))
                    vec = v[:, idx]
                    # split u,v and compute site weights
                    u = vec[:N]
                    vpart = vec[N:]
                    weights = np.abs(u)**2 + np.abs(vpart)**2
                    max_site = int(np.argmax(weights))
                    max_weight = float(weights[max_site])
                    # Majorana-like measure: how close u == conj(v)
                    num = np.linalg.norm(u - np.conjugate(vpart))
                    denom = np.linalg.norm(u) + np.linalg.norm(vpart) + 1e-12
                    maj_score = float(num / denom)

                    out_rows.append((cxy, cyx, cxz, min_abs, max_site, max_weight, maj_score))
                    if min_abs < threshold:
                        candidates.append({'cxy': float(cxy), 'cyx': float(cyx), 'cxz': float(cxz), 'min_abs': min_abs, 'site': max_site, 'weight': max_weight, 'maj_score': maj_score})

    csv_path = os.path.join('results', 'abs_scan_refined.csv')
    with open(csv_path, 'w') as f:
        f.write('cxy,cyx,cxz,min_abs,site,max_weight,maj_score\n')
        for r in out_rows:
            f.write(','.join(str(x) for x in r) + '\n')

    json_path = os.path.join('results', 'abs_scan_refined_candidates.json')
    with open(json_path, 'w') as f:
        json.dump(candidates, f, indent=2)

    print('Refined scan complete. CSV:', csv_path)
    print('Candidates:', json_path, 'count=', len(candidates))


if __name__ == '__main__':
    main()
