#!/usr/bin/env python3
"""Validate top-k candidates from refined multi-scan.

Reads results/abs_scan_refined_candidates.json, picks top-k by min_abs,
recomputes BdG for multiple N, and writes numeric diagnostics to
results/abs_validate_top.json (includes per-N min_abs, site, weight,
majorana measures and site_weights arrays saved as lists).
"""
import os
import json
import numpy as np

from scipy.linalg import eigh

from xxz_R_and_H import R_xxz, dRdu, expand_on_pauli
from verify_mzm import map_c_to_params, build_kitaev_bdg

ROOT = os.path.dirname(os.path.dirname(__file__))
RES = os.path.join(ROOT, 'results')
os.makedirs(RES, exist_ok=True)


def make_base_mapping():
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
    return expand_on_pauli(h_local)


def validate_candidate(c_entry, base_mapping, Ns=(100, 200, 400)):
    # c_entry expected to have cxy, cyx, cxz (float)
    cxy = float(c_entry.get('cxy', c_entry.get('XY', 0.0)))
    cyx = float(c_entry.get('cyx', c_entry.get('YX', 0.0)))
    cxz = float(c_entry.get('cxz', c_entry.get('XZ', 0.0)))

    mapping_mod = dict(base_mapping)
    mapping_mod['X_Y'] = mapping_mod.get('X_Y', 0.0) + complex(cxy)
    mapping_mod['Y_X'] = mapping_mod.get('Y_X', 0.0) + complex(cyx)
    mapping_mod['X_Z'] = mapping_mod.get('X_Z', 0.0) + complex(cxz)

    c_xx = mapping_mod.get('X_X', 0.0)
    c_yy = mapping_mod.get('Y_Y', 0.0)
    c_xy = mapping_mod.get('X_Y', 0.0)
    c_yx = mapping_mod.get('Y_X', 0.0)
    c_zz = mapping_mod.get('Z_Z', 0.0)
    c_z0 = mapping_mod.get('Z_I', 0.0)
    c_0z = mapping_mod.get('I_Z', 0.0)

    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)

    result = {'cxy': cxy, 'cyx': cyx, 'cxz': cxz, 'params': {'t': complex(t), 'Delta': complex(Delta), 'mu': complex(mu), 'U': complex(U)}, 'per_N': {}}

    for N in Ns:
        H = build_kitaev_bdg(N, t, Delta, mu)
        w, v = eigh(H)
        idx = int(np.argmin(np.abs(w)))
        min_abs = float(np.abs(w[idx]))
        vec = v[:, idx]
        u = vec[:N]
        vpart = vec[N:]
        weights = (np.abs(u)**2 + np.abs(vpart)**2).tolist()
        max_site = int(np.argmax(weights))
        max_weight = float(max(weights))
        num = np.linalg.norm(u - np.conjugate(vpart))
        denom = np.linalg.norm(u) + np.linalg.norm(vpart) + 1e-12
        maj_norm = float(num / denom)
        maj_sim = float(1.0 - maj_norm)

        result['per_N'][str(N)] = {'min_abs': min_abs, 'site': max_site, 'max_weight': max_weight, 'maj_norm': maj_norm, 'maj_sim': maj_sim, 'site_weights': weights}

    return result


def main(top_k=3, candidates_path=os.path.join(RES, 'abs_scan_refined_candidates.json')):
    with open(candidates_path, 'r') as f:
        cand = json.load(f)

    if not cand:
        print('No candidates found in', candidates_path)
        return

    # sort by min_abs ascending
    cand_sorted = sorted(cand, key=lambda x: x.get('min_abs', 1e9))
    top = cand_sorted[:top_k]

    base = make_base_mapping()
    results = []
    for i, c in enumerate(top):
        print(f'Validating top {i+1}:', c)
        res = validate_candidate(c, base)
        results.append(res)

    out_path = os.path.join(RES, 'abs_validate_top.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=lambda x: (x.real if isinstance(x, complex) else x))

    print('Validation complete. Results:', out_path)


if __name__ == '__main__':
    main()
