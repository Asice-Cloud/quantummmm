#!/usr/bin/env python3
"""Refined pure-mixed scan (exclude (0,0)).

Scan c_xy, c_yx in [0.08,0.22] with a fine grid, skip (0,0). Save CSV and JSON.
"""
import os
import json
import numpy as np
from scipy.linalg import eigh

from xxz_R_and_H import R_xxz, dRdu, expand_on_pauli
from verify_mzm import map_c_to_params, build_kitaev_bdg


def main():
    os.makedirs('results', exist_ok=True)
    cmin, cmax, pts = 0.08, 0.22, 31
    cxy_vals = np.linspace(cmin, cmax, pts)
    cyx_vals = np.linspace(cmin, cmax, pts)
    N = 200
    threshold = 1e-6

    out_rows = []
    candidates = []

    for cxy in cxy_vals:
        for cyx in cyx_vals:
            # skip trivial zero point (not in this range but keep check)
            if abs(cxy) < 1e-12 and abs(cyx) < 1e-12:
                continue

            c_xx = 0.0
            c_yy = 0.0
            c_xy = complex(cxy)
            c_yx = complex(cyx)
            c_zz = 0.0
            c_z0 = 0.0
            c_0z = 0.0

            t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
            H = build_kitaev_bdg(N, t, Delta, mu)
            w, v = eigh(H)
            idx = int(np.argmin(np.abs(w)))
            min_abs = float(np.abs(w[idx]))
            vec = v[:, idx]
            u = vec[:N]
            vpart = vec[N:]
            weights = np.abs(u)**2 + np.abs(vpart)**2
            max_site = int(np.argmax(weights))
            max_weight = float(max(weights))
            num = np.linalg.norm(u - np.conjugate(vpart))
            denom = np.linalg.norm(u) + np.linalg.norm(vpart) + 1e-12
            maj_norm = float(num / denom)
            maj_sim = float(1.0 - maj_norm)

            out_rows.append((float(cxy), float(cyx), min_abs, max_site, max_weight, maj_sim))
            if min_abs < threshold:
                candidates.append({'cxy': float(cxy), 'cyx': float(cyx), 'min_abs': min_abs, 'site': max_site, 'weight': max_weight, 'maj_sim': maj_sim})

    csv_path = os.path.join('results', 'abs_scan_pure_mixed_refined.csv')
    with open(csv_path, 'w') as f:
        f.write('cxy,cyx,min_abs,site,max_weight,maj_sim\n')
        for r in out_rows:
            f.write(','.join(str(x) for x in r) + '\n')

    json_path = os.path.join('results', 'abs_scan_pure_mixed_refined_candidates.json')
    with open(json_path, 'w') as f:
        json.dump(candidates, f, indent=2)

    print('Refined pure-mixed scan complete. CSV:', csv_path)
    print('Candidates:', json_path, 'count=', len(candidates))


if __name__ == '__main__':
    main()
