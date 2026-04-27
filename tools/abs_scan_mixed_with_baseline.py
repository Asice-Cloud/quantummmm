#!/usr/bin/env python3
"""Scan mixed c_xy/c_yx with small baseline c_xx=c_yy=epsilon.

Scans c_xy,c_yx in [0.08,0.22] with c_xx=c_yy=epsilon and records candidates.
"""
import os
import json
import numpy as np
from scipy.linalg import eigh

from xxz_R_and_H import R_xxz, dRdu, expand_on_pauli
from verify_mzm import map_c_to_params, build_kitaev_bdg


def main(epsilon=0.05, cmin=0.08, cmax=0.22, pts=31, N=200, threshold=1e-6):
    os.makedirs('results', exist_ok=True)
    cxy_vals = np.linspace(cmin, cmax, pts)
    cyx_vals = np.linspace(cmin, cmax, pts)

    out_rows = []
    candidates = []

    for cxy in cxy_vals:
        for cyx in cyx_vals:
            c_xx = float(epsilon)
            c_yy = float(epsilon)
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

    csv_path = os.path.join('results', 'abs_scan_mixed_with_baseline.csv')
    with open(csv_path, 'w') as f:
        f.write('cxy,cyx,min_abs,site,max_weight,maj_sim\n')
        for r in out_rows:
            f.write(','.join(str(x) for x in r) + '\n')

    json_path = os.path.join('results', 'abs_scan_mixed_with_baseline_candidates.json')
    with open(json_path, 'w') as f:
        json.dump(candidates, f, indent=2)

    print('Mixed+baseline scan complete. CSV:', csv_path)
    print('Candidates:', json_path, 'count=', len(candidates))
    if candidates:
        print('Top 5 by min_abs:')
        for c in sorted(candidates, key=lambda x: x['min_abs'])[:5]:
            print(c)


if __name__ == '__main__':
    main()
