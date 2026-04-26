#!/usr/bin/env python3
"""Scan single-site Pauli contributions for ABS signatures.

For each single-site term among X_I, I_X, Y_I, I_Y, Z_I, I_Z we add
small coefficients and map to BdG parameters using existing mapping
utilities. Results are saved to results/abs_scan_single_pauli.csv and
results/abs_scan_single_pauli_candidates.json
"""
import json
import os
import sys
import numpy as np

# Ensure repo root is on sys.path so `tools` modules can be imported when
# running this script from the repository root.
ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from tools.verify_mzm import map_c_to_params, build_kitaev_bdg

RESULTS = os.path.join(ROOT, 'results')
os.makedirs(RESULTS, exist_ok=True)

# Base two-site coefficient dictionary (start from a typical candidate from multi-scan)
base_c = {
    'XX': 0.2,
    'YY': 0.2,
    'XY': 0.0,
    'YX': 0.0,
    'ZZ': 0.0,
    '0Z': 0.0,
    'Z0': 0.0,
    '00': 0.0,
}

single_keys = ['X_I', 'I_X', 'Y_I', 'I_Y', 'Z_I', 'I_Z']
vals = np.linspace(0.01, 0.5, 10)

out_rows = []
candidates = []

N = 200

for key in single_keys:
    for v in vals:
        c = base_c.copy()
        # Only Z single-site terms map directly to the mu parameters used by
        # `map_c_to_params`. X/Y single-site Pauli terms are nonlocal under
        # Jordan-Wigner and are not handled by the simple map; mark them
        # as unsupported here.
        if key == 'Z_I':
            c['Z0'] = float(v)
        elif key == 'I_Z':
            c['0Z'] = float(v)
        else:
            out_rows.append((key, float(v), None, 'unsupported_single'))
            continue

        try:
            t, Delta, mu, U = map_c_to_params(c.get('XX',0.0), c.get('YY',0.0), c.get('XY',0.0), c.get('YX',0.0), c.get('ZZ',0.0), c.get('Z0',0.0), c.get('0Z',0.0))
        except Exception:
            out_rows.append((key, float(v), None, 'map_error'))
            continue

        H = build_kitaev_bdg(N, t, Delta, mu)
        evals = np.sort(np.abs(np.linalg.eigvals(H)))
        min_abs = float(evals[0])

        # compute localization of lowest eigenvector (from H diagonalizer returned)
        # build_kitaev_bdg returns (H, diag) where diag are eigenvalues from eigh of BdG
        # For speed we recompute eigenvectors for lowest pair
        from scipy.linalg import eigh
        w, vvecs = eigh(H)
        idx = np.argmin(np.abs(w))
        vec = vvecs[:, idx]
        # split into site weights (u,v per site => 2*N entries)
        weights = np.abs(vec)**2
        site_weights = weights.reshape((N, 2)).sum(axis=1)
        max_site = int(np.argmax(site_weights))
        max_weight = float(site_weights[max_site])

        out_rows.append((key, float(v), min_abs, max_site, max_weight))
        if min_abs < 1e-3:
            candidates.append({'key': key, 'val': float(v), 'min_abs': min_abs, 'site': max_site, 'weight': max_weight})

csv_path = os.path.join(RESULTS, 'abs_scan_single_pauli.csv')
with open(csv_path, 'w') as f:
    f.write('key,val,min_abs,site,max_weight\n')
    for r in out_rows:
        f.write(','.join(str(x) for x in r))
        f.write('\n')

json_path = os.path.join(RESULTS, 'abs_scan_single_pauli_candidates.json')
with open(json_path, 'w') as f:
    json.dump(candidates, f, indent=2)

print('Scan complete. CSV:', csv_path)
print('Candidates:', json_path, 'count=', len(candidates))
