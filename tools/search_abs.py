#!/usr/bin/env python3
"""Search for ABS-like parameter sets: embed a short paired region in metallic background
and look for low-energy states that are NOT Majorana-like and that are sensitive to small perturbations.

Saves candidates to results/abs_search.json and results/abs_search.csv
"""
import os, sys, json
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np
from numpy.linalg import eigh, norm
from tools.run_midpair_demo import build_midpair_chain


def majorana_score_from_vec(vec):
    Nloc = vec.shape[0] // 2
    u = vec[:Nloc]
    v = vec[Nloc:]
    diff = u - np.conjugate(v)
    return float(norm(diff) / (norm(u) + 1e-30))


def is_sensitive(H, vec_idx, perturb=1e-3):
    # apply small random on-site potential at interface and recompute lowest energy shift
    N2 = H.shape[0] // 2
    N = N2
    # build small local potential at center site
    P = np.zeros((2*N,2*N), dtype=complex)
    center = N // 2
    # add potential to particle sector at center
    P[center, center] = perturb
    H2 = H + P
    vals2, _ = eigh(H2)
    # compare lowest abs eigenvalues
    return float(np.min(np.abs(vals2)))


def run_search(N=100, widths=(2,4,8), deltas=(0.2,0.5,1.0), threshold_majorana=1e-3, energy_tol=1e-2):
    os.makedirs('results', exist_ok=True)
    candidates = []
    for w in widths:
        for d in deltas:
            H, t_hop, Delta_bond, mu_site = build_midpair_chain(N, t_val=1.0, Delta_center=d, center_width=w, mu=0.0)
            vals, vecs = eigh(H)
            idx = np.argsort(np.abs(vals))
            low = float(np.abs(vals[idx[0]]))
            vec0 = vecs[:, idx[0]]
            mscore = majorana_score_from_vec(vec0)
            # sensitivity test: energy after small perturbation
            val2 = is_sensitive(H, idx[0], perturb=1e-2)
            shift = abs(val2 - low)
            # candidate ABS: low energy but NOT Majorana-like (mscore large) and sensitive (shift > energy_tol)
            if low < 1e-1 and mscore > threshold_majorana and shift > energy_tol:
                candidates.append({'N': N, 'width': w, 'Delta_center': d, 'low_energy': low, 'majorana_score': mscore, 'shift': shift})
    # save
    with open('results/abs_search.json','w') as f:
        json.dump(candidates, f, indent=2)
    # CSV
    if candidates:
        with open('results/abs_search.csv','w') as f:
            f.write('N,width,Delta_center,low_energy,majorana_score,shift\n')
            for c in candidates:
                f.write(f"{c['N']},{c['width']},{c['Delta_center']},{c['low_energy']},{c['majorana_score']},{c['shift']}\n")
    print('Search complete. Found', len(candidates), 'candidates. Saved to results/abs_search.json')


if __name__ == '__main__':
    run_search()
