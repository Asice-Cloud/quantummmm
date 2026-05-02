#!/usr/bin/env python3
"""Validate candidates by computing eigenpairs and simple diagnostics.

Outputs a JSON with per-candidate metrics: min_abs_E, maj_sim, IPR, max_site, localization.
"""
import os
import json
import numpy as np
from numpy.linalg import eigh
import importlib.util

# load build_bdg_from_params from sibling module scan_all_mixtures.py
here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location('scan_all_mixtures', os.path.join(here, 'scan_all_mixtures.py'))
scan_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan_mod)
build_bdg_from_params = scan_mod.build_bdg_from_params


def maj_similarity(u, v):
    num = np.linalg.norm(u - np.conjugate(v))
    denom = np.linalg.norm(u) + np.linalg.norm(v)
    if denom == 0:
        return 0.0
    return 1.0 - num / denom


def site_ldos(u, v):
    # u/v are length-N (particle/hole)
    return np.abs(u)**2 + np.abs(v)**2


def analyze_candidate(cand, N=120):
    t = complex(cand['t'])
    Delta = complex(cand['Delta'])
    mu = complex(cand['mu'])
    H = build_bdg_from_params(N, t, Delta, mu)
    vals, vecs = eigh(H)
    idx = np.argsort(np.abs(vals))
    ev = vals[idx[0]]
    eigvec = vecs[:, idx[0]]
    half = H.shape[0] // 2
    u = eigvec[:half]
    v = eigvec[half:]
    ldos = site_ldos(u, v)
    max_site = int(np.argmax(ldos))
    max_weight = float(np.max(ldos))
    ipr = float(np.sum(ldos**2))
    maj_sim = float(maj_similarity(u, v))
    return {
        'min_abs_E': float(np.abs(ev)),
        'maj_sim': maj_sim,
        'ipr': ipr,
        'max_site': max_site,
        'max_weight': max_weight,
    }


def main():
    os.makedirs('results', exist_ok=True)
    with open('results/scan_all_mixtures_candidates.json') as f:
        cands = json.load(f)
    out = []
    for i, cand in enumerate(cands):
        try:
            metrics = analyze_candidate(cand)
            cand.update(metrics)
            out.append(cand)
        except Exception as e:
            print('Error on candidate', i, e)

    with open('results/scan_all_mixtures_validated.json', 'w') as f:
        json.dump(out, f, indent=2)
    print('Validation complete. Validated:', len(out))


if __name__ == '__main__':
    main()
