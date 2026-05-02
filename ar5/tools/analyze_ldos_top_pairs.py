#!/usr/bin/env python3
"""Analyze top-N paired LDOS to determine whether combined LDOS shows two endpoint peaks.

Outputs results/ldos_top_pairs_summary.json and prints a concise table.
"""
import os
import json
import numpy as np
import importlib.util

here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location('scan_all_mixtures', os.path.join(here, 'scan_all_mixtures.py'))
scan_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan_mod)
build_bdg_from_params = scan_mod.build_bdg_from_params


def ldos_from_vec(vec):
    half = vec.shape[0]//2
    u = vec[:half]
    v = vec[half:]
    return np.abs(u)**2 + np.abs(v)**2


def compute_pair_ldos(cand, N=120):
    t = complex(cand['t'])
    Delta = complex(cand['Delta'])
    mu = complex(cand['mu'])
    H = build_bdg_from_params(N, t, Delta, mu)
    vals, vecs = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))
    v1 = vecs[:, idx[0]]
    v2 = vecs[:, idx[1]]
    l1 = ldos_from_vec(v1)
    l2 = ldos_from_vec(v2)
    s = l1 + l2
    return s, vals[idx[0]], vals[idx[1]]


def is_two_endpoint_peaks(s, left_frac=0.2, right_frac=0.2, min_ratio=0.1):
    N = len(s)
    left_region = int(N*left_frac)
    right_region = N - int(N*right_frac)
    # find two largest peaks
    peaks = np.argsort(s)[-2:][::-1]
    p1, p2 = int(peaks[0]), int(peaks[1])
    w1, w2 = float(s[p1]), float(s[p2])
    # consider endpoint if within left_region or >= right_region
    left = (p1 < left_region) or (p2 < left_region)
    right = (p1 >= right_region) or (p2 >= right_region)
    both_ends = left and right
    # also if two peaks both on same end -> not two-endpoints
    return {
        'peak_indices': [p1, p2],
        'peak_weights': [w1, w2],
        'both_ends': both_ends,
        'left_region': left_region,
        'right_region': right_region,
    }


def main(topk=6, N=120):
    os.makedirs('results', exist_ok=True)
    data = json.load(open('results/scan_all_mixtures_validated.json'))
    data_sorted = sorted(data, key=lambda d: d.get('min_abs', 1e9))
    chosen = data_sorted[:topk]
    summary = []
    for i, cand in enumerate(chosen, start=1):
        s, e1, e2 = compute_pair_ldos(cand, N=N)
        info = is_two_endpoint_peaks(s)
        rec = {
            'rank': i,
            'combo': cand.get('combo'),
            'E1': float(e1), 'E2': float(e2),
            'peak_indices': info['peak_indices'],
            'peak_weights': info['peak_weights'],
            'both_ends': info['both_ends']
        }
        summary.append(rec)
    open('results/ldos_top_pairs_summary.json','w').write(json.dumps(summary, indent=2))
    # print concise table
    print('rank  both_ends  peaks(indices)  peak_weights  combo')
    for r in summary:
        print(f"{r['rank']:2d}    {str(r['both_ends']):5s}     {r['peak_indices']}    {np.round(r['peak_weights'],4).tolist()}    {r['combo']}")


if __name__ == '__main__':
    main()
