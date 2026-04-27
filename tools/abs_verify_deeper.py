#!/usr/bin/env python3
"""Deeper verification for ABS-like candidates:
- robustness to local endpoint perturbation
- partner-mode localization check

Reads `results/scan_all_mixtures_validated.json` and writes
`results/abs_like_deeper_report.json`.
"""
import os
import json
import numpy as np
from numpy.linalg import eigh
import importlib.util

# load build_bdg_from_params from sibling module
here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location('scan_all_mixtures', os.path.join(here, 'scan_all_mixtures.py'))
scan_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan_mod)
build_bdg_from_params = scan_mod.build_bdg_from_params


def partner_check(H, topk=2):
    vals, vecs = eigh(H)
    idx = np.argsort(np.abs(vals))[:topk]
    modes = []
    half = H.shape[0] // 2
    for i in idx:
        v = vecs[:, i]
        u = v[:half]
        w = v[half:]
        ldos = np.abs(u)**2 + np.abs(w)**2
        modes.append({'E': float(vals[i]), 'max_site': int(np.argmax(ldos)), 'max_weight': float(np.max(ldos))})
    return modes


def robustness_test(t, Delta, mu, N=120, ntrials=100, sigma=0.01):
    H0 = build_bdg_from_params(N, t, Delta, mu)
    vals0 = np.sort(np.abs(np.linalg.eigvals(H0)))
    base = float(vals0[0])
    mins = []
    for _ in range(ntrials):
        mu_pert = np.array([0.0]*N, dtype=float)
        # local endpoint perturbation at site 0 and random small bulk noise
        mu_pert[0] += np.random.normal(scale=sigma)
        mu_pert += np.random.normal(scale=sigma*0.1, size=N)
        # build H with site-dependent mu: replace diagonal -mu by -mu - mu_pert
        A = np.zeros((N,N), dtype=complex)
        B = np.zeros((N,N), dtype=complex)
        for i in range(N):
            A[i,i] = - (mu + mu_pert[i])
            if i < N-1:
                A[i,i+1] = -t
                A[i+1,i] = -t
                B[i,i+1] = Delta
                B[i+1,i] = -Delta
        top = np.concatenate((A,B), axis=1)
        bottom = np.concatenate((-B.conj(), -A.T), axis=1)
        H = np.concatenate((top,bottom), axis=0)
        vals = np.abs(np.linalg.eigvals(H))
        mins.append(float(np.min(vals)))
    mins = np.array(mins)
    return {'base_min': base, 'min_mean': float(np.mean(mins)), 'min_std': float(np.std(mins)), 'min_median': float(np.median(mins))}


def main():
    os.makedirs('results', exist_ok=True)
    data = json.load(open('results/scan_all_mixtures_validated.json'))
    N = 120
    report = []
    # select ABS-like heuristically
    abs_like = [d for d in data if d.get('maj_sim',0) < 0.5 and (d.get('max_site',-1) <=3 or d.get('max_site',-1) >= N-4)]
    for i, d in enumerate(abs_like):
        print(f'Analyzing candidate {i+1}/{len(abs_like)}: {d}')
        t = complex(d['t'])
        Delta = complex(d['Delta'])
        mu = complex(d['mu'])
        modes = partner_check(build_bdg_from_params(N, t, Delta, mu), topk=2)
        rob = robustness_test(t, Delta, mu, N=N, ntrials=100, sigma=0.01)
        report.append({'candidate': d, 'modes': modes, 'robustness': rob})
    with open('results/abs_like_deeper_report.json', 'w') as f:
        json.dump(report, f, indent=2)
    print('Deeper verification complete. Reported:', len(report))


if __name__ == '__main__':
    main()
