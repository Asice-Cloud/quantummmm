#!/usr/bin/env python3
"""Compute maj_sim, Majorana component densities and overlaps for top-N candidates.

Produces results/majorana_overlap_topN.json and a CSV summary.
"""
import os
import json
import numpy as np
import importlib.util
from math import isfinite

here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location('scan_all_mixtures', os.path.join(here, 'scan_all_mixtures.py'))
scan_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan_mod)
build_bdg_from_params = scan_mod.build_bdg_from_params


def maj_sim_from_vec(vec):
    half = vec.shape[0]//2
    u = vec[:half]
    v = vec[half:]
    num = np.linalg.norm(u - np.conj(v))
    den = np.linalg.norm(u) + np.linalg.norm(v)
    return float(1.0 - (num / den)) if den>0 else 0.0


def psi_A_B_from_vec(vec):
    half = vec.shape[0]//2
    u = vec[:half]
    v = vec[half:]
    psiA = u + np.conj(v)
    psiB = u - np.conj(v)
    return psiA, psiB


def overlap_norm(a, b):
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na==0 or nb==0:
        return 0.0
    return float(np.sum(np.abs(a)*np.abs(b)) / (na*nb))


def centroid(rho):
    tot = np.sum(rho)
    if tot<=0:
        return None
    x = np.arange(len(rho))
    return float(np.sum(x * rho) / tot)


def analyze_topN(topk=6, N=120):
    os.makedirs('results', exist_ok=True)
    data = json.load(open('results/scan_all_mixtures_validated.json'))
    data_sorted = sorted(data, key=lambda d: d.get('min_abs', 1e9))
    chosen = data_sorted[:topk]
    out = []
    csv_lines = ["rank,combo,E1,E2,maj_sim1,maj_sim2,peak1,peak2,centroidA1,centroidB1,centroidA2,centroidB2,S_A1B1,S_A2B2,S_A1_A2,S_B1_B2"]
    for i, cand in enumerate(chosen, start=1):
        t = complex(cand['t'])
        Delta = complex(cand['Delta'])
        mu = complex(cand['mu'])
        H = build_bdg_from_params(N, t, Delta, mu)
        vals, vecs = np.linalg.eigh(H)
        idx = np.argsort(np.abs(vals))
        v1 = vecs[:, idx[0]]
        v2 = vecs[:, idx[1]]
        E1 = float(vals[idx[0]]);
        E2 = float(vals[idx[1]]);
        maj1 = maj_sim_from_vec(v1)
        maj2 = maj_sim_from_vec(v2)
        # A/B components and densities
        A1, B1 = psi_A_B_from_vec(v1)
        A2, B2 = psi_A_B_from_vec(v2)
        rhoA1 = np.abs(A1)**2; rhoB1 = np.abs(B1)**2
        rhoA2 = np.abs(A2)**2; rhoB2 = np.abs(B2)**2
        centroidA1 = centroid(rhoA1); centroidB1 = centroid(rhoB1)
        centroidA2 = centroid(rhoA2); centroidB2 = centroid(rhoB2)
        # overlaps
        S_A1B1 = overlap_norm(A1, B1)
        S_A2B2 = overlap_norm(A2, B2)
        S_A1_A2 = overlap_norm(A1, A2)
        S_B1_B2 = overlap_norm(B1, B2)
        peak1 = int(np.argmax(np.abs(v1[:N])**2 + np.abs(v1[N:])**2))
        peak2 = int(np.argmax(np.abs(v2[:N])**2 + np.abs(v2[N:])**2))
        rec = {
            'rank': i,
            'combo': cand.get('combo'),
            'E1': E1, 'E2': E2,
            'maj_sim1': maj1, 'maj_sim2': maj2,
            'peak1': peak1, 'peak2': peak2,
            'centroidA1': centroidA1, 'centroidB1': centroidB1,
            'centroidA2': centroidA2, 'centroidB2': centroidB2,
            'S_A1B1': S_A1B1, 'S_A2B2': S_A2B2, 'S_A1_A2': S_A1_A2, 'S_B1_B2': S_B1_B2
        }
        out.append(rec)
        csv_lines.append(','.join([str(rec.get(k,'')) for k in ['rank','combo','E1','E2','maj_sim1','maj_sim2','peak1','peak2','centroidA1','centroidB1','centroidA2','centroidB2','S_A1B1','S_A2B2','S_A1_A2','S_B1_B2']]))
    json.dump(out, open('results/majorana_overlap_topN.json','w'), indent=2)
    open('results/majorana_overlap_topN.csv','w').write('\n'.join(csv_lines))
    print('Wrote results/majorana_overlap_topN.json and .csv')


if __name__ == '__main__':
    analyze_topN()
