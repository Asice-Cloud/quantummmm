#!/usr/bin/env python3
"""Perturb endpoint mu slightly and recompute Pfaffian to probe topology.

Reads `results/finite_chain_pfaffian_top10.json` and for each candidate
computes Pfaffian for mu +/- eps (endpoint shift) and reports signs/values.
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

from finite_chain_pfaffian import majorana_matrix_from_bdG, pfaffian


def perturb_mu_and_pf(t, Delta, mu, N=120, eps=1e-3):
    # apply small change to site 0 chemical potential
    def build_H(mu_arr):
        A = np.zeros((N, N), dtype=complex)
        B = np.zeros((N, N), dtype=complex)
        for i in range(N):
            A[i, i] = -mu_arr[i]
            if i < N - 1:
                A[i, i + 1] = -t
                A[i + 1, i] = -t
                B[i, i + 1] = Delta
                B[i + 1, i] = -Delta
        top = np.concatenate((A, B), axis=1)
        bottom = np.concatenate((-B.conj(), -A.T), axis=1)
        return np.concatenate((top, bottom), axis=0)

    mu0 = np.real(mu)
    mu_plus = np.array([mu0]*N)
    mu_minus = np.array([mu0]*N)
    mu_plus[0] += eps
    mu_minus[0] -= eps
    H_plus = build_H(mu_plus)
    H_minus = build_H(mu_minus)
    M_plus = majorana_matrix_from_bdG(H_plus)
    M_minus = majorana_matrix_from_bdG(H_minus)
    pf_plus = pfaffian(M_plus)
    pf_minus = pfaffian(M_minus)
    return {'pf_plus': complex(pf_plus), 'pf_minus': complex(pf_minus), 'sign_plus': np.sign(np.real(pf_plus)) if abs(pf_plus)!=0 else 0, 'sign_minus': np.sign(np.real(pf_minus)) if abs(pf_minus)!=0 else 0}


def main():
    os.makedirs('results', exist_ok=True)
    top = json.load(open('results/finite_chain_pfaffian_top10.json'))
    out = []
    for item in top:
        d = item['candidate']
        t = complex(d['t'])
        Delta = complex(d['Delta'])
        mu = complex(d['mu'])
        res = perturb_mu_and_pf(t, Delta, mu, N=120, eps=1e-3)
        out.append({'candidate': d, **res})
    with open('results/finite_chain_pfaffian_top10_perturbed.json','w') as f:
        json.dump(out, f, indent=2, default=str)
    print('Perturbation Pfaffian test complete for top 10')


if __name__ == '__main__':
    main()
