#!/usr/bin/env python3
"""Majorana-ness check: compute ||u - conj(v)||/||u|| for lowest mode.
Usage: python3 tools/majorana_check.py
"""
import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np
from numpy.linalg import eigh, norm
from tools.pulse_abs import build_chain_from_params

def run(N=200, t_left=1.0, Delta_left=0.5, mu_left=0.0, t_right=1.0, Delta_right=0.0, mu_right=3.0, interface_width=4):
    os.makedirs('results', exist_ok=True)
    H, t_hop, Delta_bond, mu_site = build_chain_from_params(N, t_left, Delta_left, mu_left, t_right, Delta_right, mu_right, interface_width=interface_width)
    vals, vecs = eigh(H)
    idx = np.argsort(np.abs(vals))
    vec0 = vecs[:, idx[0]]
    Nloc = vec0.shape[0]//2
    u = vec0[:Nloc]
    v = vec0[Nloc:]
    diff = u - np.conjugate(v)
    score = norm(diff) / (norm(u) + 1e-30)
    # also compute overlap u·conj(v)
    overlap = abs(np.vdot(u, np.conjugate(v)))
    with open('results/majorana_check.txt','w') as f:
        f.write(f'N={N} lowest_energy={float(vals[idx[0]]):.6e}\n')
        f.write(f'majorana_score={score:.6e}\n')
        f.write(f'abs_overlap_u_conjv={overlap:.6e}\n')
    print('Wrote results/majorana_check.txt')
    print(f'majorana_score={score:.6e} overlap={overlap:.6e}')

if __name__ == '__main__':
    run()
