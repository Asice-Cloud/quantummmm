#!/usr/bin/env python3
"""Size-scaling test: run BdG for several N and record lowest absolute energy.
Usage: python3 tools/test_size_scaling.py
"""
import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np
from numpy.linalg import eigh
from tools.pulse_abs import build_chain_from_params

def run(ns=(200,400,800), t_left=1.0, Delta_left=0.5, mu_left=0.0, t_right=1.0, Delta_right=0.0, mu_right=3.0, interface_width=4):
    os.makedirs('results', exist_ok=True)
    out = []
    for N in ns:
        H, t_hop, Delta_bond, mu_site = build_chain_from_params(N, t_left, Delta_left, mu_left, t_right, Delta_right, mu_right, interface_width=interface_width)
        vals, vecs = eigh(H)
        idx = np.argsort(np.abs(vals))
        lowvals = np.real_if_close(vals[idx[:6]])
        out.append((N, [float(v) for v in lowvals]))
        np.savetxt(f'results/scale_N{N}_vals.txt', np.sort(np.real_if_close(vals)))
        # save density for lowest mode
        vec0 = vecs[:, idx[0]]
        Nloc = vec0.shape[0]//2
        dens0 = np.abs(vec0[:Nloc])**2 + np.abs(vec0[Nloc:])**2
        np.savetxt(f'results/scale_N{N}_density0.txt', dens0)
    with open('results/scale_summary.txt','w') as f:
        for N, lows in out:
            f.write(f'N={N} lows=' + ' '.join(f'{v:.6e}' for v in lows) + '\n')
    print('Wrote results/scale_summary.txt and per-N files')

if __name__ == '__main__':
    run()
