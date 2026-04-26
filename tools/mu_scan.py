#!/usr/bin/env python3
"""Scan mu_right and record lowest abs energy vs mu.
Usage: python3 tools/mu_scan.py
"""
import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np
from numpy.linalg import eigh
from tools.pulse_abs import build_chain_from_params

def run(mu_vals=None, N=200):
    if mu_vals is None:
        mu_vals = np.linspace(2.5, 3.5, 21)
    os.makedirs('results', exist_ok=True)
    out = []
    for mu in mu_vals:
        H, _, _, _ = build_chain_from_params(N, 1.0, 0.5, 0.0, 1.0, 0.0, mu, interface_width=4)
        vals, _ = eigh(H)
        low = float(np.min(np.abs(np.real_if_close(vals))))
        out.append((mu, low))
    np.savetxt('results/mu_scan.txt', np.array(out))
    print('Wrote results/mu_scan.txt')

if __name__ == '__main__':
    run()
