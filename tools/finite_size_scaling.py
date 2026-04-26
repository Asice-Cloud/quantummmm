#!/usr/bin/env python3
"""Compute lowest BdG energy vs system size N for a chosen eta.

Usage:
  python3 tools/finite_size_scaling.py --eta 0.6 --Ns 50 100 200 --force-delta 0.3

Saves results/finite_size_eta_<eta>.json and results/finite_size_eta_<eta>.png
"""
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# ensure project root on sys.path
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from tools.pulse_abs import h_from_xxz, pauli_map_to_params, build_chain_from_params


def run_scaling(eta, Ns, force_delta=None, out_json=None, out_png=None):
    os.makedirs('results', exist_ok=True)
    rows = []
    mapping = h_from_xxz(float(eta))
    t_val, Delta_val, mu_val, U = pauli_map_to_params(mapping)
    for N in Ns:
        Delta_left = float(force_delta) if force_delta is not None else Delta_val
        H, t_hop, Delta_bond, mu_site = build_chain_from_params(int(N), t_val, Delta_left, mu_val, t_val, 0.0, mu_val, interface_width=4)
        vals, vecs = np.linalg.eigh(H)
        idx = np.argsort(np.real_if_close(vals))
        lowE = float(np.real_if_close(vals[idx[0]]))
        rows.append({'N': int(N), 'lowE': lowE})

    if out_json is None:
        out_json = f'results/finite_size_eta_{eta:.3f}.json'
    if out_png is None:
        out_png = f'results/finite_size_eta_{eta:.3f}.png'
    with open(out_json, 'w') as f:
        json.dump({'eta': eta, 'force_delta': force_delta, 'data': rows}, f, indent=2)

    Ns_arr = np.array([r['N'] for r in rows])
    Es = np.array([r['lowE'] for r in rows])
    plt.figure(figsize=(6,3))
    plt.semilogy(Ns_arr, np.abs(Es), '-o')
    plt.xlabel('N')
    plt.ylabel('|lowest E|')
    plt.title(f'Finite-size scaling eta={eta}')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png)
    print('Wrote', out_json, 'and', out_png)


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--eta', type=float, required=True)
    p.add_argument('--Ns', type=int, nargs='+', default=[50,100,150,200])
    p.add_argument('--force-delta', type=float, default=None)
    args = p.parse_args()
    run_scaling(args.eta, args.Ns, force_delta=args.force_delta)


if __name__ == '__main__':
    main()
