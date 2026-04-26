#!/usr/bin/env python3
"""Scan R_xxz parameter eta and record mapped (t, Delta, mu, U) from local h.

Saves results/eta_delta_scan.json and results/eta_delta_scan.png
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

from tools.pulse_abs import h_from_xxz, pauli_map_to_params


def run_scan(etas, out_json='results/eta_delta_scan.json', out_png='results/eta_delta_scan.png'):
    os.makedirs('results', exist_ok=True)
    rows = []
    for eta in etas:
        try:
            mapping = h_from_xxz(float(eta))
            t, Delta, mu, U = pauli_map_to_params(mapping)
            rows.append({'eta': float(eta), 't': complex(t), 'Delta': complex(Delta), 'mu': float(mu), 'U': float(U)})
        except Exception as e:
            rows.append({'eta': float(eta), 'error': str(e)})

    with open(out_json, 'w') as f:
        json.dump(rows, f, default=str, indent=2)

    # plot |Delta| vs eta for successful rows
    et = []
    dabs = []
    for r in rows:
        if 'Delta' in r:
            et.append(r['eta'])
            dabs.append(abs(complex(r['Delta'])))
    plt.figure(figsize=(6,3))
    plt.plot(et, dabs, '-o')
    plt.xlabel('eta')
    plt.ylabel('|Delta|')
    plt.title('|Delta| from R->h mapping')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png)
    print('Wrote', out_json, 'and', out_png)


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--n', type=int, default=40)
    p.add_argument('--eta-min', type=float, default=0.05)
    p.add_argument('--eta-max', type=float, default=1.5)
    args = p.parse_args()
    etas = np.linspace(args.eta_min, args.eta_max, args.n)
    run_scan(etas)


if __name__ == '__main__':
    main()
