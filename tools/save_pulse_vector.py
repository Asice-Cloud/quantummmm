#!/usr/bin/env python3
"""Save lowest BdG eigenvector for a pulse-derived S-N chain.

Creates `results/vec_lowest.npy` (by default) and optionally runs
`tools/pauli_to_majorana.py` to produce a Pauli expansion JSON.

Usage:
  python3 tools/save_pulse_vector.py --eta 0.6 --force-delta 0.3 --N 200 --run-pauli
"""
import argparse
import os
import subprocess
import numpy as np
import sys

# Ensure project root is on sys.path when script is run as
# `python3 tools/save_pulse_vector.py` so `import tools.*` works.
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from tools.pulse_abs import h_from_xxz, pauli_map_to_params, build_chain_from_params


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--eta', type=float, default=0.6)
    p.add_argument('--force-delta', type=float, default=None,
                   help='If given, override left Delta to this value')
    p.add_argument('--N', type=int, default=200)
    p.add_argument('--interface-width', type=int, default=4)
    p.add_argument('--out', default='results/vec_lowest.npy')
    p.add_argument('--run-pauli', action='store_true', help='Run pauli_to_majorana on saved vector')
    p.add_argument('--top', type=int, default=20, help='Top Pauli strings to print if running pauli script')
    args = p.parse_args()

    os.makedirs('results', exist_ok=True)

    mapping = h_from_xxz(args.eta)
    t_val, Delta_val, mu_val, U = pauli_map_to_params(mapping)

    t_right = t_val
    Delta_right = 0.0
    mu_right = mu_val
    t_left = t_val
    Delta_left = args.force_delta if args.force_delta is not None else Delta_val
    mu_left = mu_val

    N = args.N
    H, t_hop, Delta_bond, mu_site = build_chain_from_params(N, t_left, Delta_left, mu_left,
                                                          t_right, Delta_right, mu_right,
                                                          interface_width=args.interface_width)

    vals, vecs = np.linalg.eigh(H)
    idx = np.argsort(np.real_if_close(vals))
    vec = vecs[:, idx[0]]
    np.save(args.out, vec)
    print(f'saved {args.out} (N={N}) ; lowest energy = {np.real_if_close(vals[idx[0]]):.6e}')

    if args.run_pauli:
        cmd = ["python3", "tools/pauli_to_majorana.py", "--vec", args.out, "--top", str(args.top)]
        out_path = 'results/pauli_gamma_pulse.txt'
        print('Running pauli_to_majorana and saving output to', out_path)
        with open(out_path, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        print('Pauli expansion written to', out_path)


if __name__ == '__main__':
    main()
