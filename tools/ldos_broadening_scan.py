#!/usr/bin/env python3
"""Scan LDOS zero-energy broadening (eta) to mimic finite-T / experimental broadening.

Generates LDOS snapshots for a few eta values and saves PNGs to `results/`.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# ensure repo root importable
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import embed_kitaev
import tools.paper_params as P


def run(L=None, t0=None, Delta0=None, mu0=None, VD=None, qd_width=None, etas=None):
    L = P.L if L is None else L
    t0 = P.t0 if t0 is None else t0
    Delta0 = P.DELTA if Delta0 is None else Delta0
    mu0 = P.mu0 if mu0 is None else mu0
    VD = P.VD if VD is None else VD
    qd_width = P.QD_WIDTH if qd_width is None else qd_width
    etas = etas if etas is not None else [1e-2, 1e-3, 1e-4]

    snapshots = {
        'init': (1.0, 1.0, 0.0, 0.0),
        'after_step1': (0.0, 1.0, 1.0, 0.0),
        'after_step2': (1.0, 0.0, 1.0, 0.0),
        'after_step3': (1.0, 1.0, 0.0, 0.0),
    }

    os.makedirs('results', exist_ok=True)
    for eta in etas:
        for name, (g1, g2, g3, g4) in snapshots.items():
            # construct base arrays
            mu = mu0 * np.ones(L)
            t_links = t0 * np.ones(L - 1)
            Delta_links = Delta0 * np.ones(L - 1)

            # QD at left end
            for i in range(qd_width):
                mu[i] = mu[i] - VD

            # modify first two links according to g1 and g3
            t_links_mod = t_links.copy()
            t_links_mod[0] = t0 * g1
            if L > 2:
                t_links_mod[1] = t0 * g3
            Delta_mod = Delta_links.copy()
            Delta_mod[0] = Delta0 * (g1 if g1 > 0 else 0.01)
            if L > 2:
                Delta_mod[1] = Delta0 * (g3 if g3 > 0 else 0.01)

            H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
            ldos, E = embed_kitaev.compute_zero_ldos(H, eta=eta)

            plt.figure(figsize=(8,3))
            plt.plot(np.arange(L), ldos, '-o')
            plt.xlabel('site')
            plt.ylabel(f'LDOS(E=0), eta={eta:.0e}')
            plt.title(f'LDOS snapshot: {name} (g1,g2,g3,g4)={g1,g2,g3,g4}')
            out = f'results/ldos_snapshot_{name}_eta{int(-np.log10(eta))}.png'
            plt.tight_layout()
            plt.savefig(out)
            plt.close()
            print('Saved', out)


if __name__ == '__main__':
    run()
