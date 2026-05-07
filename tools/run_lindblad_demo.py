#!/usr/bin/env python3
"""Run a small Lindblad demo using the mapped two-level Hamiltonian along the tetron path.

Produces a few PNGs in `results/` showing Bloch components under decoherence.
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

from tools import tetron_path_sim as tetron
from tools import tetron_mapped_sim as mapped
from tools import paper_params as P
from tools import lindblad_twolevel as lindblad


def build_H_timegrid(T_step=200.0, n_per_step=200):
    # get time grid and step/s fractions
    tgrid, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    # use default mapping constants if available
    try:
        import numpy as _np
        d = _np.load('results/mapping_fit_ABC.npz')
        A0 = float(d.get('A0', 0.0023))
        B0 = float(d.get('B0', 0.03))
        C0 = float(d.get('C0', 0.16))
    except Exception:
        A0, B0, C0 = 0.0023, 0.03, 0.16

    H_list = []
    # Build H(t) using gates_at and theta_from_time
    for step, s in zip(step_idx, slist):
        g1, g2, g3, g4 = tetron.gates_at(int(step), float(s))
        # simple surrogate S and M from gate amplitudes for demo
        S = float(g1 + g3)
        M = float(g2 - g4)
        theta = tetron.theta_from_time(int(step), float(s))
        dx = A0 * S
        dy = B0 * np.sin(theta)
        dz = C0 * M
        H = dx * lindblad.sigma_x + dy * lindblad.sigma_y + dz * lindblad.sigma_z
        H_list.append(H)
    return tgrid, H_list


def main():
    os.makedirs('results', exist_ok=True)
    tgrid, H_list = build_H_timegrid(T_step=200.0, n_per_step=200)

    cases = [
        ('coherent', 0.0, 0.0),
        ('dephasing', 0.5, 0.0),
        ('relax', 0.0, 0.1),
        ('both', 0.3, 0.05),
    ]

    for name, g_deph, g_rel in cases:
        bloch = lindblad.run_lindblad_time_series(H_list, dt=(tgrid[1]-tgrid[0]), gamma_deph=g_deph, gamma_relax=g_rel)
        plt.figure(figsize=(6,3))
        plt.plot(tgrid, bloch[:,0], label='Bx')
        plt.plot(tgrid, bloch[:,1], label='By')
        plt.plot(tgrid, bloch[:,2], label='Bz')
        plt.legend()
        plt.xlabel('t')
        plt.ylabel('Bloch')
        out = f'results/lindblad_{name}.png'
        plt.title(f'Lindblad demo: {name} (γϕ={g_deph}, γ1={g_rel})')
        plt.tight_layout()
        plt.savefig(out)
        plt.close()
        print('Saved', out)


if __name__ == '__main__':
    main()
