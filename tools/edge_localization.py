#!/usr/bin/env python3
"""Compute and plot edge localization of the lowest BdG eigenstate.

Uses the eight-vertex effective mapping to Kitaev params, builds BdG chain,
and scans chain lengths to report E0(L) and endpoint weight.

Saves results to `results/edge_localization/` by default.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Pauli
I2 = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(a, b):
    return np.kron(a, b)


def H4_eight_vertex(u: float, delta: float) -> np.ndarray:
    # using the effective Pauli expansion from ybe223.md
    H = (
        np.cos(u) * kron(X, X)
        + 0.5 * np.sin(u) * (kron(Y, X) - kron(X, Y))
        + 0.5 * delta * (kron(Z, I2) - kron(I2, Z))
    )
    return H


def pauli_expand(H4: np.ndarray) -> dict:
    # returns h[(a,b)] with a,b in {'I','X','Y','Z'} using 1/4 Tr(P H)
    syms = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
    h = {}
    for a, A in syms.items():
        for b, B in syms.items():
            P = kron(A, B)
            h[(a, b)] = 0.25 * np.trace(P @ H4)
    return h


def map_to_kitaev_from_h(h: dict) -> tuple:
    # indices: 'X'->xx, 'Y'->yy etc.
    c_xx = h.get(('X', 'X'), 0.0)
    c_yy = h.get(('Y', 'Y'), 0.0)
    c_xy = h.get(('X', 'Y'), 0.0)
    c_yx = h.get(('Y', 'X'), 0.0)
    c_zz = h.get(('Z', 'Z'), 0.0)
    c_z0 = h.get(('Z', 'I'), 0.0)
    c_0z = h.get(('I', 'Z'), 0.0)
    t = c_xx + c_yy + 1j * (c_xy - c_yx)
    Delta = c_xx - c_yy - 1j * (c_xy + c_yx)
    mu = 4.0 * c_zz - 2.0 * (c_z0 + c_0z)
    return complex(t), complex(Delta), complex(mu)


def build_bdg_chain(L: int, t: complex, Delta: complex, mu: complex, open_boundary: bool = True) -> np.ndarray:
    # Build single-particle BdG Hamiltonian (2L x 2L)
    h = np.zeros((L, L), dtype=complex)
    Delta_mat = np.zeros((L, L), dtype=complex)
    for i in range(L):
        h[i, i] = -mu / 2.0
    for i in range(L - 1):
        h[i, i + 1] = -t / 2.0
        h[i + 1, i] = -t.conjugate() / 2.0
        Delta_mat[i, i + 1] = Delta / 2.0
        Delta_mat[i + 1, i] = -Delta / 2.0
    if not open_boundary:
        h[0, L - 1] = -t / 2.0
        h[L - 1, 0] = -t.conjugate() / 2.0
        Delta_mat[0, L - 1] = Delta / 2.0
        Delta_mat[L - 1, 0] = -Delta / 2.0
    top = np.hstack([h, Delta_mat])
    bottom = np.hstack([-Delta_mat.conjugate(), -h.T])
    H_bdg = np.vstack([top, bottom])
    H_bdg = 0.5 * (H_bdg + H_bdg.conj().T)
    return H_bdg


def edge_weights_from_eigvec(v: np.ndarray, L: int, edge_sites: int = 1) -> tuple:
    # v length 2L, first L are particle, next L are hole
    probs = np.abs(v[:L]) ** 2 + np.abs(v[L:]) ** 2
    left = float(np.sum(probs[:edge_sites]))
    right = float(np.sum(probs[-edge_sites:]))
    return left, right


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--delta', type=float, default=0.015)
    p.add_argument('--u', type=float, default=0.0, help='path parameter u to analyze')
    p.add_argument('--Ls', type=str, default='40,80,160,320', help='comma list of chain lengths')
    p.add_argument('--edge-sites', type=int, default=1)
    p.add_argument('--outdir', type=str, default='results/edge_localization')
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    Ls = [int(x) for x in args.Ls.split(',') if x.strip()]
    u = float(args.u)
    delta = float(args.delta)

    # build H4 and map
    H4 = H4_eight_vertex(u, delta)
    h = pauli_expand(H4)
    t, Delta, mu = map_to_kitaev_from_h(h)

    results = {'L': [], 'E0': [], 'left_w': [], 'right_w': []}

    for L in Ls:
        Hbdg = build_bdg_chain(L, t, Delta, mu, open_boundary=True)
        E, V = eigh(Hbdg)
        idx = int(np.argmin(np.abs(E)))
        e0 = float(np.abs(E[idx]))
        v = V[:, idx]
        left, right = edge_weights_from_eigvec(v, L, edge_sites=args.edge_sites)
        results['L'].append(L)
        results['E0'].append(e0)
        results['left_w'].append(left)
        results['right_w'].append(right)

    np.savez(outdir / f'edge_loc_u{u:.3f}_d{delta:.3g}.npz', **results)

    # plot E0 vs L
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5), dpi=150)
    ax[0].plot(results['L'], results['E0'], '-o')
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_xlabel('L')
    ax[0].set_ylabel('E0 (abs)')
    ax[0].set_title(f'E0 vs L (u={u:.3f}, δ={delta:.3g})')

    ax[1].plot(results['L'], results['left_w'], '-o', label='left')
    ax[1].plot(results['L'], results['right_w'], '-o', label='right')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('L')
    ax[1].set_ylabel('edge weight')
    ax[1].set_title('Edge weight of lowest eigenstate')
    ax[1].legend()

    fig.tight_layout()
    fig.savefig(outdir / f'edge_loc_u{u:.3f}_d{delta:.3g}.png')
    plt.close(fig)

    print('saved results to', outdir)


if __name__ == '__main__':
    main()
