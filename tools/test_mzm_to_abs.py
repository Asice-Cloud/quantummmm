#!/usr/bin/env python3
"""Test MZM -> ABS by adding Majorana-like coupling to the 4x4 H along the path.

Produces:
- results/mzm_abs/min_d_vs_epsilon.npz  (data)
- results/mzm_abs/min_d_vs_epsilon.png (plot)
- results/mzm_abs/d_u_examples.png     (example |d(u)| curves)

Usage: python3 tools/test_mzm_to_abs.py
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import math
from scipy.linalg import eigh

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Pauli
I2 = np.array([[1, 0], [0, 1]], dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(a, b):
    return np.kron(a, b)


def h4_from_path(g1, g2, g3, g4, tc=1.0):
    # same microscopic mapping used elsewhere in the repo
    H = (
        -(g1 * g2) * kron(SZ, I2)
        - (g1 * g3) * kron(SY, SX)
        + (g1 * g4) * kron(SY, SY)
        - (g2 * g3) * kron(SX, SX)
        - (g2 * g4) * kron(SX, SY)
        - (g3 * g4) * kron(I2, SZ)
    )
    return tc * H


def path_segment(k: int, u: float):
    if k == 1:
        return u, 0.0, 1.0 - u, 1.0
    if k == 2:
        return 1.0 - u, u, 0.0, 1.0
    if k == 3:
        return 0.0, 1.0 - u, u, 1.0
    raise ValueError("bad segment")


def project_to_two_level(H4: np.ndarray) -> np.ndarray:
    psi01 = np.array([0, 1, 0, 0], dtype=complex)
    psi10 = np.array([0, 0, 1, 0], dtype=complex)
    Heff = np.zeros((2, 2), dtype=complex)
    Heff[0, 0] = psi01.conj().T @ H4 @ psi01
    Heff[0, 1] = psi01.conj().T @ H4 @ psi10
    Heff[1, 0] = psi10.conj().T @ H4 @ psi01
    Heff[1, 1] = psi10.conj().T @ H4 @ psi10
    return Heff


def d_from_Heff(Heff: np.ndarray):
    d0 = 0.5 * np.trace(Heff)
    dx = 0.5 * np.trace(Heff @ SX)
    dy = 0.5 * np.trace(Heff @ SY)
    dz = 0.5 * np.trace(Heff @ SZ)
    return float(dx.real), float(dy.real), float(dz.real), float(d0.real)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=300, help="u grid points per segment")
    p.add_argument("--eps-list", type=str, default="0.0,0.02,0.05,0.1,0.2,0.5", help="comma list of coupling epsilons")
    p.add_argument("--tc", type=float, default=1.0)
    p.add_argument("--outdir", type=str, default="results/mzm_abs")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    us = np.linspace(0.0, 1.0, args.N, endpoint=False)
    # build full u array for three segments concatenated (for plotting)
    u_full = []
    H4s = []
    for k in (1, 2, 3):
        for u in us:
            u_full.append((k - 1) + u)
            g1, g2, g3, g4 = path_segment(k, float(u))
            H4s.append(h4_from_path(g1, g2, g3, g4, tc=args.tc))
    H4s = np.array(H4s)
    u_full = np.array(u_full)

    eps_list = [float(x) for x in args.eps_list.split(",") if x.strip()]

    # coupling operator chosen to map to sigma_z in logical subspace: (Z⊗I - I⊗Z)/2
    coupling_op = 0.5 * (kron(SZ, I2) - kron(I2, SZ))

    min_d_for_eps = []
    example_curves = {}

    for eps in eps_list:
        ds = []
        for H4 in H4s:
            H4_c = H4 + eps * coupling_op
            Heff = project_to_two_level(H4_c)
            dx, dy, dz, d0 = d_from_Heff(Heff)
            ds.append(np.sqrt(dx * dx + dy * dy + dz * dz))
        ds = np.array(ds)
        min_d_for_eps.append(float(np.min(ds)))
        # store example for eps == 0 and a middle epsilon for illustration
        if eps == 0.0 or eps == eps_list[len(eps_list) // 2]:
            example_curves[eps] = ds.copy()

    min_d_for_eps = np.array(min_d_for_eps)

    # save data
    np.savez(outdir / "min_d_vs_epsilon.npz", eps_list=np.array(eps_list), min_d=min_d_for_eps, u_full=u_full)

    # plot min|d| vs epsilon
    fig, ax = plt.subplots(figsize=(5, 3.5), dpi=150)
    ax.plot(eps_list, min_d_for_eps, '-o')
    ax.set_xlabel('coupling epsilon')
    ax.set_ylabel('min_u |d(u)|')
    ax.set_title('Minimum |d| along path vs Majorana coupling')
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(outdir / 'min_d_vs_epsilon.png')
    plt.close(fig)

    # plot example |d(u)| curves
    fig, ax = plt.subplots(figsize=(6, 3.5), dpi=150)
    for eps, ds in example_curves.items():
        ax.plot(u_full, ds, label=f'eps={eps}')
    ax.set_xlabel('path parameter (segment.u)')
    ax.set_ylabel('|d(u)|')
    ax.set_title('|d(u)| along path (examples)')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(outdir / 'd_u_examples.png')
    plt.close(fig)

    print('saved results to', outdir)


if __name__ == '__main__':
    main()
