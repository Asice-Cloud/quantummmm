#!/usr/bin/env python3
"""Compute Kitaev params and a simple topological indicator across (u,delta).

Maps the eight-vertex effective H to (t,Delta,mu), computes bulk gap by
sampling k, and applies the simple Kitaev criterion |mu| < 2|t|.

Saves results under results/topo_invariant/.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import math

# Pauli
I2 = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(a, b):
    return np.kron(a, b)


def H4_eight_vertex(u: float, delta: float) -> np.ndarray:
    return (
        math.cos(u) * kron(X, X)
        + 0.5 * math.sin(u) * (kron(Y, X) - kron(X, Y))
        + 0.5 * delta * (kron(Z, I2) - kron(I2, Z))
    )


def pauli_expand(H4: np.ndarray) -> dict:
    syms = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
    h = {}
    for a, A in syms.items():
        for b, B in syms.items():
            P = kron(A, B)
            h[(a, b)] = 0.25 * np.trace(P @ H4)
    return h


def map_to_kitaev_from_h(h: dict) -> tuple:
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


def bulk_gap(t: complex, Delta: complex, mu: complex, nk: int = 801) -> float:
    ks = np.linspace(0.0, np.pi, nk)
    minE = 1e9
    for k in ks:
        e = energy_k(t, Delta, mu, k)
        if e < minE:
            minE = e
    return float(minE)


def energy_k(t: complex, Delta: complex, mu: complex, k: float) -> float:
    eik = np.exp(1j * k)
    h_k = -mu / 2.0 - 0.5 * (t * eik + np.conjugate(t) * np.conjugate(eik))
    Delta_k = 0.5 * (Delta * eik - Delta * np.conjugate(eik))
    # BdG 2x2
    H = np.array([[h_k, Delta_k], [np.conjugate(Delta_k), -np.conjugate(h_k)]], dtype=complex)
    vals = np.linalg.eigvalsh(H)
    return float(np.min(np.abs(vals)))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--u-list', type=str, default='0,1.5708,3.1416,4.7124')
    p.add_argument('--delta-list', type=str, default='0,0.015,0.1')
    p.add_argument('--outdir', type=str, default='results/topo_invariant')
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    u_list = [float(x) for x in args.u_list.split(',') if x.strip()]
    d_list = [float(x) for x in args.delta_list.split(',') if x.strip()]

    records = []
    for u in u_list:
        for d in d_list:
            H4 = H4_eight_vertex(u, d)
            h = pauli_expand(H4)
            t, Delta, mu = map_to_kitaev_from_h(h)
            gap = bulk_gap(t, Delta, mu)
            topo = (abs(mu) < 2.0 * abs(t)) and (gap > 1e-12)
            records.append((u, d, complex(t), complex(Delta), complex(mu), float(abs(mu)), float(2.0 * abs(t)), float(gap), bool(topo)))
            print(f'u={u:.4g}, d={d:.4g} -> t={t:.4g}, Delta={Delta:.4g}, mu={mu:.4g}, gap={gap:.4g}, topo={topo}')

    np.savez(outdir / 'topo_invariant.npz', data=records)
    print('saved summary to', outdir / 'topo_invariant.npz')


if __name__ == '__main__':
    main()
