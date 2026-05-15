#!/usr/bin/env python3
"""Compute LDOS (energy vs position) for the effective Kitaev chain
constructed from the eight-vertex mapping used in the notes.

Saves a heatmap and a zero-energy spatial profile.
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


def build_bdg_chain(L: int, t: complex, Delta: complex, mu: complex, open_boundary: bool = True) -> np.ndarray:
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


def compute_ldos(Hbdg: np.ndarray, energies: np.ndarray, eta: float = 1e-3) -> np.ndarray:
    # LDOS at each site: sum_n (|u_n(i)|^2 + |v_n(i)|^2) * Lorentzian(E- E_n)
    E, V = eigh(Hbdg)
    # ensure ascending
    idx = np.argsort(E)
    E = E[idx]
    V = V[:, idx]
    L = Hbdg.shape[0] // 2
    nE = energies.size
    ldos = np.zeros((L, nE), dtype=float)
    lorentz_norm = eta / math.pi
    for n in range(E.size):
        En = float(E[n])
        probs = np.abs(V[:L, n]) ** 2 + np.abs(V[L:, n]) ** 2
        # Lorentzian over energies
        lor = lorentz_norm / ((energies - En) ** 2 + eta * eta)
        ldos += np.outer(probs, lor)
    return energies, ldos, E


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--u', type=float, required=True)
    p.add_argument('--delta', type=float, default=0.015)
    p.add_argument('--L', type=int, default=160)
    p.add_argument('--Emin', type=float, default=-0.02)
    p.add_argument('--Emax', type=float, default=0.02)
    p.add_argument('--nE', type=int, default=801)
    p.add_argument('--eta', type=float, default=1e-4)
    p.add_argument('--outdir', type=str, default='results/ldos')
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    H4 = H4_eight_vertex(args.u, args.delta)
    h = pauli_expand(H4)
    t, Delta, mu = map_to_kitaev_from_h(h)
    Hbdg = build_bdg_chain(args.L, t, Delta, mu, open_boundary=True)

    energies = np.linspace(args.Emin, args.Emax, args.nE)
    energies, ldos, spec = compute_ldos(Hbdg, energies, eta=args.eta)

    # save
    fname = outdir / f'ldos_u{args.u:.4g}_d{args.delta:.3g}_L{args.L}.npz'
    np.savez(fname, energies=energies, ldos=ldos, spec=spec)

    # plot heatmap
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), dpi=150)
    im = ax[0].imshow(ldos, aspect='auto', origin='lower', extent=[energies[0], energies[-1], 0, args.L])
    ax[0].set_xlabel('Energy')
    ax[0].set_ylabel('Site index')
    ax[0].set_title(f'LDOS (u={args.u}, d={args.delta}, L={args.L})')
    fig.colorbar(im, ax=ax[0], fraction=0.046)

    # zero-energy spatial profile
    # find index closest to zero energy
    idx0 = np.argmin(np.abs(energies))
    ax[1].plot(np.arange(args.L), ldos[:, idx0])
    ax[1].set_xlabel('Site index')
    ax[1].set_ylabel('LDOS(E=0)')
    ax[1].set_title('LDOS at E=0')

    fig.tight_layout()
    outpng = outdir / f'ldos_u{args.u:.4g}_d{args.delta:.3g}_L{args.L}.png'
    fig.savefig(outpng)
    plt.close(fig)

    print('saved', fname, outpng)


if __name__ == '__main__':
    main()
