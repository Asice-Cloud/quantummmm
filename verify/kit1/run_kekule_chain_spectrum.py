#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

"""Toy 1D Kitaev chain with a Kekule-like modulation of the pairing channel.

Goal:
  - Start from a uniform topological Kitaev chain (t, Delta, mu).
  - In a finite central region D, modulate the complex pairing Delta_j
    along the bonds with a 2π phase winding (Kekule-like texture).
  - Diagonalize the resulting inhomogeneous BdG Hamiltonian and inspect
    the low-energy spectrum and spatial profile of near-zero modes.

This numerically illustrates the idea from kit-new3: using a spatial
texture in a "topological channel" (here the pairing) as an alternative
to purely temporal parameter paths.
"""


def build_kekule_chain_bdg(N: int,
                            t: float,
                            mu: float,
                            Delta_bonds: np.ndarray) -> np.ndarray:
    """Build 2N x 2N BdG Hamiltonian for an open Kitaev chain with
    bond-dependent complex pairing Delta_bonds[j] on bond (j, j+1).

    Nambu basis: Psi = (c_0,...,c_{N-1}, c_0^†,...,c_{N-1}^†)^T.
    """
    assert len(Delta_bonds) == N - 1

    H = np.zeros((2 * N, 2 * N), dtype=complex)

    # onsite chemical potential: h_ii = -mu/2, h̄_ii = +mu/2
    for i in range(N):
        H[i, i] = -mu / 2.0
        H[i + N, i + N] = mu / 2.0

    # nearest-neighbor hopping and pairing with inhomogeneous Delta_j
    for j in range(N - 1):
        k = j + 1
        Delta_j = Delta_bonds[j]

        # hopping -t: particle block h_{j,k} = -t
        H[j, k] += -t
        H[k, j] += -t
        # corresponding hole block
        H[j + N, k + N] += t
        H[k + N, j + N] += t

        # pairing Delta_j on bond (j, k)
        # particle -> hole block
        H[j, k + N] += Delta_j
        H[k, j + N] += -Delta_j
        # hole -> particle block
        H[j + N, k] += -Delta_j.conjugate()
        H[k + N, j] += Delta_j.conjugate()

    return H


def make_kekule_profile(N: int,
                         Delta0: complex,
                         jL: int,
                         jR: int) -> np.ndarray:
    """Construct a bond-dependent pairing profile Delta_bonds[j].

    - For bonds j < jL and j >= jR, use uniform pairing Delta0.
    - For bonds jL <= j < jR, let Delta_j = |Delta0| * exp(i * theta_j)
      where theta_j winds from 0 to 2π across the region [jL, jR).
    """
    assert 0 <= jL < jR <= N - 1
    Delta_bonds = np.empty(N - 1, dtype=complex)

    # length of distorted region in number of bonds
    Ld = jR - jL
    for j in range(N - 1):
        if j < jL or j >= jR:
            Delta_bonds[j] = Delta0
        else:
            # position within distorted region
            x = (j - jL) / max(Ld - 1, 1)
            theta = 2.0 * np.pi * x
            Delta_bonds[j] = np.abs(Delta0) * np.exp(1j * theta)

    return Delta_bonds


def lowest_spectrum_and_modes(H: np.ndarray, n_show: int = 8,
                              tol_zero: float = 1e-3):
    """Diagonalize H, return lowest-|E| eigenvalues and indices of
    near-zero modes (|E| < tol_zero).
    """
    w, v = np.linalg.eigh(H)
    idx = np.argsort(np.abs(w))
    w_sorted = w[idx]
    v_sorted = v[:, idx]
    near_zero_mask = np.abs(w_sorted) < tol_zero
    near_zero_indices = np.where(near_zero_mask)[0]
    return w_sorted[:n_show], v_sorted, near_zero_indices


def mode_density_on_chain(v_mode: np.ndarray) -> np.ndarray:
    """Compute site-resolved density |u_j|^2 + |v_j|^2 of a BdG eigenvector.

    v_mode has length 2N in Nambu basis; we split into particle/hole
    components and sum densities per site.
    """
    L = v_mode.shape[0]
    assert L % 2 == 0
    N = L // 2
    u = v_mode[:N]
    v = v_mode[N:]
    dens = np.abs(u) ** 2 + np.abs(v) ** 2
    # normalize to 1 for readability
    norm = dens.sum()
    if norm > 0:
        dens = dens / norm
    return dens


def main():
    # Chain parameters (choose a point well inside the topological phase)
    N = 60
    t = 1.0
    mu = 0.0
    Delta0 = 0.8  # real, positive background pairing

    # Define distorted region D on bonds [jL, jR)
    jL = 20
    jR = 40
    Delta_bonds = make_kekule_profile(N, Delta0, jL, jR)

    H = build_kekule_chain_bdg(N, t, mu, Delta_bonds)

    # Inspect low-energy spectrum
    w_low, v_all, near_zero_idx = lowest_spectrum_and_modes(H, n_show=12,
                                                            tol_zero=1e-3)
    print("=== Lowest-|E| eigenvalues (sorted by |E|) ===")
    print(np.round(w_low, 6))
    print(f"Number of near-zero modes (|E| < 1e-3): {len(near_zero_idx)}")

    if len(near_zero_idx) > 0:
        # Take up to first two near-zero modes and inspect their spatial profiles
        max_modes = min(2, len(near_zero_idx))
        x_sites = np.arange(N)
        plt.figure(figsize=(6, 3 * max_modes))
        for k in range(max_modes):
            idx_mode = near_zero_idx[k]
            v_mode = v_all[:, idx_mode]
            dens = mode_density_on_chain(v_mode)
            plt.subplot(max_modes, 1, k + 1)
            plt.plot(x_sites, dens, "-o", markersize=3)
            plt.axvspan(jL, jR, color="orange", alpha=0.15,
                        label="distorted region" if k == 0 else None)
            plt.ylabel(f"mode {k} density")
        plt.xlabel("site j")
        plt.tight_layout()
        plt.savefig("kekule_chain_near_zero_modes.png", dpi=200)
        plt.close()
        print("Saved kekule_chain_near_zero_modes.png (near-zero mode profiles).")


if __name__ == "__main__":
    main()
