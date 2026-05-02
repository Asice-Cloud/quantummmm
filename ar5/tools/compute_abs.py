#!/usr/bin/env python3
"""Compute ABS examples: S-N interface and short Josephson junction.
Saves spectra and spatial densities into results/abs_sn.png and results/abs_josephson.png
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eigh

# Local builder for spatially varying BdG (N sites, bond indexed 0..N-2)
def build_bdg_spatial(N, t_hop, Delta_bond, mu_site):
    A = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)
    for i in range(N):
        A[i, i] = -mu_site[i]
        if i < N - 1:
            t = t_hop[i]
            A[i, i + 1] = -t
            A[i + 1, i] = -t
            Delta = Delta_bond[i]
            B[i, i + 1] = Delta
            B[i + 1, i] = -Delta
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    Hbdg = np.concatenate((top, bottom), axis=0)
    return Hbdg


def example_sn(N=200, t0=1.0, Delta_L=0.3, Delta_R=0.0, mu=0.0):
    # Left half is superconducting, right half normal
    mid = N // 2
    t_hop = np.full(N - 1, t0, dtype=complex)
    Delta_bond = np.zeros(N - 1, dtype=complex)
    Delta_bond[:mid - 1] = Delta_L
    Delta_bond[mid - 1:] = Delta_R
    mu_site = np.full(N, mu, dtype=complex)
    H = build_bdg_spatial(N, t_hop, Delta_bond, mu_site)
    vals, vecs = eigh(H)
    vals = np.real_if_close(vals)
    # sort by magnitude
    idx = np.argsort(np.abs(vals))
    return vals, vecs, idx


def example_josephson(N=200, t0=1.0, Delta=0.3, w=5, mu=0.0, phi=0.5):
    # Left/right superconductors with phases +/- phi/2, middle (w sites) normal
    left_end = (N - w) // 2
    right_start = left_end + w
    t_hop = np.full(N - 1, t0, dtype=complex)
    Delta_bond = np.zeros(N - 1, dtype=complex)
    # left bonds
    for i in range(0, left_end - 0):
        Delta_bond[i] = Delta * np.exp(-1j * phi / 2)
    # right bonds
    for i in range(right_start, N - 1):
        Delta_bond[i] = Delta * np.exp(1j * phi / 2)
    mu_site = np.full(N, mu, dtype=complex)
    H = build_bdg_spatial(N, t_hop, Delta_bond, mu_site)
    vals, vecs = eigh(H)
    vals = np.real_if_close(vals)
    idx = np.argsort(np.abs(vals))
    return vals, vecs, idx


def plot_sn(vals, vecs, idx, N, out_path):
    plt.figure(figsize=(6, 4))
    plt.plot(np.sort(vals), '.', markersize=2)
    plt.title('Spectrum: S-N interface')
    plt.xlabel('Level')
    plt.ylabel('Energy')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()

    # plot density of the lowest subgap mode (if any)
    # pick smallest nonzero mode
    mags = np.abs(vals)
    order = np.argsort(mags)
    for k in range(4):
        i = order[k]
        E = vals[i]
        if abs(E) < 1e-6:
            continue
        # get wavefunction (Nambu)
        psi = vecs[:, i]
        # particle components are first N, hole components are last N
        rho_particle = np.abs(psi[:N])**2
        rho_hole = np.abs(psi[N:])**2
        rho = rho_particle + rho_hole
        plt.figure(figsize=(6,3))
        plt.plot(np.arange(N), rho)
        plt.title(f'Density of mode E={E:.4e}')
        plt.xlabel('Site')
        plt.ylabel('Particle+hole density')
        plt.tight_layout()
        out_file = out_path.replace('.png', f'_mode{k}.png')
        plt.savefig(out_file)
        plt.close()


def plot_josephson(vals_phi, phi_vals, out_path):
    energies = [np.sort(v) for v in vals_phi]
    # take lowest positive magnitude eigenvalue for each phi
    lows = []
    for v in energies:
        # v is sorted, but may include negative; pick smallest abs positive
        mags = np.abs(v)
        idx = np.argmin(mags)
        lows.append(v[idx])
    plt.figure(figsize=(6,4))
    plt.plot(phi_vals, lows, '-o')
    plt.xlabel('phi')
    plt.ylabel('lowest |E|')
    plt.title('ABS vs Josephson phase')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


if __name__ == '__main__':
    os.makedirs('results', exist_ok=True)
    # run S-N example
    N = 200
    vals, vecs, idx = example_sn(N=N)
    sn_path = 'results/abs_sn_spectrum.png'
    plot_sn(vals, vecs, idx, N, sn_path)
    print(f'S-N spectrum saved to {sn_path}')

    # run Josephson scan
    phi_vals = np.linspace(0, np.pi, 13)
    vals_phi = []
    for phi in phi_vals:
        vals, vecs, idx = example_josephson(N=N, phi=phi)
        vals_phi.append(vals)
    j_path = 'results/abs_josephson_phi.png'
    plot_josephson(vals_phi, phi_vals, j_path)
    print(f'Josephson scan saved to {j_path}')
