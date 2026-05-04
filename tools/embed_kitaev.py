#!/usr/bin/env python3
"""Small BdG Kitaev-chain simulator for LDOS snapshots under gate snapshots.

This script constructs a spinless Kitaev chain BdG Hamiltonian with
site-dependent hopping/pairing and a local QD potential at the left end
to generate an ABS. It computes LDOS (zero-energy spectral weight)
for several gate snapshots and saves figures.
"""
import os
import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt


def build_bdg(mu, t_links, Delta_links):
    # mu: length L array
    L = len(mu)
    h = np.zeros((L, L), dtype=complex)
    Delta = np.zeros((L, L), dtype=complex)
    # onsite
    for i in range(L):
        h[i, i] = -mu[i]
    # hopping and pairing on links i <-> i+1
    for i in range(L - 1):
        t = t_links[i]
        d = Delta_links[i]
        h[i, i + 1] = -t
        h[i + 1, i] = -t
        # p-wave pairing: antisymmetric pairing matrix
        Delta[i, i + 1] = d
        Delta[i + 1, i] = -d
    # BdG block matrix
    top = np.hstack([h, Delta])
    bottom = np.hstack([-Delta.conj(), -h.T.conj()])
    H = np.vstack([top, bottom])
    return H


def compute_zero_ldos(H, eta=1e-2):
    # eigen-decomposition
    E, V = eigh(H)
    L2 = H.shape[0]
    L = L2 // 2
    # spectral weight at E=0 with lorentzian broadening
    weights = np.zeros(L)
    for n in range(len(E)):
        En = E[n]
        # lorentzian weight at zero energy
        lor = eta / (np.pi * (En * En + eta * eta))
        vec = V[:, n]
        # particle + hole weight at site j is |u_j|^2 + |v_j|^2
        for j in range(L):
            uj = vec[j]
            vj = vec[L + j]
            weights[j] += (np.abs(uj) ** 2 + np.abs(vj) ** 2) * lor
    return weights, E


def snapshot_ldos(L=80, t0=1.0, Delta0=0.3, mu0=0.0, VD=2.0, qd_width=3):
    # construct base arrays
    mu = mu0 * np.ones(L)
    t_links = t0 * np.ones(L - 1)
    Delta_links = Delta0 * np.ones(L - 1)

    # QD at left end: apply potential well on first qd_width sites
    for i in range(qd_width):
        mu[i] = mu[i] - VD

    # gate snapshots (g1,g2,g3,g4) -- simplified mapping to left links
    snapshots = {
        'init': (1.0, 1.0, 0.0, 0.0),
        'after_step1': (0.0, 1.0, 1.0, 0.0),
        'after_step2': (1.0, 0.0, 1.0, 0.0),
        'after_step3': (1.0, 1.0, 0.0, 0.0),
    }

    os.makedirs('results', exist_ok=True)
    for name, (g1, g2, g3, g4) in snapshots.items():
        # map gates to left-most t_links changes (simple model)
        t_links_mod = t_links.copy()
        # modify first two links according to g1 and g3
        t_links_mod[0] = t0 * g1
        if L > 2:
            t_links_mod[1] = t0 * g3
        # pairing scaled similarly
        Delta_mod = Delta_links.copy()
        Delta_mod[0] = Delta0 * (g1 if g1 > 0 else 0.01)
        if L > 2:
            Delta_mod[1] = Delta0 * (g3 if g3 > 0 else 0.01)

        H = build_bdg(mu, t_links_mod, Delta_mod)
        ldos, E = compute_zero_ldos(H, eta=1e-3)

        plt.figure(figsize=(8,3))
        plt.plot(np.arange(L), ldos, '-o')
        plt.xlabel('site')
        plt.ylabel('LDOS(E=0)')
        plt.title(f'LDOS snapshot: {name} (g1,g2,g3,g4)={g1,g2,g3,g4}')
        out = f'results/ldos_snapshot_{name}.png'
        plt.tight_layout()
        plt.savefig(out)
        plt.close()
        print('Saved', out)

        # also save lowest eigenvalues (near zero)
        ev_file = f'results/eigs_{name}.txt'
        with open(ev_file, 'w') as fh:
            for val in E[:10]:
                fh.write(f"{val:.6e}\n")
        print('Saved', ev_file)


if __name__ == '__main__':
    snapshot_ldos()
