#!/usr/bin/env python3
"""Construct a BdG chain with a central paired region and save spectrum + density.

Usage: python3 tools/run_midpair_demo.py
"""
import os
import sys
import numpy as np
from numpy.linalg import eigh

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def build_midpair_chain(N, t_val=1.0, Delta_center=0.5, center_width=20, mu=0.0):
    # bonds are N-1
    t_hop = np.zeros(N-1, dtype=complex)
    Delta_bond = np.zeros(N-1, dtype=complex)
    mu_site = np.zeros(N, dtype=float)
    # fill t uniformly
    t_hop[:] = t_val
    # central bond region
    mid = N // 2
    half = center_width // 2
    start = max(0, mid - half)
    end = min(N-1, mid + half)
    for i in range(start, end):
        Delta_bond[i] = Delta_center
    # mu all zero
    A = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)
    for i in range(N):
        A[i,i] = -mu_site[i]
        if i < N-1:
            t = t_hop[i]
            A[i,i+1] = -t
            A[i+1,i] = -t
            Delta = Delta_bond[i]
            B[i,i+1] = Delta
            B[i+1,i] = -Delta
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    H = np.concatenate((top, bottom), axis=0)
    return H, t_hop, Delta_bond, mu_site


def run():
    os.makedirs('results', exist_ok=True)
    N = 100
    H, t_hop, Delta_bond, mu_site = build_midpair_chain(N, t_val=1.0, Delta_center=0.5, center_width=20, mu=0.0)
    vals, vecs = eigh(H)
    # save spectrum plot if matplotlib present
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6,4))
        plt.plot(np.sort(vals), '.', markersize=3)
        plt.title('Mid-pair spectrum')
        plt.tight_layout()
        plt.savefig('results/xyz_midpair_spectrum.png')
        plt.close()
    except Exception:
        pass

    idx = np.argsort(np.abs(vals))
    low_idx = idx[0]
    vec0 = vecs[:, low_idx]
    Nloc = vec0.shape[0] // 2
    u = vec0[:Nloc]
    v = vec0[Nloc:]
    dens = np.abs(u)**2 + np.abs(v)**2
    np.savetxt('results/xyz_midpair_density0.txt', dens)

    # write report
    peak = int(np.argmax(dens))
    peak_val = float(dens[peak])
    right = peak
    while right < N and dens[right] >= peak_val/np.e:
        right += 1
    xi_est = (right - peak) if right < N else float('nan')

    with open('results/xyz_midpair_report.txt', 'w') as f:
        f.write(f'N={N}, center_width=20\n')
        f.write(f'center Delta=0.5, t=1.0, mu=0\n')
        f.write('\nLowest energies (abs-sorted first 6):\n')
        for i in idx[:6]:
            f.write(f'{vals[i]:.6e}\n')
        f.write(f'\npeak_idx={peak}\npeak_val={peak_val:.6e}\nxi_est={xi_est}\n')

    print('Saved results/xyz_midpair_spectrum.png, results/xyz_midpair_density0.txt, results/xyz_midpair_report.txt')


if __name__ == '__main__':
    run()
