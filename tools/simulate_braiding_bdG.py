#!/usr/bin/env python3
"""Simulate time-dependent BdG braiding-like gate sequence and compare projected U to ideal braid.

This uses a simplified mapping: three local bond couplings around the chain mid act like the
gate-controlled couplings in the supplement. The path g1,g2,g3 is cycled to emulate teleportation.
"""
import argparse
import numpy as np
from scipy.linalg import eigh, expm
from pathlib import Path


def build_bdg(L, t_base, Delta, mu_vec, bond_extras=None):
    H0 = np.zeros((L, L), dtype=complex)
    for i in range(L-1):
        tt = -t_base
        if bond_extras and i in bond_extras:
            tt += -bond_extras[i]
        H0[i, i+1] = tt
        H0[i+1, i] = tt
    for i in range(L):
        H0[i, i] += -mu_vec[i]
    D = np.zeros((L, L), dtype=complex)
    for i in range(L-1):
        dd = 0
        if bond_extras and i in bond_extras:
            dd += bond_extras[i]
        D[i, i+1] = dd
        D[i+1, i] = -dd
    top = np.hstack([H0, D])
    bottom = np.hstack([-D.conj(), -H0.T])
    return np.vstack([top, bottom])


def gate_functions(t, T):
    # three-step cyclic schedule for g1,g2,g3 over time [0,T)
    # each step lasts T/3; within each step we linearly ramp two gates
    s = t / (T / 3.0)
    step = int(np.floor(s))
    frac = s - step
    g = [0.0, 0.0, 0.0, 1.0]
    # pattern: (0,0,1)->(1,0,0)->(0,1,0)->(0,0,1)
    if step == 0:
        # (0,0,1) -> (1,0,0): ramp g1 up, ramp g3 down
        g1 = frac
        g3 = 1.0 - frac
        g = [g1, 0.0, g3, 1.0]
    elif step == 1:
        # (1,0,0) -> (0,1,0): ramp g2 up, ramp g1 down
        g2 = frac
        g1 = 1.0 - frac
        g = [g1, g2, 0.0, 1.0]
    else:
        # (0,1,0) -> (0,0,1): ramp g3 up, ramp g2 down
        g3 = frac
        g2 = 1.0 - frac
        g = [0.0, g2, g3, 1.0]
    return g


def ideal_braid_2x2():
    sy = np.array([[0, -1j],[1j, 0]], dtype=complex)
    B = np.cos(np.pi/4.0) * np.eye(2, dtype=complex) + 1j * np.sin(np.pi/4.0) * sy
    return B


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h', type=float, required=True)
    p.add_argument('--L', type=int, default=120)
    p.add_argument('--t', type=float, default=2.1)
    p.add_argument('--Delta', type=float, default=0.1)
    p.add_argument('--T', type=float, default=300.0)
    p.add_argument('--steps', type=int, default=300)
    args = p.parse_args()

    L = args.L
    mid = L//2
    mu0 = 0.0
    mu_vec = np.zeros(L) + mu0
    mu_vec[mid] += 2.0 * args.h

    dt = args.T / float(args.steps)
    dim = 2*L
    U = np.eye(dim, dtype=complex)

    # map δ components to bond extras near mid: use bonds mid-2, mid-1, mid
    bond_indices = [mid-2, mid-1, mid]
    tc = 1.0

    print(f'Simulating braiding-like path: L={L}, h={args.h}, T={args.T}, steps={args.steps}')
    for k in range(args.steps):
        tnow = (k + 0.5) * dt
        g = gate_functions(tnow, args.T)
        # δ components
        delta = [tc * (g[0]*g[3]), tc * (g[1]*g[3]), tc * (g[2]*g[3])]
        # define bond extras: add delta as extra pairing/hopping on these bonds
        bond_extras = {bond_indices[i]: delta[i] for i in range(len(bond_indices))}
        Hbdg = build_bdg(L, args.t, args.Delta, mu_vec, bond_extras)
        U = expm(-1j * Hbdg * dt) @ U

    # initial low-energy subspace from H(0)
    g0 = gate_functions(0.0, args.T)
    delta0 = [tc * (g0[0]*g0[3]), tc * (g0[1]*g0[3]), tc * (g0[2]*g0[3])]
    bond_extras0 = {bond_indices[i]: delta0[i] for i in range(len(bond_indices))}
    H0 = build_bdg(L, args.t, args.Delta, mu_vec, bond_extras0)
    evals0, evecs0 = eigh(H0)
    idx0 = np.argsort(np.abs(evals0))
    Vlow0 = evecs0[:, idx0[:2]]

    Ueff = Vlow0.conj().T @ U @ Vlow0
    B2 = ideal_braid_2x2()
    frob = np.linalg.norm(Ueff - B2)
    frob_rel = frob / np.linalg.norm(B2)
    fid_like = np.real(np.trace(Ueff.conj().T @ B2)) / 2.0

    out = Path('results')
    out.mkdir(exist_ok=True)
    np.savez(out / f'braid_sim_h{args.h:.3f}_L{L}.npz', U=U, Ueff=Ueff, Vlow0=Vlow0)
    print('Saved results to', out / f'braid_sim_h{args.h:.3f}_L{L}.npz')
    print('frob_diff=',frob,'  rel=',frob_rel,'  fid_like=',fid_like)


if __name__ == '__main__':
    main()
