#!/usr/bin/env python3
"""Demo: Baxterization-like mapping R(u)=exp(i g(u) H0) -> time-dependent H(t).

This script builds a BdG `H0` (using the existing simple builder), constructs
R(u)=expm(1j * g(u) * H0) and shows the equivalent H(t)=g'(u) u'(t) H0.
It then simulates the time evolution U = T exp(-i ∫ H dt) by discretization
and computes the projected 2x2 `U_eff` fidelity against the ideal braid.

This is a constructive demonstration: if an algebraic R(u) can be written as
an exponential of a (local) generator H0 with scalar prefactor g(u), then
letting the spectral parameter u depend on time yields a physically
interpretable time-dependent Hamiltonian H(t).
"""
import argparse
import numpy as np
from scipy.linalg import expm, eigh
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


def ideal_braid_2x2():
    sy = np.array([[0, -1j],[1j, 0]], dtype=complex)
    B = np.cos(np.pi/4.0) * np.eye(2, dtype=complex) + 1j * np.sin(np.pi/4.0) * sy
    return B


def default_u_of_t(t, T):
    # example: sweep u from 0 to u0 smoothly (sine ramp)
    u0 = np.pi/2.0
    return u0 * (0.5 - 0.5 * np.cos(np.pi * t / T))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h', type=float, required=True)
    p.add_argument('--L', type=int, default=40)
    p.add_argument('--t', type=float, default=2.1)
    p.add_argument('--Delta', type=float, default=0.1)
    p.add_argument('--T', type=float, default=10.0)
    p.add_argument('--steps', type=int, default=200)
    args = p.parse_args()

    L = args.L
    mid = L//2
    mu0 = 0.0
    mu_vec = np.zeros(L) + mu0
    mu_vec[mid] += 2.0 * args.h

    # define a base generator H0: BdG at small initial bond extras
    bond_indices = [mid-2, mid-1, mid]
    bond_extras0 = {bond_indices[0]: 0.0, bond_indices[1]: 0.0, bond_indices[2]: 0.0}
    H0 = build_bdg(L, args.t, args.Delta, mu_vec, bond_extras0)

    # we will treat R(u) = exp(i * g(u) * H0); pick g(u)=u for demo
    def g(u):
        return u

    def dg_du(u):
        return 1.0

    def u_of_t(t):
        return default_u_of_t(t, args.T)

    def du_dt(t):
        # derivative of default_u_of_t
        return (np.pi * (np.pi/2.0) * 0.5 * np.sin(np.pi * t / args.T)) / args.T if args.T!=0 else 0.0

    dt = args.T / float(args.steps)
    dim = H0.shape[0]
    U = np.eye(dim, dtype=complex)

    for k in range(args.steps):
        tnow = (k + 0.5) * dt
        u = u_of_t(tnow)
        # H(t) = dg/du * du/dt * H0
        Ht = dg_du(u) * du_dt(tnow) * H0
        U = expm(-1j * Ht * dt) @ U

    # build low-energy projector from H0
    evals0, evecs0 = eigh(H0)
    idx0 = np.argsort(np.abs(evals0))
    Vlow0 = evecs0[:, idx0[:2]]

    Ueff = Vlow0.conj().T @ U @ Vlow0
    B2 = ideal_braid_2x2()
    frob = np.linalg.norm(Ueff - B2)
    fid_like = np.real(np.trace(Ueff.conj().T @ B2)) / 2.0

    out = Path('results')
    out.mkdir(exist_ok=True)
    np.savez(out / f'baxterize_demo_h{args.h:.3f}_L{L}.npz', U_final=U, Ueff_final=Ueff, Vlow0=Vlow0)
    print('Saved results to', out / f'baxterize_demo_h{args.h:.3f}_L{L}.npz')
    print('frob_diff=',frob,'  fid_like=',fid_like)


if __name__ == '__main__':
    main()
