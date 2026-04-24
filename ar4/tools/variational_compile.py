#!/usr/bin/env python3
"""Variationally compile a short pulse sequence in the projected 2x2 subspace.

Workflow:
- Load Vlow0 from results/braid_sim_h0.150_L120.npz (must exist)
- Build a small set of local full-space BdG operators near chain mid
- Project these to 2x2 and optimize sequence coefficients in 2x2
- Lift optimized sequence to full space and compute full fidelity
"""
import argparse
import numpy as np
from scipy.linalg import expm, eigh
from scipy.optimize import minimize
from pathlib import Path


def build_bdg_ops_small(L, t_base, Delta, mid):
    # Create small set of full-space BdG operators: three bond pairing/hopping extras and on-site mu at mid
    dim = 2 * L
    ops = []
    def make_op(bond_extra_idx=None, mu_site=None, val=1.0):
        H0 = np.zeros((L, L), dtype=complex)
        for i in range(L-1):
            tt = -t_base
            if bond_extra_idx is not None and i == bond_extra_idx:
                tt += -val
            H0[i, i+1] = tt
            H0[i+1, i] = tt
        mu_vec = np.zeros(L)
        if mu_site is not None:
            mu_vec[mu_site] = val
        for i in range(L):
            H0[i, i] += -mu_vec[i]
        D = np.zeros((L, L), dtype=complex)
        for i in range(L-1):
            dd = 0
            if bond_extra_idx is not None and i == bond_extra_idx:
                dd += val
            D[i, i+1] = dd
            D[i+1, i] = -dd
        top = np.hstack([H0, D])
        bottom = np.hstack([-D.conj(), -H0.T])
        return np.vstack([top, bottom])

    # choose bonds in a window around mid and a few on-site mu nearby
    bond_idxs = list(range(max(0, mid-4), min(L-1, mid+3)))
    for bi in bond_idxs:
        ops.append((f'bond_{bi}', make_op(bond_extra_idx=bi, val=1.0)))
    # on-site mu at mid-1, mid, mid+1 (if in range)
    for s in [mid-1, mid, mid+1]:
        if 0 <= s < L:
            ops.append((f'mu_{s}', make_op(mu_site=s, val=1.0)))
    return ops


def ideal_braid_2x2():
    sy = np.array([[0, -1j],[1j, 0]], dtype=complex)
    B = np.cos(np.pi/4.0) * np.eye(2, dtype=complex) + 1j * np.sin(np.pi/4.0) * sy
    return B


def run_projected_optimization(Vlow, basis_proj, B2, k, dt):
    n_ops = len(basis_proj)
    # params vector of length k*n_ops
    x0 = np.zeros(k * n_ops)

    def unpack(x):
        return x.reshape((k, n_ops))

    def loss_fn(x):
        coeffs = unpack(x)
        U = np.eye(2, dtype=complex)
        for s in range(k):
            Hs = sum(coeffs[s, j] * basis_proj[j] for j in range(n_ops))
            U = expm(-1j * Hs * dt) @ U
        # fidelity-like
        fid_like = np.real(np.trace(U.conj().T @ B2)) / 2.0
        return 1.0 - fid_like

    res = minimize(loss_fn, x0, method='L-BFGS-B', options={'maxiter':200, 'ftol':1e-8})
    return res


def lift_and_eval(seq_coeffs, full_ops, Vlow, B2, dt):
    # compute full U_seq and project
    dim = full_ops[0].shape[0]
    k, n_ops = seq_coeffs.shape
    U = np.eye(dim, dtype=complex)
    for s in range(k):
        Hs_full = sum(seq_coeffs[s, j] * full_ops[j] for j in range(n_ops))
        U = expm(-1j * Hs_full * dt) @ U
    Ueff = Vlow.conj().T @ U @ Vlow
    B2_full = B2
    frob = np.linalg.norm(Ueff - B2_full)
    fid_like = np.real(np.trace(Ueff.conj().T @ B2_full)) / 2.0
    return frob, fid_like, Ueff


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h', type=float, required=True)
    p.add_argument('--L', type=int, default=120)
    p.add_argument('--t', type=float, default=2.1)
    p.add_argument('--Delta', type=float, default=0.1)
    p.add_argument('--k', type=int, default=8)
    p.add_argument('--T', type=float, default=10.0)
    args = p.parse_args()

    L = args.L
    mid = L//2
    # load Vlow0 from previous braid sim results
    data = np.load(Path('results') / f'braid_sim_h{args.h:.3f}_L{L}.npz')
    Vlow = data['Vlow0']

    # build small full ops and project
    full_ops_named = build_bdg_ops_small(L, args.t, args.Delta, mid)
    names = [n for n, _ in full_ops_named]
    full_ops = [M for _, M in full_ops_named]
    basis_proj = [Vlow.conj().T @ M @ Vlow for M in full_ops]

    B2 = ideal_braid_2x2()
    dt = args.T / float(args.k)

    print('Optimizing in projected 2x2 space:', len(basis_proj), 'ops, sequence length k=', args.k)
    res = run_projected_optimization(Vlow, basis_proj, B2, args.k, dt)
    print('Optimization done. success=', res.success, 'message=', res.message)

    coeffs = res.x.reshape((args.k, len(basis_proj)))
    frob, fid_like, Ueff = lift_and_eval(coeffs, full_ops, Vlow, B2, dt)

    out = Path('results')
    out.mkdir(exist_ok=True)
    np.savez(out / f'variational_compile_h{args.h:.3f}_L{L}.npz', coeffs=coeffs, frob=frob, fid_like=fid_like, Ueff=Ueff)
    print('Saved variational compile results to', out / f'variational_compile_h{args.h:.3f}_L{L}.npz')
    print('lifted frob_diff=', frob, ' fid_like=', fid_like)


if __name__ == '__main__':
    main()
