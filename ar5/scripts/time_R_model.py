#!/usr/bin/env python3
"""Build a time-dependent Hamiltonian H(t) from R(u) by converting
the R-matrix generator into a Schrödinger Hamiltonian.

Usage:
 - Import functions: build_time_H_sequence(u0,u1,eta,T,n_steps)
 - The script returns an array of times and a list of 4x4 Hermitian H(t)
   matrices (local two-site operators) ready to embed into larger models.

Conversion logic (kept explicit):
 - R(u) is the XXZ R-matrix from `xxz_R_and_H.R_xxz`.
 - Local generator (canonical) is computed as H_local(u) = -i R(u)^{-1} dR/du.
 - For a schedule u(t) with du/dt = alpha, the Schrödinger Hamiltonian
   that produces the discrete update U = R(u+du) R(u)^{-1} via
   U ≈ exp(-i H_schr dt) is H_schr(t) = - (du/dt) * H_local(u).

This file provides utility functions and a small CLI example that
writes `results/time_H_sequence.npz` when run.
"""
import argparse
import importlib.util
import numpy as np
import os


def R_xxz(u, eta):
    """Return the 4x4 trigonometric XXZ R-matrix.

    Standard six-vertex parametrization:
      a(u) = sin(u + eta)
      b(u) = sin(u)
      c = sin(eta)

    R = [[a,0,0,0],[0,b,c,0],[0,c,b,0],[0,0,0,a]]
    """
    a = np.sin(u + eta)
    b = np.sin(u)
    c = np.sin(eta)
    R = np.array([
        [a, 0.0, 0.0, 0.0],
        [0.0, b, c, 0.0],
        [0.0, c, b, 0.0],
        [0.0, 0.0, 0.0, a],
    ], dtype=complex)
    return R


def dRdu(u, eta):
    """Derivative dR/du for the above parametrization.

    a' = cos(u + eta), b' = cos(u), c' = 0
    """
    ap = np.cos(u + eta)
    bp = np.cos(u)
    cp = 0.0
    dR = np.array([
        [ap, 0.0, 0.0, 0.0],
        [0.0, bp, cp, 0.0],
        [0.0, cp, bp, 0.0],
        [0.0, 0.0, 0.0, ap],
    ], dtype=complex)
    return dR



def compute_H_local(u, eta):
    """Compute local generator H_local(u) = -i R^{-1} dR/du (4x4).

    Returns a 4x4 complex Hermitian matrix (within numerical error).
    """
    R = R_xxz(u, eta)
    dR = dRdu(u, eta)
    Rinv = np.linalg.inv(R)
    H_local = -1j * (Rinv @ dR)
    # symmetrize to reduce numerical non-hermiticity
    H_local = 0.5 * (H_local + H_local.conj().T)
    return H_local


def build_time_H_sequence(u0, u1, eta=0.6, T=1.0, n_steps=200):
    """Build time array and list of H_schr matrices.

    Parameters:
      u0,u1: start and end values of spectral parameter u
      eta: R-matrix parameter
      T: total physical time (units arbitrary)
      n_steps: number of time steps

    Returns: times (n_steps+1,), H_list list of (n_steps+1) arrays 4x4
    """
    times = np.linspace(0.0, T, n_steps + 1)
    u_vals = np.linspace(u0, u1, n_steps + 1)
    du_dt = (u1 - u0) / float(T)

    H_list = []
    for u in u_vals:
        H_local = compute_H_local(u, eta)
        # Schrödinger Hamiltonian for time evolution: H_schr = - (du/dt) * H_local
        H_schr = -du_dt * H_local
        # enforce Hermiticity
        H_schr = 0.5 * (H_schr + H_schr.conj().T)
        H_list.append(H_schr)

    return times, H_list


def save_sequence(path, times, H_list):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    arrs = {f'H{i}': H_list[i] for i in range(len(H_list))}
    np.savez(path, times=times, **arrs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build time-dependent H(t) from R(u).')
    parser.add_argument('--u0', type=float, default=-0.5)
    parser.add_argument('--u1', type=float, default=0.5)
    parser.add_argument('--eta', type=float, default=0.6)
    parser.add_argument('--T', type=float, default=1.0)
    parser.add_argument('--n_steps', type=int, default=200)
    parser.add_argument('--out', type=str, default='results/time_H_sequence.npz')
    parser.add_argument('--Rmodule', type=str, default=None,
                        help='Path to a Python file exposing a function R(u) returning 4x4 array')
    parser.add_argument('--eps', type=float, default=1e-6,
                        help='Finite-difference epsilon for numerical derivative')
    parser.add_argument('--mode', type=str, default='Rinv', choices=['Rinv'],
                        help="Generator mode: only 'Rinv' (canonical) supported: -i R^{-1} dR/du")
    parser.add_argument('--singular_tol', type=float, default=1e-8,
                        help='Threshold for |det R(u)| below which u is treated as singular')
    parser.add_argument('--u_shift', type=float, default=1e-6,
                        help='Small shift to apply to u when encountering a singular point')
    args = parser.parse_args()

    # load user-supplied R(u) if provided
    Rfunc = None
    if args.Rmodule:
        spec = importlib.util.spec_from_file_location('user_R_module', args.Rmodule)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        if not hasattr(mod, 'R'):
            raise SystemExit('Provided module must define R(u) -> 4x4 array')
        Rfunc = mod.R

    def R_wrapper(u):
        if Rfunc is not None:
            return np.array(Rfunc(u), dtype=complex)
        else:
            return R_xxz(u, args.eta)

    def dR_numeric(u, eps=args.eps):
        return (R_wrapper(u + eps) - R_wrapper(u - eps)) / (2.0 * eps)

    # choose derivative function (analytic for xxz, numeric for module or by default)
    def dR_choose(u):
        if args.Rmodule is None:
            return dRdu(u, args.eta)
        return dR_numeric(u)

    # rebind compute_H_local to use chosen derivative
    def compute_H_local_chosen(u):
        # Strict (canonical) generator: H_local = -i R^{-1} dR/du
        dR = dR_choose(u)
        R = R_wrapper(u)
        # check invertibility
        detR = np.linalg.det(R)
        if abs(detR) < 1e-12:
            raise SystemExit(f'R(u) is singular or nearly singular at u={u}; cannot compute canonical generator')
        Rinv = np.linalg.inv(R)
        H_local = -1j * (Rinv @ dR)
        H_local = 0.5 * (H_local + H_local.conj().T)
        return H_local

    def build_time_H_sequence_chosen(u0, u1, T=1.0, n_steps=200,
                                    singular_tol=1e-8, u_shift=1e-6):
        times = np.linspace(0.0, T, n_steps + 1)
        u_vals = np.linspace(u0, u1, n_steps + 1)
        du_dt = (u1 - u0) / float(T)

        H_list = []
        shifted_points = []
        for i, u in enumerate(u_vals):
            # check invertibility; if near-singular, shift u slightly
            try:
                Rm = R_wrapper(u)
                detR = np.linalg.det(Rm)
            except Exception:
                detR = 0.0
            u_eff = u
            if abs(detR) <= singular_tol:
                u_eff = u + u_shift
                shifted_points.append((i, u, u_eff))
            H_local = compute_H_local_chosen(u_eff)
            H_schr = -du_dt * H_local
            H_schr = 0.5 * (H_schr + H_schr.conj().T)
            H_list.append(H_schr)

        if shifted_points:
            print('Warning: encountered near-singular R(u) at these indices, applied small shifts:')
            for idx, u_old, u_new in shifted_points:
                print(f'  index {idx}: u={u_old} -> u_eff={u_new}')

        return times, H_list

    # Build time sequence using canonical generator by default
    times, Hs = build_time_H_sequence_chosen(args.u0, args.u1, T=args.T, n_steps=args.n_steps)
    save_sequence(args.out, times, Hs)
    print('Saved time-dependent H sequence to', args.out)
