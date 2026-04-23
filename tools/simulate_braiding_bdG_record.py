#!/usr/bin/env python3
"""Run the braiding-like BdG path but record per-step diagnostics (instantaneous gap, LDOS at mid, Ueff evolution).

Saves a .npz containing:
 - U_final, Ueff_final, Vlow0
 - times, gaps, ldos_mid (per step), Ueff_list (per step small 2x2 projected)
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
    s = t / (T / 3.0)
    step = int(np.floor(s))
    frac = s - step
    if step == 0:
        g1 = frac
        g3 = 1.0 - frac
        g = [g1, 0.0, g3, 1.0]
    elif step == 1:
        g2 = frac
        g1 = 1.0 - frac
        g = [g1, g2, 0.0, 1.0]
    else:
        g3 = frac
        g2 = 1.0 - frac
        g = [0.0, g2, g3, 1.0]
    return g


def project_lowest_two(H, num=2):
    evals, evecs = eigh(H)
    idx = np.argsort(np.abs(evals))
    return evecs[:, idx[:num]], evals[idx[:num]]


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
    p.add_argument('--T', type=float, default=600.0)
    p.add_argument('--steps', type=int, default=600)
    p.add_argument('--save_every', type=int, default=10)
    p.add_argument('--mode', choices=['gate','baxter'], default='gate',
                   help="Which evolution to run: 'gate' (bond extras) or 'baxter' (spectral-parameter H(t))")
    p.add_argument('--local_slow_center', type=float, default=None,
                   help='Center time for local slowdown (same units as T)')
    p.add_argument('--local_slow_width', type=float, default=None,
                   help='Width (sigma) of local slowdown Gaussian')
    p.add_argument('--local_slow_factor', type=float, default=None,
                   help='Max slowdown factor (>1) at center (multiplier of time denominator)')
    args = p.parse_args()

    L = args.L
    mid = L//2
    mu0 = 0.0
    mu_vec = np.zeros(L) + mu0
    mu_vec[mid] += 2.0 * args.h

    dt = args.T / float(args.steps)
    dim = 2*L
    U = np.eye(dim, dtype=complex)

    bond_indices = [mid-2, mid-1, mid]
    tc = 1.0

    times = []
    gaps = []
    ldos_mid = []
    Ueff_list = []

    # initial projector
    H0 = build_bdg(L, args.t, args.Delta, mu_vec)
    Vlow0, evals_low0 = project_lowest_two(H0)

    # prepare baxter path functions if requested
    def default_u_of_t(t, T):
        u0 = np.pi/2.0
        return u0 * (0.5 - 0.5 * np.cos(np.pi * t / T))

    def du_dt(t, T):
        u0 = np.pi/2.0
        return u0 * 0.5 * (np.pi / T) * np.sin(np.pi * t / T) if T != 0 else 0.0

    # local slowdown parameters: multiply du/dt by multiplier(t) in baxter mode
    p = args
    slow_center = getattr(p, 'local_slow_center', None)
    slow_width = getattr(p, 'local_slow_width', None)
    slow_factor = getattr(p, 'local_slow_factor', None)
    def slowdown_multiplier(t):
        if slow_center is None or slow_width is None or slow_factor is None:
            return 1.0
        # gaussian-shaped slowdown: multiplier in (1/slow_factor .. 1)
        sigma = slow_width
        x = (t - slow_center) / sigma
        # at center multiplier = 1/slow_factor, far away -> 1.0
        return 1.0 / (1.0 + (slow_factor - 1.0) * np.exp(-0.5 * x * x))

    for k in range(args.steps):
        tnow = (k + 0.5) * dt
        if args.mode == 'gate':
            g = gate_functions(tnow, args.T)
            delta = [tc * (g[0]*g[3]), tc * (g[1]*g[3]), tc * (g[2]*g[3])]
            bond_extras = {bond_indices[i]: delta[i] for i in range(len(bond_indices))}
            Hbdg = build_bdg(L, args.t, args.Delta, mu_vec, bond_extras)
        else:
            # Baxterization-like H(t) = dg/du * du/dt * H0, here g(u)=u so dg/du=1
            u = default_u_of_t(tnow, args.T)
            scale = du_dt(tnow, args.T)
            # apply local slowdown multiplier
            scale *= slowdown_multiplier(tnow)
            Hbdg = scale * H0

        # instantaneous diagnostics
        evals, evecs = eigh(Hbdg)
        evals_sorted = np.sort(np.abs(evals))
        gap = evals_sorted[2] - evals_sorted[1] if len(evals_sorted)>2 else evals_sorted[0]
        gaps.append(gap)

        # LDOS at mid site: use particle component amplitude of lowest-energy eigenstate
        idx = np.argsort(np.abs(evals))
        psi0 = evecs[:, idx[0]]
        # particle block is first L entries
        ldos = np.abs(psi0[mid])**2
        ldos_mid.append(ldos)

        # evolve
        U = expm(-1j * Hbdg * dt) @ U

        # project evolving U to initial low-energy basis every save_every
        if (k % args.save_every) == 0:
            Ueff = Vlow0.conj().T @ U @ Vlow0
            Ueff_list.append(Ueff)
            times.append(tnow)

    # final evaluation
    Ueff_final = Vlow0.conj().T @ U @ Vlow0
    B2 = ideal_braid_2x2()
    frob = np.linalg.norm(Ueff_final - B2)
    fid_like = np.real(np.trace(Ueff_final.conj().T @ B2)) / 2.0

    out = Path('results')
    out.mkdir(exist_ok=True)
    np.savez(out / f'braid_sim_record_h{args.h:.3f}_L{L}.npz', U_final=U, Ueff_final=Ueff_final,
             Vlow0=Vlow0, times=np.array(times), gaps=np.array(gaps), ldos_mid=np.array(ldos_mid),
             Ueff_list=np.array(Ueff_list))
    print('Saved results to', out / f'braid_sim_record_h{args.h:.3f}_L{L}.npz')
    print('frob_diff=',frob,' fid_like=',fid_like)


if __name__ == '__main__':
    main()
