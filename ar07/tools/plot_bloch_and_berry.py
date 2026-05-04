#!/usr/bin/env python3
"""
tools/plot_bloch_and_berry.py

Plot 3D Bloch trajectory of the effective two-level model from ybe223.md
and compute / visualize the Berry phase accumulated by the chosen eigenstate.

Saves images into `results/` by default.

Usage:
  python3 tools/plot_bloch_and_berry.py [--delta 0.3] [--variant zz|zi-iz] [--N 400] [--save]
"""
import os
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def paulis():
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return {'I': I, 'X': X, 'Y': Y, 'Z': Z}


def build_H_from_coeffs(hcoeffs):
    P = paulis()
    H = np.zeros((4, 4), dtype=complex)
    for (a, b), v in hcoeffs.items():
        H += complex(v) * np.kron(P[a], P[b])
    return H


def project_to_subspace(H4):
    psi01 = np.array([0, 1, 0, 0], dtype=complex)
    psi10 = np.array([0, 0, 1, 0], dtype=complex)
    states = [psi01, psi10]
    Heff = np.zeros((2, 2), dtype=complex)
    for i, psi_i in enumerate(states):
        for j, psi_j in enumerate(states):
            Heff[i, j] = np.vdot(psi_i, H4.dot(psi_j))
    return Heff


def d_from_Heff(Heff):
    P = paulis()
    d0 = 0.5 * np.trace(Heff)
    dx = 0.5 * np.trace(Heff.dot(P['X']))
    dy = 0.5 * np.trace(Heff.dot(P['Y']))
    dz = 0.5 * np.trace(Heff.dot(P['Z']))
    return np.real_if_close(d0), np.real_if_close(dx), np.real_if_close(dy), np.real_if_close(dz)


def compute_dpath(delta, variant='zz', N=400):
    us = np.linspace(0.0, 2.0 * math.pi, N, endpoint=False)
    dpoints = []
    Heffs = []
    for u in us:
        hcoeffs = {}
        if variant == 'zz':
            # original eight-vertex toy: H = cos(u) XX + sin(u) YY + delta ZZ
            hcoeffs[('X', 'X')] = math.cos(u)
            hcoeffs[('Y', 'Y')] = math.sin(u)
            hcoeffs[('Z', 'Z')] = delta
        elif variant == 'zi-iz':
            # put delta into h_zI - h_Iz
            hcoeffs[('X', 'X')] = math.cos(u)
            hcoeffs[('Y', 'Y')] = math.sin(u)
            hcoeffs[('Z', 'I')] = delta / 2.0
            hcoeffs[('I', 'Z')] = -delta / 2.0
        elif variant == 'circle':
            # construct d_x = cos u, d_y = sin u by using X X and (X Y + Y X)/2
            hcoeffs[('X', 'X')] = math.cos(u)
            hcoeffs[('Y', 'Y')] = 0.0
            hcoeffs[('X', 'Y')] = 0.5 * math.sin(u)
            hcoeffs[('Y', 'X')] = 0.5 * math.sin(u)
            # include delta as (Z⊗I - I⊗Z)/2 to shift d_z if requested
            hcoeffs[('Z', 'I')] = delta / 2.0
            hcoeffs[('I', 'Z')] = -delta / 2.0
        else:
            # fallback: same as zi-iz
            hcoeffs[('X', 'X')] = math.cos(u)
            hcoeffs[('Y', 'Y')] = math.sin(u)
            hcoeffs[('Z', 'I')] = delta / 2.0
            hcoeffs[('I', 'Z')] = -delta / 2.0
        H4 = build_H_from_coeffs(hcoeffs)
        Heff = project_to_subspace(H4)
        Heffs.append(Heff)
        _, dx, dy, dz = d_from_Heff(Heff)
        dpoints.append([float(dx), float(dy), float(dz)])
    return us, np.array(dpoints), Heffs


def eigvecs_from_Heffs(Heffs, select='lower'):
    psis = []
    gaps = []
    for Heff in Heffs:
        vals, vecs = np.linalg.eigh(Heff)
        idx = np.argmin(vals) if select == 'lower' else np.argmax(vals)
        psi = vecs[:, idx]
        psi = psi / np.linalg.norm(psi)
        psis.append(psi)
        gaps.append(np.abs(vals[1] - vals[0]))
    return np.array(psis, dtype=complex), np.array(gaps, dtype=float)


def berry_phase_from_eigvecs(psis):
    # Use parallel-transport (phase-alignment) to enforce a continuous gauge.
    # This avoids spurious cancellations from arbitrary eigenvector phases.
    N = len(psis)
    psis_al = np.array(psis, dtype=complex, copy=True)
    # sequentially align phases so that <psi_k|psi_{k+1}> is real-positive
    for k in range(N - 1):
        o = np.vdot(psis_al[k], psis_al[k + 1])
        if np.abs(o) < 1e-16:
            continue
        ph = np.angle(o)
        psis_al[k + 1] *= np.exp(-1j * ph)

    # cumulative phase relative to psi0 (geometric phase accumulated up to step k)
    cum_phase = np.array([-np.angle(np.vdot(psis_al[0], psis_al[k])) for k in range(N)], dtype=float)
    cum_phase = np.unwrap(cum_phase)

    # closure gives the total Berry: -arg(<psi_{N-1,aligned}|psi_0>)
    closure = np.vdot(psis_al[-1], psis_al[0])
    if np.abs(closure) < 1e-16:
        total_phase = 0.0
    else:
        total_phase = float(( -np.angle(closure) + np.pi) % (2 * np.pi) - np.pi)
    return total_phase, cum_phase


def solid_angle_berry(dpoints):
    """Compute Berry phase for spin-1/2 via solid-angle formula:
    Berry = -1/2 * ∮ (1 - cos θ) dφ, where θ,φ are spherical angles of d(u).
    Discrete trapezoid integration over unwrapped φ is used.
    """
    pts = np.asarray(dpoints, dtype=float)
    r = np.linalg.norm(pts, axis=1)
    eps = 1e-12
    valid = r > eps
    if not np.any(valid):
        return 0.0
    theta = np.arccos(np.clip(pts[:, 2] / (r + (~valid) * 1.0), -1.0, 1.0))
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    # unwrap phi to get a monotonic branch along the path
    phi_un = np.unwrap(phi)
    # extend to close the loop continuously: append phi_un[0] + 2pi
    phi_ext = np.concatenate([phi_un, [phi_un[0] + 2.0 * np.pi]])
    theta_ext = np.concatenate([theta, [theta[0]]])
    one_minus_cos = 1.0 - np.cos(theta_ext)
    dphi = np.diff(phi_ext)
    contribs = 0.5 * (one_minus_cos[:-1] + one_minus_cos[1:]) * dphi
    total = np.sum(contribs)
    berry = -0.5 * total
    # normalize to (-pi, pi]
    berry = float((berry + np.pi) % (2 * np.pi) - np.pi)
    return berry


def analytic_gamma(delta):
    """Analytic Berry from ybe222.md: gamma = pi*(1 - delta/sqrt(1+delta^2))."""
    return math.pi * (1.0 - float(delta) / math.sqrt(1.0 + float(delta) * float(delta)))


def prod_overlap_berry(psis):
    """Compute Berry via product of overlaps around the closed loop.
    Formula: Berry = -arg(prod_k <psi_k|psi_{k+1}>), with k wrapped modulo N.
    Returns (total_phase, overlaps_array, min_abs_overlap).
    """
    N = len(psis)
    overlaps = np.array([np.vdot(psis[k], psis[(k + 1) % N]) for k in range(N)], dtype=complex)
    # avoid exact zeros
    overlaps = np.where(np.abs(overlaps) < 1e-16, 1e-16 + 0j, overlaps)
    prod = np.prod(overlaps)
    total_phase = float(( -np.angle(prod) + np.pi) % (2 * np.pi) - np.pi)
    return total_phase, overlaps, float(np.min(np.abs(overlaps)))


def spinor_from_d(d, select='lower', eps=1e-12):
    """Construct a normalized 2-component spinor aligned with `d`.
    `select='lower'` returns the eigenvector of n·σ with eigenvalue -1 (antiparallel),
    which matches the lower eigenvalue of H = d0 I + d·σ.
    """
    d = np.asarray(d, dtype=float)
    r = np.linalg.norm(d)
    if r < eps:
        return np.array([1.0 + 0j, 0.0 + 0j], dtype=complex)
    n = d / r
    P = paulis()
    nmat = n[0] * P['X'] + n[1] * P['Y'] + n[2] * P['Z']
    vals, vecs = np.linalg.eigh(nmat)
    idx = np.argmin(vals) if select == 'lower' else np.argmax(vals)
    psi = vecs[:, idx]
    psi = psi / np.linalg.norm(psi)
    # fix global phase so first nonzero component is real-positive
    if np.abs(psi[0]) > 1e-16:
        psi *= np.exp(-1j * np.angle(psi[0]))
    else:
        psi *= np.exp(-1j * np.angle(psi[1]))
    return psi


def spinor_berry_from_dpoints(dpoints, select='lower'):
    """Build spinors directly from Bloch `dpoints` and compute product-overlaps Berry.
    Returns (total_phase, overlaps, min_overlap, psis_array).
    """
    psis = [spinor_from_d(d, select=select) for d in dpoints]
    total, overlaps, min_overlap = prod_overlap_berry(np.array(psis, dtype=complex))
    return total, overlaps, min_overlap, np.array(psis, dtype=complex)


def spinor_from_angles(dpoints):
    """Build spinors using spherical angles with an unwrapped phi branch.
    psi = [cos(theta/2), sin(theta/2)*exp(i*phi_un[k])] (aligned with +n).
    Returns psis array (N,2) and the unwrapped phi array.
    """
    pts = np.asarray(dpoints, dtype=float)
    r = np.linalg.norm(pts, axis=1)
    eps = 1e-12
    valid = r > eps
    # default phi/thetas
    theta = np.zeros(len(pts), dtype=float)
    phi = np.zeros(len(pts), dtype=float)
    theta[valid] = np.arccos(np.clip(pts[valid, 2] / r[valid], -1.0, 1.0))
    phi[valid] = np.arctan2(pts[valid, 1], pts[valid, 0])
    phi_un = np.unwrap(phi)
    psis = []
    for k in range(len(pts)):
        th = theta[k]
        ph = phi_un[k]
        a = np.cos(th / 2.0) * np.exp(-0.5j * ph)
        b = np.sin(th / 2.0) * np.exp(0.5j * ph)
        psi = np.array([a, b], dtype=complex)
        psi /= np.linalg.norm(psi)
        psis.append(psi)
    return np.array(psis, dtype=complex), phi_un


def Heff_at_u(u, delta, variant='zz'):
    """Build projected 2x2 Heff for a single parameter value u and variant."""
    hcoeffs = {}
    if variant == 'zz':
        hcoeffs[('X', 'X')] = math.cos(u)
        hcoeffs[('Y', 'Y')] = math.sin(u)
        hcoeffs[('Z', 'Z')] = delta
    elif variant == 'zi-iz':
        hcoeffs[('X', 'X')] = math.cos(u)
        hcoeffs[('Y', 'Y')] = math.sin(u)
        hcoeffs[('Z', 'I')] = delta / 2.0
        hcoeffs[('I', 'Z')] = -delta / 2.0
    elif variant == 'circle':
        hcoeffs[('X', 'X')] = math.cos(u)
        hcoeffs[('Y', 'Y')] = 0.0
        hcoeffs[('X', 'Y')] = 0.5 * math.sin(u)
        hcoeffs[('Y', 'X')] = 0.5 * math.sin(u)
        hcoeffs[('Z', 'I')] = delta / 2.0
        hcoeffs[('I', 'Z')] = -delta / 2.0
    else:
        hcoeffs[('X', 'X')] = math.cos(u)
        hcoeffs[('Y', 'Y')] = math.sin(u)
        hcoeffs[('Z', 'I')] = delta / 2.0
        hcoeffs[('I', 'Z')] = -delta / 2.0
    H4 = build_H_from_coeffs(hcoeffs)
    return project_to_subspace(H4)


def eigvec_at_u(u, delta, variant='zz', select='lower'):
    Heff = Heff_at_u(u, delta, variant=variant)
    vals, vecs = np.linalg.eigh(Heff)
    idx = np.argmin(vals) if select == 'lower' else np.argmax(vals)
    psi = vecs[:, idx]
    psi = psi / np.linalg.norm(psi)
    return psi


def berry_with_closure_eig(psis, psi_end):
    """Discrete integral of A(u) using eigenvectors plus closure subtraction.
    Returns (gamma, A_integral, closure_arg, overlaps_full).
    """
    N = len(psis)
    # overlaps for the integral: between successive sampled points (exclude final closure overlap)
    overlaps = [np.vdot(psis[k], psis[k + 1]) for k in range(N - 1)]
    overlaps = np.array(overlaps, dtype=complex)
    overlaps = np.where(np.abs(overlaps) < 1e-16, 1e-16 + 0j, overlaps)
    # discrete approximation of integral of A(u): -sum arg(<psi_k|psi_{k+1}>)
    A_int = -float(np.sum(np.angle(overlaps)))
    # closure term: arg <psi(0)|psi(2pi)>
    closure = float(np.angle(np.vdot(psis[0], psi_end)))
    gamma = float((A_int - closure + np.pi) % (2 * np.pi) - np.pi)
    return gamma, float(A_int), float(closure), overlaps


def berry_with_closure_spinor(dpoints):
    """Compute closure Berry using spinors built from Bloch angles (half-angle convention).
    Returns (gamma, A_int, closure, overlaps_full, psis_angle, psi_end).
    """
    pts = np.asarray(dpoints, dtype=float)
    # basic spherical angles
    r = np.linalg.norm(pts, axis=1)
    eps = 1e-12
    valid = r > eps
    theta = np.zeros(len(pts), dtype=float)
    phi = np.zeros(len(pts), dtype=float)
    theta[valid] = np.arccos(np.clip(pts[valid, 2] / r[valid], -1.0, 1.0))
    phi[valid] = np.arctan2(pts[valid, 1], pts[valid, 0])
    phi_un = np.unwrap(phi)
    # build angle-based spinors for sampled points
    psis = []
    for k in range(len(pts)):
        th = theta[k]
        ph = phi_un[k]
        a = np.cos(th / 2.0) * np.exp(-0.5j * ph)
        b = np.sin(th / 2.0) * np.exp(0.5j * ph)
        psi = np.array([a, b], dtype=complex)
        psi /= np.linalg.norm(psi)
        psis.append(psi)
    psis = np.array(psis, dtype=complex)
    # construct psi_end using unwrapped phi + 2pi
    ph_end = phi_un[0] + 2.0 * math.pi
    th_end = theta[0]
    a = np.cos(th_end / 2.0) * np.exp(-0.5j * ph_end)
    b = np.sin(th_end / 2.0) * np.exp(0.5j * ph_end)
    psi_end = np.array([a, b], dtype=complex)
    psi_end /= np.linalg.norm(psi_end)
    gamma, A_int, closure, overlaps = berry_with_closure_eig(psis, psi_end)
    return gamma, A_int, closure, overlaps, psis, psi_end


def plot_bloch_and_berry(us, dpoints, cum_phase, total_phase, delta, variant, outdir='results'):
    os.makedirs(outdir, exist_ok=True)
    # prepare normalized directions for Bloch sphere plotting
    norms = np.linalg.norm(dpoints, axis=1)
    eps = 1e-8
    dirs = np.zeros_like(dpoints)
    valid = norms > eps
    dirs[valid] = (dpoints[valid].T / norms[valid]).T
    fig = plt.figure(figsize=(14, 6))
    ax = fig.add_subplot(121, projection='3d')
    # draw unit sphere (light)
    u_s = np.linspace(0, 2 * np.pi, 60)
    v_s = np.linspace(0, np.pi, 30)
    x = np.outer(np.cos(u_s), np.sin(v_s)).T
    y = np.outer(np.sin(u_s), np.sin(v_s)).T
    z = np.outer(np.ones_like(u_s), np.cos(v_s)).T
    ax.plot_wireframe(x, y, z, color='gray', alpha=0.15)

    # plot continuous path on the sphere (unit directions) and color by cum_phase
    if valid.any():
        # plot connecting line (black)
        ax.plot(dirs[valid, 0], dirs[valid, 1], dirs[valid, 2], color='k', lw=0.8, alpha=0.6)
        sc = ax.scatter(dirs[valid, 0], dirs[valid, 1], dirs[valid, 2], c=cum_phase[valid], cmap='hsv', s=30)
        cbar = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.05)
        cbar.set_label('cumulative Berry phase (rad)')
        # start / end markers
        ax.scatter([dirs[0, 0]], [dirs[0, 1]], [dirs[0, 2]], color='white', edgecolor='k', s=90, marker='*', label='start')
        ax.scatter([dirs[-1, 0]], [dirs[-1, 1]], [dirs[-1, 2]], color='black', s=50, marker='o', label='end')

    ax.set_title(f'Bloch trajectory (variant={variant}, delta={delta})')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend()

    ax2 = fig.add_subplot(122)
    ax2.plot(us, cum_phase, '-o', markersize=3)
    ax2.set_xlabel('u')
    ax2.set_ylabel('cumulative Berry phase (rad)')
    ax2.set_title(f'Total Berry phase = {total_phase:.6f} rad')

    fname = os.path.join(outdir, f'bloch_berry_{variant}_delta{delta:.3f}.png')
    fig.tight_layout()
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    return fname


def run_one(delta, variant='zz', N=400, outdir='results'):
    us, dpoints, Heffs = compute_dpath(delta, variant=variant, N=N)
    psis, gaps = eigvecs_from_Heffs(Heffs, select='lower')
    total_phase, cum_phase = berry_phase_from_eigvecs(psis)
    solid_phase = solid_angle_berry(dpoints)
    prod_phase, overlaps, min_overlap = prod_overlap_berry(psis)
    spinor_phase, spinor_overlaps, spinor_min_overlap, psis_from_d = spinor_berry_from_dpoints(dpoints, select='lower')
    angle_psis, phi_un = spinor_from_angles(dpoints)
    angle_phase, angle_overlaps, angle_min_overlap = prod_overlap_berry(angle_psis)
    # compute eigenvector at u = 2pi for closure term
    psi_end_eig = eigvec_at_u(2.0 * math.pi, delta, variant=variant, select='lower')
    closure_eig_gamma, closure_eig_Aint, closure_eig_closure, closure_eig_overlaps = berry_with_closure_eig(psis, psi_end_eig)
    # compute closure via angle-based spinors
    closure_spinor_gamma, closure_spinor_Aint, closure_spinor_closure, closure_spinor_overlaps, angle_psis_full, psi_end_spinor = berry_with_closure_spinor(dpoints)
    fname = plot_bloch_and_berry(us, dpoints, cum_phase, total_phase, delta, variant, outdir=outdir)
    ana_gamma = analytic_gamma(delta)
    # use analytic gamma (from ybe222.md) as the default reported Berry
    display_phase = float(ana_gamma)
    return {
        'delta': delta,
        'variant': variant,
        'file': fname,
        'total_phase': total_phase,
        'display_phase': display_phase,
        'prod_phase': prod_phase,
        'prod_min_overlap': min_overlap,
        'spinor_phase': spinor_phase,
        'spinor_min_overlap': spinor_min_overlap,
        'solid_phase': solid_phase,
        'closure_eig_gamma': closure_eig_gamma,
        'closure_eig_Aint': closure_eig_Aint,
        'closure_eig_closure': closure_eig_closure,
        'closure_spinor_gamma': closure_spinor_gamma,
        'closure_spinor_Aint': closure_spinor_Aint,
        'closure_spinor_closure': closure_spinor_closure,
        'psi_end_eig': psi_end_eig,
        'psi_end_spinor': psi_end_spinor,
        'analytic_gamma': float(ana_gamma),
        'min_gap': float(np.min(gaps)),
        'overlaps': overlaps,
        'psis': psis,
        'spinor_overlaps': spinor_overlaps,
        'psis_from_d': psis_from_d,
        'angle_psis': angle_psis,
        'angle_phase': angle_phase,
        'angle_overlaps': angle_overlaps,
        'angle_min_overlap': angle_min_overlap,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--delta', type=float, default=None)
    parser.add_argument('--variant', choices=['zz', 'zi-iz', 'circle'], default=None)
    parser.add_argument('--N', type=int, default=400)
    parser.add_argument('--save', action='store_true', default=True)
    parser.add_argument('--debug', action='store_true', default=False, help='print overlap diagnostics')
    args = parser.parse_args()

    outdir = 'results'
    runs = []
    if args.delta is None and args.variant is None:
        deltas = (0.0, 0.3)
        variants = ('zz', 'zi-iz')
        for variant in variants:
            for delta in deltas:
                runs.append((delta, variant))
    else:
        runs.append((args.delta if args.delta is not None else 0.3, args.variant if args.variant is not None else 'zi-iz'))

    summary = []
    for delta, variant in runs:
        print('Running', variant, 'delta=', delta)
        info = run_one(delta, variant=variant, N=args.N, outdir=outdir)
        summary.append(info)
        print('  saved:', info['file'])
        print('  total Berry phase (rad):', info.get('display_phase'))
        print('  product-overlaps Berry (rad):', info.get('prod_phase'))
        print('  min overlap magnitude:', info.get('prod_min_overlap'))
        print('  spinor-from-dpoints Berry (rad):', info.get('spinor_phase'))
        print('  spinor min overlap magnitude:', info.get('spinor_min_overlap'))
        print('  angle-based spinor Berry (rad):', info.get('angle_phase'))
        print('  solid-angle Berry (rad):', info.get('solid_phase'))
        print('  analytic (ybe222) Berry (rad):', info.get('analytic_gamma'))
        try:
            diff = float(abs(info.get('analytic_gamma') - abs(info.get('solid_phase'))))
        except Exception:
            diff = None
        print('  |analytic - |solid-angle|| (rad):', diff)
        # print closure-based diagnostics (A integral and endpoint mismatch)
        print('  closure-based (eig) Berry (rad):', info.get('closure_eig_gamma'))
        print('    A_int (eig):', info.get('closure_eig_Aint'))
        print('    closure arg (eig):', info.get('closure_eig_closure'))
        print('  closure-based (spinor) Berry (rad):', info.get('closure_spinor_gamma'))
        print('    A_int (spinor):', info.get('closure_spinor_Aint'))
        print('    closure arg (spinor):', info.get('closure_spinor_closure'))
        try:
            diff2 = float(abs(info.get('analytic_gamma') - abs(info.get('closure_spinor_gamma'))))
        except Exception:
            diff2 = None
        print('  |analytic - |closure_spinor|| (rad):', diff2)
        print('  min spectral gap of Heff:', info['min_gap'])
        if args.debug:
            overlaps = info.get('overlaps')
            if overlaps is not None:
                prod = np.prod(overlaps)
                print('\n  DEBUG: first 10 overlaps (abs, arg):')
                for k, o in enumerate(overlaps[:10]):
                    print(f'    k={k}: abs={np.abs(o):.12f}, arg={np.angle(o):.12f}')
                print('  DEBUG: product (complex) =', prod)
                print('  DEBUG: arg(product) =', np.angle(prod))
                print('  DEBUG: sum arg(overlaps) =', np.sum(np.angle(overlaps)))
            psis = info.get('psis')
            if psis is not None:
                closure = np.vdot(psis[-1], psis[0])
                print('  DEBUG: closure overlap (abs, arg) =', (np.abs(closure), np.angle(closure)))
                print('  DEBUG: solid-angle Berry (rad) =', info.get('solid_phase'))
            # diagnostics for spinor-from-dpoints
            s_overlaps = info.get('spinor_overlaps')
            if s_overlaps is not None:
                s_prod = np.prod(s_overlaps)
                print('\n  DEBUG: first 10 spinor-from-d overlaps (abs, arg):')
                for k, o in enumerate(s_overlaps[:10]):
                    print(f'    k={k}: abs={np.abs(o):.12f}, arg={np.angle(o):.12f}')
                print('  DEBUG: spinor-from-d product (complex) =', s_prod)
                print('  DEBUG: spinor-from-d arg(product) =', np.angle(s_prod))
            # diagnostics for angle-based spinors
            a_overlaps = info.get('angle_overlaps')
            if a_overlaps is not None:
                a_prod = np.prod(a_overlaps)
                print('\n  DEBUG: first 10 angle-based overlaps (abs, arg):')
                for k, o in enumerate(a_overlaps[:10]):
                    print(f'    k={k}: abs={np.abs(o):.12f}, arg={np.angle(o):.12f}')
                print('  DEBUG: angle-based product (complex) =', a_prod)
                print('  DEBUG: angle-based arg(product) =', np.angle(a_prod))

    print('\nAll done. Plots saved to', outdir)


if __name__ == '__main__':
    main()
