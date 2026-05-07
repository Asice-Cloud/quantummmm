#!/usr/bin/env python3
"""Scan alpha (and optional delta list) to collect final SU(2) rotations
on the Bloch sphere produced by the H(u) pipeline and plot coverage.

Saves results to results/rotation_reachability_delta{delta:.3g}.npz and
plots to PNG.
"""
import os
import argparse
import numpy as np
import math
from scipy.linalg import norm

import verify_from_R as vr


def scan_for_delta(delta, N=200, nalpha=40, alphas=None, outdir='results'):
    os.makedirs(outdir, exist_ok=True)
    us, du, H4s = vr.compute_Hs_from_R(delta, N=N)
    H2_list, ds = vr.build_H2_list_from_H4s(H4s)
    dt = du
    if alphas is None:
        alphas = np.logspace(-3, 1, nalpha)
    axes = np.zeros((len(alphas), 3))
    angles = np.zeros(len(alphas))
    U_diffs = np.zeros(len(alphas))
    for i, a in enumerate(alphas):
        Hs = [a * H for H in H2_list]
        U_final, _ = vr.compute_U_from_Hlist(Hs, dt)
        axis, angle = vr.rot_axis_angle_from_U(U_final)
        axes[i] = axis
        angles[i] = angle
    # Compute coverage on sphere grid
    thetas = np.arccos(np.clip(axes[:, 2], -1.0, 1.0))  # 0..pi
    phis = np.arctan2(axes[:, 1], axes[:, 0])  # -pi..pi
    # bin grid
    ntheta = 36
    nphi = 72
    th_edges = np.linspace(0, math.pi, ntheta + 1)
    ph_edges = np.linspace(-math.pi, math.pi, nphi + 1)
    H, _, _ = np.histogram2d(thetas, phis, bins=[th_edges, ph_edges])
    occupancy = (H > 0).sum()
    total_bins = H.size
    coverage_frac = float(occupancy) / float(total_bins)

    # Save results
    out = os.path.join(outdir, f'rotation_reachability_delta{delta:.3g}.npz')
    np.savez(out, delta=delta, us=us, alphas=alphas, axes=axes, angles=angles,
             coverage_frac=coverage_frac, hist=H, th_edges=th_edges, ph_edges=ph_edges)

    # Make plots
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa

        # 3D sphere scatter
        fig = plt.figure(figsize=(5, 5), dpi=150)
        ax = fig.add_subplot(111, projection='3d')
        u = np.linspace(0, 2 * np.pi, 60)
        v = np.linspace(0, np.pi, 30)
        x = np.outer(np.cos(u), np.sin(v))
        y = np.outer(np.sin(u), np.sin(v))
        z = np.outer(np.ones_like(u), np.cos(v))
        ax.plot_surface(x, y, z, color='lightgray', alpha=0.2, linewidth=0)
        sc = ax.scatter(axes[:, 0], axes[:, 1], axes[:, 2], c=angles, cmap='viridis', s=20)
        ax.set_title(f'Rotation axes (delta={delta:.3g})')
        plt.colorbar(sc, ax=ax, fraction=0.02, pad=0.04, label='angle (rad)')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f'rotation_axes_delta{delta:.3g}.png'))
        plt.close(fig)

        # 2D projection (phi vs cos theta)
        fig, ax = plt.subplots(1, 2, figsize=(10, 4), dpi=150)
        sc2 = ax[0].scatter(phis, np.cos(thetas), c=angles, cmap='viridis', s=16)
        ax[0].set_xlabel('phi'); ax[0].set_ylabel('cos(theta)'); ax[0].set_title('phi vs cos(theta)')
        plt.colorbar(sc2, ax=ax[0], fraction=0.04, pad=0.02, label='angle (rad)')

        # histogram of angles
        ax[1].hist(angles, bins=36)
        ax[1].set_xlabel('angle (rad)'); ax[1].set_title('Histogram of rotation angles')
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f'rotation_angle_hist_delta{delta:.3g}.png'))
        plt.close(fig)
    except Exception as e:
        print('Plotting failed:', e)

    return out, coverage_frac


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--deltas', type=float, nargs='+', default=[0.015])
    p.add_argument('--N', type=int, default=200)
    p.add_argument('--nalpha', type=int, default=40)
    p.add_argument('--outdir', type=str, default='results')
    args = p.parse_args()

    results = {}
    for dlt in args.deltas:
        print('Scanning delta=', dlt)
        out, cov = scan_for_delta(dlt, N=args.N, nalpha=args.nalpha, outdir=args.outdir)
        print('Saved', out, ' coverage=', cov)
        results[dlt] = {'file': out, 'coverage': cov}
    print('All done.')


if __name__ == '__main__':
    main()
