#!/usr/bin/env python3
"""Compute fermion mode density rho_j from a BdG eigenvector, fit exponential decay
from the density peak to both edges, and save plot and summary.

Usage:
  python3 tools/plot_and_fit_density.py --vec results/vec_lowest.npy
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os


def fit_exponential(x, y, min_points=5):
    # Fit y = A * exp(-x/xi) where x>=0
    mask = (y > 0) & np.isfinite(y)
    if mask.sum() < min_points:
        return None
    lx = np.log(y[mask])
    xx = x[mask]
    # linear fit lx = logA - x/xi -> slope = -1/xi
    p = np.polyfit(xx, lx, 1)
    slope, intercept = p[0], p[1]
    xi = -1.0 / slope if slope != 0 else np.inf
    A = np.exp(intercept)
    # compute residual RMS on log scale
    pred = intercept + slope * xx
    res = np.sqrt(np.mean((lx - pred) ** 2))
    return {'A': float(A), 'xi': float(xi), 'res_rms': float(res), 'n': int(mask.sum())}


def fit_cosh_grid(x, y, p=1, xi_min=0.5, xi_max=200.0, nsteps=400, min_points=5):
    # Fit y = A * cosh(x/xi)^{-p} by grid-searching xi and solving linear A for each xi.
    if len(x) < min_points:
        return None
    # remove nonpositive
    mask = (y > 0) & np.isfinite(y)
    if mask.sum() < min_points:
        return None
    xg = x[mask]
    yg = y[mask]
    xis = np.linspace(xi_min, xi_max, nsteps)
    best = None
    for xi in xis:
        f = 1.0 / (np.cosh(xg / xi) ** p)
        # linear least squares for A: minimize ||yg - A f||^2 => A = (yg·f)/(f·f)
        denom = np.sum(f * f)
        if denom == 0:
            continue
        A = np.sum(yg * f) / denom
        pred = A * f
        res = np.sqrt(np.mean((yg - pred) ** 2))
        if (best is None) or (res < best['res']):
            best = {'xi': float(xi), 'A': float(A), 'res': float(res), 'n': int(mask.sum())}
    if best is None:
        return None
    best['res_rms'] = best.pop('res')
    best['p'] = p
    return best


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--vec', default='results/vec_lowest.npy')
    p.add_argument('--out', default='results/gamma_density.png')
    p.add_argument('--summary', default='results/gamma_density_summary.txt')
    p.add_argument('--min-weight', type=float, default=1e-8, help='min density to include in fit')
    p.add_argument('--fit', choices=['exp','cosh','both'], default='exp', help='which fit to perform')
    p.add_argument('--cosh-p', type=int, default=1, help='power p in cosh^{-p} fit')
    args = p.parse_args()

    if not os.path.exists(args.vec):
        raise SystemExit(f'vector file not found: {args.vec}')
    vec = np.load(args.vec)
    L = len(vec)
    if L % 2 != 0:
        raise SystemExit('vector length must be even (2N)')
    N = L // 2
    u = vec[:N]
    v = vec[N:]
    rho = np.abs(u)**2 + np.abs(v)**2

    # locate peak (largest rho)
    peak_idx = int(np.argmax(rho))

    # prepare left fit (from peak towards left end)
    xl = np.arange(0, peak_idx+1)[::-1]  # distances: 0,1,2,... from peak to left
    yl = rho[:peak_idx+1][::-1]
    maskl = yl >= args.min_weight
    fit_left = fit_exponential(xl[maskl], yl[maskl]) if maskl.sum() >= 5 else None

    # prepare right fit (from peak towards right end)
    xr = np.arange(0, N-peak_idx)
    yr = rho[peak_idx:]
    maskr = yr >= args.min_weight
    fit_right = fit_exponential(xr[maskr], yr[maskr]) if maskr.sum() >= 5 else None

    # Plot density and fits
    plt.figure(figsize=(8,4))
    x = np.arange(N)
    plt.plot(x, rho, '-', label=r'$\rho_j$')
    plt.semilogy(x, rho, '.')
    if args.fit in ('exp','both') and fit_left is not None:
        xs = np.arange(0, peak_idx+1)
        plt.semilogy(peak_idx - xs, fit_left['A'] * np.exp(-xs/fit_left['xi']), '--', label=f'left exp xi={fit_left["xi"]:.2f}')
    if args.fit in ('exp','both') and fit_right is not None:
        xs = np.arange(0, N-peak_idx)
        plt.semilogy(peak_idx + xs, fit_right['A'] * np.exp(-xs/fit_right['xi']), '--', label=f'right exp xi={fit_right["xi"]:.2f}')
    if args.fit in ('cosh','both'):
        # perform cosh fits (grid) with chosen p
        cosh_left = None
        cosh_right = None
        try:
            cosh_left = fit_cosh_grid(np.arange(0, peak_idx+1)[::-1], rho[:peak_idx+1][::-1], p=args.cosh_p)
            cosh_right = fit_cosh_grid(np.arange(0, N-peak_idx), rho[peak_idx:], p=args.cosh_p)
        except Exception:
            cosh_left = cosh_right = None
        if cosh_left is not None:
            xs = np.arange(0, peak_idx+1)
            plt.semilogy(peak_idx - xs, cosh_left['A'] * (1.0/np.cosh(xs/cosh_left['xi'])**args.cosh_p), '-.', label=f'left cosh p={args.cosh_p} xi={cosh_left["xi"]:.2f}')
        if cosh_right is not None:
            xs = np.arange(0, N-peak_idx)
            plt.semilogy(peak_idx + xs, cosh_right['A'] * (1.0/np.cosh(xs/cosh_right['xi'])**args.cosh_p), '-.', label=f'right cosh p={args.cosh_p} xi={cosh_right["xi"]:.2f}')
    plt.axvline(peak_idx, color='k', linestyle=':', label=f'peak={peak_idx}')
    plt.xlabel('site j')
    plt.ylabel(r'$\rho_j$')
    plt.title('Fermion mode density and exponential fits')
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out)
    plt.close()

    # write summary
    with open(args.summary, 'w') as f:
        f.write(f'N = {N}\\n')
        f.write(f'peak_idx = {peak_idx}\\n')
        f.write(f'max_rho = {float(rho[peak_idx]):.6e}\\n')
        if fit_left is not None:
            f.write('left_exp_fit: A={A:.6e}, xi={xi:.6f}, res_rms={res_rms:.6e}, n={n}\n'.format(**fit_left))
        else:
            f.write('left_exp_fit: None or insufficient points\n')
        if fit_right is not None:
            f.write('right_exp_fit: A={A:.6e}, xi={xi:.6f}, res_rms={res_rms:.6e}, n={n}\n'.format(**fit_right))
        else:
            f.write('right_exp_fit: None or insufficient points\n')
        # cosh results
        if args.fit in ('cosh','both'):
            if cosh_left is not None:
                f.write('left_cosh_fit: A={A:.6e}, xi={xi:.6f}, p={p}, res_rms={res_rms:.6e}, n={n}\n'.format(**cosh_left))
            else:
                f.write('left_cosh_fit: None or insufficient points\n')
            if cosh_right is not None:
                f.write('right_cosh_fit: A={A:.6e}, xi={xi:.6f}, p={p}, res_rms={res_rms:.6e}, n={n}\n'.format(**cosh_right))
            else:
                f.write('right_cosh_fit: None or insufficient points\n')
        

    print('Saved density plot to', args.out)
    print('Saved summary to', args.summary)


if __name__ == '__main__':
    main()
