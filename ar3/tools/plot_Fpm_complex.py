#!/usr/bin/env python3
"""Plot complex loci for F_- and F_+ given by
F_-: XY - XZ - YZ - 1 = 0  => Z = (XY - 1)/(X+Y)
F_+: XY + XZ + YZ - 1 = 0  => Z = (1 - XY)/(X+Y)

We sample X in complex plane and fix Y = exp(2i b_y) for a few b_y values.
Save heatmaps of | |Z|-1 | and contours where |Z|=1 (real b_z solutions).
"""
import numpy as np
import matplotlib.pyplot as plt
import os

outdir = 'results'
os.makedirs(outdir, exist_ok=True)

# grid
N = 800
re = np.linspace(-4, 4, N)
im = np.linspace(-4, 4, N)
Re, Im = np.meshgrid(re, im)
Xgrid = Re + 1j * Im

# choose some sample b_y values
b_y_list = [0.0, 0.3, 1.0]

for b_y in b_y_list:
    Y = np.exp(2j * b_y)
    X = Xgrid
    denom = X + Y

    # avoid division by near-zero denom
    eps = 1e-12
    denom_safe = denom.copy()
    denom_safe[np.abs(denom_safe) < eps] = eps

    Z_minus = (X * Y - 1.0) / denom_safe
    Z_plus = (1.0 - X * Y) / denom_safe

    # compute metrics
    diff_minus = np.abs(np.abs(Z_minus) - 1.0)
    diff_plus = np.abs(np.abs(Z_plus) - 1.0)

    # log scale for visualization
    log_minus = np.log10(diff_minus + 1e-16)
    log_plus = np.log10(diff_plus + 1e-16)

    # plot heatmaps and contours
    fig, axs = plt.subplots(1,2,figsize=(12,5))
    im0 = axs[0].imshow(log_minus, extent=(re[0], re[-1], im[0], im[-1]), origin='lower', cmap='inferno')
    axs[0].set_title(f'F_- | |Z|-1 | (log10)  b_y={b_y}')
    axs[0].set_xlabel('Re X'); axs[0].set_ylabel('Im X')
    fig.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(log_plus, extent=(re[0], re[-1], im[0], im[-1]), origin='lower', cmap='inferno')
    axs[1].set_title(f'F_+ | |Z|-1 | (log10)  b_y={b_y}')
    axs[1].set_xlabel('Re X'); axs[1].set_ylabel('Im X')
    fig.colorbar(im1, ax=axs[1])

    # overlay contour where |Z|-1 ~ 0 -> use increasing contour levels
    contour_levels = np.sort(np.array([1e-6, 1e-4, 1e-3]))
    CS0 = axs[0].contour(Re, Im, diff_minus, levels=contour_levels, colors=['cyan','lime','white'], linewidths=0.8)
    CS1 = axs[1].contour(Re, Im, diff_plus, levels=contour_levels, colors=['cyan','lime','white'], linewidths=0.8)

    plt.tight_layout()
    fname = os.path.join(outdir, f'Fpm_b_y_{b_y:.3f}.png')
    fig.savefig(fname, dpi=150)
    plt.close(fig)

    # also save scatter of approximate solutions where diff < tol
    tol = 1e-3
    mask_minus = diff_minus < tol
    mask_plus = diff_plus < tol

    fig2, ax2 = plt.subplots(1,2,figsize=(10,4))
    if np.any(mask_minus):
        ax2[0].scatter(Re[mask_minus].ravel(), Im[mask_minus].ravel(), s=1, color='navy')
    else:
        # fallback: extract contour segments from CS0 (support multiple matplotlib APIs)
        try:
            for coll in CS0.collections:
                for path in coll.get_paths():
                    v = path.vertices
                    ax2[0].plot(v[:,0], v[:,1], '.', markersize=1, color='navy')
        except Exception:
            # QuadContourSet may expose segments via .allsegs
            try:
                for seglist in CS0.allsegs:
                    for seg in seglist:
                        if len(seg) > 0:
                            ax2[0].plot(seg[:,0], seg[:,1], '.', markersize=1, color='navy')
            except Exception:
                pass
    ax2[0].set_title(f'Approx |Z_-|=1 (tol={tol}) b_y={b_y}')
    ax2[0].set_xlim(re[0],re[-1]); ax2[0].set_ylim(im[0],im[-1])
    ax2[0].set_xlabel('Re X'); ax2[0].set_ylabel('Im X')

    if np.any(mask_plus):
        ax2[1].scatter(Re[mask_plus].ravel(), Im[mask_plus].ravel(), s=1, color='maroon')
    else:
        try:
            for coll in CS1.collections:
                for path in coll.get_paths():
                    v = path.vertices
                    ax2[1].plot(v[:,0], v[:,1], '.', markersize=1, color='maroon')
        except Exception:
            try:
                for seglist in CS1.allsegs:
                    for seg in seglist:
                        if len(seg) > 0:
                            ax2[1].plot(seg[:,0], seg[:,1], '.', markersize=1, color='maroon')
            except Exception:
                pass
    ax2[1].set_title(f'Approx |Z_+|=1 (tol={tol}) b_y={b_y}')
    ax2[1].set_xlim(re[0],re[-1]); ax2[1].set_ylim(im[0],im[-1])
    ax2[1].set_xlabel('Re X'); ax2[1].set_ylabel('Im X')

    plt.tight_layout()
    fname2 = os.path.join(outdir, f'Fpm_scatter_b_y_{b_y:.3f}.png')
    fig2.savefig(fname2, dpi=150)
    plt.close(fig2)

print('Plots saved to', outdir)
