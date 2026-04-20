#!/usr/bin/env python3
"""3D surface plots for F_- and F_+ loci over complex X plane.
Plots log10(| |Z|-1 |) as a surface over Re(X), Im(X) for fixed b_y values.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

outdir = 'results'
os.makedirs(outdir, exist_ok=True)

# grid (moderate resolution)
N = 300
re = np.linspace(-4, 4, N)
im = np.linspace(-4, 4, N)
Re, Im = np.meshgrid(re, im)
Xgrid = Re + 1j * Im

b_y_list = list(np.linspace(0.0, 1.5, 9))  # more samples for observation

for b_y in b_y_list:
    Y = np.exp(2j * b_y)
    X = Xgrid
    denom = X + Y
    eps = 1e-12
    denom_safe = denom.copy()
    denom_safe[np.abs(denom_safe) < eps] = eps

    Z_minus = (X * Y - 1.0) / denom_safe
    Z_plus = (1.0 - X * Y) / denom_safe

    diff_minus = np.abs(np.abs(Z_minus) - 1.0)
    diff_plus = np.abs(np.abs(Z_plus) - 1.0)

    # log scale
    log_minus = np.log10(diff_minus + 1e-16)
    log_plus = np.log10(diff_plus + 1e-16)

    # clip extreme values for better visualization
    clim_min, clim_max = -8, 1
    log_minus_clip = np.clip(log_minus, clim_min, clim_max)
    log_plus_clip = np.clip(log_plus, clim_min, clim_max)

    # 3D plot for F_-
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(Re, Im, log_minus_clip, rcount=200, ccount=200, cmap='viridis', linewidth=0, antialiased=False)
    ax.set_title(f'F_- log10(| |Z|-1 |)  b_y={b_y}')
    ax.set_xlabel('Re X'); ax.set_ylabel('Im X'); ax.set_zlabel('log10(| |Z|-1 |)')
    fig.colorbar(surf, shrink=0.6, aspect=8)
    ax.view_init(elev=30, azim=135)
    plt.tight_layout()
    fname = os.path.join(outdir, f'3d_Fminus_b_y_{b_y:.3f}.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)

    # 3D plot for F_+
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(Re, Im, log_plus_clip, rcount=200, ccount=200, cmap='magma', linewidth=0, antialiased=False)
    ax.set_title(f'F_+ log10(| |Z|-1 |)  b_y={b_y}')
    ax.set_xlabel('Re X'); ax.set_ylabel('Im X'); ax.set_zlabel('log10(| |Z|-1 |)')
    fig.colorbar(surf, shrink=0.6, aspect=8)
    ax.view_init(elev=30, azim=135)
    plt.tight_layout()
    fname2 = os.path.join(outdir, f'3d_Fplus_b_y_{b_y:.3f}.png')
    fig.savefig(fname2, dpi=200)
    plt.close(fig)

print('3D plots saved to', outdir)
