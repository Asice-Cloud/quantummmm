#!/usr/bin/env python3
"""Plot eigenenergy surfaces of a two-site Hamiltonian

H = sum_{mu in {0,x,y,z}} c_{mu mu} (sigma^mu ⊗ sigma^mu)

Produces PNGs saved to `results/` for several c_zz values.
"""
import numpy as np
import os
from numpy import kron
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

outdir = 'results'
os.makedirs(outdir, exist_ok=True)

sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
si = np.eye(2, dtype=complex)

def H_matrix(c00, cxx, cyy, czz):
    H = c00 * kron(si, si) + cxx * kron(sx, sx) + cyy * kron(sy, sy) + czz * kron(sz, sz)
    return H

def compute_grid(czz, c00=0.0, rx=(-2,2), ry=(-2,2), N=201):
    cx = np.linspace(rx[0], rx[1], N)
    cy = np.linspace(ry[0], ry[1], N)
    CX, CY = np.meshgrid(cx, cy)
    E0 = np.zeros_like(CX)
    E1 = np.zeros_like(CX)
    E2 = np.zeros_like(CX)
    E3 = np.zeros_like(CX)
    for i in range(N):
        for j in range(N):
            H = H_matrix(c00, CX[i,j], CY[i,j], czz)
            evals = np.linalg.eigvalsh(H)
            E0[i,j], E1[i,j], E2[i,j], E3[i,j] = evals
    return CX, CY, E0, E1, E2, E3

def plot_surface(X, Y, Z, title, fname, cmap='viridis'):
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, rcount=200, ccount=200, cmap=cmap, linewidth=0, antialiased=False)
    ax.set_xlabel('c_xx'); ax.set_ylabel('c_yy'); ax.set_zlabel('Energy')
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(fname, dpi=200)
    plt.close(fig)

def main():
    # choose c_zz samples to inspect (more samples, wider range)
    import numpy as _np
    czz_list = list(_np.linspace(-2.0, 2.0, 13))
    for czz in czz_list:
        print('Computing grid for c_zz=', czz)
        CX, CY, E0, E1, E2, E3 = compute_grid(czz, c00=0.0, rx=(-4,4), ry=(-4,4), N=301)
        # plot lowest two bands E0, E1 and highest two E2, E3
        fname = os.path.join(outdir, f'energy_E0_czz_{czz:.3f}.png')
        plot_surface(CX, CY, E0, f'Lowest eigenvalue E0 (c_zz={czz})', fname, cmap='viridis')
        fname = os.path.join(outdir, f'energy_E1_czz_{czz:.3f}.png')
        plot_surface(CX, CY, E1, f'Second eigenvalue E1 (c_zz={czz})', fname, cmap='plasma')
        fname = os.path.join(outdir, f'energy_E2_czz_{czz:.3f}.png')
        plot_surface(CX, CY, E2, f'Third eigenvalue E2 (c_zz={czz})', fname, cmap='magma')
        fname = os.path.join(outdir, f'energy_E3_czz_{czz:.3f}.png')
        plot_surface(CX, CY, E3, f'Highest eigenvalue E3 (c_zz={czz})', fname, cmap='cividis')
        print('Saved energy plots for c_zz=', czz)

if __name__ == '__main__':
    main()
