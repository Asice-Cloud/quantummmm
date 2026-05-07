#!/usr/bin/env python3
"""Quantify non-Abelian Berry curvature / holonomy for the eight-vertex model.

Produces a (u,delta) grid of the projected subspace basis W(u,delta), computes
matrix-valued Berry connections A_u, A_delta, curvature F_{u,delta}, Wilson
plaquette holonomies and gauge-invariant diagnostics.

Saves results to `results/nonabelian_grid.npz` and heatmaps in `results/`.
"""
import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh


def R_eight_vertex(u, delta):
    c = math.cos(u)
    s = math.sin(u)
    return np.array([
        [c, 0, 0, 1j * delta * s],
        [0, s, c, 0],
        [0, c, s, 0],
        [1j * delta * s, 0, 0, c]
    ], dtype=complex)


def compute_H4_at(u, delta, du):
    # periodic in u -> use wrap-around for u-steps
    Rp = R_eight_vertex(u + du, delta)
    Rm = R_eight_vertex(u - du, delta)
    R0 = R_eight_vertex(u, delta)
    dR = (Rp - Rm) / (2.0 * du)
    try:
        Rinv = np.linalg.inv(R0)
    except Exception:
        Rinv = np.linalg.pinv(R0)
    H = 1j * dR.dot(Rinv)
    H = 0.5 * (H + H.conj().T)
    return H


def project_to_subspace(H4):
    # Project to instantaneous low-energy m-dimensional subspace of H4.
    # Default behavior: return orthonormal columns spanning the m states
    # with smallest absolute eigenvalues (near-zero modes), m=2 for our use.
    w, v = eigh(H4)
    # pick eigenvectors with smallest |E|
    idx = np.argsort(np.abs(w))[:2]
    Vsmall = v[:, idx]
    # re-orthonormalize (numerical)
    Q, _ = np.linalg.qr(Vsmall)
    return Q[:, :2]


# Note: previous implementation embedded a 2x2 effective H into the 4-dim
# |01>,|10> subspace. We now construct W directly from H4 instantaneous
# eigenvectors (see project_to_subspace).


def align_grid_W(Wgrid, nu, nd, passes=2):
    # Wgrid shape (nu, nd, 4, m)
    for _ in range(passes):
        # align along u for each delta
        for j in range(nd):
            for i in range(1, nu):
                Wprev = Wgrid[i - 1, j]
                Wcur = Wgrid[i, j]
                M = Wprev.conj().T @ Wcur
                U, s, Vh = np.linalg.svd(M)
                R = Vh.conj().T @ U.conj().T
                Wgrid[i, j] = Wcur @ R
        # align along delta for each u
        for i in range(nu):
            for j in range(1, nd):
                Wprev = Wgrid[i, j - 1]
                Wcur = Wgrid[i, j]
                M = Wprev.conj().T @ Wcur
                U, s, Vh = np.linalg.svd(M)
                R = Vh.conj().T @ U.conj().T
                Wgrid[i, j] = Wcur @ R
    return Wgrid


def compute_connections(Wgrid, us, deltas):
    nu, nd, _, m = Wgrid.shape
    du = us[1] - us[0]
    dd = deltas[1] - deltas[0]
    A_u = np.zeros((nu, nd, m, m), dtype=complex)
    A_d = np.zeros((nu, nd, m, m), dtype=complex)
    # use periodic wrap in u direction
    for i in range(nu):
        ip = (i + 1) % nu
        im = (i - 1) % nu
        for j in range(nd):
            dW_du = (Wgrid[ip, j] - Wgrid[im, j]) / (2.0 * du)
            A_u[i, j] = 1j * (Wgrid[i, j].conj().T @ dW_du)
    for j in range(nd):
        for i in range(nu):
            if 0 < j < nd - 1:
                jp = j + 1
                jm = j - 1
                dW_dd = (Wgrid[i, jp] - Wgrid[i, jm]) / (2.0 * dd)
            elif j == 0:
                dW_dd = (Wgrid[i, 1] - Wgrid[i, 0]) / dd
            else:
                dW_dd = (Wgrid[i, -1] - Wgrid[i, -2]) / dd
            A_d[i, j] = 1j * (Wgrid[i, j].conj().T @ dW_dd)
    return A_u, A_d


def compute_curvature_and_wilson(Wgrid, A_u, A_d, us, deltas):
    nu, nd, _, m = Wgrid.shape
    du = us[1] - us[0]
    dd = deltas[1] - deltas[0]
    Fnorm = np.zeros((nu, nd))
    Udev = np.zeros((nu - 1, nd - 1))
    eigvars = {}
    for i in range(1, nu - 1):
        for j in range(1, nd - 1):
            # derivatives of A
            A_d_ip = (A_d[i + 1, j] - A_d[i - 1, j]) / (2.0 * du)
            A_u_jp = (A_u[i, j + 1] - A_u[i, j - 1]) / (2.0 * dd)
            comm = A_u[i, j] @ A_d[i, j] - A_d[i, j] @ A_u[i, j]
            F = A_d_ip - A_u_jp + comm
            Fnorm[i, j] = np.linalg.norm(F, ord='fro')
    # Wilson plaquette (discrete holonomy) and deviation
    for i in range(nu - 1):
        for j in range(nd - 1):
            M1 = Wgrid[i, j].conj().T @ Wgrid[i + 1, j]
            M2 = Wgrid[i + 1, j].conj().T @ Wgrid[i + 1, j + 1]
            M3 = Wgrid[i + 1, j + 1].conj().T @ Wgrid[i, j + 1]
            M4 = Wgrid[i, j + 1].conj().T @ Wgrid[i, j]
            Ucell = M1 @ M2 @ M3 @ M4
            area = du * dd
            Udev[i, j] = np.linalg.norm(Ucell - np.eye(m), ord='fro') / (area + 1e-30)
            eigs = np.linalg.eigvals(Ucell)
            eigvars[(i, j)] = eigs
    return Fnorm, Udev, eigvars


def main():
    os.makedirs('results', exist_ok=True)
    Nu = 200
    Nd = 21
    us = np.linspace(0.0, 2.0 * math.pi, Nu, endpoint=False)
    deltas = np.linspace(0.0, 0.3, Nd)
    du = us[1] - us[0]
    m = 2

    Wgrid = np.zeros((Nu, Nd, 4, m), dtype=complex)
    print('Computing subspace bases on grid...')
    for j, delta in enumerate(deltas):
        for i, u in enumerate(us):
            H4 = compute_H4_at(u, delta, du)
            W = project_to_subspace(H4)
            Wgrid[i, j] = W

    print('Aligning gauge across grid...')
    Wgrid = align_grid_W(Wgrid, Nu, Nd, passes=2)

    print('Computing connections...')
    A_u, A_d = compute_connections(Wgrid, us, deltas)

    print('Computing curvature and Wilson diagnostics...')
    Fnorm, Udev, eigvars = compute_curvature_and_wilson(Wgrid, A_u, A_d, us, deltas)

    # integrated curvature norm
    total = np.sum(Fnorm**2) * (du * (deltas[1] - deltas[0]))
    np.savez('results/nonabelian_grid.npz', us=us, deltas=deltas, Fnorm=Fnorm, Udev=Udev, total=total)
    print('Saved results/nonabelian_grid.npz, total integrated ||F||^2 =', total)

    plt.figure(figsize=(6, 3))
    plt.subplot(1, 2, 1)
    plt.title('||F||_F')
    plt.imshow(Fnorm.T, origin='lower', aspect='auto', extent=[0, 2 * math.pi, deltas[0], deltas[-1]])
    plt.xlabel('u')
    plt.ylabel('delta')
    plt.colorbar()
    plt.subplot(1, 2, 2)
    plt.title('Wilson dev / area')
    plt.imshow(Udev.T, origin='lower', aspect='auto', extent=[0, 2 * math.pi, deltas[0], deltas[-2]])
    plt.xlabel('u')
    plt.colorbar()
    plt.tight_layout()
    plt.savefig('results/nonabelian_heatmaps.png', dpi=200)
    plt.close()

    print('Saved results/nonabelian_heatmaps.png')


if __name__ == '__main__':
    main()
