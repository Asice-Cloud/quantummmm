#!/usr/bin/env python3
"""Compute zero-mode eigenvectors for Kitaev BdG chain and estimate localization length.
Usage: run as script; parameters are set in main block.
"""
import numpy as np
from numpy.linalg import eigh
import matplotlib.pyplot as plt
from matplotlib import patches

def map_c_to_params(c_xx, c_yy, c_xy=0.0, c_yx=0.0, c_zz=0.0, c_z0=0.0, c_0z=0.0):
    t = c_xx + c_yy + 1j*(c_xy - c_yx)
    Delta = c_xx - c_yy - 1j*(c_xy + c_yx)
    U = 4.0 * c_zz
    mu = 4.0*c_zz - 2.0*(c_z0 + c_0z)
    return t, Delta, mu, U

def build_kitaev_bdg(N, t, Delta, mu):
    A = np.zeros((N,N), dtype=complex)
    B = np.zeros((N,N), dtype=complex)
    for i in range(N):
        A[i,i] = -mu
        if i+1 < N:
            A[i,i+1] = -t
            A[i+1,i] = -t
            B[i,i+1] = Delta
            B[i+1,i] = -Delta
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    Hbdg = np.concatenate((top, bottom), axis=0)
    return Hbdg

def compute_site_density(psi, N):
    # psi is length 2N (u; v), site density = |u_j|^2 + |v_j|^2
    u = psi[:N]
    v = psi[N:]
    dens = np.abs(u)**2 + np.abs(v)**2
    s = np.sum(dens)
    if s>0:
        dens = dens / s
    return dens

def fit_exponential_from_edge(dens, max_sites=20):
    # fit log(dens) = a + b * x  => decay length xi = -1/b
    # take first max_sites where dens>threshold
    threshold = 1e-12
    idx = np.where(dens>threshold)[0]
    if idx.size==0:
        return None, None
    # take contiguous from 0 up to max_sites
    end = min(max_sites, len(dens))
    xs = np.arange(end)
    ys = dens[:end]
    mask = ys>threshold
    if np.sum(mask) < 3:
        return None, None
    xs = xs[mask]
    ys = ys[mask]
    logy = np.log(ys)
    # linear fit
    A = np.vstack([xs, np.ones_like(xs)]).T
    b, a = np.linalg.lstsq(A, logy, rcond=None)[0]
    xi = -1.0/b if abs(b)>1e-16 else None
    return xi, (a,b)


def fit_exponential_from_peak(dens, max_sites=30, threshold=1e-12):
    """Find peak and fit exponential decay outward from peak on both sides.
    Returns best fit (xi,(a,b),peak_index,side,points_used) or None if fit fails.
    side = 'right' means fit to increasing index from peak; 'left' means decreasing index.
    """
    N = len(dens)
    p = int(np.argmax(dens))

    def fit_segment(start, direction):
        xs = []
        ys = []
        for i in range(max_sites):
            idx = start + direction * i
            if idx < 0 or idx >= N:
                break
            val = dens[idx]
            if val <= threshold:
                break
            xs.append(i)
            ys.append(val)
        if len(xs) < 3:
            return None
        logy = np.log(ys)
        A = np.vstack([xs, np.ones_like(xs)]).T
        b, a = np.linalg.lstsq(A, logy, rcond=None)[0]
        xi = -1.0 / b if abs(b) > 1e-16 else None
        # compute residual norm
        fit = a + b * np.array(xs)
        res = np.linalg.norm(logy - fit)
        return xi, (a, b), res, len(xs)

    right = fit_segment(p, 1)
    left = fit_segment(p, -1)

    candidates = []
    if right is not None:
        xi, ab, res, pts = right
        candidates.append((xi, ab, p, 'right', res, pts))
    if left is not None:
        xi, ab, res, pts = left
        candidates.append((xi, ab, p, 'left', res, pts))

    if not candidates:
        # fallback: fixed-N points from peak outward (even if below threshold)
        N = len(dens)
        p = int(np.argmax(dens))
        eps = 1e-300
        def fit_fixed(start, direction, npoints=20):
            xs = []
            ys = []
            for i in range(npoints):
                idx = start + direction * i
                if idx < 0 or idx >= N:
                    break
                xs.append(i)
                ys.append(dens[idx] + eps)
            if len(xs) < 3:
                return None
            logy = np.log(ys)
            A = np.vstack([xs, np.ones_like(xs)]).T
            b, a = np.linalg.lstsq(A, logy, rcond=None)[0]
            xi = -1.0 / b if abs(b) > 1e-16 else None
            fit = a + b * np.array(xs)
            res = np.linalg.norm(logy - fit)
            return xi, (a, b), res, len(xs)

        right = fit_fixed(p, 1, npoints=min(max_sites, N-p))
        left = fit_fixed(p, -1, npoints=min(max_sites, p+1))
        candidates = []
        if right is not None:
            xi, ab, res, pts = right
            candidates.append((xi, ab, p, 'right_fixed', res, pts))
        if left is not None:
            xi, ab, res, pts = left
            candidates.append((xi, ab, p, 'left_fixed', res, pts))
        if not candidates:
            return None
    # choose candidate with more points, then smaller residual
    candidates.sort(key=lambda x: (-x[5], x[4]))
    best = candidates[0]
    return best[0], best[1], best[2], best[3], best[5]

if __name__ == '__main__':
    # parameters (example matching verify_mzm)
    c_xx = 1.0
    c_yy = 0.0
    c_xy = 0.0
    c_yx = 0.0
    c_zz = 0.0
    c_z0 = 0.0
    c_0z = 0.0

    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    N = 120
    Hbdg = build_kitaev_bdg(N, t, Delta, mu)
    eigs, vecs = eigh(Hbdg)
    eigs = np.real_if_close(eigs)
    idx_sorted = np.argsort(np.abs(eigs))

    # find near-zero modes
    tol = 1e-6
    zero_indices = [i for i in idx_sorted if abs(eigs[i])<tol]
    print(f'Found {len(zero_indices)} near-zero modes (tol={tol}).')
    if len(zero_indices)==0:
        print('No zero modes found; adjust parameters.')
        exit(1)

    # take first two zero modes (if present)
    saved_modes = []
    for k,i in enumerate(zero_indices[:2]):
        psi = vecs[:, i]
        dens = compute_site_density(psi, N)
        # robust fit from peak
        fit_res = fit_exponential_from_peak(dens, max_sites=40, threshold=1e-14)
        if fit_res is None:
            xi = None
            a = b = None
            peak_idx = np.argmax(dens)
            side = 'unknown'
            pts_used = 0
        else:
            xi, (a, b), peak_idx, side, pts_used = fit_res

        print(f'Zero mode {k+1} (eigenvalue={eigs[i]:.3e}), peak_index={peak_idx}, fitted_side={side}, points_used={pts_used}')
        print('site, density (around peak, +/-5):')
        start = max(0, peak_idx-5)
        end = min(N, peak_idx+6)
        for j in range(start, end):
            print(f'{j:3d}: {dens[j]:.4e}')
        if xi is None:
            print('Could not fit exponential decay for this mode (insufficient data).')
        else:
            print(f'Estimated localization length xi ~ {xi:.3f} sites (fit slope={b:.3e}).')
        print('')
        saved_modes.append((peak_idx, dens, xi, side))
    # also print symmetry: right edge densities
    # compute densities for mode that localizes on right by reversing order
    # (if only one mode found, skip second)
    # Plot densities for found zero modes with schematic chain below
    if saved_modes:
        fig = plt.figure(figsize=(10,4))
        gs = fig.add_gridspec(2, 1, height_ratios=[3, 0.6], hspace=0.12)
        ax = fig.add_subplot(gs[0])
        x = np.arange(N)
        # plot each mode: smooth a little for cartoon-like curve
        for idx, (peak_idx, dens, xi, side) in enumerate(saved_modes):
            # small smoothing for visual clarity
            window = max(1, min(7, N//50))
            kernel = np.ones(window) / window
            smooth = np.convolve(dens, kernel, mode='same')
            label = f'mode {idx+1} peak={peak_idx} xi={xi:.3f}' if xi is not None else f'mode {idx+1} peak={peak_idx}'
            ax.plot(x, smooth, '-', lw=1.5, alpha=0.9, label=label)
            ax.plot(x, dens, '.', markersize=3, alpha=0.6)
        ax.set_xlim(0, N-1)
        ax.set_xlabel('site')
        ax.set_ylabel('density')
        ax.set_title('Zero-mode site density')
        ax.legend(loc='upper right')
        ax.grid(True, linestyle='--', alpha=0.4)

        # schematic chain row
        ax2 = fig.add_subplot(gs[1])
        ax2.set_xlim(-0.5, N-0.5)
        ax2.set_ylim(-1, 1)
        ax2.axis('off')
        # draw chain bar
        bar = patches.Rectangle((0, -0.25), N-1, 0.5, facecolor='#f5e6c8', edgecolor='none')
        ax2.add_patch(bar)
        # draw endpoint circles
        r = 0.6
        left_circle = patches.Circle((0, 0), radius=r, facecolor='#ff9f1c', edgecolor='k')
        right_circle = patches.Circle((N-1, 0), radius=r, facecolor='#ff9f1c', edgecolor='k')
        ax2.add_patch(left_circle)
        ax2.add_patch(right_circle)
        # small caption above bar
        ax2.text((N-1)/2, 0.6, 'Density of states', ha='center', va='bottom', fontsize=10)

        outpng = 'tools/zero_mode_profile_cartoon.png'
        plt.tight_layout()
        plt.savefig(outpng, dpi=300)
        print(f'Zero-mode cartoon density plot saved to {outpng}')
