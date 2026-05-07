#!/usr/bin/env python3
"""Verify braid geometric vs dynamical decomposition for mapped two-level model.

Produces `results/verify_braid_geom_dyn.pkl` and prints summary.

Workflow:
- build two-level H(t) from tetron path (uses existing mapping)
- compute exact time-ordered U(T)
- compute adiabatic decomposition U = U_geom * U_dyn (approx.)
- find scaling alpha s.t. dynamical phases become equal (U_dyn ~ global phase * I)
"""
import os
import sys
import pickle
import numpy as np
from scipy.linalg import expm, eigh

# repo root
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import generate_comparison_panels as g
import tools.paper_params as P

# Pauli matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA = [sigma_x, sigma_y, sigma_z]


def compute_U_from_Hlist(H_list, dt):
    U = np.eye(2, dtype=complex)
    Ulist = []
    for H in H_list:
        U = expm(-1j * H * dt) @ U
        Ulist.append(U.copy())
    return U, Ulist


def rot_axis_angle_from_U(U):
    R = np.zeros((3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            R[i, j] = 0.5 * np.trace(SIGMA[i] @ U @ SIGMA[j] @ U.conj().T).real
    # clamp to real
    R = R.real
    tr = np.trace(R)
    # numerical clamp for argument
    val = (tr - 1.0) / 2.0
    val = min(1.0, max(-1.0, val))
    angle = float(np.arccos(val))
    vals, vecs = np.linalg.eig(R)
    idx = int(np.argmin(np.abs(vals - 1.0)))
    axis = vecs[:, idx].real
    if np.linalg.norm(axis) < 1e-12:
        axis = np.array([1.0, 0.0, 0.0])
    else:
        axis = axis / np.linalg.norm(axis)
    return axis, angle


def eigenbasis_series(H_list):
    N = len(H_list)
    e_vals = np.zeros((N, 2))
    e_vecs = np.zeros((N, 2, 2), dtype=complex)
    overlaps_products = np.ones(2, dtype=complex)
    prev_v = None
    for i, H in enumerate(H_list):
        w, v = eigh(H)
        # v columns are eigenvectors
        if i == 0:
            e_vals[0] = w
            e_vecs[0] = v
            prev_v = v.copy()
        else:
            # match ordering by maximizing overlaps
            ov = np.abs(prev_v.conj().T @ v)
            score1 = ov[0, 0] + ov[1, 1]
            score2 = ov[0, 1] + ov[1, 0]
            if score1 >= score2:
                perm = [0, 1]
            else:
                perm = [1, 0]
            vperm = v[:, perm]
            # accumulate raw overlaps for geometric phase
            for k in range(2):
                overlaps_products[k] *= np.vdot(prev_v[:, k], vperm[:, k])
            # fix phase of vperm for continuity (make overlap real positive)
            for k in range(2):
                ovk = np.vdot(prev_v[:, k], vperm[:, k])
                if abs(ovk) > 1e-16:
                    vperm[:, k] *= (ovk / abs(ovk)).conjugate()
            e_vals[i] = w[perm]
            e_vecs[i] = vperm
            prev_v = vperm.copy()
    return e_vals, e_vecs, overlaps_products


def adiabatic_decomposition(H_list, dt):
    e_vals, e_vecs, overlaps_products = eigenbasis_series(H_list)
    # dynamical phases: theta_n = integral of eigenvalues
    theta = np.sum(e_vals, axis=0) * dt
    # geometric gamma from overlaps product: use -angle(product)
    gamma = -np.angle(overlaps_products)
    V0 = e_vecs[0]
    VT = e_vecs[-1]
    U_dyn = V0 @ np.diag(np.exp(-1j * theta)) @ V0.conj().T
    U_geom = VT @ np.diag(np.exp(1j * gamma)) @ V0.conj().T
    return U_geom, U_dyn, theta, gamma


def compute_bloch_trajectory(Ulist, psi0):
    traj = np.zeros((len(Ulist), 3), dtype=float)
    for i, U in enumerate(Ulist):
        psi = U @ psi0
        traj[i, 0] = np.real(np.vdot(psi, sigma_x @ psi))
        traj[i, 1] = np.real(np.vdot(psi, sigma_y @ psi))
        traj[i, 2] = np.real(np.vdot(psi, sigma_z @ psi))
    return traj


def plot_bloch_sphere(traj, out, title=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(4, 4), dpi=200)
    ax = fig.add_subplot(111, projection='3d')
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, color='lightgray', alpha=0.25, linewidth=0)
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], '-o', markersize=2, lw=1)
    ax.scatter([traj[0, 0]], [traj[0, 1]], [traj[0, 2]], color='green', s=40, label='start')
    ax.scatter([traj[-1, 0]], [traj[-1, 1]], [traj[-1, 2]], color='red', s=40, label='end')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    if title:
        ax.set_title(title)
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def alpha_scan_and_plot(H_list, dt, alphas, alpha_target=None, out_prefix='results/alpha_scan'):
    axes = np.zeros((len(alphas), 3), dtype=float)
    angles = np.zeros(len(alphas), dtype=float)
    for i, a in enumerate(alphas):
        Hs = [a * H for H in H_list]
        U_final, _ = compute_U_from_Hlist(Hs, dt)
        axis, angle = rot_axis_angle_from_U(U_final)
        axes[i, :] = axis
        angles[i] = angle
        if (i + 1) % max(1, len(alphas)//10) == 0:
            print(f'alpha_scan progress {i+1}/{len(alphas)}')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 1, figsize=(6, 6), dpi=200)
    ax[0].semilogx(alphas, angles, '-o')
    ax[0].set_ylabel('rotation angle (rad)')
    ax[0].set_xlabel('alpha')
    if alpha_target is not None:
        ax[0].axvline(alpha_target, color='C1', linestyle='--', label='alpha_target')
        ax[0].legend()

    ax[1].semilogx(alphas, axes[:, 0], label='axis x')
    ax[1].semilogx(alphas, axes[:, 1], label='axis y')
    ax[1].semilogx(alphas, axes[:, 2], label='axis z')
    ax[1].set_xlabel('alpha')
    ax[1].set_ylabel('axis component')
    ax[1].legend()
    plt.tight_layout()
    out = f'{out_prefix}.png'
    fig.savefig(out)
    plt.close(fig)
    print('Saved', out)
    return alphas, axes, angles


def plot_U_comparison(U_final, U_geom, U_dyn, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    prod = U_geom @ U_dyn
    diff = np.abs(U_final - prod)
    fig, axs = plt.subplots(1, 3, figsize=(9, 3), dpi=200)
    im0 = axs[0].imshow(np.abs(U_final), cmap='viridis')
    axs[0].set_title('abs(U_final)')
    fig.colorbar(im0, ax=axs[0])
    im1 = axs[1].imshow(np.abs(prod), cmap='viridis')
    axs[1].set_title('abs(U_geom*U_dyn)')
    fig.colorbar(im1, ax=axs[1])
    im2 = axs[2].imshow(diff, cmap='magma')
    axs[2].set_title('abs(diff)')
    fig.colorbar(im2, ax=axs[2])
    plt.suptitle('U vs U_geom*U_dyn')
    plt.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print('Saved', out)


def main():
    os.makedirs('results', exist_ok=True)
    # load mapping parameters if present
    mpf = os.path.join('results', 'mapping_fit.npz')
    if os.path.exists(mpf):
        dd = np.load(mpf)
        if 'A0_fit' in dd.files:
            A0 = float(dd['A0_fit']); B0 = float(dd['B0_fit']); C0 = float(dd['C0_fit']); ts = float(dd['ts_fit'])
        else:
            res = dd['res'] if 'res' in dd.files else dd
            A0, B0, C0, ts = [float(x) for x in res[:4]]
    else:
        A0 = 0.00216; B0 = 0.03209; C0 = 0.16; ts = 0.7576

    print('Using A0,B0,C0,ts =', A0, B0, C0, ts)

    # Build path (use Fig.3 parameters as braid-like path)
    T = P.FIG3_T
    Nper = P.FIG3_N_PER_STEP
    tlist, E_bdg, S, M, Theta, step_idx, slist, dt = g.compute_bdg_minE(T_step=T, n_per_step=Nper, mode='minE')
    H_list = g.build_two_level_H_list(A0, B0, C0, ts, step_idx, slist)

    # base eigenvalue integral
    epos = np.array([np.max(np.linalg.eigvalsh(H)) for H in H_list])
    I0 = float(epos.sum() * dt)
    print('Base integral I0 =', I0)
    if I0 <= 0:
        print('Warning: base integral non-positive; cannot scale to target')

    alpha_target = np.pi / I0 if I0 > 0 else np.nan
    print('Alpha target (m=1) =', alpha_target)

    results = {}
    for alpha in [1.0, alpha_target]:
        if not np.isfinite(alpha):
            continue
        Hs = [alpha * H for H in H_list]
        U_final, Ulist = compute_U_from_Hlist(Hs, dt)
        axis, angle = rot_axis_angle_from_U(U_final)
        U_geom, U_dyn, theta, gamma = adiabatic_decomposition(Hs, dt)
        prod = U_geom @ U_dyn
        res_norm = float(np.linalg.norm(U_final - prod))
        evals = np.linalg.eigvals(U_dyn)
        diff_eval = float(np.max(np.abs(evals - evals[0])))
        print(f'alpha={alpha:.6g}: axis={axis}, angle(rad)={angle:.6g}, ||U-Ug*Ud||={res_norm:.3e}, Udyn_eigs={evals}, maxdiff={diff_eval:.3e}')
        results[f'alpha_{alpha:.6g}'] = {
            'axis': axis,
            'angle': angle,
            'res_norm': res_norm,
            'U_final': U_final,
            'U_geom': U_geom,
            'U_dyn': U_dyn,
            'evals_U_dyn': evals,
            'diff_eval': diff_eval,
            'theta': theta,
            'gamma': gamma,
            'I0': I0,
            'alpha_target': alpha_target,
        }

    # --- plots and report ---
    # recompute Ulist for each alpha entry to get trajectories and comparisons
    for alpha_key in list(results.keys()):
        try:
            kval = float(alpha_key.split('_')[1])
        except Exception:
            continue
        Hs = [kval * H for H in H_list]
        U_final, Ulist = compute_U_from_Hlist(Hs, dt)
        # initial state: ground eigenvector of H_list[0]
        w0, v0 = eigh(H_list[0])
        psi0 = v0[:, 0]
        traj = compute_bloch_trajectory(Ulist, psi0)
        outb = f'results/bloch_traj_alpha{kval:.6g}.png'
        plot_bloch_sphere(traj, outb, title=f'Bloch traj alpha={kval:.6g}')
        # comparison plot U vs product
        U_geom = results[alpha_key]['U_geom']
        U_dyn = results[alpha_key]['U_dyn']
        outc = f'results/compare_U_prod_alpha{kval:.6g}.png'
        plot_U_comparison(U_final, U_geom, U_dyn, outc)

    # alpha scan (log-spaced)
    alphas = np.logspace(-2, 0.5, 80)
    _alphas, _axes_arr, _angles = alpha_scan_and_plot(H_list, dt, alphas, alpha_target=alpha_target, out_prefix='results/alpha_scan')

    with open('results/verify_braid_geom_dyn.pkl', 'wb') as f:
        pickle.dump(results, f)
    print('Saved results/verify_braid_geom_dyn.pkl')


if __name__ == '__main__':
    main()
