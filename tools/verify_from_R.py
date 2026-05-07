#!/usr/bin/env python3
"""Re-derive H(u) from R(u) (eight-vertex), project to low-energy subspace,
construct H_eff = d·σ, simulate U(t), decompose U = U_geom * U_dyn,
scan alpha and save diagnostics/plots.

Saves PNGs and a results .npz per delta.
"""
import os
import sys
import math
import pickle
import numpy as np
from scipy.linalg import expm, eigh

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# Pauli matrices (2x2)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_I = np.eye(2, dtype=complex)

# 4x4 tensor basis helpers
def kron(a, b):
    return np.kron(a, b)


def R_eight_vertex(u, delta):
    c = math.cos(u)
    s = math.sin(u)
    return np.array([
        [c, 0, 0, 1j * delta * s],
        [0, s, c, 0],
        [0, c, s, 0],
        [1j * delta * s, 0, 0, c]
    ], dtype=complex)


def compute_Hs_from_R(delta, N=600):
    us = np.linspace(0.0, 2.0 * np.pi, N)
    Rs = [R_eight_vertex(u, delta) for u in us]
    du = us[1] - us[0]
    H4s = []
    for i, u in enumerate(us):
        if 0 < i < N - 1:
            dR = (Rs[i + 1] - Rs[i - 1]) / (2.0 * du)
        elif i == 0:
            dR = (Rs[1] - Rs[0]) / du
        else:
            dR = (Rs[-1] - Rs[-2]) / du
        try:
            Rinv = np.linalg.inv(Rs[i])
        except Exception:
            Rinv = np.linalg.pinv(Rs[i])
        H = 1j * dR.dot(Rinv)
        # make Hermitian numeric
        H = 0.5 * (H + H.conj().T)
        H4s.append(H)
    return us, du, H4s


def project_to_two_level(H4):
    # basis order |00,01,10,11> index 0..3
    # low-energy subspace chosen as |01> (idx=1), |10> (idx=2)
    H2 = np.zeros((2, 2), dtype=complex)
    H2[0, 0] = H4[1, 1]
    H2[0, 1] = H4[1, 2]
    H2[1, 0] = H4[2, 1]
    H2[1, 1] = H4[2, 2]
    return H2


def decompose_H2_to_d(H2):
    # H2 = d0 I + d·σ
    d0 = 0.5 * np.trace(H2)
    dx = 0.5 * np.trace(H2 @ sigma_x)
    dy = 0.5 * np.trace(H2 @ sigma_y)
    dz = 0.5 * np.trace(H2 @ sigma_z)
    dvec = np.array([dx, dy, dz], dtype=complex)
    return float(d0), dvec.real


def build_H2_list_from_H4s(H4s):
    H2_list = []
    ds = []
    for H4 in H4s:
        H2 = project_to_two_level(H4)
        d0, dvec = decompose_H2_to_d(H2)
        H2_list.append(dvec[0] * sigma_x + dvec[1] * sigma_y + dvec[2] * sigma_z)
        ds.append(dvec)
    return H2_list, np.array(ds)


def compute_U_from_Hlist(H_list, dt):
    U = np.eye(2, dtype=complex)
    Ulist = []
    for H in H_list:
        U = expm(-1j * H * dt) @ U
        Ulist.append(U.copy())
    return U, Ulist


def eigenbasis_series(H_list):
    N = len(H_list)
    e_vals = np.zeros((N, 2))
    e_vecs = np.zeros((N, 2, 2), dtype=complex)
    overlaps_products = np.ones(2, dtype=complex)
    prev_v = None
    for i, H in enumerate(H_list):
        w, v = eigh(H)
        if i == 0:
            e_vals[0] = w
            e_vecs[0] = v
            prev_v = v.copy()
        else:
            ov = np.abs(prev_v.conj().T @ v)
            score1 = ov[0, 0] + ov[1, 1]
            score2 = ov[0, 1] + ov[1, 0]
            perm = [0, 1] if score1 >= score2 else [1, 0]
            vperm = v[:, perm]
            for k in range(2):
                overlaps_products[k] *= np.vdot(prev_v[:, k], vperm[:, k])
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
    theta = np.sum(e_vals, axis=0) * dt
    gamma = -np.angle(overlaps_products)
    V0 = e_vecs[0]
    VT = e_vecs[-1]
    U_dyn = V0 @ np.diag(np.exp(-1j * theta)) @ V0.conj().T
    U_geom = VT @ np.diag(np.exp(1j * gamma)) @ V0.conj().T
    return U_geom, U_dyn, theta, gamma


def rot_axis_angle_from_U(U):
    # convert SU(2) to rotation axis and angle
    # Remove global phase
    det = np.linalg.det(U)
    U0 = U / (det ** 0.5)
    # rotation vector via trace
    tr = np.trace(U0)
    angle = math.acos(max(-1.0, min(1.0, (tr.real - 1.0) / 2.0)))
    # extract axis from antihermitian part
    ax = np.array([np.real(U0[1, 0] + U0[0, 1]), np.imag(U0[1, 0] - U0[0, 1]), np.real(U0[0, 0] - U0[1, 1])])
    if np.linalg.norm(ax) < 1e-12:
        ax = np.array([1.0, 0.0, 0.0])
    else:
        ax = ax / np.linalg.norm(ax)
    return ax, float(angle)


def compute_bloch_trajectory(Ulist, psi0):
    traj = np.zeros((len(Ulist), 3), dtype=float)
    for i, U in enumerate(Ulist):
        psi = U @ psi0
        traj[i, 0] = float(np.real(np.vdot(psi, sigma_x @ psi)))
        traj[i, 1] = float(np.real(np.vdot(psi, sigma_y @ psi)))
        traj[i, 2] = float(np.real(np.vdot(psi, sigma_z @ psi)))
    return traj


def plot_and_save_traj(traj, out, title=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa
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


def alpha_scan(H_list, dt, alpha_vals, out_prefix):
    axes = np.zeros((len(alpha_vals), 3), dtype=float)
    angles = np.zeros(len(alpha_vals), dtype=float)
    for i, a in enumerate(alpha_vals):
        Hs = [a * H for H in H_list]
        U_final, _ = compute_U_from_Hlist(Hs, dt)
        axis, angle = rot_axis_angle_from_U(U_final)
        axes[i] = axis
        angles[i] = angle
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(6, 6), dpi=200)
    ax[0].semilogx(alpha_vals, angles, '-o')
    ax[0].set_ylabel('rotation angle (rad)')
    ax[0].set_xlabel('alpha')
    ax[1].semilogx(alpha_vals, axes[:, 0], label='axis x')
    ax[1].semilogx(alpha_vals, axes[:, 1], label='axis y')
    ax[1].semilogx(alpha_vals, axes[:, 2], label='axis z')
    ax[1].set_xlabel('alpha')
    ax[1].set_ylabel('axis component')
    ax[1].legend()
    plt.tight_layout()
    out = out_prefix + '.png'
    fig.savefig(out)
    plt.close(fig)
    return axes, angles


def run_for_delta(delta, N=600):
    os.makedirs('results', exist_ok=True)
    us, du, H4s = compute_Hs_from_R(delta, N=N)
    H2_list, ds = build_H2_list_from_H4s(H4s)
    Evals = np.array([np.linalg.norm(d) for d in ds])
    dt = du
    I0 = float(np.sum(Evals) * dt)
    alpha_target = math.pi / I0 if I0 > 0 else np.nan

    report = {'delta': delta, 'I0': I0, 'alpha_target': alpha_target}

    for alpha in [1.0, alpha_target]:
        if not np.isfinite(alpha):
            continue
        Hs = [alpha * H for H in H2_list]
        U_final, Ulist = compute_U_from_Hlist(Hs, dt)
        axis, angle = rot_axis_angle_from_U(U_final)
        U_geom, U_dyn, theta, gamma = adiabatic_decomposition(Hs, dt)
        prod = U_geom @ U_dyn
        res_norm = float(np.linalg.norm(U_final - prod))
        evals_dyn = np.linalg.eigvals(U_dyn)
        report[f'alpha_{alpha:.6g}'] = {
            'axis': axis,
            'angle': angle,
            'res_norm': res_norm,
            'evals_U_dyn': evals_dyn,
            'theta': theta,
            'gamma': gamma,
        }

        # initial state: ground eigenvector of first H
        w0, v0 = eigh(Hs[0])
        psi0 = v0[:, 0]
        traj = compute_bloch_trajectory(Ulist, psi0)
        outb = f'results/bloch_delta{delta:.3g}_alpha{alpha:.6g}.png'
        plot_and_save_traj(traj, outb, title=f'delta={delta:.3g}, alpha={alpha:.6g}')
        outc = f'results/compareU_delta{delta:.3g}_alpha{alpha:.6g}.png'
        plot_U_comparison(U_final, U_geom, U_dyn, outc)

    # alpha scan
    alphas = np.logspace(-3, 0.5, 80)
    axes, angles = alpha_scan(H2_list, dt, alphas, out_prefix=f'results/alpha_scan_delta{delta:.3g}')
    report['alpha_scan'] = {'alphas': alphas, 'axes': axes, 'angles': angles}

    # save numeric results
    np.savez(f'results/deriveR_delta{delta:.3g}.npz', us=us, ds=ds, E=Evals, I0=I0, alpha_target=alpha_target)
    with open(f'results/deriveR_report_delta{delta:.3g}.pkl', 'wb') as f:
        pickle.dump(report, f)
    print('Saved results for delta=', delta)
    return report


def main():
    deltas = [0.0, 0.1, 0.3]
    all_reports = {}
    for d in deltas:
        print('Running delta=', d)
        r = run_for_delta(d, N=600)
        all_reports[d] = r
    with open('results/deriveR_all_reports.pkl', 'wb') as f:
        pickle.dump(all_reports, f)
    print('All done. Reports saved to results/')


if __name__ == '__main__':
    main()
