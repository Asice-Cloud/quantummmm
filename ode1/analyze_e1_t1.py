#!/usr/bin/env python3
"""
E1-t1 Interaction Analysis Tool
================================
Systematic study of ABS effects in Majorana braiding via SO(5) RK4.
Handles all parameter regimes including E1=0 and t1=0.

Produces six diagnostic figures:
  Fig 1: Bloch vector components vs time
  Fig 2: Effective omega (rotation generator) components
  Fig 3: Accumulated rotation angle phi(t)
  Fig 4: Bloch sphere 3D trajectories
  Fig 5: Final Bloch vector vs lg(t1/E1) sweep
  Fig 6: Multi-tau fidelity comparison (with --tau-list)

Usage:
  python analyze_e1_t1.py
  python analyze_e1_t1.py --E1 0.01 --tau 5.0
  python analyze_e1_t1.py --tau-list 1,5,20,50,100
"""
import numpy as np; pi = np.pi
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

# ═══════════════════════════════════════════════════════════════
# Physical parameters
# ═══════════════════════════════════════════════════════════════
tc = 0.3   # gate coupling amplitude (meV)
E0 = 0.3   # QD level amplitude (meV)

# Majorana indices (0-based): gamma_0=gamma_1, gamma_1=gamma_2,
#   gamma_2=gamma_3, gamma_3=gamma_a, gamma_4=gamma_b

# ═══════════════════════════════════════════════════════════════
# Gating functions
# ═══════════════════════════════════════════════════════════════
def fp(t, tau): return 0.5 * (1.0 + np.cos(pi * t / tau))
def fm(t, tau): return 0.5 * (1.0 - np.cos(pi * t / tau))

# ═══════════════════════════════════════════════════════════════
# SO(5) antisymmetric generator A(t)  --  dR/dt = A @ R
# ═══════════════════════════════════════════════════════════════
# H_EM = iEd gamma_a gamma_b + iE1 gamma_1 gamma_2 + i|t2| gamma_a gamma_2
#        - i|t1| gamma_b gamma_1 - i|t3| gamma_a gamma_3
# Convention: A[i,j] = 2 * (coefficient of i gamma_i gamma_j in H)

def A_step1(t, tau, E1, t1c):
    """G1 off, t2 opens, Ed -> 0."""
    A = np.zeros((5, 5))
    A[0, 1] = 2 * E1;                        A[1, 0] = -2 * E1
    A[1, 3] = -2 * tc * fm(t, tau);          A[3, 1] =  2 * tc * fm(t, tau)
    A[0, 4] =  2 * t1c * fm(t, tau);         A[4, 0] = -2 * t1c * fm(t, tau)
    A[3, 4] =  2 * E0 * fp(t, tau);          A[4, 3] = -2 * E0 * fp(t, tau)
    return A

def A_step2(t, tau, E1, t1c):
    """G2 off, G1 on, t3 opens, t2 closes."""
    A = np.zeros((5, 5))
    A[0, 1] = 2 * E1;                        A[1, 0] = -2 * E1
    A[1, 3] = -2 * tc * fp(t, tau);          A[3, 1] =  2 * tc * fp(t, tau)
    A[2, 3] =  2 * tc * fm(t, tau);          A[3, 2] = -2 * tc * fm(t, tau)
    A[0, 4] =  2 * t1c * fp(t, tau);         A[4, 0] = -2 * t1c * fp(t, tau)
    return A

def A_step3(t, tau, E1, t1c):
    """G2 on, t3 closes, Ed recovers."""
    A = np.zeros((5, 5))
    A[0, 1] = 2 * E1;                        A[1, 0] = -2 * E1
    A[2, 3] =  2 * tc * fp(t, tau);          A[3, 2] = -2 * tc * fp(t, tau)
    A[3, 4] =  2 * E0 * fm(t, tau);          A[4, 3] = -2 * E0 * fm(t, tau)
    return A

_A_FUNCS = [A_step1, A_step2, A_step3]

# ═══════════════════════════════════════════════════════════════
# SO(5) RK4 solver
# ═══════════════════════════════════════════════════════════════
def so5_protocol(tau, E1, t1c, n_per_step=600, record_traj=False):
    """Propagate SO(5) R through 3-step braid. Returns final R + optional traj."""
    dt = tau / n_per_step
    R = np.eye(5)

    if record_traj:
        n_total = 3 * n_per_step + 1
        t_hist = np.zeros(n_total)
        R_traj = np.zeros((n_total, 5, 5))
        t_hist[0] = 0.0; R_traj[0] = R.copy()
        idx = 1

    for step_idx, A_fn in enumerate(_A_FUNCS):
        for s in range(n_per_step):
            t = s * dt
            A0 = A_fn(t, tau, E1, t1c)
            k1 = A0 @ R
            A_half = A_fn(t + 0.5*dt, tau, E1, t1c)
            k2 = A_half @ (R + 0.5*dt*k1)
            k3 = A_half @ (R + 0.5*dt*k2)
            A_full = A_fn(t + dt, tau, E1, t1c)
            k4 = A_full @ (R + dt*k3)
            R = R + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

            if record_traj:
                t_hist[idx] = step_idx*tau + (s+1)*dt
                R_traj[idx] = R.copy()
                idx += 1

    if record_traj:
        return R, t_hist, R_traj
    return R

# ═══════════════════════════════════════════════════════════════
# Observables from SO(5) matrix
# ═══════════════════════════════════════════════════════════════
# sigma_x <-> i gamma_2 gamma_3 (indices 1,2)
# sigma_y <-> i gamma_3 gamma_1 (indices 2,0)
# sigma_z <-> i gamma_1 gamma_2 (indices 0,1)

def R_to_bloch_vec(R):
    """Bloch vector from {gamma_1,gamma_2,gamma_3} SO(3) sub-block."""
    vx = np.real(1j * (R[1, 2] - R[2, 1]))
    vy = np.real(R[2, 0] + R[0, 2])
    vz = np.real(R[0, 0] - R[1, 1])
    v = np.array([vx, vy, vz])
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else np.zeros(3)

def R_to_axis_angle(R):
    """Extract rotation axis n (3-vector) and angle phi (radians) from SO(3) sub-block."""
    R3 = R[:3, :3]
    cos_phi = np.clip((np.trace(R3) - 1.0) / 2.0, -1.0, 1.0)
    phi = np.arccos(cos_phi)
    if phi < 1e-12:
        return np.array([0.0, 0.0, 1.0]), 0.0
    n = np.array([R3[2, 1] - R3[1, 2],
                  R3[0, 2] - R3[2, 0],
                  R3[1, 0] - R3[0, 1]])
    n_norm = np.linalg.norm(n)
    return (n / n_norm, phi) if n_norm > 1e-12 else (np.array([0.,0.,1.]), phi)

def fidelity_double_braid(R):
    """Fidelity |<psi_1^-|U(6tau)|psi_1^+>|^2 after double braid."""
    ov = 0.5 * (R[0, 0] + 1j*R[1, 0] + 1j*R[0, 1] - R[1, 1])
    return np.abs(ov)**2

def extract_omega(R, A_full):
    """Instantaneous angular velocity omega from SO(5) generator.
    omega = projection of R^T A R onto the {gamma_1,gamma_2,gamma_3} so(3) subspace.
    wx <-> gamma_2 gamma_3 (idx 1,2), wy <-> gamma_3 gamma_1 (idx 2,0),
    wz <-> gamma_1 gamma_2 (idx 0,1)."""
    Omega_full = R.T @ A_full @ R
    return np.array([Omega_full[1, 2], Omega_full[2, 0], Omega_full[0, 1]])

# ═══════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════
def run_and_plot(tau, E1_list, t1_list, labels, out_prefix="e1t1_analysis"):
    n_cases = len(E1_list)
    colors = plt.cm.tab10(np.linspace(0, 1, max(n_cases, 3)))

    # ── Solve all cases ──
    all_data = []
    for i in range(n_cases):
        label_short = labels[i].split('(')[0].strip()
        print(f"  [{i+1}/{n_cases}] {label_short}")
        R_final, t_hist, R_traj = so5_protocol(tau, E1_list[i], t1_list[i], record_traj=True)

        n_pts = len(t_hist)
        bloch_traj = np.zeros((n_pts, 3))
        omega_traj = np.zeros((n_pts, 3))
        phi_traj = np.zeros(n_pts)
        for j in range(n_pts):
            bloch_traj[j] = R_to_bloch_vec(R_traj[j])
            if j < n_pts - 1:
                t = t_hist[j]
                step_idx = min(int(t / tau), 2)
                A_fn = _A_FUNCS[step_idx]
                t_in_step = t % tau if t > 0 else 0.0
                A_full = A_fn(t_in_step, tau, E1_list[i], t1_list[i])
                omega_traj[j] = extract_omega(R_traj[j], A_full)
            _, phi_j = R_to_axis_angle(R_traj[j])
            phi_traj[j] = phi_j

        axis_final, angle_final = R_to_axis_angle(R_final)
        fid = fidelity_double_braid(R_final @ R_final)

        all_data.append({
            't': t_hist, 'bloch': bloch_traj, 'omega': omega_traj,
            'phi': phi_traj, 'v_final': bloch_traj[-1],
            'axis_final': axis_final, 'angle_final': angle_final,
            'fidelity_double': fid
        })

    print("  All cases solved.\n")

    # ═══ Fig 1: Bloch components vs time ═══
    fig1, axes1 = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    comp_lbl = [r'$v_x$ ($\sigma_x$, i$\gamma_2\gamma_3$)',
                r'$v_y$ ($\sigma_y$, i$\gamma_3\gamma_1$)',
                r'$v_z$ ($\sigma_z$, i$\gamma_1\gamma_2$)']
    for ci in range(3):
        ax = axes1[ci]
        for i in range(n_cases):
            ax.plot(all_data[i]['t'], all_data[i]['bloch'][:, ci],
                    color=colors[i], lw=1.2, label=labels[i])
        for s in [tau, 2*tau]:
            ax.axvline(x=s, color='gray', ls='--', lw=0.6, alpha=0.5)
        ax.set_ylabel(comp_lbl[ci], fontsize=10)
        ax.axhline(y=0, color='black', lw=0.3); ax.grid(True, alpha=0.2)
        if ci == 0: ax.legend(fontsize=7, ncol=2)
    axes1[-1].set_xlabel('t')
    fig1.suptitle(f'Bloch Vector Components  (tau={tau})', fontsize=14, fontweight='bold')
    plt.tight_layout(); fig1.savefig(f'{out_prefix}_fig1_bloch_components.png', dpi=200)
    plt.close(fig1); print("  Fig 1: Bloch components saved.")

    # ═══ Fig 2: Omega components ═══
    fig2, axes2 = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    omega_lbl = [r'$\omega_x$ ($\sigma_x$ rotation)',
                 r'$\omega_y$ ($\sigma_y$ rotation)',
                 r'$\omega_z$ ($\sigma_z$ rotation)']
    for ci in range(3):
        ax = axes2[ci]
        for i in range(n_cases):
            skip = max(1, len(all_data[i]['t']) // 500)
            ax.plot(all_data[i]['t'][::skip], all_data[i]['omega'][::skip, ci],
                    color=colors[i], lw=1.2, label=labels[i])
        for s in [tau, 2*tau]:
            ax.axvline(x=s, color='gray', ls='--', lw=0.6, alpha=0.5)
        ax.set_ylabel(omega_lbl[ci], fontsize=10)
        ax.axhline(y=0, color='black', lw=0.3); ax.grid(True, alpha=0.2)
        if ci == 0: ax.legend(fontsize=7, ncol=2)
    for i in range(n_cases):
        d = all_data[i]
        net_wz = np.trapezoid(d['omega'][:, 2], d['t'])
        axes2[2].text(0.98, 0.95 - i*0.08, f'{labels[i]}: net wz = {net_wz:.4f}',
                      transform=axes2[2].transAxes, ha='right', fontsize=7.5,
                      color=colors[i], va='top')
    axes2[-1].set_xlabel('t')
    fig2.suptitle(r'$\omega(t)$: Rotation Generator Components  (tau=' + f'{tau})',
                  fontsize=14, fontweight='bold')
    plt.tight_layout(); fig2.savefig(f'{out_prefix}_fig2_omega_components.png', dpi=200)
    plt.close(fig2); print("  Fig 2: Omega components saved.")

    # ═══ Fig 3: Rotation angle phi(t) ═══
    fig3, ax3 = plt.subplots(figsize=(10, 5))
    for i in range(n_cases):
        ax3.plot(all_data[i]['t'], all_data[i]['phi'],
                 color=colors[i], lw=1.5, label=labels[i])
    for s in [tau, 2*tau]:
        ax3.axvline(x=s, color='gray', ls='--', lw=0.6, alpha=0.5)
    ax3.axhline(y=pi/2, color='green', ls=':', lw=0.8, alpha=0.5,
                label=r'$\pi/2$ (geometric braid)')
    ax3.set_xlabel('t'); ax3.set_ylabel(r'$\phi(t)$ (rad)')
    ax3.set_title(f'Accumulated Rotation Angle  (tau={tau})', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=8); ax3.grid(True, alpha=0.2)
    plt.tight_layout(); fig3.savefig(f'{out_prefix}_fig3_phi_evolution.png', dpi=200)
    plt.close(fig3); print("  Fig 3: Phi evolution saved.")

    # ═══ Fig 4: Bloch sphere 3D ═══
    fig4 = plt.figure(figsize=(16, 7))
    for spi, (va, ts) in enumerate([((25, -50), '3D View'), ((5, -90), 'Side (xz plane)')]):
        ax = fig4.add_subplot(1, 2, spi+1, projection='3d')
        u, v = np.mgrid[0:2*pi:40j, 0:pi:20j]
        ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v),
                          color='gray', alpha=0.06, linewidth=0.2)
        ph = np.linspace(0, 2*pi, 200)
        ax.plot(np.cos(ph), np.sin(ph), 0*ph, 'gray', alpha=0.2, ls='--', lw=0.5)
        for i in range(n_cases):
            d = all_data[i]
            skip = max(1, len(d['t']) // 600)
            tr = d['bloch'][::skip]
            ax.plot(tr[:,0], tr[:,1], tr[:,2], color=colors[i], lw=1.2, label=labels[i])
            ax.scatter(*tr[0], color=colors[i], s=30, marker='o',
                       edgecolors='white', linewidths=0.5, zorder=5)
            ax.scatter(*tr[-1], color=colors[i], s=60, marker='*',
                       edgecolors='white', linewidths=0.5, zorder=5)
        ax.set_xlabel(r'$\sigma_x$'); ax.set_ylabel(r'$\sigma_y$'); ax.set_zlabel(r'$\sigma_z$')
        ax.set_xlim(-1.1,1.1); ax.set_ylim(-1.1,1.1); ax.set_zlim(-1.1,1.1)
        ax.view_init(*va); ax.set_title(f'{ts}  (o start, * end)', fontsize=11)
        if spi == 0: ax.legend(fontsize=7, loc='upper left', bbox_to_anchor=(1.02,1))
    fig4.suptitle(f'Bloch Sphere Trajectories  (tau={tau})', fontsize=14, fontweight='bold')
    plt.tight_layout(); fig4.savefig(f'{out_prefix}_fig4_bloch_3d.png', dpi=200)
    plt.close(fig4); print("  Fig 4: Bloch 3D saved.")

    # ═══ Fig 5: Sweep final state vs lg(t1/E1) ═══
    E1_sw = 0.01; n_sw = 50
    lg_r = np.linspace(-2, 2, n_sw)
    t1_sw = E1_sw * 10**lg_r
    vx_a = np.zeros(n_sw); vy_a = np.zeros(n_sw); vz_a = np.zeros(n_sw)
    phi_a = np.zeros(n_sw); fid_a = np.zeros(n_sw)
    print(f"  Sweeping lg(t1/E1): E1={E1_sw}, {n_sw} points...")
    for j in range(n_sw):
        R = so5_protocol(tau, E1_sw, t1_sw[j])
        v = R_to_bloch_vec(R); vx_a[j], vy_a[j], vz_a[j] = v
        _, phi_a[j] = R_to_axis_angle(R)
        fid_a[j] = fidelity_double_braid(R @ R)
    print("  Sweep done.\n")

    fig5, axes5 = plt.subplots(2, 2, figsize=(14, 10))
    ax = axes5[0,0]
    ax.plot(lg_r, vx_a, 'o-', color='#2196F3', ms=3, lw=1, label=r'$v_x$')
    ax.plot(lg_r, vy_a, 's-', color='#FF5722', ms=3, lw=1, label=r'$v_y$')
    ax.plot(lg_r, vz_a, '^-', color='#4CAF50', ms=3, lw=1, label=r'$v_z$')
    ax.set_xlabel(r'$\lg(t_1/E_1)$'); ax.set_ylabel('Bloch component')
    ax.set_title('Final Bloch Vector'); ax.legend(fontsize=8)
    ax.axhline(y=0, color='gray', lw=0.4); ax.grid(True, alpha=0.3)

    ax = axes5[0,1]
    ax.plot(lg_r, phi_a, 'o-', color='#9C27B0', ms=3, lw=1)
    ax.axhline(y=pi/2, color='green', ls=':', lw=0.8, alpha=0.5, label=r'$\pi/2$')
    ax.set_xlabel(r'$\lg(t_1/E_1)$'); ax.set_ylabel(r'$\phi$ (rad)')
    ax.set_title('Rotation Angle'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    ax = axes5[1,0]
    ax.plot(lg_r, fid_a, 'o-', color='#E91E63', ms=3, lw=1)
    ax.axhline(y=1.0, color='green', ls=':', lw=0.8, alpha=0.5, label='Perfect braid')
    ax.set_xlabel(r'$\lg(t_1/E_1)$'); ax.set_ylabel('Fidelity')
    ax.set_title('Double-Braid Fidelity'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    ax = axes5[1,1]; ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    ax.text(0.5, 0.80, 'Regime interpretation:', ha='center', fontsize=12, fontweight='bold', transform=ax.transAxes)
    ax.text(0.5, 0.62, r'$t_1 \ll E_1$: E1-dominated, axis near $\sigma_z$', ha='center', fontsize=10, transform=ax.transAxes, color='#2196F3')
    ax.text(0.5, 0.48, r'$t_1 \approx E_1$: competing axes, rich interference', ha='center', fontsize=10, transform=ax.transAxes, color='#9C27B0')
    ax.text(0.5, 0.34, r'$t_1 \gg E_1$: t1-dominated, fast $\sigma_y$ precession', ha='center', fontsize=10, transform=ax.transAxes, color='#FF5722')
    ax.text(0.5, 0.18, f'E1={E1_sw} meV, tau={tau}, tc=E0={tc} meV', ha='center', fontsize=9, transform=ax.transAxes, color='gray')

    fig5.suptitle(f'Parametric Sweep vs t1/E1  (E1={E1_sw}, tau={tau})', fontsize=14, fontweight='bold')
    plt.tight_layout(); fig5.savefig(f'{out_prefix}_fig5_sweep.png', dpi=200)
    plt.close(fig5); print("  Fig 5: Sweep saved.")

    # ═══ Summary table ═══
    print("=" * 80)
    print(f"{'Case':35s} {'axis_x':>7s} {'axis_y':>7s} {'axis_z':>7s} {'angle':>7s} {'fidelity':>9s}")
    print("-" * 80)
    for i in range(n_cases):
        d = all_data[i]; ax_v = d['axis_final']
        print(f"{labels[i]:35s} {ax_v[0]:+7.4f} {ax_v[1]:+7.4f} {ax_v[2]:+7.4f} "
              f"{d['angle_final']:7.4f} {d['fidelity_double']:9.6f}")
    print("=" * 80)
    print(f"\nPlots: {out_prefix}_fig*.png")
    return all_data

# ═══════════════════════════════════════════════════════════════
# Multi-tau sweep
# ═══════════════════════════════════════════════════════════════
def multitau_sweep(tau_vals, E1v, out_prefix):
    n_sw = 30
    lg_grid = np.linspace(-2, 2, n_sw)
    F_map = np.zeros((len(tau_vals), n_sw))
    print(f"\n  Multi-tau: {len(tau_vals)}x{n_sw} grid...")
    for ti, t_val in enumerate(tau_vals):
        for j in range(n_sw):
            R = so5_protocol(t_val, E1v, E1v * 10**lg_grid[j])
            F_map[ti, j] = fidelity_double_braid(R @ R)
        print(f"    tau={t_val:.1f} done")
    print("  Heatmap done.\n")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: line plots
    for ti, t_val in enumerate(tau_vals):
        color = plt.cm.viridis(ti / max(1, len(tau_vals)-1))
        fid_line = np.array([fidelity_double_braid(
            so5_protocol(t_val, E1v, E1v * 10**lg) @ so5_protocol(t_val, E1v, E1v * 10**lg))
            for lg in lg_grid])
        axes[0,0].plot(lg_grid, fid_line, '-', color=color, lw=1.5, label=f'tau={t_val}')
    axes[0,0].set_xlabel(r'$\lg(t_1/E_1)$'); axes[0,0].set_ylabel('Fidelity')
    axes[0,0].set_title('Fidelity vs t1/E1 (all tau)'); axes[0,0].legend(fontsize=7); axes[0,0].grid(True, alpha=0.3)

    # Panel 2: heatmap
    im = axes[0,1].pcolormesh(lg_grid, tau_vals, F_map, shading='auto', cmap='viridis')
    axes[0,1].set_xlabel(r'$\lg(t_1/E_1)$'); axes[0,1].set_ylabel(r'$\tau$')
    axes[0,1].set_title('Fidelity Heatmap')
    plt.colorbar(im, ax=axes[0,1], label='Fidelity')

    # Panel 3: contours (interpolated)
    try:
        from scipy.ndimage import zoom
        Fz = zoom(F_map, 3, order=1)
        tau_z = np.linspace(tau_vals[0], tau_vals[-1], len(tau_vals)*3)
        lg_z = np.linspace(lg_grid[0], lg_grid[-1], n_sw*3)
        axes[1,0].contourf(lg_z, tau_z, Fz, levels=12, cmap='viridis')
        cs = axes[1,0].contour(lg_z, tau_z, Fz, levels=[0.3,0.5,0.7,0.9],
                               colors='white', linewidths=0.5)
        axes[1,0].clabel(cs, inline=True, fontsize=7, fmt='%.1f')
    except Exception:
        axes[1,0].text(0.5, 0.5, 'Insufficient data', ha='center', transform=axes[1,0].transAxes)
    axes[1,0].set_xlabel(r'$\lg(t_1/E_1)$'); axes[1,0].set_ylabel(r'$\tau$')
    axes[1,0].set_title('Fidelity Contours (interpolated)')

    # Panel 4: annotation
    ax = axes[1,1]; ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    ax.text(0.5, 0.85, 'Multi-tau Analysis', ha='center', fontsize=13, fontweight='bold', transform=ax.transAxes)
    ax.text(0.5, 0.68, 'Small tau: non-adiabatic, fidelity low', ha='center', fontsize=10, transform=ax.transAxes, color='#2196F3')
    ax.text(0.5, 0.52, 'Large tau: adiabatic, geometric braid dominates', ha='center', fontsize=10, transform=ax.transAxes, color='#4CAF50')
    ax.text(0.5, 0.36, 'Diagonal stripes = E1-t1 interference', ha='center', fontsize=10, transform=ax.transAxes, color='#FF5722')
    ax.text(0.5, 0.18, f'E1={E1v} meV, tc=E0={tc} meV', ha='center', fontsize=9, transform=ax.transAxes, color='gray')

    fig.suptitle(f'Multi-tau Fidelity Landscape  (E1={E1v} meV)', fontsize=14, fontweight='bold')
    plt.tight_layout(); fig.savefig(f'{out_prefix}_fig6_multitau.png', dpi=200)
    plt.close(fig); print("  Fig 6: Multi-tau saved.")

# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='E1-t1 analysis via SO(5) RK4')
    parser.add_argument('--tau', type=float, default=50.0)
    parser.add_argument('--E1', type=float, default=0.01)
    parser.add_argument('--prefix', type=str, default='e1t1_analysis')
    parser.add_argument('--tau-list', type=str, default=None,
                        help='Comma-separated tau values, e.g. 1,5,20,50')
    args = parser.parse_args()

    tau = args.tau; E1v = args.E1

    cases = [
        (0.0, 0.0,     'E1=0, t1=0 (pure MZM)'),
        (0.0, 0.01,    'E1=0, t1=0.01 (E1=0, ABS phase only)'),
        (E1v, 0.0,     f'E1={E1v}, t1=0 (pure E1)'),
        (E1v, E1v,     f'E1=t1={E1v} (equal coupling)'),
        (E1v, E1v*0.1, f'E1={E1v}, t1={E1v*0.1:.4f} (t1 << E1)'),
        (E1v, E1v*10,  f'E1={E1v}, t1={E1v*10:.2f} (t1 >> E1)'),
    ]

    print(f"\n{'='*80}")
    print(f"E1-t1 Interaction Analysis  |  SO(5) RK4  |  tau={tau}  |  tc=E0={tc}")
    print(f"{'='*80}\n")

    run_and_plot(tau, [c[0] for c in cases], [c[1] for c in cases],
                 [c[2] for c in cases], out_prefix=args.prefix)

    if args.tau_list:
        tvs = [float(x) for x in args.tau_list.split(',')]
        multitau_sweep(tvs, E1v, args.prefix)

    print("\nDone.")
