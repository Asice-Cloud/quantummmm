#!/usr/bin/env python3
"""
E1 从 0 增大：和乐群从 SO(2) → SU(2) 的过渡
================================================
理论问题：E1=0 时 A=D → ∫ω_z=0 → 和乐限制在 SO(2)（xy 平面旋转）
         E1≠0 时 A≠D → ∫ω_z≠0 → 和乐扩展为 SU(2)

关键问题：这是突变还是渐变？σ_z 分量如何随 E1 增长？

诊断量：
  n̂_z = 旋转轴在 σ_z 方向的分量（=0 为纯 SO(2)，≠0 为 SU(2)）
  φ   = 总旋转角
  ∫ω_z dt = 动力学 σ_z 相位累积

用法：
  python e1_so2_to_su2.py
"""
import numpy as np; pi = np.pi
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

tc=0.3; E0=0.3

def fp(t,tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t,tau): return 0.5*(1-np.cos(pi*t/tau))

def A1(t,tau,e,t1c):
    A=np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
    A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
    A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau); return A
def A2(t,tau,e,t1c):
    A=np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3]=-2*tc*fp(t,tau); A[3,1]=2*tc*fp(t,tau)
    A[2,3]=2*tc*fm(t,tau); A[3,2]=-2*tc*fm(t,tau)
    A[0,4]=2*t1c*fp(t,tau); A[4,0]=-2*t1c*fp(t,tau); return A
def A3(t,tau,e,t1c):
    A=np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[2,3]=2*tc*fp(t,tau); A[3,2]=-2*tc*fp(t,tau)
    A[3,4]=2*E0*fm(t,tau); A[4,3]=-2*E0*fm(t,tau); return A

_A_FUNCS = [A1, A2, A3]

def so5_propagate_with_traj(tau, e, t1c, n=600):
    """SO(5) propagation with omega(t) trajectory recording."""
    dt = tau / n
    R = np.eye(5)
    n_total = 3*n + 1
    t_hist = np.zeros(n_total)
    omega_hist = np.zeros((n_total, 3))
    t_hist[0] = 0.0; omega_hist[0] = [0,0,0]
    idx = 1

    for step_idx, A_fn in enumerate(_A_FUNCS):
        for s in range(n):
            t = s * dt
            A0 = A_fn(t, tau, e, t1c)
            k1 = A0 @ R
            A_h = A_fn(t+0.5*dt, tau, e, t1c)
            k2 = A_h @ (R+0.5*dt*k1)
            k3 = A_h @ (R+0.5*dt*k2)
            A_f = A_fn(t+dt, tau, e, t1c)
            k4 = A_f @ (R+dt*k3)
            R = R + dt/6*(k1+2*k2+2*k3+k4)

            # Extract omega
            Omega = R.T @ A_f @ R
            t_hist[idx] = step_idx*tau + (s+1)*dt
            omega_hist[idx] = [Omega[1,2], Omega[2,0], Omega[0,1]]
            idx += 1

    return R, t_hist, omega_hist

def axis_angle(R):
    """从 SO(5) 的 {γ1,γ2,γ3} 子块提取旋转轴和角"""
    R3 = R[:3,:3]
    cos_phi = np.clip((np.trace(R3)-1)/2, -1, 1)
    phi = np.arccos(cos_phi)
    if phi < 1e-12: return np.array([0.,0.,1.]), 0.
    n = np.array([R3[2,1]-R3[1,2], R3[0,2]-R3[2,0], R3[1,0]-R3[0,1]])
    nn = np.linalg.norm(n)
    return (n/nn, phi) if nn > 1e-12 else (np.array([0.,0.,1.]), phi)

# ═══════════════════════════════════════════════════════════════
# 主扫参：E1 from 0 to Emax
# ═══════════════════════════════════════════════════════════════
def sweep_E1(tau=50.0, t1c=0.01, Emax=0.05, n_pts=40):
    """扫 E1 ∈ [0, Emax]，提取 n̂_z, φ, ∫ω_z dt"""
    E1_vals = np.linspace(0, Emax, n_pts)
    # 在 E1=0 附近加密
    E1_vals = np.sort(np.concatenate([
        np.linspace(0, 0.005, 15), np.linspace(0.005, Emax, n_pts-15)
    ]))
    E1_vals = np.unique(E1_vals)
    n_pts = len(E1_vals)

    nz_arr = np.zeros(n_pts)
    phi_arr = np.zeros(n_pts)
    wz_int_arr = np.zeros(n_pts)
    fid_arr = np.zeros(n_pts)

    print(f"  Sweeping E1: 0 → {Emax}, {n_pts} points, tau={tau}, t1={t1c}")
    for i, e in enumerate(E1_vals):
        R, t_hist, omega_hist = so5_propagate_with_traj(tau, e, t1c)
        n_hat, phi = axis_angle(R)
        # n_hat from axis_angle: need to ensure consistent sign
        # We want n_hat_z as a signed quantity
        nz_arr[i] = n_hat[2]
        phi_arr[i] = phi
        wz_int_arr[i] = np.trapezoid(omega_hist[:,2], t_hist)
        # Fidelity after double braid
        Rd = R @ R
        ov = 0.5*(Rd[0,0]+1j*Rd[1,0]+1j*Rd[0,1]-Rd[1,1])
        fid_arr[i] = np.abs(ov)**2
        if (i+1) % 10 == 0: print(f"    E1={e:.4f}: nz={n_hat[2]:+.4f}, phi={phi:.4f}, fid={fid_arr[i]:.4f}")

    return E1_vals, nz_arr, phi_arr, wz_int_arr, fid_arr

# ═══════════════════════════════════════════════════════════════
# 画图
# ═══════════════════════════════════════════════════════════════
print("="*60)
print("E1: SO(2) → SU(2) Transition")
print("="*60)

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# 多组 tau 对比
tau_list = [100, 50, 20, 10]
colors = plt.cm.viridis(np.linspace(0, 1, len(tau_list)))

for ti, tau in enumerate(tau_list):
    print(f"\n  tau = {tau}")
    E1v, nz, phi, wz_int, fid = sweep_E1(tau=tau, t1c=0.01, Emax=0.05)
    c = colors[ti]; label = f'τ={tau}'

    # Panel 1: n̂_z vs E1 — 核心诊断量
    ax = axes[0,0]
    ax.plot(E1v, np.abs(nz), 'o-', color=c, ms=3, lw=1.5, label=label)
    ax.set_xlabel('E₁ (meV)'); ax.set_ylabel(r'$|\hat{n}_z|$')
    ax.set_title(r'Rotation axis $\sigma_z$ component (=0 for SO(2), ≠0 for SU(2))', fontsize=10)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # Panel 2: ∫ω_z dt vs E1
    ax = axes[0,1]
    ax.plot(E1v, np.abs(wz_int), 'o-', color=c, ms=3, lw=1.5, label=label)
    ax.set_xlabel('E₁ (meV)'); ax.set_ylabel(r'$|\int \omega_z\,dt|$')
    ax.set_title(r'Net $\sigma_z$ phase accumulation', fontsize=10)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # Panel 3: φ vs E1
    ax = axes[0,2]
    ax.plot(E1v, phi, 'o-', color=c, ms=3, lw=1.5, label=label)
    ax.axhline(y=pi/2, color='green', ls=':', lw=0.8, alpha=0.5)
    ax.set_xlabel('E₁ (meV)'); ax.set_ylabel(r'$\phi$ (rad)')
    ax.set_title('Rotation angle (π/2 = pure geometric)', fontsize=10)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # Panel 4: fidelity vs E1
    ax = axes[1,0]
    ax.plot(E1v, fid, 'o-', color=c, ms=3, lw=1.5, label=label)
    ax.set_xlabel('E₁ (meV)'); ax.set_ylabel('Fidelity')
    ax.set_title('Double-braid fidelity', fontsize=10)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # Panel 5: |n̂_z| vs φ (parametric plot)
    ax = axes[1,1]
    ax.plot(np.abs(nz), phi, 'o-', color=c, ms=3, lw=1.5, label=label)
    ax.axhline(y=pi/2, color='green', ls=':', lw=0.8, alpha=0.5)
    ax.set_xlabel(r'$|\hat{n}_z|$'); ax.set_ylabel(r'$\phi$ (rad)')
    ax.set_title(r'$\phi$ vs $|\hat{n}_z|$ (SU(2) distance from SO(2))', fontsize=10)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # Panel 6: log-log to see scaling at small E1
    ax = axes[1,2]
    mask = E1v > 0
    ax.loglog(E1v[mask], np.abs(nz[mask]), 'o-', color=c, ms=3, lw=1.5, label=label)
    ax.set_xlabel('E₁ (meV)'); ax.set_ylabel(r'$|\hat{n}_z|$')
    ax.set_title(r'Log-log: $|\hat{n}_z| \propto E_1^p$', fontsize=10)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

# Annotate the transition
axes[0,0].text(0.98, 0.95, 'SO(2): n̂_z = 0\n  ↓\nSU(2): n̂_z ≠ 0',
               transform=axes[0,0].transAxes, ha='right', va='top',
               fontsize=9, bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

fig.suptitle('Holonomy Group Transition: SO(2) → SU(2) as E₁ increases',
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig('e1_so2_to_su2.png', dpi=200)
plt.close(fig)
print("\n  Saved: e1_so2_to_su2.png")

# ═══════════════════════════════════════════════════════════════
# Figure 2: omega_z(t) for selected E1 values
# ═══════════════════════════════════════════════════════════════
print("\n  Plotting omega_z(t) at selected E1...")
tau_demo = 50.0
E1_demo = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05]

fig2, axes2 = plt.subplots(2, 1, figsize=(14, 8))

for e_val in E1_demo:
    R, t_hist, omega_hist = so5_propagate_with_traj(tau_demo, e_val, 0.01)
    wz_cum = np.cumsum(omega_hist[:,2]) * (t_hist[1]-t_hist[0])

    axes2[0].plot(t_hist, omega_hist[:,2], lw=1, alpha=0.8, label=f'E₁={e_val}')
    axes2[1].plot(t_hist, wz_cum, lw=1.2, alpha=0.8, label=f'E₁={e_val}')

for ax in axes2:
    for s in [tau_demo, 2*tau_demo]:
        ax.axvline(x=s, color='gray', ls='--', lw=0.5, alpha=0.4)
    ax.axhline(y=0, color='black', lw=0.3)

axes2[0].set_ylabel(r'$\omega_z(t)$'); axes2[0].set_title('Instantaneous σ_z rotation rate')
axes2[0].legend(fontsize=7, ncol=3)
axes2[1].set_xlabel('t'); axes2[1].set_ylabel(r'$\int_0^t \omega_z\,dt$')
axes2[1].set_title('Cumulative σ_z phase')
axes2[1].legend(fontsize=7, ncol=3)

# Mark A=D symmetry point
axes2[1].text(0.98, 0.05, 'E₁=0: ∫ω_z=0 (A=D)\nE₁≠0: ∫ω_z≠0 (A≠D)',
              transform=axes2[1].transAxes, ha='right', va='bottom',
              fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

fig2.suptitle(r'$\omega_z(t)$ Evolution: How E₁ Breaks the A=D Symmetry',
              fontsize=14, fontweight='bold')
plt.tight_layout()
fig2.savefig('e1_omega_z_evolution.png', dpi=200)
plt.close(fig2)
print("  Saved: e1_omega_z_evolution.png")

print("\nDone.")
