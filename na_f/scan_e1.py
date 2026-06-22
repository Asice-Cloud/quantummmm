#!/usr/bin/env python3
"""E₁ 扫描：研究 E₁ 从 0 增大时 Bloch 矢量的变化"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3

def fp(t, tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1-np.cos(pi*t/tau))

def step_matrix(t, tau, E1, t1c, step):
    A = np.zeros((5,5))
    A[0,1]=2*E1; A[1,0]=-2*E1
    if step == 1:
        A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
        A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
        A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau)
    elif step == 2:
        A[1,3]=-2*tc*fp(t,tau); A[3,1]=2*tc*fp(t,tau)
        A[2,3]=2*tc*fm(t,tau); A[3,2]=-2*tc*fm(t,tau)
        A[0,4]=2*t1c*fp(t,tau); A[4,0]=-2*t1c*fp(t,tau)
    else:
        A[2,3]=2*tc*fp(t,tau); A[3,2]=-2*tc*fp(t,tau)
        A[3,4]=2*E0*fm(t,tau); A[4,3]=-2*E0*fm(t,tau)
    return A

def propagate_step(tau_c, E1, t1c, step):
    """自适应步数：确保每 E₁ 周期至少 20 步"""
    omega_max = max(2*E1, 2*tc, 2*E0, 2*abs(t1c)) + 1e-6
    period = 2*pi/omega_max
    n_pts = max(300, int(tau_c/period * 20))
    dt = tau_c/n_pts
    R = np.eye(5)
    for s in range(n_pts):
        t = s*dt
        A = step_matrix(t, tau_c, E1, t1c, step)
        # RK4
        k1 = A @ R
        k2 = step_matrix(t+0.5*dt, tau_c, E1, t1c, step) @ (R+0.5*dt*k1)
        k3 = step_matrix(t+0.5*dt, tau_c, E1, t1c, step) @ (R+0.5*dt*k2)
        k4 = step_matrix(t+dt, tau_c, E1, t1c, step) @ (R+dt*k3)
        R = R + dt/6*(k1+2*k2+2*k3+k4)
    return R

def full_braid(tau_p, E1, t1c):
    tau_c = tau_p*100
    R = np.eye(5)
    for step in [1,2,3]:
        R = propagate_step(tau_c, E1, t1c, step)
    return R

def R_to_bloch_axis(R):
    """从 MZM 子块提取 Bloch 矢量（n̂ 方向，未归一化）"""
    x = 1j*(R[0,1]-R[1,0]); y = R[0,1]+R[1,0]; z = R[0,0]-R[1,1]
    return np.real(np.array([x, y, z]))

def fidelity(R):
    """双次 swap 保真度 |⟨ψ₁⁻|U|ψ₁⁺⟩|²"""
    return abs(0.5*(R[0,0] + 1j*R[1,0] + 1j*R[0,1] - R[1,1]))**2

# ═══════════════════════════════════════════════════════════
# 扫描参数
E1_vals = np.concatenate([
    [0.0],
    np.logspace(-3, np.log10(0.3), 20),  # 0.001 → 0.3
    np.linspace(0.35, 0.5, 5),           # 0.35 → 0.5
])
tau_vals = [10, 30]              # 非绝热 vs 中间
t1_vals  = [0.0, 0.01, 0.05]     # 不同 t₁ 强度

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('E₁ scan: Bloch vector components $(v_x, v_y, v_z)$ vs $E_1$ (meV)', fontsize=13)

for row, tau in enumerate(tau_vals):
    for col, t1 in enumerate(t1_vals):
        ax = axes[row, col]
        vx_list, vy_list, vz_list, fid_list = [], [], [], []
        
        for E1 in E1_vals:
            R = full_braid(tau, E1, t1)
            v_raw = R_to_bloch_axis(R)
            nrm = np.linalg.norm(v_raw)
            if nrm > 1e-12:
                vx_list.append(v_raw[0]/nrm)
                vy_list.append(v_raw[1]/nrm)
                vz_list.append(v_raw[2]/nrm)
            else:
                vx_list.append(0); vy_list.append(0); vz_list.append(0)
            fid_list.append(fidelity(R))
        
        # 画 v_x, v_y, v_z
        ax.plot(E1_vals, vx_list, 'o-', ms=3, lw=1.2, label='$v_x$', alpha=0.8)
        ax.plot(E1_vals, vy_list, 's-', ms=3, lw=1.2, label='$v_y$ (σ_z flux)', alpha=0.8)
        ax.plot(E1_vals, vz_list, '^-', ms=3, lw=1.2, label='$v_z$', alpha=0.8)
        
        # 标记 t_c = 0.3
        ax.axvline(x=0.3, color='gray', ls='--', lw=0.8, alpha=0.5)
        ax.text(0.3, ax.get_ylim()[1]*0.95, '$t_c$', ha='center', fontsize=8, color='gray')
        
        ax.set_xscale('symlog', linthresh=0.005)
        ax.set_xlabel('$E_1$ (meV)')
        ax.set_ylabel('Bloch component')
        ax.set_title(f'τ={tau}, t₁={t1}')
        ax.legend(fontsize=7)
        ax.set_ylim(-1.1, 1.1)
        ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('na_f/fig_e1_scan_bloch.png', dpi=200)
print("✓ na_f/fig_e1_scan_bloch.png")

# ═══════════════════════════════════════════════════════════
# 第二张图：fidelity 和 v_y 的焦点图
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))

for ax, tau in zip(axes2, [15, 30]):
    for t1, style, color in zip([0.0, 0.01, 0.05],
                                 ['-', '--', ':'],
                                 ['#1f77b4', '#ff7f0e', '#2ca02c']):
        vy_list, fid_list = [], []
        for E1 in E1_vals:
            R = full_braid(tau, E1, t1)
            v_raw = R_to_bloch_axis(R)
            nrm = np.linalg.norm(v_raw)
            vy = v_raw[1]/nrm if nrm > 1e-12 else 0
            vy_list.append(vy)
            fid_list.append(fidelity(R))
        ax.plot(E1_vals, vy_list, style, lw=1.5, color=color,
                label=f't₁={t1}')
    ax.axvline(x=0.3, color='gray', ls='--', lw=0.8)
    ax.set_xscale('symlog', linthresh=0.005)
    ax.set_xlabel('$E_1$ (meV)')
    ax.set_ylabel('$v_y$ (σ_z flux indicator)')
    ax.set_title(f'$v_y$ vs $E_1$, τ={tau}')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('na_f/fig_e1_scan_vy.png', dpi=200)
print("✓ na_f/fig_e1_scan_vy.png")

# ═══════════════════════════════════════════════════════════
# 第三张：标记关键区域
print("\n=== Key values ===")
for tau in [10, 30]:
    print(f"\nτ={tau}:")
    print(f"  E₁=0: ", end="")
    R0 = full_braid(tau, 0.0, 0.01)
    v0 = R_to_bloch_axis(R0); v0 /= np.linalg.norm(v0)
    print(f"v=({v0[0]:.4f}, {v0[1]:.4f}, {v0[2]:.4f}), fid={fidelity(R0):.4f}")
    
    for E1_label, E1_val in [("0.01", 0.01), ("0.1", 0.1), ("0.3", 0.3), ("0.5", 0.5)]:
        R = full_braid(tau, E1_val, 0.01)
        v = R_to_bloch_axis(R); v /= np.linalg.norm(v)
        print(f"  E₁={E1_label}: v=({v[0]:.4f}, {v[1]:.4f}, {v[2]:.4f}), fid={fidelity(R):.4f}")
