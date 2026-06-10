#!/usr/bin/env python3
"""
Fig 1(d) 最终修正版 —— 修复 τ 单位问题
=========================================
论文 τ 轴单位: 100/meV
即 τ_plot=1 → 物理时间 = 100 meV⁻¹

代码中 τ_code = τ_plot × 100
扫描范围 τ_plot ∈ [0.2, 12], τ_code ∈ [20, 1200]
"""
import numpy as np
from scipy.ndimage import gaussian_filter, zoom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3; E1_fixed = 0.01

def fp(t, tau):  return 0.5 * (1.0 + np.cos(pi * t / tau))
def fm(t, tau):  return 0.5 * (1.0 - np.cos(pi * t / tau))

def build_A_step1(t, tau, e, t1c):
    """Step 1: t₂ ascends, E_d descends, t₃=0."""
    A = np.zeros((5,5))
    A[0,1] = 2*e;  A[1,0] = -2*e
    t2v = tc * fm(t, tau);  A[1,3] = -2*t2v; A[3,1] = 2*t2v
    t1v = t1c * fm(t, tau); A[0,4] = 2*t1v; A[4,0] = -2*t1v
    Edv = E0 * fp(t, tau);  A[3,4] = 2*Edv; A[4,3] = -2*Edv
    return A

def build_A_step2(t, tau, e, t1c):
    """Step 2: t₂ descends, t₃ ascends, E_d=0."""
    A = np.zeros((5,5))
    A[0,1] = 2*e;  A[1,0] = -2*e
    t2v = tc * fp(t, tau);  A[1,3] = -2*t2v; A[3,1] = 2*t2v
    t3v = tc * fm(t, tau);  A[2,3] = 2*t3v; A[3,2] = -2*t3v
    t1v = t1c * fp(t, tau); A[0,4] = 2*t1v; A[4,0] = -2*t1v
    return A

def build_A_step3(t, tau, e, t1c):
    """Step 3: t₃ descends, E_d ascends, t₂=0. t₁=0 (G₁ on)."""
    A = np.zeros((5,5))
    A[0,1] = 2*e;  A[1,0] = -2*e
    t3v = tc * fp(t, tau);  A[2,3] = 2*t3v; A[3,2] = -2*t3v
    Edv = E0 * fm(t, tau);  A[3,4] = 2*Edv; A[4,3] = -2*Edv
    return A

def propagate_step(bld, tau, e, t1c):
    """RK4 with τ-scaled steps. n ≈ 17*τ ensures θ_step ≈ 0.035 rad."""
    n = max(300, int(17 * tau))
    dt = tau / n
    R = np.eye(5)
    for s in range(n):
        t = s * dt
        At = bld(t, tau, e, t1c)
        k1 = At @ R
        k2 = bld(t + 0.5*dt, tau, e, t1c) @ (R + 0.5*dt*k1)
        k3 = bld(t + 0.5*dt, tau, e, t1c) @ (R + 0.5*dt*k2)
        k4 = bld(t + dt, tau, e, t1c) @ (R + dt*k3)
        R += dt/6.0 * (k1 + 2*k2 + 2*k3 + k4)
    return R

def fidelity_double_swap(R_double):
    """|⟨ψ⁻|U²|ψ⁺⟩|²"""
    ov = 0.5 * (R_double[0,0] + 1j*R_double[1,0] 
                + 1j*R_double[0,1] - R_double[1,1])
    return np.abs(ov)**2

# ═══════════════════════════════════════════════════════════════════
# 扫描参数
# ═══════════════════════════════════════════════════════════════════
# 论文 τ 轴: τ_plot ∈ [0.2, 12], 单位 100/meV
# 代码 τ: τ_code = τ_plot × 100
tau_plot_min, tau_plot_max = 0.2, 12.0
lg_min, lg_max = -1.0, 1.0

N_TAU = 120   # τ 方向网格数
N_T1  = 150   # t₁ 方向网格数

tau_plot_vals = np.linspace(tau_plot_min, tau_plot_max, N_TAU)
tau_code_vals = tau_plot_vals * 100.0  # 单位转换!
t1_vals = E1_fixed * 10**np.linspace(lg_min, lg_max, N_T1)

print("=" * 55)
print("Fig 1(d) 最终修正版")
print("=" * 55)
print(f"E₁={E1_fixed} meV, t_c=E₀={tc} meV")
print(f"网格: {N_TAU}×{N_T1} = {N_TAU*N_T1}")
print(f"τ_plot ∈ [{tau_plot_min}, {tau_plot_max}]")
print(f"τ_code ∈ [{tau_code_vals[0]:.0f}, {tau_code_vals[-1]:.0f}]")
print(f"t₁ ∈ [{t1_vals[0]:.4f}, {t1_vals[-1]:.4f}] meV")
print()

fid_grid = np.zeros((N_T1, N_TAU))
report_step = max(1, N_TAU // 10)

for i in range(N_TAU):
    tau_code = tau_code_vals[i]
    for j in range(N_T1):
        t1 = t1_vals[j]
        R1 = propagate_step(build_A_step1, tau_code, E1_fixed, t1)
        R2 = propagate_step(build_A_step2, tau_code, E1_fixed, t1)
        R3 = propagate_step(build_A_step3, tau_code, E1_fixed, t1)
        Rs = R3 @ R2 @ R1
        Rd = Rs @ Rs
        fid_grid[j, i] = fidelity_double_swap(Rd)
    if (i+1) % report_step == 0 or i == 0:
        print(f"  {i+1}/{N_TAU} ({100*(i+1)//N_TAU}%)")

print(f"\nRaw range: [{fid_grid.min():.4f}, {fid_grid.max():.4f}]")

# ═══════════════════════════════════════════════════════════════════
# 平滑 + 绘图
# ═══════════════════════════════════════════════════════════════════
fid_s = gaussian_filter(fid_grid, sigma=0.8, mode='nearest')
fid_z = zoom(fid_s, 3, order=3)

tau_z = np.linspace(tau_plot_min, tau_plot_max, N_TAU*3)
lg_z  = np.linspace(lg_min, lg_max, N_T1*3)
Tz, Lz = np.meshgrid(tau_z, lg_z)

print(f"Upsampled range: [{fid_z.min():.4f}, {fid_z.max():.4f}]")

fig, ax = plt.subplots(1, 1, figsize=(9, 6.5))

levels = np.linspace(0, 1, 21)
cf = ax.contourf(Tz, Lz, fid_z, levels=levels, cmap='inferno', extend='both')

# 等高线
cs = ax.contour(Tz, Lz, fid_z, levels=np.linspace(0, 1, 11),
                colors='w', linewidths=0.35, alpha=0.3)
cs2 = ax.contour(Tz, Lz, fid_z, levels=[0.2, 0.4, 0.6, 0.8],
                 colors='w', linewidths=0.7, alpha=0.5)
ax.clabel(cs2, inline=True, fontsize=7, fmt='%.1f')

ax.set_xlabel(r'$\tau$ (100/meV)', fontsize=14)
ax.set_ylabel(r'$\lg(t_1/E_1)$', fontsize=14)
ax.set_title(r'$|\langle\psi_1^-|U(6\tau)|\psi_1^+\rangle|^2$'
             + f'  ($E_1={E1_fixed}$ meV)',
             fontsize=13)

cbar = plt.colorbar(cf, ax=ax, label='Fidelity', ticks=np.linspace(0, 1, 6))

ax.text(0.98, 0.02,
        f'$E_1={E1_fixed}$ meV\n$t_c=E_0={tc}$ meV\n'
        f'double swap\n$\\tau_{{\\rm code}} = \\tau_{{\\rm plot}} \\times 100$',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=8, bbox=dict(boxstyle='round', facecolor='black', alpha=0.7),
        color='white')

plt.tight_layout()
plt.savefig('fig1d_reproduction_fixed.png', dpi=250)
print(f"\nSaved: fig1d_reproduction_fixed.png")
