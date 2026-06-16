#!/usr/bin/env python3
"""
Fig 1(d) 最终版 — Bug 已修复，高效扫描
=========================================
Bug 修复: 变量 t 被覆写导致 Ed 计算错误
优化: 自适应 RK4 步数, 80×100 网格 + Gaussian 平滑
"""
import numpy as np
from scipy.ndimage import gaussian_filter, zoom
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3; E1_fixed = 0.01  # 论文图注值，高分辨率验证

def fp(t, tau):  return 0.5 * (1.0 + np.cos(pi * t / tau))
def fm(t, tau):  return 0.5 * (1.0 - np.cos(pi * t / tau))

# ── 修正版: 使用独立变量名避免覆写 ──
def b1(t, tau, e, t1c):  # Step 1
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    t2v = tc*fm(t, tau); A[1,3] = -2*t2v; A[3,1] = 2*t2v
    t1v = t1c*fm(t, tau); A[0,4] = 2*t1v; A[4,0] = -2*t1v
    Edv = E0*fp(t, tau); A[3,4] = 2*Edv; A[4,3] = -2*Edv; return A

def b2(t, tau, e, t1c):  # Step 2
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    t2v = tc*fp(t, tau); A[1,3] = -2*t2v; A[3,1] = 2*t2v
    t3v = tc*fm(t, tau); A[2,3] = 2*t3v; A[3,2] = -2*t3v
    t1v = t1c*fp(t, tau); A[0,4] = 2*t1v; A[4,0] = -2*t1v; return A

def b3(t, tau, e, t1c):  # Step 3 (t1=0, G1 on)
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    t3v = tc*fp(t, tau); A[2,3] = 2*t3v; A[3,2] = -2*t3v
    Edv = E0*fm(t, tau); A[3,4] = 2*Edv; A[4,3] = -2*Edv; return A

def prop(bld, tau, e, t1c):
    """n = max(500, 2τ) 确保 θ_step ≈ 0.15 rad (已验证收敛)"""
    n = max(500, int(2 * tau))
    dt = tau / n; R = np.eye(5)
    for s in range(n):
        t = s * dt
        k1 = bld(t,tau,e,t1c) @ R
        k2 = bld(t+0.5*dt,tau,e,t1c) @ (R+0.5*dt*k1)
        k3 = bld(t+0.5*dt,tau,e,t1c) @ (R+0.5*dt*k2)
        k4 = bld(t+dt,tau,e,t1c) @ (R+dt*k3)
        R += dt/6.0 * (k1 + 2*k2 + 2*k3 + k4)
    return R

def fid(R):
    ov = 0.5*(R[0,0]+1j*R[1,0]+1j*R[0,1]-R[1,1])
    return np.abs(ov)**2

# ═══════════════════════════════════════════════════════════
# 扫描
# ═══════════════════════════════════════════════════════════
N_TAU, N_T1 = 120, 80  # 高分辨率 τ 扫，捕捉窄带
tau_p = np.linspace(0.2, 12.0, N_TAU)
tau_c = tau_p * 100.0  # 单位转换: τ(100/meV) → τ(meV⁻¹)
t1_v  = E1_fixed * 10**np.linspace(-1.0, 1.0, N_T1)

print(f"Fig 1(d) 修正版 | E₁={E1_fixed} | 网格={N_TAU}×{N_T1}")
print(f"  τ_code ∈ [{tau_c[0]:.0f},{tau_c[-1]:.0f}]")

F = np.zeros((N_T1, N_TAU))
for i in range(N_TAU):
    for j in range(N_T1):
        R1=prop(b1,tau_c[i],E1_fixed,t1_v[j])
        R2=prop(b2,tau_c[i],E1_fixed,t1_v[j])
        R3=prop(b3,tau_c[i],E1_fixed,t1_v[j])
        F[j,i] = fid((R3@R2@R1)@(R3@R2@R1))
    if (i+1)%10==0: print(f"  {i+1}/{N_TAU}")

print(f"Raw range: [{F.min():.4f}, {F.max():.4f}]")

# ── 轻平滑 + 超采样 ──
F_s = gaussian_filter(F, sigma=0.3, mode='nearest')  # 减轻平滑，保留窄带
F_z = zoom(F_s, 2, order=3)  # 2x 超采样

tau_z = np.linspace(0.2, 12.0, N_TAU*2)
lg_z  = np.linspace(-1.0, 1.0, N_T1*2)

# ── 绘图 ──
from matplotlib.colors import LinearSegmentedColormap
fig, ax = plt.subplots(figsize=(8.5, 6))
levels = np.linspace(0, 1, 13)  # 少层级 = 更接近论文
paper_c = LinearSegmentedColormap.from_list('paper',
    ['#0d0887','#46039f','#7201a8','#9711a1','#c94d71','#d76e56',
     '#de8d3e','#e8ab31','#f0c92b','#fae724'], N=256)
cf = ax.contourf(*np.meshgrid(tau_z, lg_z), F_z, levels=levels,
                 cmap=paper_c, extend='both')
cs = ax.contour(*np.meshgrid(tau_z, lg_z), F_z,
                levels=np.arange(0.2, 1.0, 0.2),
                colors='white', linewidths=0.4, alpha=0.5)
ax.clabel(cs, inline=True, fontsize=7, fmt='%.1f')

ax.set_xlabel(r'$\tau$ (100/meV)', fontsize=13)
ax.set_ylabel(r'$\lg(t_1/E_1)$', fontsize=13)
ax.set_title(r'$|\langle\psi_1^-|U(6\tau)|\psi_1^+\rangle|^2$'
             + f'  ($E_1={E1_fixed}$ meV)', fontsize=12)
plt.colorbar(cf, ax=ax, label='Fidelity', ticks=np.linspace(0, 1, 6))
ax.text(0.98, 0.02,
        f'$E_1={E1_fixed}$ meV\n$t_c=E_0={tc}$ meV\ndouble swap',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=8, bbox=dict(boxstyle='round', facecolor='black', alpha=0.7),
        color='white')

plt.tight_layout()
plt.savefig('fig1d_E1_0_01_hires.png', dpi=200)
print(f"\n✓ Saved: fig1d_E1_0_01_hires.png")
