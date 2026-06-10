#!/usr/bin/env python3
"""
Fig 1(d) 最终版 —— 高分辨率平滑 contour 图
============================================
使用 RK4-fixed-step so(5) 传播 + 高斯平滑 + contourf.

论文 Fig 1(d): 双次 swap braiding 后 fidelity = |⟨ψ₁⁻|U(6τ)|ψ₁⁺⟩|²
参数: E₁=0.01 meV, t_c=E₀=0.3 meV
"""
import numpy as np
from scipy.ndimage import gaussian_filter, zoom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3; E1_fixed = 0.01
def fp(t, tau):  return 0.5 * (1.0 + np.cos(pi * t / tau))
def fm(t, tau):  return 0.5 * (1.0 - np.cos(pi * t / tau))

def build_A_step1(t, tau, e, t1):
    A = np.zeros((5,5))
    A[0,1]=2*e; A[1,0]=-2*e
    t2=tc*fm(t,tau); A[1,3]=-2*t2; A[3,1]=2*t2
    t=t1*fm(t,tau); A[0,4]=2*t; A[4,0]=-2*t
    Ed=E0*fp(t,tau); A[3,4]=2*Ed; A[4,3]=-2*Ed; return A

def build_A_step2(t, tau, e, t1):
    A = np.zeros((5,5))
    A[0,1]=2*e; A[1,0]=-2*e
    t2=tc*fp(t,tau); A[1,3]=-2*t2; A[3,1]=2*t2
    t3=tc*fm(t,tau); A[2,3]=2*t3; A[3,2]=-2*t3
    t=t1*fp(t,tau); A[0,4]=2*t; A[4,0]=-2*t; return A

def build_A_step3(t, tau, e, t1):
    A = np.zeros((5,5))
    A[0,1]=2*e; A[1,0]=-2*e
    t3=tc*fp(t,tau); A[2,3]=2*t3; A[3,2]=-2*t3
    Ed=E0*fm(t,tau); A[3,4]=2*Ed; A[4,3]=-2*Ed
    t=t1; A[0,4]=2*t; A[4,0]=-2*t; return A

def propagate_rk4(bld, tau, e, t1):
    n = max(150, int(200 * tau / 12.0))
    dt = tau / n; R = np.eye(5)
    for s in range(n):
        t = s * dt
        k1 = bld(t, tau, e, t1) @ R
        k2 = bld(t+0.5*dt, tau, e, t1) @ (R + 0.5*dt*k1)
        k3 = bld(t+0.5*dt, tau, e, t1) @ (R + 0.5*dt*k2)
        k4 = bld(t+dt, tau, e, t1) @ (R + dt*k3)
        R += dt/6.0 * (k1 + 2*k2 + 2*k3 + k4)
    return R

def fidelity(R_total):
    """|⟨ψ⁻|U|ψ⁺⟩|² 其中 |ψ⁺⟩=(γ₁+iγ₂)/√2, |ψ⁻⟩=(γ₁-iγ₂)/√2"""
    ov = 0.5 * (R_total[0,0] + 1j*R_total[1,0] + 1j*R_total[0,1] - R_total[1,1])
    return np.abs(ov)**2

print("=" * 55)
print("Fig 1(d) 最终版 — 扫描")
print("=" * 55)
print(f"E₁={E1_fixed} meV, t_c=E₀={tc} meV")

n_tau, n_t1 = 200, 250
tau_vals = np.linspace(0.2, 12.0, n_tau)
t1_vals  = E1_fixed * 10**np.linspace(-1.0, 1.0, n_t1)

fid = np.zeros((n_t1, n_tau))
step_rpt = max(1, n_tau // 10)

for i in range(n_tau):
    tau = tau_vals[i]
    for j in range(n_t1):
        t1 = t1_vals[j]
        R1 = propagate_rk4(build_A_step1, tau, E1_fixed, t1)
        R2 = propagate_rk4(build_A_step2, tau, E1_fixed, t1)
        R3 = propagate_rk4(build_A_step3, tau, E1_fixed, t1)
        Rs = R3 @ R2 @ R1
        Rd = Rs @ Rs
        fid[j, i] = fidelity(Rd)
    if (i+1) % step_rpt == 0:
        print(f"  {i+1}/{n_tau} ({100*(i+1)//n_tau}%)")

print(f"Raw range: [{fid.min():.4f}, {fid.max():.4f}]")

# ── 平滑 + 超采样 ──
fid_s = gaussian_filter(fid, sigma=1.0, mode='nearest')
fid_z = zoom(fid_s, 3, order=3)

tau_z = np.linspace(0.2, 12.0, n_tau*3)
lg_z  = np.linspace(-1.0, 1.0, n_t1*3)
Tz, Lz = np.meshgrid(tau_z, lg_z)

print(f"Upsampled: {fid_z.shape}, range: [{fid_z.min():.4f}, {fid_z.max():.4f}]")

# ── 绘图 ──
fig, ax = plt.subplots(1, 1, figsize=(9, 6.5))

levels = np.linspace(0, 1, 21)
cf = ax.contourf(Tz, Lz, fid_z, levels=levels, cmap='inferno', extend='both')

# 细等高线
cs = ax.contour(Tz, Lz, fid_z, levels=np.linspace(0, 1, 11),
                colors='w', linewidths=0.35, alpha=0.3)
# 主要等高线标注
cs2 = ax.contour(Tz, Lz, fid_z, levels=[0.2, 0.4, 0.6],
                 colors='w', linewidths=0.7, alpha=0.5)
ax.clabel(cs2, inline=True, fontsize=7, fmt='%.1f')

ax.set_xlabel(r'$\tau$ (100/meV)', fontsize=14)
ax.set_ylabel(r'$\lg(t_1/E_1)$', fontsize=14)
ax.set_title(r'$|\langle\psi_1^-|U(6\tau)|\psi_1^+\rangle|^2$'
             + f'  ($E_1={E1_fixed}$ meV, $t_c=E_0={tc}$ meV)',
             fontsize=12)

cbar = plt.colorbar(cf, ax=ax, label='Fidelity', ticks=np.linspace(0, 1, 6))

ax.text(0.98, 0.02,
        f'$E_1={E1_fixed}$ meV\n'
        f'$t_c=E_0={tc}$ meV\n'
        f'double swap\n'
        f'grid: {n_tau}$\\times${n_t1}',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=8, bbox=dict(boxstyle='round', facecolor='black', alpha=0.7),
        color='white')

plt.tight_layout()
plt.savefig('fig1d_reproduction_v4.png', dpi=250)
print(f"\nSaved: fig1d_reproduction_v4.png")
