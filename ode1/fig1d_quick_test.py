#!/usr/bin/env python3
"""
Fig 1(d) 快速等高线 — E₁=0.004, OLD step2 (t₁保持)
=====================================================
低分辨率网格 50×60，快速验证全局匹配
"""
import numpy as np
from scipy.ndimage import gaussian_filter, zoom
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3; E1_fixed = 0.004

def fp(t, tau): return 0.5*(1.0 + np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1.0 - np.cos(pi*t/tau))

def b1(t, tau, e, t1c):
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
    A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
    A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau)
    return A

def b2(t, tau, e, t1c):  # OLD: t₁ present
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3]=-2*tc*fp(t,tau); A[3,1]=2*tc*fp(t,tau)
    A[2,3]=2*tc*fm(t,tau); A[3,2]=-2*tc*fm(t,tau)
    A[0,4]=2*t1c*fp(t,tau); A[4,0]=-2*t1c*fp(t,tau)
    return A

def b3(t, tau, e, t1c):
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[2,3]=2*tc*fp(t,tau); A[3,2]=-2*tc*fp(t,tau)
    A[3,4]=2*E0*fm(t,tau); A[4,3]=-2*E0*fm(t,tau)
    return A

def prop(bld, tau, e, t1c):
    n = max(500, int(2*tau))
    dt = tau/n; R = np.eye(5)
    for s in range(n):
        t = s*dt
        k1 = bld(t,tau,e,t1c) @ R
        k2 = bld(t+0.5*dt,tau,e,t1c) @ (R+0.5*dt*k1)
        k3 = bld(t+0.5*dt,tau,e,t1c) @ (R+0.5*dt*k2)
        k4 = bld(t+dt,tau,e,t1c) @ (R+dt*k3)
        R += dt/6.0*(k1+2*k2+2*k3+k4)
    return R

def fid(R):
    ov = 0.5*(R[0,0]+1j*R[1,0]+1j*R[0,1]-R[1,1])
    return np.abs(ov)**2

# ═══════════════════════════════════════════════════════════
N_TAU, N_T1 = 50, 60  # 低分辨率快速验证
tau_p = np.linspace(0.2, 12.0, N_TAU)
tau_c = tau_p * 100.0
t1_v  = E1_fixed * 10**np.linspace(-1.0, 1.0, N_T1)

print(f"Fig 1(d) 快速等高线 | E₁={E1_fixed} | 网格={N_TAU}×{N_T1}")
print(f"  τ ∈ [0.2, 12], t₁/E₁ ∈ [0.1, 10]")
print(f"  Step2: t₁ present (OLD), Step3: t₁=0")

F = np.zeros((N_T1, N_TAU))
for i in range(N_TAU):
    for j in range(N_T1):
        R1=prop(b1,tau_c[i],E1_fixed,t1_v[j])
        R2=prop(b2,tau_c[i],E1_fixed,t1_v[j])
        R3=prop(b3,tau_c[i],E1_fixed,t1_v[j])
        F[j,i] = fid((R3@R2@R1)@(R3@R2@R1))
    if (i+1)%10==0: print(f"  τ {i+1}/{N_TAU} ({tau_p[i]:.1f}): fid range [{F[:,i].min():.4f}, {F[:,i].max():.4f}]")

print(f"\nRaw range: [{F.min():.4f}, {F.max():.4f}]")

# Smoothing
F_s = gaussian_filter(F, sigma=0.8, mode='nearest')
F_z = zoom(F_s, 3, order=3)
tau_z = np.linspace(0.2, 12.0, N_TAU*3)
lg_z  = np.linspace(-1.0, 1.0, N_T1*3)

# Plot
from matplotlib.colors import LinearSegmentedColormap
fig, ax = plt.subplots(figsize=(8.5, 6))
levels = np.linspace(0, 1, 13)
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
             + f'  ($E_1={E1_fixed}$ meV, step2 t₁≠0)', fontsize=12)
plt.colorbar(cf, ax=ax, label='Fidelity', ticks=np.linspace(0, 1, 6))
ax.text(0.98, 0.02,
        f'$E_1={E1_fixed}$ meV\n$t_c=E_0={tc}$ meV\nstep2 t₁≠0, step3 t₁=0',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=8, bbox=dict(boxstyle='round', facecolor='black', alpha=0.7),
        color='white')
plt.tight_layout()
plt.savefig('fig1d_E1_0_004_oldstep2.png', dpi=200)
print(f"\n✓ Saved: fig1d_E1_0_004_oldstep2.png")

# Also find the main peak
idx = np.unravel_index(np.argmax(F_z), F_z.shape)
print(f"Main peak: τ={tau_z[idx[1]]:.2f}, lg(t₁/E₁)={lg_z[idx[0]]:.2f}, fid={F_z[idx]:.4f}")
