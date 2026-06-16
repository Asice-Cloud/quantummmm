#!/usr/bin/env python3
"""
对比测试：单编织 vs 双编织 @ E₁=0.01
=====================================
目标：判断论文 Fig 1(d) 到底是单编织还是双编织
"""
import numpy as np
from scipy.ndimage import gaussian_filter, zoom
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3; E1_fixed = 0.01

def fp(t, tau): return 0.5*(1.0 + np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1.0 - np.cos(pi*t/tau))

def b1(t, tau, e, t1c):
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
    A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
    A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau)
    return A

def b2(t, tau, e, t1c):
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
N_TAU = 200
tau_p = np.linspace(0.2, 12.0, N_TAU)
tau_c = tau_p * 100.0
t1_list = [E1_fixed * 0.1, E1_fixed * 1.0, E1_fixed * 10.0]  # lg=-1, 0, 1

fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

for idx, t1_val in enumerate(t1_list):
    F_single = np.zeros(N_TAU)
    F_double = np.zeros(N_TAU)
    
    for i in range(N_TAU):
        R1=prop(b1,tau_c[i],E1_fixed,t1_val)
        R2=prop(b2,tau_c[i],E1_fixed,t1_val)
        R3=prop(b3,tau_c[i],E1_fixed,t1_val)
        U3 = R3@R2@R1
        F_single[i] = fid(U3)          # 单编织 U(3τ)
        F_double[i] = fid(U3 @ U3)     # 双编织 U(6τ)
    
    # Count oscillation peaks (>0.3 threshold)
    def count_peaks(F):
        peaks = 0
        for i in range(1, len(F)-1):
            if F[i] > F[i-1] and F[i] > F[i+1] and F[i] > 0.3:
                peaks += 1
        return peaks
    
    n_s = count_peaks(F_single)
    n_d = count_peaks(F_double)
    
    ax = axes[idx]
    ax.plot(tau_p, F_single, 'b-', label=f'SINGLE 3τ ({n_s} peaks)', alpha=0.7, lw=1.5)
    ax.plot(tau_p, F_double, 'r--', label=f'DOUBLE 6τ ({n_d} peaks)', alpha=0.8, lw=1.5)
    ax.set_ylabel('Fidelity', fontsize=11)
    ax.set_title(f'E₁={E1_fixed}, t₁/E₁={t1_val/E1_fixed:.0f} (lg={np.log10(t1_val/E1_fixed):.0f})', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)

axes[-1].set_xlabel(r'$\tau$ (100/meV)', fontsize=12)
plt.tight_layout()
plt.savefig('compare_single_vs_double.png', dpi=150)

# ═══════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════
print(f"E₁ = {E1_fixed} meV")
print(f"Expected oscillations (from Φ=6E₁·1200 for double):")
print(f"  Φ_total = {6*E1_fixed*1200:.1f} rad")
print(f"  Double braid cycles = {6*E1_fixed*1200/(2*pi):.1f}")
print(f"  Single braid cycles = {3*E1_fixed*1200/(2*pi):.1f}")
print(f"  Paper shows ~4-5 cycles")
print(f"\n✓ Saved: compare_single_vs_double.png")
