#!/usr/bin/env python3
"""
诊断测试: Step 2 中 t₁=0 (G1 开), 线扫描 τ 看峰值位置
"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3

def fp(t, tau): return 0.5 * (1.0 + np.cos(pi * t / tau))
def fm(t, tau): return 0.5 * (1.0 - np.cos(pi * t / tau))

def b1(t, tau, e, t1c):
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3] = -2*tc*fm(t,tau); A[3,1] = 2*tc*fm(t,tau)
    A[0,4] = 2*t1c*fm(t,tau); A[4,0] = -2*t1c*fm(t,tau)
    A[3,4] = 2*E0*fp(t,tau); A[4,3] = -2*E0*fp(t,tau)
    return A

def b2_no_t1(t, tau, e, t1c):  # ★ t₁=0 全程 (G1 开)
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3] = -2*tc*fp(t,tau); A[3,1] = 2*tc*fp(t,tau)
    A[2,3] = 2*tc*fm(t,tau); A[3,2] = -2*tc*fm(t,tau)
    # t1 = 0 — G1 is ON
    return A

def b2_old(t, tau, e, t1c):  # 旧版: t₁ 在 step2 也渐变
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3] = -2*tc*fp(t,tau); A[3,1] = 2*tc*fp(t,tau)
    A[2,3] = 2*tc*fm(t,tau); A[3,2] = -2*tc*fm(t,tau)
    A[0,4] = 2*t1c*fp(t,tau); A[4,0] = -2*t1c*fp(t,tau)
    return A

def b3(t, tau, e, t1c):
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[2,3] = 2*tc*fp(t,tau); A[3,2] = -2*tc*fp(t,tau)
    A[3,4] = 2*E0*fm(t,tau); A[4,3] = -2*E0*fm(t,tau)
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
        R += dt/6.*(k1+2*k2+2*k3+k4)
    return R

def fid(R):
    return np.abs(0.5*(R[0,0]+1j*R[1,0]+1j*R[0,1]-R[1,1]))**2

# ═══ 线扫描: 固定 lg(t₁/E₁)=0, 扫描 τ ═══
tau_p = np.linspace(0.2, 12.0, 200)
tau_c = tau_p * 100.0

E1_vals = [0.001, 0.004, 0.005, 0.01]
versions = [('t₁=0 in Step2', b2_no_t1), ('旧版 t₁≠0', b2_old)]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for idx, E1 in enumerate(E1_vals):
    ax = axes[idx]
    for label, b2_func in versions:
        F = np.zeros(len(tau_p))
        for i in range(len(tau_p)):
            t1_val = E1  # lg=0 → t₁=E₁
            R1 = prop(b1, tau_c[i], E1, t1_val)
            R2 = prop(b2_func, tau_c[i], E1, t1_val)
            R3 = prop(b3, tau_c[i], E1, t1_val)
            F[i] = fid((R3@R2@R1)@(R3@R2@R1))
        ax.plot(tau_p, F, lw=1.5, label=label, alpha=0.8)
    
    ax.set_title(f'E₁={E1} meV,  lg(t₁/E₁)=0')
    ax.set_xlabel('τ (100/meV)'); ax.set_ylabel('Fidelity')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)

plt.suptitle('Step 2 t₁ 诊断: 旧版 vs t₁=0', fontsize=14)
plt.tight_layout()
plt.savefig('debug_step2_t1_line.png', dpi=150)
print("✓ Saved: debug_step2_t1_line.png")

# ═══ 也打印峰值位置 ═══
print("\n峰值位置 (τ):")
print(f"{'E₁':>8} {'旧版 τ_peak':>12} {'t₁=0 τ_peak':>12} {'差':>8}")
for E1 in E1_vals:
    peaks = {}
    for label, b2_func in versions:
        F = np.zeros(len(tau_p))
        for i in range(len(tau_p)):
            R1 = prop(b1, tau_c[i], E1, E1)
            R2 = prop(b2_func, tau_c[i], E1, E1)
            R3 = prop(b3, tau_c[i], E1, E1)
            F[i] = fid((R3@R2@R1)@(R3@R2@R1))
        # 找峰值
        peak_idx = np.argmax(F[10:]) + 10  # 跳过 τ < 1
        peaks[label] = tau_p[peak_idx]
    diff = peaks['t₁=0 in Step2'] - peaks['旧版 t₁≠0']
    print(f"{E1:8.3f} {peaks['旧版 t₁≠0']:12.2f} {peaks['t₁=0 in Step2']:12.2f} {diff:8.2f}")
