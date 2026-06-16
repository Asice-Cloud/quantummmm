#!/usr/bin/env python3
"""
快速线扫描：测试 Step2 t₁=0 对峰值位置的影响
=============================================
只扫 τ（1D），固定 t₁/E₁=1 (lg=0)，大幅加速
对比: 旧方案(step2 t₁≠0) vs 新方案(step2 t₁=0)
"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3

def fp(t, tau): return 0.5*(1.0 + np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1.0 - np.cos(pi*t/tau))

def b1_old(t, tau, e, t1c):  # Step 1 (same for both)
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    t2v = tc*fm(t,tau); A[1,3]=-2*t2v; A[3,1]=2*t2v
    t1v = t1c*fm(t,tau); A[0,4]=2*t1v; A[4,0]=-2*t1v
    Edv = E0*fp(t,tau); A[3,4]=2*Edv; A[4,3]=-2*Edv; return A

def b2_old(t, tau, e, t1c):  # OLD: t₁ present in step 2
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    t2v = tc*fp(t,tau); A[1,3]=-2*t2v; A[3,1]=2*t2v
    t3v = tc*fm(t,tau); A[2,3]=2*t3v; A[3,2]=-2*t3v
    t1v = t1c*fp(t,tau); A[0,4]=2*t1v; A[4,0]=-2*t1v; return A

def b2_new(t, tau, e, t1c):  # NEW: t₁=0 in step 2 (G1 on)
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    t2v = tc*fp(t,tau); A[1,3]=-2*t2v; A[3,1]=2*t2v
    t3v = tc*fm(t,tau); A[2,3]=2*t3v; A[3,2]=-2*t3v
    # t₁=0 — G1 on suppresses coupling
    return A

def b3(t, tau, e, t1c):  # Step 3 (t₁=0)
    A = np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    t3v = tc*fp(t,tau); A[2,3]=2*t3v; A[3,2]=-2*t3v
    Edv = E0*fm(t,tau); A[3,4]=2*Edv; A[4,3]=-2*Edv; return A

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

def scan_line(e1, tau_arr, t1_val, use_new_step2):
    """1D τ scan, returns fidelity array"""
    b2 = b2_new if use_new_step2 else b2_old
    F = np.zeros(len(tau_arr))
    for i, tp in enumerate(tau_arr):
        tc_phys = tp*100.0
        R1 = prop(b1_old, tc_phys, e1, t1_val)
        R2 = prop(b2, tc_phys, e1, t1_val)
        R3 = prop(b3, tc_phys, e1, t1_val)
        F[i] = fid((R3@R2@R1) @ (R3@R2@R1))
    return F

# ═══════════════════════════════════════════════════════════
N_TAU = 200  # 密集扫 τ
tau_p = np.linspace(0.2, 12.0, N_TAU)
E1_list = [0.004, 0.005, 0.006]
# t₁/E₁ = 1 → t₁ = E₁ (lg=0, paper's mid-range)
t1_by_E1 = 1.0

print("="*60)
print("快速线扫描: 测试 Step2 t₁=0 对峰值位置的影响")
print("固定 t₁/E₁ = 1.0, N_τ =", N_TAU)
print("="*60)

fig, axes = plt.subplots(len(E1_list), 1, figsize=(10, 3*len(E1_list)), sharex=True)
if len(E1_list) == 1: axes = [axes]

for idx, e1 in enumerate(E1_list):
    t1 = e1 * t1_by_E1
    print(f"\nE₁ = {e1:.4f} meV, t₁ = {t1:.4f} meV")
    
    F_old = scan_line(e1, tau_p, t1, use_new_step2=False)
    F_new = scan_line(e1, tau_p, t1, use_new_step2=True)
    
    # Find peak positions
    def find_peaks(F, tau_p):
        peaks = []
        for i in range(1, len(F)-1):
            if F[i] > F[i-1] and F[i] > F[i+1] and F[i] > 0.3:
                peaks.append((tau_p[i], F[i]))
        return peaks
    
    p_old = find_peaks(F_old, tau_p)
    p_new = find_peaks(F_new, tau_p)
    
    print(f"  OLD (step2 t₁≠0): peaks at τ = {[f'{p[0]:.2f}' for p in p_old[:6]]}")
    print(f"  NEW (step2 t₁=0): peaks at τ = {[f'{p[0]:.2f}' for p in p_new[:6]]}")
    if p_old and p_new:
        shift = p_new[0][0] - p_old[0][0]
        print(f"  → Main peak shift: {shift:+.2f}τ")
    
    ax = axes[idx]
    ax.plot(tau_p, F_old, 'b-', label='OLD: step2 t₁≠0', alpha=0.7, lw=1.5)
    ax.plot(tau_p, F_new, 'r--', label='NEW: step2 t₁=0', alpha=0.8, lw=1.5)
    ax.set_ylabel('Fidelity', fontsize=11)
    ax.set_title(f'E₁ = {e1} meV, t₁/E₁ = {t1_by_E1}', fontsize=12)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)
    
    # Mark peaks
    for pi_tau, pi_fid in p_old[:3]:
        ax.axvline(pi_tau, color='blue', alpha=0.2, ls=':')
    for pi_tau, pi_fid in p_new[:3]:
        ax.axvline(pi_tau, color='red', alpha=0.2, ls=':')

axes[-1].set_xlabel(r'$\tau$ (100/meV)', fontsize=12)
plt.tight_layout()
plt.savefig('test_step2_t1_zero_linescan.png', dpi=150)
print(f"\n✓ Saved: test_step2_t1_zero_linescan.png")
print("\nIf NEW peak at τ≈5.8 → step2 t₁=0 was the main fix.")
print("If peak still at τ≈7 → need to explore tc/E0 ratio next.")
