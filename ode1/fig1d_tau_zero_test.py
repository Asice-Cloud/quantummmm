#!/usr/bin/env python3
"""τ 从 0 开始的高密度线扫：检查 τ<2 区域的窄带"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3; E1 = 0.01

def fp(t, tau): return 0.5*(1+np.cos(pi*t/tau)) if tau>0 else 0.0
def fm(t, tau): return 0.5*(1-np.cos(pi*t/tau)) if tau>0 else 0.0

def b1(t, tau, t1c):
    A = np.zeros((5,5)); A[0,1]=2*E1; A[1,0]=-2*E1
    if tau > 1e-10:
        A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
        A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
        A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau)
    return A

def b2(t, tau, t1c):
    A = np.zeros((5,5)); A[0,1]=2*E1; A[1,0]=-2*E1
    if tau > 1e-10:
        A[1,3]=-2*tc*fp(t,tau); A[3,1]=2*tc*fp(t,tau)
        A[2,3]=2*tc*fm(t,tau); A[3,2]=-2*tc*fm(t,tau)
        A[0,4]=2*t1c*fp(t,tau); A[4,0]=-2*t1c*fp(t,tau)
    return A

def b3(t, tau):
    A = np.zeros((5,5)); A[0,1]=2*E1; A[1,0]=-2*E1
    if tau > 1e-10:
        A[2,3]=2*tc*fp(t,tau); A[3,2]=-2*tc*fp(t,tau)
        A[3,4]=2*E0*fm(t,tau); A[4,3]=-2*E0*fm(t,tau)
    return A

def prop(bld, tau, t1c=None):
    """propagate: tau=0 → R=I (no evolution)"""
    if tau < 1e-10:
        return np.eye(5)
    n = max(500, int(2*tau))
    dt = tau/n; R = np.eye(5)
    for s in range(n):
        t = s*dt
        if t1c is not None:
            A = bld(t,tau,t1c)
        else:
            A = bld(t,tau)
        k1 = A @ R
        if t1c is not None:
            k2 = bld(t+0.5*dt,tau,t1c) @ (R+0.5*dt*k1)
            k3 = bld(t+0.5*dt,tau,t1c) @ (R+0.5*dt*k2)
            k4 = bld(t+dt,tau,t1c) @ (R+dt*k3)
        else:
            k2 = bld(t+0.5*dt,tau) @ (R+0.5*dt*k1)
            k3 = bld(t+0.5*dt,tau) @ (R+0.5*dt*k2)
            k4 = bld(t+dt,tau) @ (R+dt*k3)
        R += dt/6*(k1+2*k2+2*k3+k4)
    return R

def fid(R):
    ov = 0.5*(R[0,0]+1j*R[1,0]+1j*R[0,1]-R[1,1])
    return np.abs(ov)**2

# ═══════════════════════════════════════════════════════════
# τ from 0 to 3, dense scan
N_TAU = 400
tau_p = np.linspace(0.0, 3.0, N_TAU)  # start from exactly 0
tau_c = tau_p * 100.0

t1_vals = [E1*0.1, E1*1.0, E1*10.0]
labels = ['lg=-1', 'lg=0', 'lg=1']

fig, ax = plt.subplots(figsize=(10, 5))
colors = ['#2196F3', '#FF5722', '#4CAF50']

for t1, label, color in zip(t1_vals, labels, colors):
    F = np.zeros(N_TAU)
    for i in range(N_TAU):
        if tau_c[i] < 1e-10:
            F[i] = 0.0  # τ=0: no evolution, |ψ₁⁺⟩⊥|ψ₁⁻⟩
        else:
            R1 = prop(b1, tau_c[i], t1)
            R2 = prop(b2, tau_c[i], t1)
            R3 = prop(b3, tau_c[i])
            U3 = R3 @ R2 @ R1
            F[i] = fid(U3 @ U3)
    
    ax.plot(tau_p, F, color=color, lw=1.2, label=f't₁/E₁={t1/E1:.0f} ({label})')

ax.set_xlabel(r'$\tau$ (100/meV)', fontsize=13)
ax.set_ylabel('Fidelity', fontsize=13)
ax.set_title(f'E₁={E1} meV, τ∈[0,3], {N_TAU} points', fontsize=12)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 1.05)

# Zoom inset: τ∈[0, 0.5]
axins = ax.inset_axes([0.55, 0.5, 0.4, 0.4])
for t1, label, color in zip(t1_vals, labels, colors):
    F = np.zeros(N_TAU)
    for i in range(N_TAU):
        if tau_c[i] < 1e-10:
            F[i] = 0.0
        else:
            R1 = prop(b1, tau_c[i], t1)
            R2 = prop(b2, tau_c[i], t1)
            R3 = prop(b3, tau_c[i])
            U3 = R3 @ R2 @ R1
            F[i] = fid(U3 @ U3)
    axins.plot(tau_p, F, color=color, lw=1.0)
axins.set_xlim(0, 0.5)
axins.set_ylim(0, 1.05)
axins.grid(True, alpha=0.3)
axins.set_title('τ∈[0, 0.5] zoom', fontsize=9)
ax.indicate_inset_zoom(axins, edgecolor='gray')

plt.tight_layout()
plt.savefig('fig1d_tau_from_zero.png', dpi=200)
print("✓ fig1d_tau_from_zero.png")
print(f"  E₁={E1}, τ∈[0,3], {N_TAU} points")
print("  τ=0: fid=0 (no evolution, |ψ₁⁺⟩⊥|ψ₁⁻⟩)")
