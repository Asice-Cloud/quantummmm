#!/usr/bin/env python3
"""E₁ 扫描 —— 轻量版，只算终点"""
import numpy as np, matplotlib; matplotlib.use('Agg')
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
    omega = max(2*abs(E1), 0.6)  # dominant frequency
    n = max(200, int(tau_c * omega / pi * 10))  # ~10 steps per half-period
    dt = tau_c/n; R = np.eye(5)
    for s in range(n):
        t = s*dt; A = step_matrix(t, tau_c, E1, t1c, step)
        k1=dt*A@R; k2=dt*step_matrix(t+0.5*dt,tau_c,E1,t1c,step)@(R+0.5*k1)
        k3=dt*step_matrix(t+0.5*dt,tau_c,E1,t1c,step)@(R+0.5*k2)
        k4=dt*step_matrix(t+dt,tau_c,E1,t1c,step)@(R+k3)
        R += (k1+2*k2+2*k3+k4)/6
        # re-orthogonalize occasionally
        if s % 50 == 0:
            U,_,Vt = np.linalg.svd(R); R = U@Vt
    return R

def full_braid(tau_p, E1, t1c):
    tau_c = tau_p*100; R = np.eye(5)
    for step in [1,2,3]:
        R = propagate_step(tau_c, E1, t1c, step) @ R
    return R

def get_bloch(R):
    x=1j*(R[0,1]-R[1,0]); y=R[0,1]+R[1,0]; z=R[0,0]-R[1,1]
    v=np.real([x,y,z]); n=np.linalg.norm(v)
    return v/n if n>1e-12 else np.zeros(3)

# ═══════════════════════════════════════════════════════════
E1s = np.concatenate([[0.0], np.logspace(-3, np.log10(0.3), 20), np.linspace(0.35,0.5,5)])
t1s = [0.0, 0.01, 0.05]
taus = [5, 15]

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
for ri, tau in enumerate(taus):
    for ci, t1 in enumerate(t1s):
        ax = axes[ri, ci]; vv = np.zeros((len(E1s), 3))
        for i, E1 in enumerate(E1s):
            R = full_braid(tau, E1, t1)
            vv[i] = get_bloch(R)
            if i % 5 == 0: print(f"  τ={tau} t₁={t1} E₁={E1:.4f} done")
        ax.plot(E1s, vv[:,0], 'o-', ms=3, label='$v_x$')
        ax.plot(E1s, vv[:,1], 's-', ms=3, label='$v_y$ (σz flux)')
        ax.plot(E1s, vv[:,2], '^-', ms=3, label='$v_z$')
        ax.axvline(0.3, color='gray', ls='--', lw=0.8)
        ax.set_xscale('symlog', linthresh=0.005)
        ax.set_xlabel('E₁ (meV)'); ax.set_ylabel('Bloch')
        ax.set_title(f'τ={tau}, t₁={t1}'); ax.legend(fontsize=7)
        ax.set_ylim(-1.1, 1.1); ax.grid(True, alpha=0.3)

plt.tight_layout(); plt.savefig('na_f/fig_e1_scan_v2.png', dpi=200)
print("\n✓ na_f/fig_e1_scan_v2.png")
