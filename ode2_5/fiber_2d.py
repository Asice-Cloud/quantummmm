#!/usr/bin/env python3
"""Sp(2) 纤维分析 2D：扫描 E₁ × t₁"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3

def fp(t, tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1-np.cos(pi*t/tau))

def quat_mul(p, q):
    a,b,c,d = p; e,f,g,h = q
    return np.array([a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g, a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e])

def quat_norm(q): return np.sqrt(q[0]**2+q[1]**2+q[2]**2+q[3]**2)
def quat_conj(q): return np.array([q[0],-q[1],-q[2],-q[3]])
def quat_inv(q): return quat_conj(q)/quat_norm(q)**2

def get_ABCD(t, tau_c, E1, t1):
    t2 = tc*(fm(t,tau_c) if t<tau_c else (fp(t,tau_c) if t<2*tau_c else 0))
    t3 = tc*(0 if t<tau_c else (fm(t,tau_c) if t<2*tau_c else fp(t,tau_c)))
    Ed = E0*(fp(t,tau_c) if t<tau_c else (0 if t<2*tau_c else fm(t,tau_c)))
    t1v= t1*(fm(t,tau_c) if t<tau_c else (fp(t,tau_c) if t<2*tau_c else 0))
    return (np.array([0,(E1+t3)/2,t2/2,0]),
            np.array([0,(-E1+t3)/2,t2/2,0]),
            np.array([abs(t1v)/2,0,0,Ed/2]),
            np.array([-abs(t1v)/2,0,0,Ed/2]))

def propagate_sp2(tau_p, E1, t1, nps=300):
    tau_c = tau_p*100
    X=np.array([1.,0,0,0]); Y=np.array([0.,0,0,0])
    Z=np.array([0.,0,0,0]); W=np.array([1.,0,0,0])
    for step in range(3):
        t0=step*tau_c; dt=tau_c/nps
        for s in range(nps):
            t=t0+s*dt; A,B,C,D=get_ABCD(t,tau_c,E1,t1)
            def f(X_,Z_,Y_,W_, A_,B_,C_,D_):
                return (quat_mul(A_,X_)+quat_mul(B_,Z_),
                        quat_mul(C_,X_)+quat_mul(D_,Z_),
                        quat_mul(A_,Y_)+quat_mul(B_,W_),
                        quat_mul(C_,Y_)+quat_mul(D_,W_))
            k1X,k1Z,k1Y,k1W=f(X,Z,Y,W, A,B,C,D)
            k1X*=dt;k1Z*=dt;k1Y*=dt;k1W*=dt
            A2,B2,C2,D2=get_ABCD(t+0.5*dt,tau_c,E1,t1)
            k2X,k2Z,k2Y,k2W=f(X+0.5*k1X,Z+0.5*k1Z,Y+0.5*k1Y,W+0.5*k1W,A2,B2,C2,D2)
            k2X*=dt;k2Z*=dt;k2Y*=dt;k2W*=dt
            k3X,k3Z,k3Y,k3W=f(X+0.5*k2X,Z+0.5*k2Z,Y+0.5*k2Y,W+0.5*k2W,A2,B2,C2,D2)
            k3X*=dt;k3Z*=dt;k3Y*=dt;k3W*=dt
            A4,B4,C4,D4=get_ABCD(t+dt,tau_c,E1,t1)
            k4X,k4Z,k4Y,k4W=f(X+k3X,Z+k3Z,Y+k3Y,W+k3W,A4,B4,C4,D4)
            k4X*=dt;k4Z*=dt;k4Y*=dt;k4W*=dt
            X+=(k1X+2*k2X+2*k3X+k4X)/6; Z+=(k1Z+2*k2Z+2*k3Z+k4Z)/6
            Y+=(k1Y+2*k2Y+2*k3Y+k4Y)/6; W+=(k1W+2*k2W+2*k3W+k4W)/6
        # re-orthogonalize
        c0=np.hstack([X,Z]); c1=np.hstack([Y,W])
        c0/=np.sqrt(np.sum(c0**2)); c1-=np.sum(c0*c1)*c0; c1/=np.sqrt(np.sum(c1**2))
        X=c0[:4];Z=c0[4:];Y=c1[:4];W=c1[4:]
    return X,Y,Z,W

def decompose(X,Y,Z,W):
    q=quat_mul(Z,quat_inv(X))
    nX=quat_norm(X); nW=quat_norm(W)
    u=X/nX if nX>1e-12 else np.array([1.,0,0,0])
    v=W/nW if nW>1e-12 else np.array([1.,0,0,0])
    th_u=2*np.arccos(np.clip(u[0],-1,1)); th_v=2*np.arccos(np.clip(v[0],-1,1))
    return q,th_u,th_v,quat_norm(q)

# ═══════════════════════════════════════════════════════════
E1s = np.concatenate([[0.0], np.logspace(-3, np.log10(0.5), 20)])
t1s = np.concatenate([[0.0], np.logspace(-3, np.log10(0.1), 15)])
tau = 15

grid_u = np.zeros((len(E1s), len(t1s)))
grid_v = np.zeros((len(E1s), len(t1s)))
grid_q = np.zeros((len(E1s), len(t1s)))

for i, E1 in enumerate(E1s):
    for j, t1 in enumerate(t1s):
        X,Y,Z,W = propagate_sp2(tau, E1, t1, nps=200)
        q,th_u,th_v,qn = decompose(X,Y,Z,W)
        grid_u[i,j] = th_u/np.pi
        grid_v[i,j] = th_v/np.pi
        grid_q[i,j] = qn
    print(f"E₁={E1:.4f} done ({i+1}/{len(E1s)})")

# ═══════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
TT, EE = np.meshgrid(t1s, E1s)

for ax, data, title in zip(axes,
    [grid_u, grid_v, grid_q],
    [r'$\theta_u / \pi$ (MZM fiber)',
     r'$\theta_v / \pi$ (ancilla fiber)',
     r'$|q|$ (base space)']):
    im = ax.pcolormesh(TT, EE, data, shading='auto', cmap='RdYlBu_r')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('t₁ (meV)'); ax.set_ylabel('E₁ (meV)')
    ax.set_title(title)
    plt.colorbar(im, ax=ax)

plt.suptitle(f'Sp(2) fiber-base decomposition, τ={tau}', fontsize=13)
plt.tight_layout()
plt.savefig('ode2_5/fig_fiber_2d.png', dpi=200)
print("\n✓ ode2_5/fig_fiber_2d.png")
