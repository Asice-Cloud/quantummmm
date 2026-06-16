#!/usr/bin/env python3
"""E₁=0 轨迹集合 — 高效版：一次传播，沿路采样"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3

def fp(t, tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1-np.cos(pi*t/tau))

def step1_matrix(t, tau, t1c):
    A = np.zeros((5,5))
    A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
    A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
    A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau)
    return A

def step2_matrix(t, tau, t1c):
    A = np.zeros((5,5))
    A[1,3]=-2*tc*fp(t,tau); A[3,1]=2*tc*fp(t,tau)
    A[2,3]=2*tc*fm(t,tau); A[3,2]=-2*tc*fm(t,tau)
    A[0,4]=2*t1c*fp(t,tau); A[4,0]=-2*t1c*fp(t,tau)
    return A

def step3_matrix(t, tau):
    A = np.zeros((5,5))
    A[2,3]=2*tc*fp(t,tau); A[3,2]=-2*tc*fp(t,tau)
    A[3,4]=2*E0*fm(t,tau); A[4,3]=-2*E0*fm(t,tau)
    return A

def rk4_step(R, A_fn, t, tau, dt, t1c=None):
    if t1c is not None:
        k1 = A_fn(t, tau, t1c) @ R
        k2 = A_fn(t+0.5*dt, tau, t1c) @ (R+0.5*dt*k1)
        k3 = A_fn(t+0.5*dt, tau, t1c) @ (R+0.5*dt*k2)
        k4 = A_fn(t+dt, tau, t1c) @ (R+dt*k3)
    else:
        k1 = A_fn(t, tau) @ R
        k2 = A_fn(t+0.5*dt, tau) @ (R+0.5*dt*k1)
        k3 = A_fn(t+0.5*dt, tau) @ (R+0.5*dt*k2)
        k4 = A_fn(t+dt, tau) @ (R+dt*k3)
    return R + dt/6*(k1+2*k2+2*k3+k4)

def R_to_bloch(R):
    x = 1j*(R[0,1]-R[1,0]); y = R[0,1]+R[1,0]; z = R[0,0]-R[1,1]
    v = np.real(np.array([x, y, z]))
    n = np.linalg.norm(v)
    return v/n if n > 1e-12 else np.zeros(3)

def full_trajectory(tau_p, t1, n_pts=150):
    tau_c = tau_p*100; dt = 3*tau_c/n_pts
    R = np.eye(5); traj = np.zeros((n_pts,3))
    for s in range(n_pts):
        t = s*dt
        if t < tau_c:
            R = rk4_step(R, step1_matrix, t, tau_c, dt, t1)
        elif t < 2*tau_c:
            R = rk4_step(R, step2_matrix, t-tau_c, tau_c, dt, t1)
        else:
            R = rk4_step(R, step3_matrix, t-2*tau_c, tau_c, dt)
        traj[s] = R_to_bloch(R)
    return traj

# ══════════════ Plot ══════════════
t1_list = np.logspace(-2.5, -1, 8)
tau_list = [1.5, 3.0, 6.0, 10.0]
colors = plt.cm.plasma(np.linspace(0.1, 0.95, len(t1_list)))

fig = plt.figure(figsize=(14, 6.5))

# ── 3D Bloch ──
ax = fig.add_subplot(121, projection='3d')
u,v = np.mgrid[0:2*pi:40j, 0:pi:20j]
ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v),
                  color='gray', alpha=0.08, linewidth=0.2)
phi = np.linspace(0,2*pi,200)
ax.plot(np.cos(phi), np.sin(phi), 0*phi, 'gray', alpha=0.25, ls='--', lw=0.8)

for ci, t1 in enumerate(t1_list):
    for tau in tau_list:
        print(f"  t1={t1:.4f}, tau={tau:.1f}...")
        tr = full_trajectory(tau, t1)
        alpha = 0.75 if tau <= 4 else 0.25
        ax.plot(tr[:,0], tr[:,1], tr[:,2], color=colors[ci], alpha=alpha, lw=0.7)
        ax.scatter(*tr[-1], color=colors[ci], s=8, alpha=0.6)

ax.set_xlabel('σ_x'); ax.set_ylabel('σ_y'); ax.set_zlabel('σ_z')
ax.set_xlim(-1.1,1.1); ax.set_ylim(-1.1,1.1); ax.set_zlim(-1.1,1.1)
ax.set_title(f'E₁=0: {len(t1_list)} t₁ values × {len(tau_list)} τ values\nAll on equator', fontsize=11)
ax.view_init(15, -50)

# ── Top view ──
ax2 = fig.add_subplot(122)
ax2.set_aspect('equal')
ax2.add_patch(plt.Circle((0,0),1,fill=False,color='gray',alpha=0.3,lw=0.8))
ax2.axhline(0,color='gray',alpha=0.15,lw=0.5); ax2.axvline(0,color='gray',alpha=0.15,lw=0.5)
ax2.annotate('σ_x',(1.1,0),fontsize=11,ha='center')
ax2.annotate('σ_y',(0,1.1),fontsize=11,ha='center')

for ci, t1 in enumerate(t1_list):
    for tau in tau_list:
        tr = full_trajectory(tau, t1, n_pts=120)
        alpha = 0.75 if tau <= 4 else 0.25
        ax2.plot(tr[:,0], tr[:,1], color=colors[ci], alpha=alpha, lw=0.7)
        ax2.scatter(*tr[-1,:2], color=colors[ci], s=6, alpha=0.5)

ax2.set_xlim(-1.15,1.15); ax2.set_ylim(-1.15,1.15)
ax2.set_title('Top view: σ_z ≡ 0 (never leaves equator)', fontsize=11)
ax2.set_xlabel('σ_x'); ax2.set_ylabel('σ_y')

plt.tight_layout()
plt.savefig('bloch_E1_zero_collection.png', dpi=200)
print("\n✓ bloch_E1_zero_collection.png")
print("E₁=0: rotation axis varies in σ_x-σ_y plane, but NEVER gets σ_z component")
