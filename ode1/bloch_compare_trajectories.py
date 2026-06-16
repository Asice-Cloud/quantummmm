#!/usr/bin/env python3
"""E₁=0 vs E₁≠0 轨迹对比图 — 统一画法"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
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

def R_to_bloch(R):
    x = 1j*(R[0,1]-R[1,0]); y = R[0,1]+R[1,0]; z = R[0,0]-R[1,1]
    v = np.real(np.array([x, y, z]))
    n = np.linalg.norm(v)
    return v/n if n > 1e-12 else np.zeros(3)

def trajectory(tau_p, E1, t1c, n_pts=200):
    tau_c = tau_p*100; dt_each = tau_c/n_pts
    R = np.eye(5); traj = np.zeros((3*n_pts, 3))
    for step in [1,2,3]:
        for s in range(n_pts):
            t = s*dt_each
            A = step_matrix(t, tau_c, E1, t1c, step)
            k1 = A @ R
            k2 = step_matrix(t+0.5*dt_each, tau_c, E1, t1c, step) @ (R+0.5*dt_each*k1)
            k3 = step_matrix(t+0.5*dt_each, tau_c, E1, t1c, step) @ (R+0.5*dt_each*k2)
            k4 = step_matrix(t+dt_each, tau_c, E1, t1c, step) @ (R+dt_each*k3)
            R = R + dt_each/6*(k1+2*k2+2*k3+k4)
            traj[(step-1)*n_pts+s] = R_to_bloch(R)
    return traj

# ═══════════════════════════════════════════════════════════
fig = plt.figure(figsize=(14, 6.5))

# 画球面的辅助函数
def draw_sphere(ax):
    u,v = np.mgrid[0:2*pi:40j, 0:pi:20j]
    ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v),
                      color='gray', alpha=0.06, linewidth=0.2)
    phi = np.linspace(0,2*pi,200)
    ax.plot(np.cos(phi), np.sin(phi), 0*phi, 'gray', alpha=0.25, ls='--', lw=0.6)
    ax.set_xlim(-1.1,1.1); ax.set_ylim(-1.1,1.1); ax.set_zlim(-1.1,1.1)
    ax.set_xlabel('σ_x'); ax.set_ylabel('σ_y'); ax.set_zlabel('σ_z')

# ── 左：E₁=0 ──
ax1 = fig.add_subplot(121, projection='3d')
draw_sphere(ax1)
ax1.set_title('E₁ = 0: Trajectories on equator', fontsize=12)

t1_list = np.logspace(-2.3, -1, 6)
tau_list = [2.0, 5.0, 8.0]
colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(t1_list)*len(tau_list)))

ci = 0
for t1 in t1_list:
    for tau in tau_list:
        print(f"  E₁=0: t1={t1:.4f}, tau={tau:.1f}...")
        tr = trajectory(tau, 0.0, t1)
        ax1.plot(tr[:,0], tr[:,1], tr[:,2], color=colors[ci], lw=0.9)
        ax1.scatter(*tr[0], color=colors[ci], s=20, marker='o', edgecolors='white', linewidths=0.3)
        ax1.scatter(*tr[-1], color=colors[ci], s=35, marker='*', edgecolors='white', linewidths=0.3)
        ci += 1
ax1.view_init(20, -50)

# ── 右：E₁≠0 ──
ax2 = fig.add_subplot(122, projection='3d')
draw_sphere(ax2)
ax2.set_title('E₁ ≠ 0: Trajectories leave equator', fontsize=12)

params_nonzero = [
    (0.01, 0.01, 5.0), (0.01, 0.03, 4.0), (0.005, 0.005, 8.0),
    (0.015, 0.005, 6.0), (0.02, 0.01, 3.0), (0.008, 0.015, 7.0),
]
colors2 = plt.cm.Oranges(np.linspace(0.3, 0.95, len(params_nonzero)))

for ci, (E1, t1, tau) in enumerate(params_nonzero):
    print(f"  E₁≠0: E1={E1:.4f}, t1={t1:.4f}, tau={tau:.1f}...")
    tr = trajectory(tau, E1, t1)
    ax2.plot(tr[:,0], tr[:,1], tr[:,2], color=colors2[ci], lw=0.9)
    ax2.scatter(*tr[0], color=colors2[ci], s=20, marker='o', edgecolors='white', linewidths=0.3)
    ax2.scatter(*tr[-1], color=colors2[ci], s=35, marker='*', edgecolors='white', linewidths=0.3)
ax2.view_init(20, -50)

plt.tight_layout()
plt.savefig('bloch_compare_trajectories.png', dpi=200)
print("\n✓ bloch_compare_trajectories.png")
print("  Left (blue):  E₁=0 — all trajectories on equator")
print("  Right (orange): E₁≠0 — trajectories cover full sphere")
print("  ○ = start, ★ = end")
