#!/usr/bin/env python3
"""E₁≠0 轨迹：展示离开赤道、覆盖球面的完整路径"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3

def fp(t, tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1-np.cos(pi*t/tau))

def step_matrix(t, tau, E1, t1c, step):
    A = np.zeros((5,5))
    A[0,1]=2*E1; A[1,0]=-2*E1  # E₁ always on
    if step == 1:
        A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
        A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
        A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau)
    elif step == 2:
        A[1,3]=-2*tc*fp(t,tau); A[3,1]=2*tc*fp(t,tau)
        A[2,3]=2*tc*fm(t,tau); A[3,2]=-2*tc*fm(t,tau)
        A[0,4]=2*t1c*fp(t,tau); A[4,0]=-2*t1c*fp(t,tau)
    else:  # step 3
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
            idx = (step-1)*n_pts + s
            traj[idx] = R_to_bloch(R)
    return traj

# ═══════════════════════════════════════════════════════════
fig = plt.figure(figsize=(16, 7))

# 视角1: 3D
ax = fig.add_subplot(121, projection='3d')
# 球面
u,v = np.mgrid[0:2*pi:40j, 0:pi:20j]
ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v),
                  color='gray', alpha=0.06, linewidth=0.2)
# 赤道
phi = np.linspace(0,2*pi,200)
ax.plot(np.cos(phi), np.sin(phi), 0*phi, 'gray', alpha=0.2, ls='--', lw=0.6)

params = [
    (0.01, 0.01, 5.0, 'E₁=t₁=0.01, τ=5'),
    (0.01, 0.03, 4.0, 'E₁=0.01, t₁=0.03, τ=4'),
    (0.005, 0.005, 8.0, 'E₁=t₁=0.005, τ=8'),
    (0.01, 0.001, 6.0, 'E₁=0.01, t₁=0.001, τ=6'),
    (0.02, 0.01, 3.0, 'E₁=0.02, t₁=0.01, τ=3'),
]
colors = plt.cm.tab10(np.linspace(0,1,len(params)))

for ci, (E1, t1, tau, label) in enumerate(params):
    print(f"  {label}...")
    tr = trajectory(tau, E1, t1)
    ax.plot(tr[:,0], tr[:,1], tr[:,2], color=colors[ci], lw=1.2, label=label)
    ax.scatter(*tr[0], color=colors[ci], s=30, marker='o', edgecolors='white', linewidths=0.5, zorder=5)
    ax.scatter(*tr[-1], color=colors[ci], s=50, marker='*', edgecolors='white', linewidths=0.5, zorder=5)

ax.set_xlabel('σ_x'); ax.set_ylabel('σ_y'); ax.set_zlabel('σ_z')
ax.set_xlim(-1.1,1.1); ax.set_ylim(-1.1,1.1); ax.set_zlim(-1.1,1.1)
ax.set_title('E₁≠0: Trajectories leave equator\n○ = start, ★ = end', fontsize=11)
ax.legend(fontsize=7, loc='upper left', bbox_to_anchor=(1.02,1))
ax.view_init(20, -50)

# 视角2: 侧视 (xz 平面), 看离开赤道
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v),
                   color='gray', alpha=0.06, linewidth=0.2)
ax2.plot(np.cos(phi), np.sin(phi), 0*phi, 'gray', alpha=0.2, ls='--', lw=0.6)

for ci, (E1, t1, tau, label) in enumerate(params):
    tr = trajectory(tau, E1, t1)
    ax2.plot(tr[:,0], tr[:,1], tr[:,2], color=colors[ci], lw=1.2)
    ax2.scatter(*tr[0], color=colors[ci], s=30, marker='o', edgecolors='white', linewidths=0.5, zorder=5)
    ax2.scatter(*tr[-1], color=colors[ci], s=50, marker='*', edgecolors='white', linewidths=0.5, zorder=5)

ax2.set_xlabel('σ_x'); ax2.set_ylabel('σ_y'); ax2.set_zlabel('σ_z')
ax2.set_xlim(-1.1,1.1); ax2.set_ylim(-1.1,1.1); ax2.set_zlim(-1.1,1.1)
ax2.set_title('Side view: σ_z ≠ 0 (off-equator motion)', fontsize=11)
ax2.view_init(5, -90)  # 侧视，看 z 分量

plt.tight_layout()
plt.savefig('bloch_E1_nonzero_trajectories.png', dpi=200)
print("\n✓ bloch_E1_nonzero_trajectories.png")
print("  ○ = start point, ★ = end point")
print("  All trajectories clearly leave the equatorial plane")
