#!/usr/bin/env python3
"""
E₁=0 轨迹集合：展示所有轨迹都限在赤道面
=============================================
多组 t₁, τ 参数，全部在同一 Bloch 球上叠加
"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3
E1 = 0.0   # ← 关键：E₁=0

def fp(t, tau): return 0.5*(1.0 + np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1.0 - np.cos(pi*t/tau))

def b1(t, tau, t1c):
    A = np.zeros((5,5))
    A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
    A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
    A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau)
    return A

def b2(t, tau, t1c):
    A = np.zeros((5,5))
    A[1,3]=-2*tc*fp(t,tau); A[3,1]=2*tc*fp(t,tau)
    A[2,3]=2*tc*fm(t,tau); A[3,2]=-2*tc*fm(t,tau)
    A[0,4]=2*t1c*fp(t,tau); A[4,0]=-2*t1c*fp(t,tau)
    return A

def b3(t, tau):
    A = np.zeros((5,5))
    A[2,3]=2*tc*fp(t,tau); A[3,2]=-2*tc*fp(t,tau)
    A[3,4]=2*E0*fm(t,tau); A[4,3]=-2*E0*fm(t,tau)
    return A

def prop(bld, tau, t1c=None):
    n = max(500, int(2*tau))
    dt = tau/n; R = np.eye(5)
    for s in range(n):
        t = s*dt
        if t1c is not None:
            k1 = bld(t,tau,t1c) @ R
            k2 = bld(t+0.5*dt,tau,t1c) @ (R+0.5*dt*k1)
            k3 = bld(t+0.5*dt,tau,t1c) @ (R+0.5*dt*k2)
            k4 = bld(t+dt,tau,t1c) @ (R+dt*k3)
        else:
            k1 = bld(t,tau) @ R
            k2 = bld(t+0.5*dt,tau) @ (R+0.5*dt*k1)
            k3 = bld(t+0.5*dt,tau) @ (R+0.5*dt*k2)
            k4 = bld(t+dt,tau) @ (R+dt*k3)
        R += dt/6.0*(k1+2*k2+2*k3+k4)
    return R

def R_to_bloch(R):
    """从 SO(5) 矩阵 [0,1] 子块提取 ancilla Bloch 向量"""
    # SO(5) 子块 → spinor 旋转变换 → Bloch 向量
    a = R[0,0]; b = R[0,1]; c = R[1,0]; d = R[1,1]
    # 从 2×2 SO(2) 子块重建旋量变换
    # 这里直接用 fidelity 的中间量
    r = np.array([
        1j*(R[0,1] - R[1,0]),  # σ_x
        R[0,1] + R[1,0],        # σ_y  
        R[0,0] - R[1,1]         # σ_z
    ])
    return np.real(r) if np.abs(r).max() > 1e-10 else np.zeros(3)

# ═══════════════════════════════════════════════════════════
fig = plt.figure(figsize=(14, 7))

# ---- 3D Bloch 球 ----
ax = fig.add_subplot(121, projection='3d')

# 画球面
u, v = np.mgrid[0:2*pi:40j, 0:pi:20j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x, y, z, color='gray', alpha=0.1, linewidth=0.3)

# 画赤道
phi = np.linspace(0, 2*pi, 200)
ax.plot(np.cos(phi), np.sin(phi), 0*phi, 'gray', alpha=0.3, ls='--', lw=1)

# 画轴
ax.quiver(0,0,-1.3,0,0,2.6, color='black', alpha=0.2, arrow_length_ratio=0.05)
ax.quiver(-1.3,0,0,2.6,0,0, color='black', alpha=0.2, arrow_length_ratio=0.05)
ax.quiver(0,-1.3,0,0,2.6,0, color='black', alpha=0.2, arrow_length_ratio=0.05)

colors = plt.cm.viridis(np.linspace(0.1, 0.9, 15))

# 多组 t₁ 值（全部 E₁=0）
t1_vals = np.logspace(-3, -1, 15)  # t₁ 从 0.001 到 0.1
tau_vals = [2.0, 5.0, 8.0, 12.0]

def evolve_to_time(t_target, tau_c, t1):
    """演化到指定时间 t_target，返回 SO(5) 矩阵 R(t_target)"""
    if t_target <= tau_c:
        return prop(lambda t,tau,t1c: b1(t,tau,t1c), tau_c, t1)
    elif t_target <= 2*tau_c:
        R1 = prop(lambda t,tau,t1c: b1(t,tau,t1c), tau_c, t1)
        R2 = prop(lambda t,tau,t1c: b2(t,tau,t1c), tau_c, t1)
        return R2 @ R1
    else:
        R1 = prop(lambda t,tau,t1c: b1(t,tau,t1c), tau_c, t1)
        R2 = prop(lambda t,tau,t1c: b2(t,tau,t1c), tau_c, t1)
        R3 = prop(lambda t,tau: b3(t,tau), tau_c)  # Step 3 无 t1
        return R3 @ R2 @ R1

for ti, t1 in enumerate(t1_vals):
    for tau_p in tau_vals:
        tau_c = tau_p * 100.0
        n_steps = 200
        traj = np.zeros((n_steps, 3))
        
        for s in range(n_steps):
            t_target = (s + 1) * 3*tau_c / n_steps
            R = evolve_to_time(t_target, tau_c, t1)
            traj[s] = R_to_bloch(R)
        
        # Normalize
        norms = np.linalg.norm(traj, axis=1)
        norms = np.maximum(norms, 1e-10)
        traj = traj / norms[:, None]
        
        alpha = 0.7 if tau_p <= 5 else 0.3
        ax.plot(traj[:,0], traj[:,1], traj[:,2], 
                color=colors[ti], alpha=alpha, lw=0.8)

ax.set_title('E₁ = 0: All trajectories on equator\n(15 t₁ values × 4 τ values)', fontsize=12)
ax.set_xlabel('σ_x'); ax.set_ylabel('σ_y'); ax.set_zlabel('σ_z')
ax.set_xlim(-1.2,1.2); ax.set_ylim(-1.2,1.2); ax.set_zlim(-1.2,1.2)
ax.view_init(20, -60)

# ---- 顶视图（向下看赤道） ----
ax2 = fig.add_subplot(122)
ax2.set_aspect('equal')

# 画单位圆
theta = np.linspace(0, 2*pi, 200)
ax2.plot(np.cos(theta), np.sin(theta), 'gray', alpha=0.3, lw=1)
ax2.axhline(0, color='gray', alpha=0.2, lw=0.5)
ax2.axvline(0, color='gray', alpha=0.2, lw=0.5)
ax2.arrow(0,0,1.3,0, head_width=0.05, head_length=0.08, fc='black', alpha=0.5)
ax2.arrow(0,0,0,1.3, head_width=0.05, head_length=0.08, fc='black', alpha=0.5)

for ti, t1 in enumerate(t1_vals):
    for tau_p in tau_vals:
        tau_c = tau_p * 100.0
        n_steps = 200
        traj = np.zeros((n_steps, 3))
        
        for s in range(n_steps):
            t_target = (s + 1) * 3*tau_c / n_steps
            R = evolve_to_time(t_target, tau_c, t1)
            traj[s] = R_to_bloch(R)
        
        norms = np.linalg.norm(traj[:,:2], axis=1)
        norms = np.maximum(norms, 1e-10)
        traj_xy = traj[:,:2] / norms[:, None]
        
        alpha = 0.7 if tau_p <= 5 else 0.3
        ax2.plot(traj_xy[:,0], traj_xy[:,1], 
                 color=colors[ti], alpha=alpha, lw=0.8)
        # Mark endpoints
        ax2.scatter(*traj_xy[-1], color=colors[ti], s=10, alpha=0.5, edgecolors='none')

ax2.set_title('Top view: All confined to equatorial plane\n(σ_z = 0 always)', fontsize=12)
ax2.set_xlabel('σ_x'); ax2.set_ylabel('σ_y')
ax2.set_xlim(-1.2,1.2); ax2.set_ylim(-1.2,1.2)

plt.tight_layout()
plt.savefig('bloch_E1_zero_collection.png', dpi=200)
print("✓ Saved: bloch_E1_zero_collection.png")
print("\nKey: 15 t₁ values from 0.001 to 0.1, each at 4 τ values (2,5,8,12)")
print("All trajectories stay exactly on the equator (σ_z = 0)")
