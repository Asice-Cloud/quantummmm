#!/usr/bin/env python3
"""Sp(2) 纤维分析：研究 ancilla 内部相位 h ∈ U(2) 随 E₁ 的变化"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi; tc = 0.3; E0 = 0.3

def fp(t, tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1-np.cos(pi*t/tau))

def quat_mul(p, q):
    """四元数乘法 p·q, 输入为 (scalar, i, j, k)"""
    a,b,c,d = p; e,f,g,h = q
    return np.array([
        a*e - b*f - c*g - d*h,
        a*f + b*e + c*h - d*g,
        a*g - b*h + c*e + d*f,
        a*h + b*g - c*f + d*e
    ])

def quat_norm(q):
    return np.sqrt(q[0]**2 + q[1]**2 + q[2]**2 + q[3]**2)

def quat_conj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

def quat_inv(q):
    n2 = quat_norm(q)**2
    return quat_conj(q) / n2

# 四元数块 A,B,C,D (从 report.md §3.4)
def get_ABCD(t, tau, E1_val, t1_val):
    t2 = tc * (fm(t,tau) if t < tau else (fp(t,tau) if t < 2*tau else 0))
    t3 = tc * (0 if t < tau else (fm(t,tau) if t < 2*tau else fp(t,tau)))
    Ed = E0 * (fp(t,tau) if t < tau else (0 if t < 2*tau else fm(t,tau)))
    t1 = t1_val * (fm(t,tau) if t < tau else (fp(t,tau) if t < 2*tau else 0))
    
    A = np.array([0, (E1_val+t3)/2, t2/2, 0])         # (0, i, j, k)
    D = np.array([0, (-E1_val+t3)/2, t2/2, 0])
    B = np.array([abs(t1)/2, 0, 0, Ed/2])               # (1, i, j, k) — 标量 + k
    C = np.array([-abs(t1)/2, 0, 0, Ed/2])
    return A, B, C, D

def propagate_sp2(tau_p, E1_val, t1_val, n_per_step=500):
    """传播 U ∈ Sp(2) 三步协议，返回 U(3τ)"""
    tau_c = tau_p * 100
    # U = [X Y; Z W], 每个是四元数 (4-vector)
    X = np.array([1.,0,0,0]); Y = np.array([0.,0,0,0])
    Z = np.array([0.,0,0,0]); W = np.array([1.,0,0,0])
    
    for step in range(3):
        t0 = step * tau_c
        for s in range(n_per_step):
            t = t0 + s * tau_c / n_per_step
            dt = tau_c / n_per_step
            A,B,C,D = get_ABCD(t, tau_c, E1_val, t1_val)
            
            # RK4 on 16 real variables
            def f(X_,Z_,Y_,W_, A_,B_,C_,D_):
                dX = quat_mul(A_, X_) + quat_mul(B_, Z_)
                dZ = quat_mul(C_, X_) + quat_mul(D_, Z_)
                dY = quat_mul(A_, Y_) + quat_mul(B_, W_)
                dW = quat_mul(C_, Y_) + quat_mul(D_, W_)
                return dX, dZ, dY, dW
            
            k1X,k1Z,k1Y,k1W = f(X,Z,Y,W, A,B,C,D)
            k1X*=dt; k1Z*=dt; k1Y*=dt; k1W*=dt
            
            A2,B2,C2,D2 = get_ABCD(t+0.5*dt, tau_c, E1_val, t1_val)
            k2X,k2Z,k2Y,k2W = f(X+0.5*k1X, Z+0.5*k1Z, Y+0.5*k1Y, W+0.5*k1W, A2,B2,C2,D2)
            k2X*=dt; k2Z*=dt; k2Y*=dt; k2W*=dt
            
            k3X,k3Z,k3Y,k3W = f(X+0.5*k2X, Z+0.5*k2Z, Y+0.5*k2Y, W+0.5*k2W, A2,B2,C2,D2)
            k3X*=dt; k3Z*=dt; k3Y*=dt; k3W*=dt
            
            A4,B4,C4,D4 = get_ABCD(t+dt, tau_c, E1_val, t1_val)
            k4X,k4Z,k4Y,k4W = f(X+k3X, Z+k3Z, Y+k3Y, W+k3W, A4,B4,C4,D4)
            k4X*=dt; k4Z*=dt; k4Y*=dt; k4W*=dt
            
            X += (k1X+2*k2X+2*k3X+k4X)/6
            Z += (k1Z+2*k2Z+2*k3Z+k4Z)/6
            Y += (k1Y+2*k2Y+2*k3Y+k4Y)/6
            W += (k1W+2*k2W+2*k3W+k4W)/6
        
        # re-orthogonalize: project U back to Sp(2)
        # Gram-Schmidt on quaternion columns
        c0 = np.hstack([X, Z])
        c1 = np.hstack([Y, W])
        n0 = np.sqrt(np.sum(c0**2))
        c0 /= n0
        # project out
        dot = np.sum(c0 * c1)
        c1 -= dot * c0
        n1 = np.sqrt(np.sum(c1**2))
        c1 /= n1
        X = c0[:4]; Z = c0[4:]
        Y = c1[:4]; W = c1[4:]
    
    return X, Y, Z, W

def decompose_fiber(X, Y, Z, W):
    """从 U = [X Y; Z W] 分解出 q (基) 和 u,v (纤维)"""
    q = quat_mul(Z, quat_inv(X))
    nX = quat_norm(X)
    nW = quat_norm(W)
    u = X / nX if nX > 1e-12 else np.array([1.,0,0,0])
    v = W / nW if nW > 1e-12 else np.array([1.,0,0,0])
    return q, u, v, nX, nW

# ═══════════════════════════════════════════════════════════
# 扫描 E₁
E1_vals = np.concatenate([[0.0], np.logspace(-3, np.log10(0.3), 15), [0.35, 0.4, 0.45, 0.5]])
t1_fixed = 0.01
tau_fixed = 10

results = []
for E1 in E1_vals:
    X,Y,Z,W = propagate_sp2(tau_fixed, E1, t1_fixed)
    q, u, v, nX, nW = decompose_fiber(X, Y, Z, W)
    # u, v 是单位四元数 → 提取旋转角和轴
    theta_u = 2*np.arccos(np.clip(u[0], -1, 1))
    axis_u = u[1:4] / np.sin(theta_u/2) if np.sin(theta_u/2) > 1e-10 else np.zeros(3)
    theta_v = 2*np.arccos(np.clip(v[0], -1, 1))
    axis_v = v[1:4] / np.sin(theta_v/2) if np.sin(theta_v/2) > 1e-10 else np.zeros(3)
    
    results.append({
        'E1': E1,
        'q': q, 'u': u, 'v': v,
        'theta_u': theta_u, 'theta_v': theta_v,
        'axis_u': axis_u, 'axis_v': axis_v,
        '|X|': nX, '|W|': nW
    })
    print(f"  E₁={E1:.4f}: |X|={nX:.4f}, |W|={nW:.4f}, θ_u={theta_u:.4f}, θ_v={theta_v:.4f}")

# ═══════════════════════════════════════════════════════════
# 作图
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
E1a = np.array([r['E1'] for r in results])

# 图1: 纤维旋转角
ax = axes[0,0]
theta_u = np.array([r['theta_u'] for r in results])
theta_v = np.array([r['theta_v'] for r in results])
ax.plot(E1a, theta_u/np.pi, 'o-', ms=4, label=r'$\theta_u / \pi$ (MZM fiber)')
ax.plot(E1a, theta_v/np.pi, 's-', ms=4, label=r'$\theta_v / \pi$ (ancilla fiber)')
ax.set_xscale('symlog', linthresh=0.005)
ax.set_xlabel('E₁ (meV)'); ax.set_ylabel('Fiber rotation angle / π')
ax.set_title(f'Fiber holonomy vs E₁ (τ={tau_fixed}, t₁={t1_fixed})')
ax.legend(); ax.grid(True, alpha=0.3)

# 图2: 纤维旋转轴分量
ax = axes[0,1]
axis_u_arr = np.array([r['axis_u'] for r in results])
axis_v_arr = np.array([r['axis_v'] for r in results])
for i, label in enumerate(['i','j','k']):
    ax.plot(E1a, axis_u_arr[:,i], 'o-', ms=3, alpha=0.6, label=f'u-axis {label}')
    ax.plot(E1a, axis_v_arr[:,i], 's--', ms=3, alpha=0.6, label=f'v-axis {label}')
ax.set_xscale('symlog', linthresh=0.005)
ax.set_xlabel('E₁ (meV)'); ax.set_ylabel('Fiber rotation axis')
ax.set_title('Fiber rotation axes')
ax.legend(fontsize=6, ncol=2); ax.grid(True, alpha=0.3)

# 图3: 基空间 |q|
ax = axes[1,0]
q_norms = np.array([quat_norm(r['q']) for r in results])
ax.plot(E1a, q_norms, 'o-', ms=4)
ax.set_xscale('symlog', linthresh=0.005)
ax.set_xlabel('E₁ (meV)'); ax.set_ylabel('|q|')
ax.set_title('Base space: |q| = |Z|/|X|')
ax.grid(True, alpha=0.3)

# 图4: 纤维-基耦合强度：|X| vs |W|
ax = axes[1,1]
absX = np.array([r['|X|'] for r in results])
absW = np.array([r['|W|'] for r in results])
ax.plot(E1a, absX, 'o-', ms=4, label='|X| (MZM amplitude)')
ax.plot(E1a, absW, 's-', ms=4, label='|W| (ancilla amplitude)')
ax.set_xscale('symlog', linthresh=0.005)
ax.set_xlabel('E₁ (meV)'); ax.set_ylabel('Amplitude')
ax.set_title('Subspace amplitudes')
ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('ode2_5/fig_fiber_analysis.png', dpi=200)
print("\n✓ ode2_5/fig_fiber_analysis.png")

# ═══════════════════════════════════════════════════════════
# 打印关键对比
print("\n=== E₁=0 vs E₁=0.1 key comparison ===")
for E1_target in [0.0, 0.01, 0.1, 0.3]:
    idx = np.argmin(np.abs(E1a - E1_target))
    r = results[idx]
    print(f"E₁={r['E1']:.4f}: θ_u={r['theta_u']:.4f} ({r['theta_u']/np.pi:.3f}π), "
          f"θ_v={r['theta_v']:.4f} ({r['theta_v']/np.pi:.3f}π), |q|={quat_norm(r['q']):.4f}")
