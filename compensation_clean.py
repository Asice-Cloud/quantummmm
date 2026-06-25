#!/usr/bin/env python3
"""补偿搜索 —— 高精度版。匹配 U(Γ̄)·U(Γ) 本征值与理想"""
import numpy as np

pi = np.pi; tc = 0.3; E0 = 0.3
np.random.seed(42)

def quat_mul(p, q):
    a,b,c,d = p; e,f,g,h = q
    return np.array([a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                     a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e])

def quat_to_mat2x2_c(q):
    a,b,c,d = q
    return np.array([[a+1j*b, c+1j*d], [-c+1j*d, a-1j*b]])

def propagate_one_step(X, Z, Y, W, tau_c, E1, t1, nps, step):
    """RK4, 单步传播"""
    dt = tau_c / nps
    for _ in range(nps):
        t = _ * dt
        fr = 0.5 * (1.0 - np.cos(pi * t / tau_c))
        ff = 0.5 * (1.0 + np.cos(pi * t / tau_c))
        if step == 1:
            t2, t3, Ed, t1v = tc * fr, 0.0, E0 * ff, t1 * fr
        elif step == 2:
            t2, t3, Ed, t1v = tc * ff, tc * fr, 0.0, t1 * ff
        else:
            t2, t3, Ed, t1v = 0.0, tc * ff, E0 * fr, 0.0
        A = np.array([0, (E1 + t3) / 2, t2 / 2, 0])
        D = np.array([0, (-E1 + t3) / 2, t2 / 2, 0])
        B = np.array([abs(t1v) / 2, 0, 0, Ed / 2])
        C = np.array([-abs(t1v) / 2, 0, 0, Ed / 2])

        def f(X_, Z_, Y_, W_):
            return (quat_mul(A, X_) + quat_mul(B, Z_),
                    quat_mul(C, X_) + quat_mul(D, Z_),
                    quat_mul(A, Y_) + quat_mul(B, W_),
                    quat_mul(C, Y_) + quat_mul(D, W_))
        k1X, k1Z, k1Y, k1W = f(X, Z, Y, W)
        k1X *= dt; k1Z *= dt; k1Y *= dt; k1W *= dt
        k2X, k2Z, k2Y, k2W = f(X + 0.5*k1X, Z + 0.5*k1Z, Y + 0.5*k1Y, W + 0.5*k1W)
        k2X *= dt; k2Z *= dt; k2Y *= dt; k2W *= dt
        k3X, k3Z, k3Y, k3W = f(X + 0.5*k2X, Z + 0.5*k2Z, Y + 0.5*k2Y, W + 0.5*k2W)
        k3X *= dt; k3Z *= dt; k3Y *= dt; k3W *= dt
        k4X, k4Z, k4Y, k4W = f(X + k3X, Z + k3Z, Y + k3Y, W + k3W)
        k4X *= dt; k4Z *= dt; k4Y *= dt; k4W *= dt
        X += (k1X + 2*k2X + 2*k3X + k4X) / 6
        Z += (k1Z + 2*k2Z + 2*k3Z + k4Z) / 6
        Y += (k1Y + 2*k2Y + 2*k3Y + k4Y) / 6
        W += (k1W + 2*k2W + 2*k3W + k4W) / 6
    # 重正交化
    c0 = np.hstack([X, Z]); c1 = np.hstack([Y, W])
    n0 = np.sqrt(np.sum(c0**2)); c0 /= n0
    c1 -= np.sum(c0 * c1) * c0; c1 /= np.sqrt(np.sum(c1**2))
    return c0[:4], c0[4:], c1[:4], c1[4:]

def propagate_protocol(tau, E1, t1, nps=2000):
    """标准 3 步协议"""
    tau_c = tau * 100
    X = np.array([1., 0, 0, 0]); Y = np.array([0., 0, 0, 0])
    Z = np.array([0., 0, 0, 0]); W = np.array([1., 0, 0, 0])
    for s in [1, 2, 3]:
        X, Z, Y, W = propagate_one_step(X, Z, Y, W, tau_c, E1, t1, nps, s)
    return X, Y, Z, W

def propagate_general(tau1, tau2, tau3, t1a, t1b, t1c, E1, nps=2000):
    """三步独立 τ 和 t₁"""
    taus = [tau1 * 100, tau2 * 100, tau3 * 100]
    t1s = [t1a, t1b, t1c]
    X = np.array([1., 0, 0, 0]); Y = np.array([0., 0, 0, 0])
    Z = np.array([0., 0, 0, 0]); W = np.array([1., 0, 0, 0])
    for s in [1, 2, 3]:
        X, Z, Y, W = propagate_one_step(X, Z, Y, W, taus[s-1], E1, t1s[s-1], nps, s)
    return X, Y, Z, W

def mul_U(X1, Y1, Z1, W1, X2, Y2, Z2, W2):
    X = quat_mul(X1, X2) + quat_mul(Y1, Z2)
    Y = quat_mul(X1, Y2) + quat_mul(Y1, W2)
    Z = quat_mul(Z1, X2) + quat_mul(W1, Z2)
    W = quat_mul(Z1, Y2) + quat_mul(W1, W2)
    return X, Y, Z, W

def eigenvalues(X, Y, Z, W):
    U4 = np.zeros((4, 4), dtype=complex)
    U4[:2, :2] = quat_to_mat2x2_c(X); U4[:2, 2:] = quat_to_mat2x2_c(Y)
    U4[2:, :2] = quat_to_mat2x2_c(Z); U4[2:, 2:] = quat_to_mat2x2_c(W)
    return np.sort(np.abs(np.angle(np.linalg.eigvals(U4))))

# ═══════════════════════════════════════════════════════════
NPS = 2000  # dt ≈ 0.75 for τ=15

# 理想编织
print("Computing ideal (E₁=0, t₁=0, τ=50)...")
Xid, Yid, Zid, Wid = propagate_protocol(50, 0.0, 0.0, nps=NPS)
ev_ideal = eigenvalues(Xid, Yid, Zid, Wid)
print(f"  ev = {ev_ideal/pi}π")

# ─── 3τ 搜索 ───
print("\n=== 3τ (单次交换) ===")
E1 = 0.1; tau_m = 15; t1_m = 0.01
Xm, Ym, Zm, Wm = propagate_protocol(tau_m, E1, t1_m, nps=NPS)
ev_main = eigenvalues(Xm, Ym, Zm, Wm)
d_main = np.linalg.norm(ev_main - ev_ideal)
print(f"  无补偿: dist={d_main:.4f}, ev={ev_main/pi}π")

best_d = 1e9; best_p = None
for i in range(500):
    tau_v = np.random.uniform(3, 20, 3)
    t1_v = 10**np.random.uniform(np.log10(0.003), np.log10(0.06), 3)
    Xb, Yb, Zb, Wb = propagate_general(*tau_v, *t1_v, E1, nps=NPS)
    Xc, Yc, Zc, Wc = mul_U(Xb, Yb, Zb, Wb, Xm, Ym, Zm, Wm)
    ev = eigenvalues(Xc, Yc, Zc, Wc)
    d = np.linalg.norm(ev - ev_ideal)
    if d < best_d: best_d = d; best_p = (tau_v, t1_v)
    if i % 100 == 0: print(f"  {i}/500: best={best_d:.4f}")

print(f"\n  最优: dist={best_d:.4f} (改善 {d_main/best_d:.1f}x)")
print(f"  τ̄=({best_p[0][0]:.1f}, {best_p[0][1]:.1f}, {best_p[0][2]:.1f})")
print(f"  t̄₁=({best_p[1][0]:.4f}, {best_p[1][1]:.4f}, {best_p[1][2]:.4f})")

# ─── 6τ 搜索 ───
print("\n=== 6τ (双次交换) ===")
Xid2, Yid2, Zid2, Wid2 = mul_U(Xid, Yid, Zid, Wid, Xid, Yid, Zid, Wid)
ev_ideal2 = eigenvalues(Xid2, Yid2, Zid2, Wid2)
print(f"  理想: ev={ev_ideal2/pi}π")

Xm2, Ym2, Zm2, Wm2 = mul_U(Xm, Ym, Zm, Wm, Xm, Ym, Zm, Wm)
ev_main2 = eigenvalues(Xm2, Ym2, Zm2, Wm2)
d_main2 = np.linalg.norm(ev_main2 - ev_ideal2)
print(f"  无补偿: dist={d_main2:.4f}")

best_d2 = 1e9; best_p2 = None
for i in range(500):
    tau_v = np.random.uniform(3, 20, 3)
    t1_v = 10**np.random.uniform(np.log10(0.003), np.log10(0.06), 3)
    Xb, Yb, Zb, Wb = propagate_general(*tau_v, *t1_v, E1, nps=NPS)
    Xc, Yc, Zc, Wc = mul_U(Xb, Yb, Zb, Wb, Xm2, Ym2, Zm2, Wm2)
    ev = eigenvalues(Xc, Yc, Zc, Wc)
    d = np.linalg.norm(ev - ev_ideal2)
    if d < best_d2: best_d2 = d; best_p2 = (tau_v, t1_v)
    if i % 100 == 0: print(f"  {i}/500: best={best_d2:.4f}")

print(f"\n  最优: dist={best_d2:.4f} (改善 {d_main2/best_d2:.1f}x)")
print(f"  τ̄=({best_p2[0][0]:.1f}, {best_p2[0][1]:.1f}, {best_p2[0][2]:.1f})")
print(f"  t̄₁=({best_p2[1][0]:.4f}, {best_p2[1][1]:.4f}, {best_p2[1][2]:.4f})")
