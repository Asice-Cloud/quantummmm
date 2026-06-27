#!/usr/bin/env python3
"""双次交换 (6τ) 的补偿搜索"""
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

def propagate_step(X,Z,Y,W, tau_c, E1, t1, nps, step):
    dt = tau_c/nps
    for s in range(nps):
        t = s*dt; fr=0.5*(1-np.cos(pi*t/tau_c)); ff=0.5*(1+np.cos(pi*t/tau_c))
        if step==1: t2=tc*fr; t3=0; Ed=E0*ff; t1v=t1*fr
        elif step==2: t2=tc*ff; t3=tc*fr; Ed=0; t1v=t1*ff
        else: t2=0; t3=tc*ff; Ed=E0*fr; t1v=0
        A=np.array([0,(E1+t3)/2,t2/2,0]); D=np.array([0,(-E1+t3)/2,t2/2,0])
        B=np.array([abs(t1v)/2,0,0,Ed/2]); C=np.array([-abs(t1v)/2,0,0,Ed/2])
        k1X=dt*(quat_mul(A,X)+quat_mul(B,Z)); k1Z=dt*(quat_mul(C,X)+quat_mul(D,Z))
        k1Y=dt*(quat_mul(A,Y)+quat_mul(B,W)); k1W=dt*(quat_mul(C,Y)+quat_mul(D,W))
        Xh=X+0.5*k1X;Zh=Z+0.5*k1Z;Yh=Y+0.5*k1Y;Wh=W+0.5*k1W
        k2X=dt*(quat_mul(A,Xh)+quat_mul(B,Zh)); k2Z=dt*(quat_mul(C,Xh)+quat_mul(D,Zh))
        k2Y=dt*(quat_mul(A,Yh)+quat_mul(B,Wh)); k2W=dt*(quat_mul(C,Yh)+quat_mul(D,Wh))
        Xh=X+0.5*k2X;Zh=Z+0.5*k2Z;Yh=Y+0.5*k2Y;Wh=W+0.5*k2W
        k3X=dt*(quat_mul(A,Xh)+quat_mul(B,Zh)); k3Z=dt*(quat_mul(C,Xh)+quat_mul(D,Zh))
        k3Y=dt*(quat_mul(A,Yh)+quat_mul(B,Wh)); k3W=dt*(quat_mul(C,Yh)+quat_mul(D,Wh))
        Xh=X+k3X;Zh=Z+k3Z;Yh=Y+k3Y;Wh=W+k3W
        k4X=dt*(quat_mul(A,Xh)+quat_mul(B,Zh)); k4Z=dt*(quat_mul(C,Xh)+quat_mul(D,Zh))
        k4Y=dt*(quat_mul(A,Yh)+quat_mul(B,Wh)); k4W=dt*(quat_mul(C,Yh)+quat_mul(D,Wh))
        X+=(k1X+2*k2X+2*k3X+k4X)/6; Z+=(k1Z+2*k2Z+2*k3Z+k4Z)/6
        Y+=(k1Y+2*k2Y+2*k3Y+k4Y)/6; W+=(k1W+2*k2W+2*k3W+k4W)/6
    return X,Z,Y,W

def propagate_Nstep(tau_p, E1, t1, n_rounds=2, nps=150):
    """n_rounds 轮三段协议 (n_rounds=1 是单次, n_rounds=2 是双次)"""
    tau_c=tau_p*100
    X=np.array([1.,0,0,0]);Y=np.array([0.,0,0,0])
    Z=np.array([0.,0,0,0]);W=np.array([1.,0,0,0])
    for _ in range(n_rounds):
        for s in [1,2,3]:
            X,Z,Y,W=propagate_step(X,Z,Y,W,tau_c,E1,t1,nps,s)
    return X,Y,Z,W

def safe_propagate_general(tau1,tau2,tau3,t1a,t1b,t1c,E1,nps=None):
    taus=[tau1*100,tau2*100,tau3*100];t1s=[t1a,t1b,t1c]
    if nps is None: nps=max(80,int(max(taus)/5))
    X=np.array([1.,0,0,0]);Y=np.array([0.,0,0,0])
    Z=np.array([0.,0,0,0]);W=np.array([1.,0,0,0])
    for s in [1,2,3]:
        X,Z,Y,W=propagate_step(X,Z,Y,W,taus[s-1],E1,t1s[s-1],min(nps,200),s)
        if np.any(~np.isfinite(X)): return None
    return X,Y,Z,W

def mul_U(X1,Y1,Z1,W1,X2,Y2,Z2,W2):
    X=quat_mul(X1,X2)+quat_mul(Y1,Z2);Y=quat_mul(X1,Y2)+quat_mul(Y1,W2)
    Z=quat_mul(Z1,X2)+quat_mul(W1,Z2);W=quat_mul(Z1,Y2)+quat_mul(W1,W2)
    return X,Y,Z,W

def evals(X,Y,Z,W):
    U4=np.zeros((4,4),dtype=complex)
    U4[:2,:2]=quat_to_mat2x2_c(X);U4[:2,2:]=quat_to_mat2x2_c(Y)
    U4[2:,:2]=quat_to_mat2x2_c(Z);U4[2:,2:]=quat_to_mat2x2_c(W)
    return np.sort(np.abs(np.angle(np.linalg.eigvals(U4))))

# ═══════════════════════════════════════════════════════════
# 理想双次交换 (E₁=0, t₁=0, τ=50, 6τ)
# ═══════════════════════════════════════════════════════════
Xid,Yid,Zid,Wid = propagate_Nstep(50, 0.0, 0.0, n_rounds=2, nps=300)
ev_ideal = evals(Xid,Yid,Zid,Wid)
print(f"Ideal (6τ, E₁=0): ev={ev_ideal/pi}π")

# 主协议 (E₁≠0, τ=15, 6τ)
E1_target = 0.1; tau_main = 15; t1_main = 0.01
Xm,Ym,Zm,Wm = propagate_Nstep(tau_main, E1_target, t1_main, n_rounds=2, nps=200)
ev_main = evals(Xm,Ym,Zm,Wm)
d_main = np.linalg.norm(ev_main - ev_ideal)
print(f"Main (6τ, E₁={E1_target}): ev={ev_main/pi}π, dist={d_main:.4f}")

# 搜索补偿 Γ̄（3步，独立参数）
best_d = 1e9; best_params = None
for i in range(2000):
    tau_vals = np.random.uniform(3, 18, 3)
    t1_vals = 10**np.random.uniform(np.log10(0.003), np.log10(0.06), 3)
    Xb,Yb,Zb,Wb = safe_propagate_general(*tau_vals, *t1_vals, E1_target, nps=80)
    if Xb is None: continue
    Xc,Yc,Zc,Wc = mul_U(Xb,Yb,Zb,Wb, Xm,Ym,Zm,Wm)
    ev = evals(Xc,Yc,Zc,Wc)
    d = np.linalg.norm(ev - ev_ideal)
    if d < best_d: best_d = d; best_params = (tau_vals, t1_vals)
    if i % 400 == 0:
        print(f"  {i}/2000, best dist={best_d:.4f}")

print(f"\nBest: dist={best_d:.4f} (was {d_main:.4f}, {d_main/best_d:.1f}x better)")
print(f"  τ̄=({best_params[0][0]:.1f},{best_params[0][1]:.1f},{best_params[0][2]:.1f})")
print(f"  t̄₁=({best_params[1][0]:.4f},{best_params[1][1]:.4f},{best_params[1][2]:.4f})")

Xb,Yb,Zb,Wb = safe_propagate_general(*best_params[0], *best_params[1], E1_target, nps=200)
Xc,Yc,Zc,Wc = mul_U(Xb,Yb,Zb,Wb, Xm,Ym,Zm,Wm)
ev_best = evals(Xc,Yc,Zc,Wc)
print(f"  Composite ev: {ev_best/pi}π")
