#!/usr/bin/env python3
"""扫描 E₁ 范围，检验补偿方案是否普遍有效"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
        t = s*dt; fr = 0.5*(1-np.cos(pi*t/tau_c)); ff = 0.5*(1+np.cos(pi*t/tau_c))
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

def propagate_3step(tau_p, E1, t1, nps=150):
    tau_c=tau_p*100; X=np.array([1.,0,0,0]);Y=np.array([0.,0,0,0])
    Z=np.array([0.,0,0,0]);W=np.array([1.,0,0,0])
    for s in [1,2,3]: X,Z,Y,W=propagate_step(X,Z,Y,W,tau_c,E1,t1,nps,s)
    return X,Y,Z,W

def safe_propagate_general(tau1,tau2,tau3, t1a,t1b,t1c, E1, nps=None):
    taus=[tau1*100,tau2*100,tau3*100]; t1s=[t1a,t1b,t1c]
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

# 理想编织
Xid,Yid,Zid,Wid = propagate_3step(50, 0.0, 0.0, nps=300)
ev_ideal = evals(Xid,Yid,Zid,Wid)
print(f"Ideal ev: {ev_ideal/pi}π")

# ═══════════════════════════════════════════════════════════
# 扫描 E₁/t₁ 范围
# ═══════════════════════════════════════════════════════════
tau_main = 15; t1_main = 0.01
E1_ratios = np.logspace(-1, 1, 15)  # E₁/t₁: 0.1 → 10
E1_vals = E1_ratios * t1_main
dist_no = []; dist_opt = []

for idx, E1 in enumerate(E1_vals):
    # 主协议
    Xm,Ym,Zm,Wm = propagate_3step(tau_main, E1, t1_main, nps=200)
    ev_m = evals(Xm,Ym,Zm,Wm)
    d_no = np.linalg.norm(ev_m - ev_ideal)
    dist_no.append(d_no)
    
    # 搜索补偿（轻量版）
    best_d = 1e9
    for i in range(300):
        tau_vals = np.random.uniform(3, 20, 3)
        t1_vals = 10**np.random.uniform(np.log10(0.005), np.log10(0.08), 3)
        Xb,Yb,Zb,Wb = safe_propagate_general(*tau_vals, *t1_vals, E1, nps=60)
        if Xb is None: continue
        Xc,Yc,Zc,Wc = mul_U(Xb,Yb,Zb,Wb, Xm,Ym,Zm,Wm)
        ev = evals(Xc,Yc,Zc,Wc)
        d = np.linalg.norm(ev - ev_ideal)
        if d < best_d: best_d = d
    dist_opt.append(best_d)
    print(f"E₁/t₁={E1/t1_main:.2f} (E₁={E1:.4f}): no-comp={d_no:.4f}, opt={best_d:.4f}")

# ═══════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 6))
ax.semilogx(E1_ratios, dist_no, 'o-', ms=6, lw=2, label='No compensation', color='#d62728')
ax.semilogx(E1_ratios, dist_opt, 's-', ms=6, lw=2, label='With compensation', color='#2ca02c')
ax.axhline(y=0.05, color='gray', ls='--', lw=0.8, alpha=0.5)
ax.set_xlabel('E₁ / t₁'); ax.set_ylabel('Eigenvalue distance from ideal')
ax.set_title(f'Compensation effectiveness vs E₁/t₁ (t₁={t1_main}, τ={tau_main})')
ax.legend(fontsize=11); ax.grid(True, alpha=0.3)
ax.set_ylim(0, None)
plt.tight_layout()
plt.savefig('compensation_range.png', dpi=200)
print("\n✓ compensation_range.png")
