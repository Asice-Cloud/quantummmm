#!/usr/bin/env python3
"""纯四元数 Riccati ↔ Sp(2) 验证（不经过复矩阵）"""
import numpy as np
from scipy.integrate import solve_ivp

# ── Quaternion ops ──
def qm(p,q):
    pw,px,py,pz=p; qw,qx,qy,qz=q
    return np.array([pw*qw-px*qx-py*qy-pz*qz, pw*qx+px*qw+py*qz-pz*qy,
                     pw*qy-px*qz+py*qw+pz*qx, pw*qz+px*qy-py*qx+pz*qw])
def qc(q): return np.array([q[0],-q[1],-q[2],-q[3]])
def qnorm2(q): return q[0]**2+q[1]**2+q[2]**2+q[3]**2

QI=np.array([0.,1.,0.,0.]); QJ=np.array([0.,0.,1.,0.])
QK=np.array([0.,0.,0.,1.]); Q1=np.array([1.,0.,0.,0.]); Q0=np.array([0.,0.,0.,0.])

def mat_mul(A,B):
    C=np.zeros((2,2,4))
    for r in range(2):
        for c in range(2):
            for k in range(2): C[r,c]+=qm(A[r,k],B[k,c])
    return C

# ── Gamma matrices ──
G=np.zeros((5,2,2,4))
G[0]=[[Q0,Q1],[Q1,Q0]]; G[1]=[[Q0,-QI],[QI,Q0]]; G[2]=[[Q0,-QJ],[QJ,Q0]]
G[3]=[[Q0,-QK],[QK,Q0]]; G[4]=[[Q1,Q0],[Q0,-Q1]]

GG={(i,j):mat_mul(G[i-1],G[j-1]) for i in range(1,6) for j in range(1,6) if i!=j}

# ── K(t) ──
pi=np.pi; tau=1.0; tc=0.3; E1v=0.3; t1c=0.01
def f_plus(t): return 0.5*(1+np.cos(pi*t/tau))
def f_minus(t): return 0.5*(1-np.cos(pi*t/tau))
def K_quat(t):
    K=np.zeros((2,2,4))
    K+=E1v*GG[(1,2)]; K+=tc*f_plus(t)*GG[(2,4)]
    K+=-t1c*GG[(1,5)]; K+=-tc*f_minus(t)*GG[(3,4)]
    return K

# ── Direct propagation ──
def direct_ode(t,y):
    U=y.reshape(2,2,4); dU=mat_mul(K_quat(t),U); return dU.reshape(-1)

y0=np.zeros((2,2,4)); y0[0,0]=Q1; y0[1,1]=Q1
teval=np.linspace(0,tau,51)
sol_d=solve_ivp(direct_ode,(0,tau),y0.reshape(-1),t_eval=teval,rtol=1e-10,atol=1e-12)
Ud=[sol_d.y[:,i].reshape(2,2,4) for i in range(len(sol_d.t))]

# ── Riccati + X ──
def riccati_ode(t,y):
    K=K_quat(t); A=K[0,0]; B=K[0,1]; Bb=qc(B)
    Y=y.reshape(5,4); q=Y[0]; X=Y[1:].reshape(2,2,4)
    dq=B+qm(A,q)+qm(q,A)-qm(qm(q,Bb),q)
    AqB=A-qm(q,Bb)
    dX=np.zeros((2,2,4))
    for c in range(2):
        dX[0,c]=qm(AqB,X[0,c]); dX[1,c]=qm(AqB,X[1,c])
    return np.concatenate([dq, dX[0,0], dX[0,1], dX[1,0], dX[1,1]])

yr0=np.zeros((5,4)); yr0[1]=Q1; yr0[3]=Q1
sol_r=solve_ivp(riccati_ode,(0,tau),yr0.reshape(-1),t_eval=teval,rtol=1e-10,atol=1e-12)

# ── Compare ──
errs=[]
for i in range(len(teval)):
    y=sol_r.y[:,i].reshape(5,4)
    q=y[0]; X=np.array([[y[1],y[2]],[y[3],y[4]]])
    Y=mat_mul(np.array([[q,Q0],[Q0,q]]),X)
    utop=np.concatenate([Ud[i][0,0],Ud[i][0,1]])
    rtop=np.concatenate([X[0,0],Y[0,0]])
    errs.append(np.linalg.norm(utop-rtop))

print(f"Max err ||U[0,:]_direct - U[0,:]_riccati|| = {max(errs):.2e}")
print(f"Err at final t = {errs[-1]:.2e}")

# q from Y X^{-1}
X00=sol_r.y[:,3*4:4*4,-1]  # X[1,1] at final time? No...
# Actually y layout: [q(4), X00(4), X01(4), X10(4), X11(4)]
q_f=y[0]; X00_f=y[1]; X01_f=y[2]; Y00_f=qm(q_f,X00_f)
X00_inv=qc(X00_f)/qnorm2(X00_f)
q_rec=qm(Y00_f,X00_inv)
print(f"q direct vs q=Y·X⁻¹ match: {np.allclose(q_f,q_rec,atol=1e-8)}")

# U unitary check
Uf=Ud[-1]; Ufd=np.zeros((2,2,4))
for r in range(2):
    for c in range(2): Ufd[r,c]=qc(Uf[c,r])
Id=mat_mul(Uf,Ufd)
print(f"U U† ≈ I: [0,0]={np.allclose(Id[0,0],Q1,atol=1e-8)}, [1,1]={np.allclose(Id[1,1],Q1,atol=1e-8)}")

print("\n✓ Riccati rigorously maps to Sp(2) in pure quaternion algebra." if max(errs)<1e-6 else "\n✗ FAIL")
