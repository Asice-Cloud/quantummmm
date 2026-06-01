#!/usr/bin/env python3
"""Test: does q = Y X^{-1} satisfy qdot = B + A q + q A - q Bbar q ?"""
import numpy as np
from scipy.integrate import solve_ivp

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

G=np.zeros((5,2,2,4))
G[0]=[[Q0,Q1],[Q1,Q0]]; G[1]=[[Q0,-QI],[QI,Q0]]; G[2]=[[Q0,-QJ],[QJ,Q0]]
G[3]=[[Q0,-QK],[QK,Q0]]; G[4]=[[Q1,Q0],[Q0,-Q1]]
GG={(i,j):mat_mul(G[i-1],G[j-1]) for i in range(1,6) for j in range(1,6) if i!=j}

pi=np.pi; tau=1.0; tc=0.3; E1v=0.3; t1c=0.01
def fp(t): return 0.5*(1+np.cos(pi*t/tau))
def fm(t): return 0.5*(1-np.cos(pi*t/tau))

def K_quat(t):
    K=np.zeros((2,2,4))
    K+=E1v*GG[(1,2)]; K+=tc*fp(t)*GG[(2,4)]
    K+=-t1c*GG[(1,5)]; K+=-tc*fm(t)*GG[(3,4)]
    return K

# 1. Direct propagation
def direct_ode(t,y):
    U=y.reshape(2,2,4); dU=mat_mul(K_quat(t),U); return dU.reshape(-1)
y0=np.zeros((2,2,4)); y0[0,0]=Q1; y0[1,1]=Q1
teval=np.linspace(0,tau,101)
sol_d=solve_ivp(direct_ode,(0,tau),y0.reshape(-1),t_eval=teval,rtol=1e-10,atol=1e-12)

# Extract q_exact(t) = Y(t) X(t)^{-1} (top row of U)
q_exact=[]
for i in range(len(teval)):
    U=sol_d.y[:,i].reshape(2,2,4)
    X=U[0,0]; Y=U[0,1]; X_inv=qc(X)/qnorm2(X)
    q_exact.append(qm(Y,X_inv))

# 2. Riccati ODE: qdot = B + A q + q A - q Bbar q
def riccati_ode(t,q):
    K=K_quat(t); A=K[0,0]; B=K[0,1]; Bb=qc(B)
    return B + qm(A,q) + qm(q,A) - qm(qm(q,Bb),q)

sol_r=solve_ivp(riccati_ode,(0,tau),Q0,t_eval=teval,rtol=1e-10,atol=1e-12)

errs=[np.linalg.norm(q_exact[i]-sol_r.y[:,i]) for i in range(len(teval))]
print(f"Riccati ODE: max err={max(errs):.2e}, final err={errs[-1]:.2e}")
print(f"q_exact = {q_exact[-1]}")
print(f"q_ricc  = {sol_r.y[:,-1]}")
print(f"MATCH: {np.allclose(q_exact[-1], sol_r.y[:,-1], atol=1e-8)}")

# 3. Without nonlinear term
def linear_ode(t,q):
    K=K_quat(t); A=K[0,0]; B=K[0,1]
    return B + qm(A,q) + qm(q,A)

sol_l=solve_ivp(linear_ode,(0,tau),Q0,t_eval=teval,rtol=1e-10,atol=1e-12)
err_l=[np.linalg.norm(q_exact[i]-sol_l.y[:,i]) for i in range(len(teval))]
print(f"Linear ODE: max err={max(err_l):.2e}")
