#!/usr/bin/env python3
"""Fig 1(d) — Sp(2) 直接传播 + 固定步长 RK4 + 进度条"""
import numpy as np; pi=np.pi; from time import time
from scipy.ndimage import gaussian_filter, zoom
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

tc=E0=0.3; E1_fixed=0.005

# ── 四元数 ops ──
def qm(p,q):
    pw,px,py,pz=p; qw,qx,qy,qz=q
    return np.array([pw*qw-px*qx-py*qy-pz*qz, pw*qx+px*qw+py*qz-pz*qy,
                     pw*qy-px*qz+py*qw+pz*qx, pw*qz+px*qy-py*qx+pz*qw])
QI=np.array([0.,1.,0.,0.]);QJ=np.array([0.,0.,1.,0.]);QK=np.array([0.,0.,0.,1.])
Q1=np.array([1.,0.,0.,0.]);Q0=np.array([0.,0.,0.,0.])
def qc(q): return np.array([q[0],-q[1],-q[2],-q[3]])
def mm(A,B):
    return np.array([[qm(A[0,0],B[0,0])+qm(A[0,1],B[1,0]), qm(A[0,0],B[0,1])+qm(A[0,1],B[1,1])],
                     [qm(A[1,0],B[0,0])+qm(A[1,1],B[1,0]), qm(A[1,0],B[0,1])+qm(A[1,1],B[1,1])]])
def mtrace(M): return M[0,0,0]+M[1,1,0]

G=np.zeros((5,2,2,4)); G[0]=[[Q0,Q1],[Q1,Q0]]; G[1]=[[Q0,-QI],[QI,Q0]]
G[2]=[[Q0,-QJ],[QJ,Q0]]; G[3]=[[Q0,-QK],[QK,Q0]]; G[4]=[[Q1,Q0],[Q0,-Q1]]

def fp(t,tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t,tau): return 0.5*(1-np.cos(pi*t/tau))

# ── K(t) (fixed: no 0.5, signs corrected) ──
G01=mm(G[0],G[1]); G31=mm(G[3],G[1]); G04=mm(G[0],G[4]); G34=mm(G[3],G[4]); G23=mm(G[2],G[3])
def K1(t,tau,e,t1c):
    K=np.zeros((2,2,4))
    K[:]=e*G01[:]+tc*fm(t,tau)*G31[:]+t1c*fm(t,tau)*G04[:]+E0*fp(t,tau)*G34[:]
    return K
def K2(t,tau,e,t1c):
    K=np.zeros((2,2,4))
    K[:]=e*G01[:]+tc*fp(t,tau)*G31[:]+tc*fm(t,tau)*G23[:]+t1c*fp(t,tau)*G04[:]
    return K
def K3(t,tau,e,t1c):
    K=np.zeros((2,2,4))
    K[:]=e*G01[:]+tc*fp(t,tau)*G23[:]+E0*fm(t,tau)*G34[:]
    return K

# ── Sp(2) RK4 ──
def sp2_rk4_step(K_fn,tau,e,t1c,U0):
    n=max(500,int(2*tau)); dt=tau/n; U=U0.copy()
    for s in range(n):
        t=s*dt
        Kt=K_fn(t,tau,e,t1c)
        k1=mm(Kt,U)
        k2=mm(K_fn(t+0.5*dt,tau,e,t1c), U+0.5*dt*k1)
        k3=mm(K_fn(t+0.5*dt,tau,e,t1c), U+0.5*dt*k2)
        k4=mm(K_fn(t+dt,tau,e,t1c), U+dt*k3)
        U[:]=U[:]+dt/6*(k1[:]+2*k2[:]+2*k3[:]+k4[:])
    return U

def sp2_to_so5(U):
    Uct=np.array([[qc(U[0,0]),qc(U[1,0])],[qc(U[0,1]),qc(U[1,1])]])
    R=np.zeros((5,5))
    for i in range(5):
        GiU=mm(G[i],U)
        for j in range(5):R[i,j]=0.5*mtrace(mm(GiU,mm(G[j],Uct)))
    return R

def fid(R): ov=0.5*(R[0,0]+1j*R[1,0]+1j*R[0,1]-R[1,1]); return np.abs(ov)**2

# ═══════════ scan ═══════════
N_TAU,N_T1=60,80
tau_p=np.linspace(0.2,12,N_TAU); tau_c=tau_p*100
t1_v=E1_fixed*10**np.linspace(-1,1,N_T1)
U0=np.zeros((2,2,4));U0[0,0]=Q1; U0[1,1]=Q1

print(f'Sp(2) RK4 | E₁={E1_fixed} | {N_TAU}×{N_T1}={N_TAU*N_T1}')
F=np.zeros((N_T1,N_TAU))
t_start=time()
for i in range(N_TAU):
    t0=time()
    for j in range(N_T1):
        U1=sp2_rk4_step(K1,tau_c[i],E1_fixed,t1_v[j],U0)
        U2=sp2_rk4_step(K2,tau_c[i],E1_fixed,t1_v[j],U1)
        U3=sp2_rk4_step(K3,tau_c[i],E1_fixed,t1_v[j],U2)
        R=sp2_to_so5(U3); Rd=R@R; F[j,i]=fid(Rd)
    dt_i=time()-t0
    eta=dt_i*(N_TAU-i-1)
    print(f'  {i+1}/{N_TAU}  ({dt_i:.1f}s  ETA {eta/60:.0f}min)')
print(f'Total: {(time()-t_start)/60:.1f}min  Range: [{F.min():.4f},{F.max():.4f}]')

# ═══════════ plot ═══════════
F_s=gaussian_filter(F,sigma=0.8,mode='nearest');F_z=zoom(F_s,3,order=3)
tau_z=np.linspace(0.2,12,N_TAU*3);lg_z=np.linspace(-1,1,N_T1*3)
fig,ax=plt.subplots(figsize=(8.5,6))
levels=np.linspace(0,1,13)
cmap=LinearSegmentedColormap.from_list('p',['#0d0887','#46039f','#7201a8','#9711a1','#c94d71','#d76e56','#de8d3e','#e8ab31','#f0c92b','#fae724'],N=256)
ax.contourf(*np.meshgrid(tau_z,lg_z),F_z,levels=levels,cmap=cmap,extend='both')
cs=ax.contour(*np.meshgrid(tau_z,lg_z),F_z,levels=np.arange(0.2,1,0.2),colors='white',linewidths=0.4,alpha=0.5)
ax.clabel(cs,inline=True,fontsize=7,fmt='%.1f')
ax.set_xlabel(r'$\tau$ (100/meV)',fontsize=13)
ax.set_ylabel(r'$\lg(t_1/E_1)$',fontsize=13)
ax.set_title(f'Sp(2): $|\\langle\\psi_1^-|U(6\\tau)|\\psi_1^+\\rangle|^2$  $E_1={E1_fixed}$ meV',fontsize=12)
plt.colorbar(ax.contourf(*np.meshgrid(tau_z,lg_z),F_z,levels=levels,cmap=cmap,extend='both'),ax=ax,label='Fidelity',ticks=np.linspace(0,1,6))
plt.tight_layout(); plt.savefig('fig1d_riccati.png',dpi=200)
print('✓ fig1d_riccati.png')
