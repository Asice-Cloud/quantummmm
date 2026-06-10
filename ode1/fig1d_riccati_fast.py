#!/usr/bin/env python3
"""Fig 1(d) Riccati ODE 快速版 —— 固定步长 RK4, 8 分量"""
import numpy as np; pi=np.pi
from scipy.ndimage import gaussian_filter, zoom
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

tc=E0=0.3; E1_fixed=0.005

# ── 四元数 ops (行内化提速) ──
def qm(p,q):
    pw,px,py,pz=p; qw,qx,qy,qz=q
    return np.array([pw*qw-px*qx-py*qy-pz*qz, pw*qx+px*qw+py*qz-pz*qy,
                     pw*qy-px*qz+py*qw+pz*qx, pw*qz+px*qy-py*qx+pz*qw])
QI=np.array([0.,1.,0.,0.]);QJ=np.array([0.,0.,1.,0.]);QK=np.array([0.,0.,0.,1.])
Q1=np.array([1.,0.,0.,0.]);Q0=np.array([0.,0.,0.,0.])
def qc(q): return np.array([q[0],-q[1],-q[2],-q[3]])
def qn2(q): return q[0]**2+q[1]**2+q[2]**2+q[3]**2
def mm(A,B):
    C=np.zeros((2,2,4))
    for r in range(2):
        for c in range(2):
            for k in range(2):C[r,c]+=qm(A[r,k],B[k,c])
    return C
def mtrace(M): return M[0,0,0]+M[1,1,0]

G=np.zeros((5,2,2,4)); G[0]=[[Q0,Q1],[Q1,Q0]]; G[1]=[[Q0,-QI],[QI,Q0]]
G[2]=[[Q0,-QJ],[QJ,Q0]]; G[3]=[[Q0,-QK],[QK,Q0]]; G[4]=[[Q1,Q0],[Q0,-Q1]]

def fp(t,tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t,tau): return 0.5*(1-np.cos(pi*t/tau))

# ── K(t) (fixed!) ──
def K1(t,tau,e,t1c):
    K=np.zeros((2,2,4))
    K+=e*mm(G[0],G[1])+tc*fm(t,tau)*mm(G[3],G[1])+t1c*fm(t,tau)*mm(G[0],G[4])+E0*fp(t,tau)*mm(G[3],G[4])
    return K
def K2(t,tau,e,t1c):
    K=np.zeros((2,2,4))
    K+=e*mm(G[0],G[1])+tc*fp(t,tau)*mm(G[3],G[1])+tc*fm(t,tau)*mm(G[2],G[3])+t1c*fp(t,tau)*mm(G[0],G[4])
    return K
def K3(t,tau,e,t1c):
    K=np.zeros((2,2,4))
    K+=e*mm(G[0],G[1])+tc*fp(t,tau)*mm(G[2],G[3])+E0*fm(t,tau)*mm(G[3],G[4])
    return K

# ── Riccati ODE RK4 ──
def riccati_protocol(tau,e,t1c):
    """Riccati ODE q̇=C+Dq-qA-qBq + Ẋ=(A+Bq)X, 8分量, 固定步长 RK4"""
    def rhs(t,y,K_fn):
        Kt=K_fn(t,tau,e,t1c); A=Kt[0,0]; B=Kt[0,1]; C=Kt[1,0]; D=Kt[1,1]
        q=y[:4]; X=y[4:8]
        dq=C+qm(D,q)-qm(q,A)-qm(qm(q,B),q)
        dX=qm(A+qm(B,q),X)
        return np.concatenate([dq,dX])
    
    n=max(500,int(2*tau)); dt=tau/n
    y=np.concatenate([Q0,Q1])  # q=0, X=1
    for K_fn in [K1,K2,K3]:
        for _ in range(n):
            t=_*dt
            k1=rhs(t,y,K_fn); k2=rhs(t+0.5*dt,y+0.5*dt*k1,K_fn)
            k3=rhs(t+0.5*dt,y+0.5*dt*k2,K_fn); k4=rhs(t+dt,y+dt*k3,K_fn)
            y+=dt/6*(k1+2*k2+2*k3+k4)
    
    q=y[:4]; X=y[4:8]; Z=qm(q,X)
    # Sp(2)→SO(5) (reconstruct full U)
    Xc=qc(X); Zc=qc(Z); nrm=np.sqrt(qn2(X)+qn2(Z)); Xn=X/nrm; Zn=Z/nrm
    Yn=np.array([-Zn[0],-Zn[1],-Zn[2],-Zn[3]]); Wn=np.array([Xn[0],Xn[1],Xn[2],Xn[3]])
    U=np.array([[Xn,Yn],[Zn,Wn]])
    Uct=np.array([[qc(U[0,0]),qc(U[1,0])],[qc(U[0,1]),qc(U[1,1])]])
    R=np.zeros((5,5))
    for i in range(5):
        GiU=mm(G[i],U)
        for j in range(5):R[i,j]=0.5*mtrace(mm(GiU,mm(G[j],Uct)))
    return R

def fid(R): ov=0.5*(R[0,0]+1j*R[1,0]+1j*R[0,1]-R[1,1]); return np.abs(ov)**2

# ═══════════ 扫描 ═══════════
N_TAU,N_T1=60,80
tau_p=np.linspace(0.2,12,N_TAU); tau_c=tau_p*100
t1_v=E1_fixed*10**np.linspace(-1,1,N_T1)
print(f'Riccati RK4 | E₁={E1_fixed} | {N_TAU}×{N_T1}')
F=np.zeros((N_T1,N_TAU))
for i in range(N_TAU):
    for j in range(N_T1):
        R1=riccati_protocol(tau_c[i],E1_fixed,t1_v[j]); Rd=R1@R1; F[j,i]=fid(Rd)
    if (i+1)%10==0: print(f'  {i+1}/{N_TAU}')
print(f'Range: [{F.min():.4f},{F.max():.4f}]')

# ═══════════ 绘图 ═══════════
F_s=gaussian_filter(F,sigma=0.8,mode='nearest'); F_z=zoom(F_s,3,order=3)
tau_z=np.linspace(0.2,12,N_TAU*3); lg_z=np.linspace(-1,1,N_T1*3)
fig,ax=plt.subplots(figsize=(8.5,6))
levels=np.linspace(0,1,13)
paper_c=LinearSegmentedColormap.from_list('p',['#0d0887','#46039f','#7201a8','#9711a1','#c94d71','#d76e56','#de8d3e','#e8ab31','#f0c92b','#fae724'],N=256)
cf=ax.contourf(*np.meshgrid(tau_z,lg_z),F_z,levels=levels,cmap=paper_c,extend='both')
cs=ax.contour(*np.meshgrid(tau_z,lg_z),F_z,levels=np.arange(0.2,1,0.2),colors='white',linewidths=0.4,alpha=0.5)
ax.clabel(cs,inline=True,fontsize=7,fmt='%.1f')
ax.set_xlabel(r'$\tau$ (100/meV)',fontsize=13)
ax.set_ylabel(r'$\lg(t_1/E_1)$',fontsize=13)
ax.set_title(r'Riccati: $|\langle\psi_1^-|U(6\tau)|\psi_1^+\rangle|^2$  $E_1$='+f'{E1_fixed} meV',fontsize=12)
plt.colorbar(cf,ax=ax,label='Fidelity',ticks=np.linspace(0,1,6))
ax.text(0.98,0.02,f'E₁={E1_fixed} meV\ntc=E₀={tc} meV\nRiccati ODE',transform=ax.transAxes,ha='right',va='bottom',fontsize=8,bbox=dict(boxstyle='round',facecolor='black',alpha=0.7),color='white')
plt.tight_layout(); plt.savefig('fig1d_riccati.png',dpi=200)
print('\n✓ fig1d_riccati.png')
