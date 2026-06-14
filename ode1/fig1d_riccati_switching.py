#!/usr/bin/env python3
"""
Fig 1(d) — Riccati ODE with chart switching.
q = Z X⁻¹ evolves via Riccati; when |q| > 1 switches to p = q⁻¹.
Reconstructs U via symplectic condition: W = X̄/|X|², Y = 0.
"""
import numpy as np; pi=np.pi; from time import time
from scipy.ndimage import gaussian_filter
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

tc=E0=0.3; E1_fixed=0.005

# ── Quaternion ops ──
qm=lambda p,q: np.array([p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3],
                         p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2],
                         p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1],
                         p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0]])
qc=lambda q: np.array([q[0],-q[1],-q[2],-q[3]])
qn2=lambda q: q[0]**2+q[1]**2+q[2]**2+q[3]**2
qinv=lambda q: qc(q)/qn2(q)
QI=np.array([0.,1.,0.,0.]);QJ=np.array([0.,0.,1.,0.]);QK=np.array([0.,0.,0.,1.])
Q1=np.array([1.,0.,0.,0.]);Q0=np.array([0.,0.,0.,0.])

def mm(A,B):
    return np.array([[qm(A[0,0],B[0,0])+qm(A[0,1],B[1,0]),
                      qm(A[0,0],B[0,1])+qm(A[0,1],B[1,1])],
                     [qm(A[1,0],B[0,0])+qm(A[1,1],B[1,0]),
                      qm(A[1,0],B[0,1])+qm(A[1,1],B[1,1])]])

G=np.zeros((5,2,2,4)); G[0]=[[Q0,Q1],[Q1,Q0]]; G[1]=[[Q0,-QI],[QI,Q0]]
G[2]=[[Q0,-QJ],[QJ,Q0]]; G[3]=[[Q0,-QK],[QK,Q0]]; G[4]=[[Q1,Q0],[Q0,-Q1]]
G01=mm(G[0],G[1]);G31=mm(G[3],G[1]);G04=mm(G[0],G[4]);G34=mm(G[3],G[4]);G23=mm(G[2],G[3])

def fp(t,tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t,tau): return 0.5*(1-np.cos(pi*t/tau))
K1=lambda t,tau,e,t1c: e*G01+tc*fm(t,tau)*G31+t1c*fm(t,tau)*G04+E0*fp(t,tau)*G34
K2=lambda t,tau,e,t1c: e*G01+tc*fp(t,tau)*G31+tc*fm(t,tau)*G23+t1c*fp(t,tau)*G04
K3=lambda t,tau,e,t1c: e*G01+tc*fp(t,tau)*G23+E0*fm(t,tau)*G34

# ── Riccati step with switching ──
def riccati_step(Kfn,tau,e,t1c):
    n=max(500,int(2*tau));dt=tau/n;q=Q0.copy();X=Q1.copy()
    for _ in range(n):
        t=_*dt;Kt=Kfn(t,tau,e,t1c);A=Kt[0,0];B=Kt[0,1];C=Kt[1,0];D=Kt[1,1]
        if qn2(q)<=1.0:
            y=np.concatenate([q,X])
            def r(y): qq,xx=y[:4],y[4:8];return np.concatenate([C+qm(D,qq)-qm(qq,A)-qm(qm(qq,B),qq),qm(A+qm(B,qq),xx)])
            k1=r(y);k2=r(y+0.5*dt*k1);k3=r(y+0.5*dt*k2);k4=r(y+dt*k3)
            y+=dt/6*(k1+2*k2+2*k3+k4);q,X=y[:4],y[4:8]
        else:
            p=qinv(q);Z=qm(q,X)
            y=np.concatenate([p,Z])
            def rp(y): pp,zz=y[:4],y[4:8];return np.concatenate([-B-qm(pp,D)+qm(A,pp)+qm(qm(pp,C),pp),qm(D+qm(C,pp),zz)])
            k1=rp(y);k2=rp(y+0.5*dt*k1);k3=rp(y+0.5*dt*k2);k4=rp(y+dt*k3)
            y+=dt/6*(k1+2*k2+2*k3+k4);p,Z=y[:4],y[4:8]
            q=qinv(p);X=qm(p,Z)
    return q,X

def riccati_protocol(tau,e,t1c):
    """One swap via Riccati. Returns Sp(2) matrix U."""
    q1,X1=riccati_step(K1,tau,e,t1c);Z1=qm(q1,X1);W1=qc(X1)/qn2(X1)
    q2,X2=riccati_step(K2,tau,e,t1c);Z2=qm(q2,X2);W2=qc(X2)/qn2(X2)
    q3,X3=riccati_step(K3,tau,e,t1c);Z3=qm(q3,X3);W3=qc(X3)/qn2(X3)
    U1=np.array([[X1,Q0],[Z1,W1]]);U2=np.array([[X2,Q0],[Z2,W2]]);U3=np.array([[X3,Q0],[Z3,W3]])
    return mm(mm(U3,U2),U1)

def sp2_to_so5(U):
    Uct=np.array([[qc(U[0,0]),qc(U[1,0])],[qc(U[0,1]),qc(U[1,1])]])
    R=np.zeros((5,5))
    for i in range(5):
        GiU=mm(G[i],U)
        for j in range(5):R[i,j]=0.5*(mm(GiU,mm(G[j],Uct))[0,0,0]+mm(GiU,mm(G[j],Uct))[1,1,0])
    return R
fid=lambda R: abs(0.5*(R[0,0]+1j*R[1,0]+1j*R[0,1]-R[1,1]))**2

# ── Sanity check ──
print('Sanity (E1=0, t1=0)...',end=' ')
U1=riccati_protocol(500,0,0);Ud=mm(U1,U1)
print(f'fidelity={fid(sp2_to_so5(Ud)):.8f}')

# ── Quick verify vs Sp(2) reference ──
def sp2_ref(tau,e,t1c):
    def step(Kfn):
        n=max(500,int(2*tau));dt=tau/n;U=np.zeros((2,2,4));U[0,0]=Q1;U[1,1]=Q1
        for s in range(n):
            t=s*dt;Kt=Kfn(t,tau,e,t1c)
            k1=mm(Kt,U);k2=mm(Kfn(t+0.5*dt,tau,e,t1c),U+0.5*dt*k1)
            k3=mm(Kfn(t+0.5*dt,tau,e,t1c),U+0.5*dt*k2);k4=mm(Kfn(t+dt,tau,e,t1c),U+dt*k3)
            U[:]=U[:]+dt/6*(k1[:]+2*k2[:]+2*k3[:]+k4[:])
        return U
    U1=step(K1);U2=step(K2);U3=step(K3);return mm(mm(U3,U2),U1)

print('Riccati vs Sp(2) ref:')
for tp in [0.5,2,5,8]:
    tau=tp*100;U_r=riccati_protocol(tau,E1_fixed,0);U_s=sp2_ref(tau,E1_fixed,0)
    Ud_r=mm(U_r,U_r);Ud_s=mm(U_s,U_s)
    fr=fid(sp2_to_so5(Ud_r));fs=fid(sp2_to_so5(Ud_s))
    print(f'  tau={tp:.0f}  Riccati={fr:.6f}  Sp2={fs:.6f}  diff={abs(fr-fs):.2e}')

# ── Sweep (quick) ──
N_TAU,N_T1=60,20
tv=np.linspace(0.2,12,N_TAU);t1v=E1_fixed*10**np.linspace(-1,1,N_T1)
grid=np.zeros((N_T1,N_TAU))
print(f'\n{N_T1}×{N_TAU} sweep...');t0=time()
for j,tt in enumerate(t1v):
    for i,tau_p in enumerate(tv):
        tau=tau_p*100
        U=riccati_protocol(tau,E1_fixed,tt);Ud=mm(U,U)
        grid[j,i]=fid(sp2_to_so5(Ud))
    if (j+1)%5==0:print(f'  {j+1}/{N_T1} {time()-t0:.1f}s')
print(f'Done: {time()-t0:.1f}s')

# ── Smooth + Plot ──
grid_s=gaussian_filter(grid,sigma=0.6,mode='nearest')
fig,ax=plt.subplots(figsize=(8,6))
[Tg,Lg]=np.meshgrid(tv,np.linspace(-1,1,N_T1))
colors=['#0d0887','#46039f','#7201a8','#9711a1','#b4298c','#c94d71','#d76e56','#de8d3e','#e8ab31','#f0c92b','#fae724']
cmap=LinearSegmentedColormap.from_list('paper',colors,N=256)
ax.contourf(Tg,Lg,grid_s,levels=np.linspace(0,1,25),cmap=cmap,extend='both')
cs=ax.contour(Tg,Lg,grid_s,levels=[0.1,0.3,0.5,0.7,0.9],colors='white',linewidths=0.5,alpha=0.6)
ax.clabel(cs,inline=True,fontsize=7,fmt='%.1f')
plt.colorbar(ax.collections[0],ax=ax,label=r'$|\langle\psi_{1-}(6\tau)|\psi_{1+}(0)\rangle|$')
ax.set_xlabel(r'$\tau$ (100/meV)',fontsize=13);ax.set_ylabel(r'$\lg(t_1/E_1)$',fontsize=13)
ax.set_title(f'Riccati+Switching Fig 1(d): $E_1$={E1_fixed} meV',fontsize=12)
plt.tight_layout();plt.savefig('fig1d_riccati_switch.png',dpi=200)
print('Saved: fig1d_riccati_switch.png')
