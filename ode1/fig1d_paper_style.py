#!/usr/bin/env python3
"""
Fig 1(d) — 论文风格版
- colormap 仿论文 (blue-green-yellow)
- 更多等高线层级
- 论文匹配标签
"""
import numpy as np
from scipy.ndimage import gaussian_filter, zoom
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

pi=np.pi;tc=0.3;E0=0.3;E1_fixed=0.01

def fp(t,tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t,tau): return 0.5*(1-np.cos(pi*t/tau))

def b1(t,tau,e,t1c):
    A=np.zeros((5,5));A[0,1]=2*e;A[1,0]=-2*e
    t2v=tc*fm(t,tau);A[1,3]=-2*t2v;A[3,1]=2*t2v
    t1v=t1c*fm(t,tau);A[0,4]=2*t1v;A[4,0]=-2*t1v
    Edv=E0*fp(t,tau);A[3,4]=2*Edv;A[4,3]=-2*Edv;return A
def b2(t,tau,e,t1c):
    A=np.zeros((5,5));A[0,1]=2*e;A[1,0]=-2*e
    t2v=tc*fp(t,tau);A[1,3]=-2*t2v;A[3,1]=2*t2v
    t3v=tc*fm(t,tau);A[2,3]=2*t3v;A[3,2]=-2*t3v
    t1v=t1c*fp(t,tau);A[0,4]=2*t1v;A[4,0]=-2*t1v;return A
def b3(t,tau,e,t1c):
    A=np.zeros((5,5));A[0,1]=2*e;A[1,0]=-2*e
    t3v=tc*fp(t,tau);A[2,3]=2*t3v;A[3,2]=-2*t3v
    Edv=E0*fm(t,tau);A[3,4]=2*Edv;A[4,3]=-2*Edv;return A

def prop(bld,tau,e,t1c):
    n=max(500,int(2*tau));dt=tau/n;R=np.eye(5)
    for s in range(n):
        t=s*dt
        k1=bld(t,tau,e,t1c)@R;k2=bld(t+0.5*dt,tau,e,t1c)@(R+0.5*dt*k1)
        k3=bld(t+0.5*dt,tau,e,t1c)@(R+0.5*dt*k2);k4=bld(t+dt,tau,e,t1c)@(R+dt*k3)
        R+=dt/6*(k1+2*k2+2*k3+k4)
    return R

def fid(R):
    ov=0.5*(R[0,0]+1j*R[1,0]+1j*R[0,1]-R[1,1]);return np.abs(ov)**2

N_TAU,N_T1=80,100
tau_p=np.linspace(0.2,12.0,N_TAU);tau_c=tau_p*100
t1_v=E1_fixed*10**np.linspace(-1,1,N_T1)
F=np.zeros((N_T1,N_TAU))
for i in range(N_TAU):
    for j in range(N_T1):
        R1=prop(b1,tau_c[i],E1_fixed,t1_v[j])
        R2=prop(b2,tau_c[i],E1_fixed,t1_v[j])
        R3=prop(b3,tau_c[i],E1_fixed,t1_v[j])
        F[j,i]=fid((R3@R2@R1)@(R3@R2@R1))
    if (i+1)%20==0:print(f'  {i+1}/{N_TAU}')

# 平滑 + 超采样
F_s=gaussian_filter(F,sigma=0.6,mode='nearest')
F_z=zoom(F_s,4,order=3)
tau_z=np.linspace(0.2,12.0,N_TAU*4);lg_z=np.linspace(-1,1,N_T1*4)

# ── 论文风格 colormap ──
colors=['#0d0887','#46039f','#7201a8','#9711a1','#b4298c','#c94d71','#d76e56','#de8d3e','#e8ab31','#f0c92b','#fae724']
cmap=LinearSegmentedColormap.from_list('paper',colors,N=256)

fig,ax=plt.subplots(figsize=(8,6))
Tz,Lz=np.meshgrid(tau_z,lg_z)

levels=np.linspace(0,1,25)
cf=ax.contourf(Tz,Lz,F_z,levels=levels,cmap=cmap,extend='both')

# 白色细等高线
cs=ax.contour(Tz,Lz,F_z,levels=np.linspace(0.1,0.9,9),
              colors='white',linewidths=0.4,alpha=0.5)
ax.clabel(cs,inline=True,fontsize=7,fmt='%.1f')

ax.set_xlabel(r'$\tau$ (100/meV)',fontsize=13)
ax.set_ylabel(r'$\lg(t_1/E_1)$',fontsize=13)
ax.set_title(r'(d) $|\langle\psi_1^-(6\tau)|\psi_1^+(0)\rangle|^2$',fontsize=13)

cbar=plt.colorbar(cf,ax=ax,ticks=np.linspace(0,1,6))
cbar.ax.tick_params(labelsize=10)

ax.text(0.02,0.98,f'$E_1={E1_fixed}$ meV',transform=ax.transAxes,
        fontsize=9,va='top',bbox=dict(boxstyle='round',facecolor='wheat',alpha=0.8))

plt.tight_layout()
plt.savefig('fig1d_final.png',dpi=250)
print(f'\nDone. Range: [{F.min():.4f},{F.max():.4f}]')
