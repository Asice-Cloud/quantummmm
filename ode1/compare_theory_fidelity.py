#!/usr/bin/env python3
"""
Fidelity: Geometric Theory vs Full SO(5)
=========================================
The geometric fiber-bundle theory IS the SO(5) dynamics, reformulated.
This script validates the theory by comparing:

  1. Full SO(5) fidelity map (exact theory = Fig 1d)
  2. E1=0 analytic formula (from curvature/holonomy analysis)
  3. E1≠0 effect (A=D symmetry breaking → full SU(2) holonomy)

Honest answer: the E1=0 analytic formula works because curvature in the
sigma_z direction vanishes (A=D symmetry). The E1≠0 case has full non-abelian
curvature and requires numerical SO(5) propagation (= the theory's numerics).

Usage:
  python compare_theory_fidelity.py --E1 0.01
  python compare_theory_fidelity.py --E1 0.01 --full
"""
import numpy as np; pi = np.pi
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter, zoom
import argparse

tc=0.3; E0=0.3

def fp(t,tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t,tau): return 0.5*(1-np.cos(pi*t/tau))

def A1(t,tau,e,t1c):
    A=np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
    A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
    A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau); return A
def A2(t,tau,e,t1c):
    A=np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3]=-2*tc*fp(t,tau); A[3,1]=2*tc*fp(t,tau)
    A[2,3]=2*tc*fm(t,tau); A[3,2]=-2*tc*fm(t,tau)
    A[0,4]=2*t1c*fp(t,tau); A[4,0]=-2*t1c*fp(t,tau); return A
def A3(t,tau,e,t1c):
    A=np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[2,3]=2*tc*fp(t,tau); A[3,2]=-2*tc*fp(t,tau)
    A[3,4]=2*E0*fm(t,tau); A[4,3]=-2*E0*fm(t,tau); return A

def so5_step(A_fn,tau,e,t1c,n=500):
    dt=tau/n; R=np.eye(5)
    for _ in range(n):
        t=_*dt; k1=A_fn(t,tau,e,t1c)@R; k2=A_fn(t+0.5*dt,tau,e,t1c)@(R+0.5*dt*k1)
        k3=A_fn(t+0.5*dt,tau,e,t1c)@(R+0.5*dt*k2); k4=A_fn(t+dt,tau,e,t1c)@(R+dt*k3)
        R+=dt/6*(k1+2*k2+2*k3+k4)
    return R

def so5_fidelity(tau,e,t1c,n=500):
    Rs=so5_step(A3,tau,e,t1c,n)@so5_step(A2,tau,e,t1c,n)@so5_step(A1,tau,e,t1c,n)
    Rd=Rs@Rs; ov=0.5*(Rd[0,0]+1j*Rd[1,0]+1j*Rd[0,1]-Rd[1,1])
    return np.abs(ov)**2

def analytic_e1zero_fid(tau, t1c, k=1.644):
    """E1=0 analytic: U=exp(-i*phi*nhat*sigma), phi=sqrt((pi/2)^2+(k*t1c*tau)^2)"""
    Phi_D=k*t1c*tau; phi=np.sqrt((pi/2)**2+Phi_D**2)
    alpha=np.arctan2(Phi_D,pi/2)
    nx,ny=np.cos(alpha),np.sin(alpha); th=2*phi; c,s=np.cos(th),np.sin(th)
    Ud=np.array([[c,-s*(ny+1j*nx)],[-s*(ny-1j*nx),c]],dtype=complex)
    S=np.array([[1,1],[1j,-1j]])/np.sqrt(2); Sd=S.conj().T
    Rr=S@Ud@Sd; ov=0.5*(Rr[0,0]+1j*Rr[1,0]+1j*Rr[0,1]-Rr[1,1])
    return np.abs(ov)**2

def run(E1v=0.01,full=False):
    NT,NL=(120,80) if full else (30,20)
    tp=np.linspace(0.2,12,NT); tc_=tp*100
    tv=E1v*10**np.linspace(-1,1,NL)
    ss=0.3 if full else 0.8; zf=2 if full else 3

    print(f"\n{'='*60}"); print(f"Geometric Theory Fidelity  |  E1={E1v}  |  {NT}x{NL}")
    print(f"{'='*60}")

    # SO(5) at E1
    print("\n[1/3] SO(5) E1={}...".format(E1v))
    F5=np.zeros((NL,NT))
    for i in range(NT):
        for j in range(NL): F5[j,i]=so5_fidelity(tc_[i],E1v,tv[j])
        if (i+1)%max(1,NT//4)==0: print(f"  {i+1}/{NT}")

    # SO(5) at E1=0
    print("\n[2/3] SO(5) E1=0...")
    F0=np.zeros((NL,NT))
    for i in range(NT):
        for j in range(NL): F0[j,i]=so5_fidelity(tc_[i],0.0,tv[j])
        if (i+1)%max(1,NT//4)==0: print(f"  {i+1}/{NT}")

    # Analytic E1=0
    k_cal=np.sqrt(1.7726**2-(pi/2)**2)/(0.01*50)
    print(f"\n[3/3] Analytic E1=0 (k={k_cal:.4f})...")
    Fa=np.zeros((NL,NT))
    for i in range(NT):
        for j in range(NL): Fa[j,i]=analytic_e1zero_fid(tc_[i],tv[j],k_cal)
        if (i+1)%max(1,NT//4)==0: print(f"  {i+1}/{NT}")

    # Prep for plotting
    def p(D): return zoom(gaussian_filter(D,sigma=ss,mode='nearest'),zf,order=3)
    F5z=p(F5); F0z=p(F0); Faz=p(Fa)
    tz=np.linspace(0.2,12,NT*zf); lz=np.linspace(-1,1,NL*zf); Tv,Lv=np.meshgrid(tz,lz)
    pc=LinearSegmentedColormap.from_list('p',['#0d0887','#46039f','#7201a8',
        '#9711a1','#c94d71','#d76e56','#de8d3e','#e8ab31','#f0c92b','#fae724'],N=256)

    # ═══ FIGURE 1: 3-panel comparison ═══
    fig,axes=plt.subplots(2,3,figsize=(19,11))
    for ci,(ax,D,t) in enumerate([
        (axes[0,0],F5z,f'SO(5) E1={E1v} (= Fig 1d)'),
        (axes[0,1],F0z,'SO(5) E1=0'),
        (axes[0,2],Faz,'Analytic E1=0 (curvature formula)')]):
        lv=np.linspace(0,1,13); cf=ax.contourf(Tv,Lv,D,levels=lv,cmap=pc,extend='both')
        cs=ax.contour(Tv,Lv,D,levels=np.arange(0.2,1,0.2),colors='white',linewidths=0.4,alpha=0.5)
        ax.clabel(cs,inline=True,fontsize=7,fmt='%.1f')
        ax.set_xlabel(r'$\tau$ (100/meV)'); ax.set_ylabel(r'$\lg(t_1/E_1)$'); ax.set_title(t,fontsize=11)
    plt.colorbar(cf,ax=axes[0,:],label='Fidelity',shrink=0.6,pad=0.02)

    # Differences
    for ci,(ax,D,t) in enumerate([
        (axes[1,0],F0z-Faz,'|Analytic - SO(5) E1=0| (formula error)'),
        (axes[1,1],F5z-F0z,f'|E1={E1v} - E1=0| (E1 effect)')]):
        vm=max(abs(D.min()),abs(D.max()),0.05)
        im=ax.pcolormesh(Tv,Lv,D,cmap=plt.cm.RdBu_r,vmin=-vm,vmax=vm,shading='auto')
        plt.colorbar(im,ax=ax,label='Delta',shrink=0.8)
        ax.set_xlabel(r'$\tau$ (100/meV)'); ax.set_ylabel(r'$\lg(t_1/E_1)$'); ax.set_title(t,fontsize=11)

    # Explanation panel
    ax=axes[1,2]; ax.axis('off')
    lines=[
        ('GEOMETRIC THEORY VALIDATION',True),
        ('',False),
        ('Left: SO(5) = fiber bundle theory (exact).',False),
        ('Center: SO(5) E1=0 = reference.',False),
        ('Right: Analytic formula from curvature.',False),
        ('',False),
        (r'$\phi = \sqrt{(\pi/2)^2 + (k t_1 \tau)^2}$',False),
        (r'$U = e^{-i\phi(\cos\alpha\sigma_x+\sin\alpha\sigma_y)}$',False),
        (f'k = {k_cal:.4f}',False),
        ('',False),
        (f'Max |analytic-SO(5)| = {abs(F0z-Faz).max():.4f}',True),
        (f'Max |E1 effect| = {abs(F5z-F0z).max():.4f}',True),
        ('',False),
        ('E1=0: A=D symmetry → flat σ_z curvature',False),
        ('→ analytic SO(2) holonomy in xy plane.',False),
        ('E1≠0: A≠D → full non-abelian SU(2) curvature',False),
        ('→ requires numerical SO(5) (the theory).',False),
    ]
    for ei,(text,bold) in enumerate(lines):
        ax.text(0.05,0.95-ei*0.043,text,ha='left',fontsize=11 if bold else 9,
                transform=ax.transAxes,fontweight='bold' if bold else 'normal',
                family='monospace' if not bold else 'sans-serif')

    fig.suptitle('Fidelity: Geometric Fiber-Bundle Theory Predictions',fontsize=15,fontweight='bold')
    plt.tight_layout(); fig.savefig('theory_comparison.png',dpi=200); plt.close(fig)
    print("  Saved: theory_comparison.png")

    # ═══ FIGURE 2: Cross-section slices ═══
    fig2,axes2=plt.subplots(2,2,figsize=(13,9))
    for ci,(ax,title,getter,horiz) in enumerate([
        (axes2[0,0],r'$\lg(t_1/E_1)=0$',lambda:NL//2,True),
        (axes2[0,1],r'$\lg(t_1/E_1)=-1$',lambda:0,True),
        (axes2[1,0],r'$\tau=5$',lambda:np.argmin(abs(tp-5)),False),
        (axes2[1,1],r'$\tau=2$',lambda:np.argmin(abs(tp-2)),False)]):
        idx=getter()
        if horiz:
            ax.plot(tp,F5[idx],'b-',lw=1.8,label=f'SO(5) E1={E1v}')
            ax.plot(tp,F0[idx],'g-',lw=1.5,label='SO(5) E1=0')
            ax.plot(tp,Fa[idx],'r--',lw=1.5,label='Analytic E1=0')
            ax.set_xlabel(r'$\tau$ (100/meV)')
        else:
            x=np.linspace(-1,1,NL)
            ax.plot(x,F5[:,idx],'b-',lw=1.8,label=f'SO(5) E1={E1v}')
            ax.plot(x,F0[:,idx],'g-',lw=1.5,label='SO(5) E1=0')
            ax.plot(x,Fa[:,idx],'r--',lw=1.5,label='Analytic E1=0')
            ax.set_xlabel(r'$\lg(t_1/E_1)$')
        ax.set_ylabel('Fidelity'); ax.set_title(title,fontsize=11)
        ax.legend(fontsize=8); ax.grid(True,alpha=0.3); ax.set_ylim(-0.05,1.05)
    fig2.suptitle('Cross-Section Slices',fontsize=14,fontweight='bold')
    plt.tight_layout(); fig2.savefig('theory_comparison_slices.png',dpi=200); plt.close(fig2)
    print("  Saved: theory_comparison_slices.png")

    # Summary
    print(f"\n{'='*60}"); print("Summary")
    print(f"{'='*60}")
    print(f"  E1=0 analytic: phi=sqrt((pi/2)^2+(k*t1*tau)^2), k={k_cal:.4f}")
    print(f"  Max |analytic - SO(5) E1=0| = {abs(F0-Fa).max():.4f}")
    print(f"  Max |E1 effect| = {abs(F5-F0).max():.4f}")
    print(f"  SO(5) E1={E1v} range: [{F5.min():.4f}, {F5.max():.4f}]")
    print(f"{'='*60}")

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--full',action='store_true')
    p.add_argument('--E1',type=float,default=0.01)
    a=p.parse_args(); run(E1v=a.E1,full=a.full); print("\nDone.")
