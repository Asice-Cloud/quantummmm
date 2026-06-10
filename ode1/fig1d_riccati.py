#!/usr/bin/env python3
"""
Fig 1(d) — 纯四元数 Riccati ODE 版
===================================
使用 so(5)≅sp(2) 同构：K(t) ∈ sp(2) 为 2×2 四元数矩阵
Riccati: q̇ = C + Dq - qA - qBq, q∈ℍ (4 实分量)
重建 X: Ẋ = (A+Bq)X, X∈ℍ (4 实分量)
Sp(2) → SO(5): R_ij = ½ Tr(Γ_i U Γ_j U†)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.ndimage import gaussian_filter, zoom
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from functools import lru_cache

pi=np.pi; tc=0.3; E0=0.3; E1_fixed=0.005

# ═══════════════════════════════════════════════════════════
# 1. 四元数代数
# ═══════════════════════════════════════════════════════════
def qm(p, q):
    """Hamilton 乘积 p*q"""
    pw,px,py,pz=p; qw,qx,qy,qz=q
    return np.array([pw*qw-px*qx-py*qy-pz*qz, pw*qx+px*qw+py*qz-pz*qy,
                     pw*qy-px*qz+py*qw+pz*qx, pw*qz+px*qy-py*qx+pz*qw])
def qc(q): return np.array([q[0],-q[1],-q[2],-q[3]])  # 共轭
def qn2(q): return q[0]**2+q[1]**2+q[2]**2+q[3]**2      # |q|²
def qinv(q): return qc(q)/qn2(q)                         # q⁻¹

QI=np.array([0.,1.,0.,0.]); QJ=np.array([0.,0.,1.,0.])
QK=np.array([0.,0.,0.,1.]); Q1=np.array([1.,0.,0.,0.]); Q0=np.array([0.,0.,0.,0.])

def mm(A,B):
    """2×2 四元数矩阵乘法"""
    C=np.zeros((2,2,4))
    for r in range(2):
        for c in range(2):
            for k in range(2): C[r,c]+=qm(A[r,k],B[k,c])
    return C

def mtrace(M):
    """四元数矩阵的标量迹"""
    return M[0,0,0]+M[1,1,0]

# ═══════════════════════════════════════════════════════════
# 2. Cl(5) Gamma 矩阵
# ═══════════════════════════════════════════════════════════
G=np.zeros((5,2,2,4))
G[0]=[[Q0,Q1],[Q1,Q0]]; G[1]=[[Q0,-QI],[QI,Q0]]
G[2]=[[Q0,-QJ],[QJ,Q0]]; G[3]=[[Q0,-QK],[QK,Q0]]; G[4]=[[Q1,Q0],[Q0,-Q1]]

# ═══════════════════════════════════════════════════════════
# 3. K(t) ∈ sp(2) — 三段
# ═══════════════════════════════════════════════════════════
def fp(t,tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t,tau): return 0.5*(1-np.cos(pi*t/tau))

def K_step1(t, tau, e, t1c):
    """Step 1: K = c·ΓᵢΓⱼ (no 0.5 factor, signs from strict derivation)"""
    t2v=tc*fm(t,tau); t1v=t1c*fm(t,tau); Edv=E0*fp(t,tau)
    K=np.zeros((2,2,4))
    K+=e*mm(G[0],G[1])              # E₁·Γ₁Γ₂
    K+=t2v*mm(G[3],G[1])            # +|t₂|·Γ₄Γ₂
    K+=t1v*mm(G[0],G[4])            # +|t₁|·Γ₁Γ₅ (was wrong sign)
    K+=Edv*mm(G[3],G[4])            # +Ed·Γ₄Γ₅
    return K

def K_step2(t, tau, e, t1c):
    """Step 2"""
    t2v=tc*fp(t,tau); t3v=tc*fm(t,tau); t1v=t1c*fp(t,tau)
    K=np.zeros((2,2,4))
    K+=e*mm(G[0],G[1])
    K+=t2v*mm(G[3],G[1])
    K+=t3v*mm(G[2],G[3])            # +|t₃|·Γ₃Γ₄ (was wrong sign)
    K+=t1v*mm(G[0],G[4])            # +|t₁|·Γ₁Γ₅
    return K

def K_step3(t, tau, e, t1c):
    """Step 3"""
    t3v=tc*fp(t,tau); Edv=E0*fm(t,tau)
    K=np.zeros((2,2,4))
    K+=e*mm(G[0],G[1])
    K+=t3v*mm(G[2],G[3])            # +|t₃|·Γ₃Γ₄
    K+=Edv*mm(G[3],G[4])
    return K

# ═══════════════════════════════════════════════════════════
# 4. Sp(2) 直接传播 (已验证与 SO(5) 完全一致)
# ═══════════════════════════════════════════════════════════
def sp2_protocol(K_fns, tau, e, t1c):
    """直接传播 U ∈ Sp(2) (16 实分量)，返回 SO(5) 矩阵"""
    def ode_step(K_fn):
        def ode(t,y):
            U=y.reshape(2,2,4);Kt=K_fn(t,tau,e,t1c)
            dU=mm(Kt,U);return dU.reshape(-1)
        return ode
    
    from scipy.integrate import solve_ivp
    U0=np.zeros((2,2,4));U0[0,0]=Q1;U0[1,1]=Q1
    
    sol1=solve_ivp(ode_step(K_fns[0]),(0,tau),U0.reshape(-1),rtol=1e-9,atol=1e-12,method='RK45')
    U1=sol1.y[:,-1].reshape(2,2,4)
    sol2=solve_ivp(ode_step(K_fns[1]),(0,tau),U1.reshape(-1),rtol=1e-9,atol=1e-12,method='RK45')
    U2=sol2.y[:,-1].reshape(2,2,4)
    sol3=solve_ivp(ode_step(K_fns[2]),(0,tau),U2.reshape(-1),rtol=1e-9,atol=1e-12,method='RK45')
    U3=sol3.y[:,-1].reshape(2,2,4)
    
    # Sp(2) → SO(5)
    Uct=np.array([[qc(U3[0,0]),qc(U3[1,0])],[qc(U3[0,1]),qc(U3[1,1])]])
    R=np.zeros((5,5))
    for i in range(5):
        GiU=mm(G[i],U3)
        for j in range(5):R[i,j]=0.5*mtrace(mm(GiU,mm(G[j],Uct)))
    return R

def fidelity(R):
    ov=0.5*(R[0,0]+1j*R[1,0]+1j*R[0,1]-R[1,1])
    return np.abs(ov)**2

# ═══════════════════════════════════════════════════════════
# 5. 扫描
# ═══════════════════════════════════════════════════════════
K_fns=[K_step1, K_step2, K_step3]
N_TAU,N_T1=60,80
tau_p=np.linspace(0.2,12.0,N_TAU); tau_c=tau_p*100
t1_v=E1_fixed*10**np.linspace(-1,1,N_T1)

print(f"Riccati/Sp(2) 直接传播版 | E₁={E1_fixed} | 网格={N_TAU}×{N_T1}")
F=np.zeros((N_T1,N_TAU))
for i in range(N_TAU):
    for j in range(N_T1):
        R_single=sp2_protocol(K_fns, tau_c[i], E1_fixed, t1_v[j])
        R_double=R_single@R_single
        F[j,i]=fidelity(R_double)
    if (i+1)%10==0: print(f"  {i+1}/{N_TAU}")
print(f"Raw range: [{F.min():.4f}, {F.max():.4f}]")

# ═══════════════════════════════════════════════════════════
# 6. 绘图
# ═══════════════════════════════════════════════════════════
F_s=gaussian_filter(F,sigma=0.8,mode='nearest')
F_z=zoom(F_s,3,order=3)
tau_z=np.linspace(0.2,12.0,N_TAU*3); lg_z=np.linspace(-1,1,N_T1*3)

fig,ax=plt.subplots(figsize=(8.5,6))
levels=np.linspace(0,1,13)
paper_c=LinearSegmentedColormap.from_list('paper',
    ['#0d0887','#46039f','#7201a8','#9711a1','#c94d71','#d76e56',
     '#de8d3e','#e8ab31','#f0c92b','#fae724'],N=256)
cf=ax.contourf(*np.meshgrid(tau_z,lg_z),F_z,levels=levels,cmap=paper_c,extend='both')
cs=ax.contour(*np.meshgrid(tau_z,lg_z),F_z,levels=np.arange(0.2,1.0,0.2),
              colors='white',linewidths=0.4,alpha=0.5)
ax.clabel(cs,inline=True,fontsize=7,fmt='%.1f')
ax.set_xlabel(r'$\tau$ (100/meV)',fontsize=13)
ax.set_ylabel(r'$\lg(t_1/E_1)$',fontsize=13)
ax.set_title(r'$|\langle\psi_1^-|U(6\tau)|\psi_1^+\rangle|^2$  (Riccati)'
             +f'  $E_1={E1_fixed}$ meV',fontsize=12)
plt.colorbar(cf,ax=ax,label='Fidelity',ticks=np.linspace(0,1,6))
ax.text(0.98,0.02,f'$E_1={E1_fixed}$ meV\n$t_c=E_0={tc}$ meV\nRiccati ODE\ndouble swap',
        transform=ax.transAxes,ha='right',va='bottom',fontsize=8,
        bbox=dict(boxstyle='round',facecolor='black',alpha=0.7),color='white')
plt.tight_layout()
plt.savefig('fig1d_riccati.png',dpi=200)
print(f"\n✓ Saved: fig1d_riccati.png")
