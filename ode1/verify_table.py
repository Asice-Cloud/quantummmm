#!/usr/bin/env python3
"""验证 solution.md §11 表格中的数据"""
import numpy as np
from scipy.integrate import solve_ivp

pi = np.pi
t_c = 0.3
E0 = 0.3

def fp(t, tau): return 0.5*(1.0+np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1.0-np.cos(pi*t/tau))

def build_A_step1(t, tau, E1v, t1v):
    A = np.zeros((5,5))
    A[0,1]=2*E1v; A[1,0]=-2*E1v
    t2=t_c*fm(t,tau); A[1,3]=-2*t2; A[3,1]=2*t2
    t1=t1v*fm(t,tau); A[0,4]=2*t1; A[4,0]=-2*t1
    Ed=E0*fp(t,tau); A[3,4]=2*Ed; A[4,3]=-2*Ed
    return A

def build_A_step2(t, tau, E1v, t1v):
    A = np.zeros((5,5))
    A[0,1]=2*E1v; A[1,0]=-2*E1v
    t2=t_c*fp(t,tau); A[1,3]=-2*t2; A[3,1]=2*t2
    t3=t_c*fm(t,tau); A[2,3]=2*t3; A[3,2]=-2*t3
    t1=t1v*fp(t,tau); A[0,4]=2*t1; A[4,0]=-2*t1
    return A

def build_A_step3(t, tau, E1v, t1v):
    A = np.zeros((5,5))
    A[0,1]=2*E1v; A[1,0]=-2*E1v
    t3=t_c*fp(t,tau); A[2,3]=2*t3; A[3,2]=-2*t3
    Ed=E0*fm(t,tau); A[3,4]=2*Ed; A[4,3]=-2*Ed
    return A

def propagate_step(build_fn, tau, E1v, t1v):
    def rhs(t,y):
        A=build_fn(t,tau,E1v,t1v); R=y.reshape(5,5); return (A@R).reshape(-1)
    y0=np.eye(5).reshape(-1)
    sol=solve_ivp(rhs,(0,tau),y0,rtol=1e-10,atol=1e-14,method='RK45')
    return sol.y[:,-1].reshape(5,5)

def run_protocol(tau, E1v, t1v):
    R1=propagate_step(build_A_step1,tau,E1v,t1v)
    R2=propagate_step(build_A_step2,tau,E1v,t1v)
    R3=propagate_step(build_A_step3,tau,E1v,t1v)
    return R3@R2@R1

def rotation_angle(R):
    """SO(3) rotation angle from trace"""
    tr = np.trace(R[:3,:3])
    phi = np.arccos(np.clip((tr-1)/2, -1, 1))
    return phi

print("=== 验证 solution.md §11 表格 ===\n")

params = [
    ("E₁=0, t₁=0",      0.0,   0.0),
    ("E₁=0, t₁=0.01",   0.0,   0.01),
    ("E₁=0.01, t₁=0.005", 0.01, 0.005),
    ("E₁=0.3, t₁=0.01", 0.3,   0.01),
]

tau_vals = np.linspace(10, 100, 10)

for label, E1v, t1v in params:
    phis = []
    print(f"--- {label} ---")
    for tau in tau_vals:
        R = run_protocol(tau, E1v, t1v)
        phi = rotation_angle(R)
        phis.append(phi)
        if tau in [10, 20, 50, 100]:
            print(f"  τ={tau:.0f}: φ={phi:.4f}")
    
    phis = np.array(phis)
    unique = len(set(np.round(phis, 4)))
    print(f"  φ range: [{phis.min():.2f}, {phis.max():.2f}]")
    print(f"  unique values (rounded to 4dp): {unique}/10")
    print()

# 也验证 paper expected: φ=π/2 for pure MZM
print("=== 对照: 纯 MZM 期望 φ=π/2≈1.5708 ===")
for tau in [20, 50, 100]:
    R = run_protocol(tau, 0.0, 0.0)
    phi = rotation_angle(R)
    print(f"  τ={tau}: φ={phi:.4f}, diff from π/2 = {abs(phi-pi/2):.2e}")
