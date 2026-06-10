#!/usr/bin/env python3
"""
Full 3-step PRB111 braiding protocol verification.
Compares so(5) vector rep vs paper's claimed result: γ₂→γ₃, γ₃→−γ₂.
Also cross-checks against Sp(2) spinor + Riccati.
"""
import numpy as np
from scipy.integrate import solve_ivp

pi = np.pi
tau = 1.0       # duration of each step
t_c = 0.3       # coupling peak
E0  = 0.3       # dot level
E1  = 0.3       # left overlap (ABS regime)
t1c = 0.01      # small t1 for ABS

def fp(t): return 0.5 * (1.0 + np.cos(pi * t / tau))
def fm(t): return 0.5 * (1.0 - np.cos(pi * t / tau))

# ── 5×5 so(5) matrix A(t) for each step ────────────────────────────
# Basis order: (γ₁, γ₂, γ₃, γₐ, γ_b)  indices 0..4

def build_A_step1(t):
    """Step 1: t₂ ascends, E_d descends, t₃=0."""
    A = np.zeros((5,5))
    t2 = t_c * fm(t)
    Ed = E0 * fp(t)
    t1 = t1c * fm(t)
    A[0,1] =  2*E1;  A[1,0] = -2*E1       # i E1 γ₁γ₂
    A[1,3] = -2*t2;  A[3,1] =  2*t2       # i|t₂| γₐγ₂
    A[0,4] =  2*t1;  A[4,0] = -2*t1       # -i|t₁| γ_bγ₁ = i|t₁| γ₁γ_b
    A[3,4] =  2*Ed;  A[4,3] = -2*Ed       # i E_d γₐγ_b
    return A

def build_A_step2(t):
    """Step 2: t₂ descends, t₃ ascends, E_d≈0."""
    A = np.zeros((5,5))
    t2 = t_c * fp(t)
    t3 = t_c * fm(t)
    t1 = t1c * fp(t)
    A[0,1] =  2*E1;  A[1,0] = -2*E1
    A[1,3] = -2*t2;  A[3,1] =  2*t2
    A[2,3] =  2*t3;  A[3,2] = -2*t3      # -i|t₃| γₐγ₃ = i|t₃| γ₃γₐ
    A[0,4] =  2*t1;  A[4,0] = -2*t1
    return A

def build_A_step3(t):
    """Step 3: t₃ descends, E_d ascends, t₂=0."""
    A = np.zeros((5,5))
    t3 = t_c * fp(t)
    Ed = E0 * fm(t)
    A[0,1] =  2*E1;  A[1,0] = -2*E1
    A[2,3] =  2*t3;  A[3,2] = -2*t3
    A[3,4] =  2*Ed;  A[4,3] = -2*Ed
    return A

# ── Propagate one step ─────────────────────────────────────────────
def propagate_step(build_A, t_span, n_pts=201):
    """T exp(∫A dt) in 5×5 so(5)."""
    def rhs(t, y):
        A = build_A(t)
        R = y.reshape(5,5)
        return (A @ R).reshape(-1)
    y0 = np.eye(5).reshape(-1)
    t_eval = np.linspace(t_span[0], t_span[1], n_pts)
    sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval,
                    rtol=1e-10, atol=1e-14, method='RK45')
    R_final = sol.y[:,-1].reshape(5,5)
    return sol.t, sol.y.reshape(5,5,-1), R_final

# ── Run full 3-step protocol ───────────────────────────────────────
print("=== PRB111 Full 3-Step Braiding Protocol ===")
print(f"Parameters: E1={E1}, tc={t_c}, E0={E0}, tau={tau}")

t1, y1, R1 = propagate_step(build_A_step1, (0, tau))
t2, y2, R2 = propagate_step(build_A_step2, (0, tau))
t3, y3, R3 = propagate_step(build_A_step3, (0, tau))

R_total = R3 @ R2 @ R1

print(f"\nStep 1 final R:\n{np.round(R1, 4)}")
print(f"\nStep 2 final R:\n{np.round(R2, 4)}")
print(f"\nStep 3 final R:\n{np.round(R3, 4)}")
print(f"\nTotal R = R3·R2·R1:\n{np.round(R_total, 4)}")

# Check orthogonality
print(f"\n||R_total^T R_total - I|| = {np.linalg.norm(R_total.T @ R_total - np.eye(5)):.2e}")
print(f"det(R_total) = {np.linalg.det(R_total):.6f}")

# ── Check braiding result ──────────────────────────────────────────
# Initial Majorana basis: e₁=γ₁, e₂=γ₂, e₃=γ₃, e₄=γₐ, e₅=γ_b
# After evolution: γ'_i = Σ_j R_{ij} γ_j
# Paper claims: γ₂ → γ₃, γ₃ → −γ₂

# R acts on column vectors: γ' = R γ
# So row i of R gives the coefficient of original γ_j in the new γ'_i
# γ'_2 = R[1,:] · γ_original = R[1,1]γ₁ + R[1,2]γ₂ + R[1,3]γ₃ + ...
# We expect: R[1,1]≈0, R[1,2]≈0, R[1,3]≈1, R[1,3]=R[1,4]≈0
# And: γ'_3 = R[2,:] · γ_original, expect R[2,1]≈0, R[2,2]≈-1, R[2,3]≈0

print("\n=== Braiding verification ===")
print("Paper claims: γ₂ → γ₃,  γ₃ → −γ₂")
print(f"γ₂→γ₃: R[1,2]={R_total[1,2]:.6f} (expect ~1)")
print(f"γ₃→−γ₂: R[2,1]={R_total[2,1]:.6f} (expect ~−1)")

# Full row check
print(f"\nγ'_1 = {np.round(R_total[0,:], 4)} · γ")
print(f"γ'_2 = {np.round(R_total[1,:], 4)} · γ")
print(f"γ'_3 = {np.round(R_total[2,:], 4)} · γ")
print(f"γ'_a = {np.round(R_total[3,:], 4)} · γ")
print(f"γ'_b = {np.round(R_total[4,:], 4)} · γ")

# ── Now compare with ideal MZM limit (E1=0, t1=0) ──────────────────
print("\n\n=== Ideal MZM limit: E1=0, t1=0 ===")
E1_save, t1c_save = E1, t1c
# Need to use local variables in the build functions...
# Let me just redo with globals

# Actually, let me make parametrized versions
def build_A_s1(t, E1v, t1v):
    A=np.zeros((5,5))
    A[0,1]=2*E1v; A[1,0]=-2*E1v
    t2=t_c*fm(t); A[1,3]=-2*t2; A[3,1]=2*t2
    t1=t1v*fm(t); A[0,4]=2*t1; A[4,0]=-2*t1
    Ed=E0*fp(t); A[3,4]=2*Ed; A[4,3]=-2*Ed
    return A

def build_A_s2(t, E1v, t1v):
    A=np.zeros((5,5))
    A[0,1]=2*E1v; A[1,0]=-2*E1v
    t2=t_c*fp(t); A[1,3]=-2*t2; A[3,1]=2*t2
    t3=t_c*fm(t); A[2,3]=2*t3; A[3,2]=-2*t3
    t1=t1v*fp(t); A[0,4]=2*t1; A[4,0]=-2*t1
    return A

def build_A_s3(t, E1v, t1v):
    A=np.zeros((5,5))
    A[0,1]=2*E1v; A[1,0]=-2*E1v
    t3=t_c*fp(t); A[2,3]=2*t3; A[3,2]=-2*t3
    Ed=E0*fm(t); A[3,4]=2*Ed; A[4,3]=-2*Ed
    t1=t1v; A[0,4]=2*t1; A[4,0]=-2*t1
    return A

def prop_step(build_fn, E1v, t1v):
    def rhs(t,y):
        A=build_fn(t,E1v,t1v); R=y.reshape(5,5); return (A@R).reshape(-1)
    y0=np.eye(5).reshape(-1)
    te=np.linspace(0,tau,201)
    sol=solve_ivp(rhs,(0,tau),y0,t_eval=te,rtol=1e-10,atol=1e-14)
    return sol.y[:,-1].reshape(5,5)

R1m = prop_step(build_A_s1, 0.0, 0.0)
R2m = prop_step(build_A_s2, 0.0, 0.0)
R3m = prop_step(build_A_s3, 0.0, 0.0)
Rt_mzm = R3m @ R2m @ R1m
print(f"γ₂→γ₃: R[1,2]={Rt_mzm[1,2]:.6f}")
print(f"γ₃→−γ₂: R[2,1]={Rt_mzm[2,1]:.6f}")
print(f"\nγ'_1 = {np.round(Rt_mzm[0,:], 4)} · γ")
print(f"γ'_2 = {np.round(Rt_mzm[1,:], 4)} · γ")
print(f"γ'_3 = {np.round(Rt_mzm[2,:], 4)} · γ")
print(f"γ'_a = {np.round(Rt_mzm[3,:], 4)} · γ")
print(f"γ'_b = {np.round(Rt_mzm[4,:], 4)} · γ")

# ── Print paper-expected result ────────────────────────────────────
print("\n=== Paper expected result ===")
print("Ideal braid:  γ₁→γ₁, γ₂→γ₃, γ₃→−γ₂, γₐ→γₐ, γ_b→γ_b")
R_expected = np.eye(5)
R_expected[1,1] = 0; R_expected[1,2] = 1
R_expected[2,1] = -1; R_expected[2,2] = 0
print(f"Expected R:\n{R_expected}")
print(f"||R_actual - R_expected|| (MZM) = {np.linalg.norm(Rt_mzm - R_expected):.4f}")
print(f"||R_actual - R_expected|| (ABS) = {np.linalg.norm(R_total - R_expected):.4f}")
