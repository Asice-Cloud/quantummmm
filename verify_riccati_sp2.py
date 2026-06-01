#!/usr/bin/env python3
"""
Definitive verification: does the quaternion Riccati ODE
exactly reproduce the Sp(2) time-ordered exponential
for the PRB111 Hamiltonian?

Uses PURE quaternion algebra — no complex matrix conversion.
"""
import numpy as np

# ── Quaternion algebra (4-vectors [w, x, y, z]) ──────────────────
def qmul(p, q):
    """Hamilton product p*q."""
    pw, px, py, pz = p
    qw, qx, qy, qz = q
    return np.array([
        pw*qw - px*qx - py*qy - pz*qz,
        pw*qx + px*qw + py*qz - pz*qy,
        pw*qy - px*qz + py*qw + pz*qx,
        pw*qz + px*qy - py*qx + pz*qw
    ])

def qconj(q):
    """Quaternion conjugate q̄."""
    return np.array([q[0], -q[1], -q[2], -q[3]])

def qnorm2(q):
    """|q|^2."""
    return q[0]**2 + q[1]**2 + q[2]**2 + q[3]**2

def qinv(q):
    """q⁻¹ = q̄/|q|²."""
    n2 = qnorm2(q)
    return qconj(q) / n2

# Basis
Q1  = np.array([1., 0., 0., 0.])    # 1
QI  = np.array([0., 1., 0., 0.])    # i
QJ  = np.array([0., 0., 1., 0.])    # j
QK  = np.array([0., 0., 0., 1.])    # k
Q0  = np.array([0., 0., 0., 0.])    # 0

# ── 2×2 quaternion matrix algebra ────────────────────────────────
def mmul(A, B):
    """Multiply two 2×2 quaternion matrices."""
    C = np.zeros((2, 2, 4))
    for r in range(2):
        for c in range(2):
            for k in range(2):
                C[r, c] += qmul(A[r, k], B[k, c])
    return C

def mconjT(M):
    """Conjugate-transpose of 2×2 quaternion matrix."""
    return np.array([[qconj(M[0,0]), qconj(M[1,0])],
                     [qconj(M[0,1]), qconj(M[1,1])]])

# ── Cl(5) Gamma matrices in quaternion 2×2 rep ───────────────────
# {Γ_A, Γ_B} = 2δ_AB I, all Hermitian
G = np.zeros((5, 2, 2, 4))
G[0] = np.array([[Q0, Q1], [Q1, Q0]])    # Γ₁: off-diag 1
G[1] = np.array([[Q0, -QI],[QI, Q0]])    # Γ₂: off-diag -i
G[2] = np.array([[Q0, -QJ],[QJ, Q0]])    # Γ₃: off-diag -j
G[3] = np.array([[Q0, -QK],[QK, Q0]])    # Γ₄: off-diag -k  (γₐ)
G[4] = np.array([[Q1, Q0], [Q0, -Q1]])   # Γ₅: diag (γ_b)

# Verify Clifford algebra
for i in range(5):
    for j in range(5):
        anticom = mmul(G[i], G[j])
        for r in range(2):
            anticom[r,r] += (mmul(G[j], G[i]))[r,r]
        for r in range(2):
            for c in range(2):
                if r == c:
                    expected = 2.0 * (1.0 if i==j else 0.0)
                    # diagonal should be scalar
                    pass
print("Clifford check: G[0]² =", np.allclose(
    (mmul(G[0], G[0]))[0,0], [4.,0.,0.,0.]))  # should be I (but we have off-diag I)

# Actually Γ₁² = [[Q1,Q0],[Q0,Q1]]² = [[Q1·Q1 + Q0·Q1, ...]]
# Let me verify properly with trace
def mtrace(M):
    """Trace of 2×2 quaternion matrix (quaternion-valued)."""
    return M[0,0] + M[1,1]

for i in range(5):
    sq = mmul(G[i], G[i])
    tr = mtrace(sq)
    print(f"  Γ{i+1}² trace = {tr}  (expect [4,0,0,0])")
    # off-diagonal of Γ² should be zero
    print(f"  Γ{i+1}² off-diag = {sq[0,1]}")

# ── so(5) generators: Σ_ij = (1/4)[Γ_i, Γ_j] ─────────────────────
def commutator(A, B):
    return 0.25 * (mmul(A, B) + (-1)*mmul(B, A))

# Actually [Γ, Γ] = ΓΓ - ΓΓ
def comm(A, B):
    AB = mmul(A, B)
    BA = mmul(B, A)
    return np.array([[[AB[r,c][k] - BA[r,c][k] for k in range(4)]
                       for c in range(2)] for r in range(2)])

def gen(i, j):
    """Σ_ij = (1/4)[Γ_i, Γ_j]."""
    return 0.25 * comm(G[i-1], G[j-1])

# ── PRB111 Hamiltonian in Sp(2) quaternion form ──────────────────
pi = np.pi
tau = 1.0
t_c = 0.3
E1_val = 0.3
t1_c = 0.01

def fp(t): return 0.5 * (1.0 + np.cos(pi * t / tau))
def fm(t): return 0.5 * (1.0 - np.cos(pi * t / tau))

def K_quat(t):
    """K(t) = Σ h_ij(t) Σ_ij ∈ sp(2) as 2×2 quaternion matrix.
    
    H_EM = E1·X₁ + |t2|·X₅ - |t1|·X₇ - |t3|·X₆ + Ed·X₁₀
    where X₁=iγ₁γ₂, X₅=iγₐγ₂, X₇=iγ₁γ_b, X₆=iγₐγ₃, X₁₀=iγₐγ_b
    In the so(5) rep: Σ_ij.
    
    Mapping (from derivation):
    γ₁=Γ₁, γ₂=Γ₂, γ₃=Γ₃, γₐ=Γ₄, γ_b=Γ₅
    X₁ = iγ₁γ₂ → Σ₁₂
    X₅ = iγₐγ₂ → Σ₂₄  (iγₐγ₂ = -iγ₂γₐ)
    X₇ = iγ₁γ_b → Σ₁₅  (in H it's -i|t1|γ_bγ₁ = i|t1|γ₁γ_b)
    X₆ = iγₐγ₃ → Σ₃₄  (in H it's -i|t3|γₐγ₃)
    X₁₀ = iγₐγ_b → Σ₄₅
    """
    K = np.zeros((2, 2, 4))
    # E1·X₁ = E1·Σ₁₂
    K += E1_val * gen(1, 2)
    # |t2|·X₅ = |t2|·Σ₂₄
    K += t_c * fp(t) * gen(2, 4)
    # -|t1|·X₇ = -|t1|·Σ₁₅    Wait, in H it's -i|t1|γ_bγ₁
    # But H = +iE1γ₁γ₂ + i|t2|γₐγ₂ - i|t1|γ_bγ₁ - i|t3|γₐγ₃ + iEdγₐγ_b
    # So X₇ = iγ₁γ_b, and coefficient is -|t1|
    # The Hamiltonian is: H = E1·X₁ + |t2|·X₅ + (-|t1|)·X₇ + (-|t3|)·X₆ + Ed·X₁₀
    # where X₅ = iγₐγ₂, X₇ = iγ₁γ_b, X₆ = iγₐγ₃ (wait, need to be careful)
    
    # Let me be very explicit. The paper says:
    # H = iEd γₐγ_b + iE1 γ₁γ₂ + i|t2| γₐγ₂ - i|t1| γ_bγ₁ - i|t3| γₐγ₃
    
    # Rewrite with order: H = Σ c_ij · i·γ_i·γ_j
    # iE1 γ₁γ₂ = E1 · (iγ₁γ₂) → coefficient of Σ₁₂ = E1
    # i|t2| γₐγ₂ = i|t2| γ₂γₐ · (-1) = -i|t2| γ₂γₐ 
    # Actually: iγₐγ₂ uses indices (γₐ is 4, γ₂ is 2 → Σ₂₄ or Σ₄₂?)
    # Σ_ij = (1/4)[Γ_i, Γ_j], and iγ_iγ_j corresponds to (1/2)Γ_iΓ_j = Σ_ij
    
    # Actually let me use: X_ij = i γ_i γ_j maps to what in the spinor rep?
    # In the Clifford algebra: Γ_i Γ_j (no i factor)
    # The Hamiltonian in spinor rep: -i H acts as -i Σ c_ij (iγ_iγ_j)
    # = Σ c_ij γ_i γ_j = Σ c_ij Γ_i Γ_j (in the chosen Clifford rep)
    
    # So K(t) = Σ c_ij(t) · (1/2) Γ_i Γ_j
    # where the coefficient for iγ_iγ_j is c_ij
    
    # Let me rebuild with this:
    # coeffs of iγ_iγ_j:
    # E1: iγ₁γ₂ → (1,2)
    # |t2|: iγₐγ₂ → (4,2) = -(2,4) sign? iγₐγ₂ = -iγ₂γₐ so it's antisymmetric
    # -|t1|: -iγ_bγ₁ = iγ₁γ_b → (1,5)
    # -|t3|: -iγₐγ₃ → (3,4) i.e. iγ₃γₐ
    # Ed: iγₐγ_b → (4,5)
    
    # So K = E1·(iγ₁γ₂) + |t2|·(iγₐγ₂) + (-|t1|)·(iγ₁γ_b) + (-|t3|)·(iγₐγ₃) + Ed·(iγₐγ_b)
    # which is E1·(1/2)Γ₁Γ₂ + |t2|·(1/2)Γ₄Γ₂ + (-|t1|)·(1/2)Γ₁Γ₅ + (-|t3|)·(1/2)Γ₃Γ₄ + Ed·(1/2)Γ₄Γ₅
    
    GG12 = mmul(G[0], G[1])     # Γ₁Γ₂
    GG42 = mmul(G[3], G[1])     # Γ₄Γ₂ (ΓₐΓ₂)
    GG15 = mmul(G[0], G[4])     # Γ₁Γ₅ (Γ₁Γ_b)
    GG34 = mmul(G[2], G[3])     # Γ₃Γ₄ (not Γ₄Γ₃!)
    GG45 = mmul(G[3], G[4])     # Γ₄Γ₅ (ΓₐΓ_b)
    
    K = 0.5 * (E1_val * GG12 
               + t_c * fp(t) * GG42    # iγₐγ₂
               + (-t1_c) * GG15       # iγ₁γ_b
               + (-t_c * fm(t)) * GG34  # -iγₐγ₃ = iγ₃γₐ
               + 0.0 * GG45)           # Ed=0 in step 2
    
    return K

# ── Verify sp(2) structure ────────────────────────────────────────
print("=== Sp(2) structure verification ===")
K0 = K_quat(0.0)
A0 = K0[0, 0]
B0 = K0[0, 1]
C0 = K0[1, 1]

print(f"K(0) = [[A, B], [C, D]]")
print(f"  A = {A0}")
print(f"  B = {B0}")
print(f"  C = {K0[1,0]}")
print(f"  D = {C0}")
print(f"  A pure imag? {np.abs(A0[0]) < 1e-12}  (A[0]={A0[0]:.2e})")
print(f"  D pure imag? {np.abs(C0[0]) < 1e-12}  (D[0]={C0[0]:.2e})")
print(f"  C = -B̄?  {np.allclose(K0[1,0], -qconj(B0), atol=1e-12)}")
print(f"  D = -A?  {np.allclose(C0, np.array([-A0[0],-A0[1],-A0[2],-A0[3]]), atol=1e-12)}")

# ── Direct propagation ────────────────────────────────────────────
from scipy.integrate import solve_ivp

def direct_ode(t, y):
    U = y.reshape(2, 2, 4)
    dU = mmul(K_quat(t), U)
    return dU.reshape(-1)

# ── Riccati ODE ───────────────────────────────────────────────────
# General block form: K = [[A, B], [C, D]]
#   Ẋ = AX + BZ,  Ż = CX + DZ
# Define q = Z X⁻¹:
#   q̇ = C + Dq - qA - qBq
# (This is the general Riccati — does NOT assume sp(2) constraints)
def riccati_ode(t, q):
    Kt = K_quat(t)
    At = Kt[0, 0]    # A
    Bt = Kt[0, 1]    # B
    Ct = Kt[1, 0]    # C
    Dt = Kt[1, 1]    # D
    # q̇ = C + Dq - qA - qBq
    dq = Ct + qmul(Dt, q) - qmul(q, At) - qmul(qmul(q, Bt), q)
    return dq

# ── Run comparison ────────────────────────────────────────────────
dt_small = 0.01
t_eval = np.linspace(0, dt_small, 21)
y0_direct = np.zeros((2, 2, 4))
y0_direct[0, 0] = Q1
y0_direct[1, 1] = Q1

sol_dir = solve_ivp(direct_ode, (0, dt_small), y0_direct.reshape(-1),
                    t_eval=t_eval, rtol=1e-12, atol=1e-14)

# For Riccati, initial q(0) = Z(0) X(0)⁻¹ = 0 · I⁻¹ = 0
q0 = Q0.copy()
sol_ric = solve_ivp(riccati_ode, (0, dt_small), q0,
                    t_eval=t_eval, rtol=1e-12, atol=1e-14)

print("\n=== Riccati vs Direct comparison ===")
max_err = 0.0
for i, ti in enumerate(t_eval):
    Ud = sol_dir.y[:, i].reshape(2, 2, 4)
    Xd, Zd = Ud[0, 0], Ud[1, 0]
    q_direct = qmul(Zd, qinv(Xd))
    q_ricc = sol_ric.y[:, i]
    err = np.linalg.norm(q_direct - q_ricc)
    max_err = max(max_err, err)
    if i % 5 == 0:
        print(f"  t={ti:.3f}: err={err:.2e}")

print(f"\nMax error: {max_err:.2e}")
if max_err < 1e-10:
    print("✓ Riccati ODE exactly reproduces direct Sp(2) propagation!")
else:
    print("✗ Discrepancy detected — check sign conventions")

# ── Also check: is ZX⁻¹ the "right" Riccati variable? ────────────
# Try q2 = YW⁻¹ as well
print("\n=== Alternative: q2 = Y W⁻¹ ===")
for i, ti in enumerate(t_eval):
    Ud = sol_dir.y[:, i].reshape(2, 2, 4)
    Yd, Wd = Ud[0, 1], Ud[1, 1]
    if qnorm2(Wd) > 1e-14:
        q2 = qmul(Yd, qinv(Wd))
        if i % 10 == 0:
            print(f"  t={ti:.3f}: q2={q2}")
    else:
        print(f"  t={ti:.3f}: W singular")

# ── Also scan over longer times ───────────────────────────────────
print("\n=== Longer propagation (up to 0.2) ===")
te_long = np.linspace(0, 0.2, 21)
sol_dl = solve_ivp(direct_ode, (0, 0.2), y0_direct.reshape(-1),
                   t_eval=te_long, rtol=1e-12, atol=1e-14)
sol_rl = solve_ivp(riccati_ode, (0, 0.2), q0,
                   t_eval=te_long, rtol=1e-12, atol=1e-14)

max_err_long = 0.0
for i, ti in enumerate(te_long):
    Ud = sol_dl.y[:, i].reshape(2, 2, 4)
    Xd, Zd = Ud[0, 0], Ud[1, 0]
    q_direct = qmul(Zd, qinv(Xd))
    q_ricc = sol_rl.y[:, i]
    err = np.linalg.norm(q_direct - q_ricc)
    max_err_long = max(max_err_long, err)

print(f"Max error over [0,0.2]: {max_err_long:.2e}")
if max_err_long < 1e-10:
    print("✓ Riccati holds over extended time!")
else:
    print("✗ Longer propagation shows discrepancy")

# ── Save summary ──────────────────────────────────────────────────
print("\n=== Summary ===")
print("General Riccati ODE: q̇ = C + Dq - qA - qBq")
print(f"where q = Z X⁻¹, U = [[X,Y],[Z,W]], U̇ = [[A,B],[C,D]] U")
print(f"If K ∈ sp(2) then C = -B̄, D = -A and this reduces to:")
print(f"  q̇ = -B̄ - Aq - qA - qBq")
