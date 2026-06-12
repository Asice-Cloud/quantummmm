#!/usr/bin/env python3
"""
Reproduce PRB 111, 205411 (2025) Figure 1(d)
============================================
Effective Majorana Hamiltonian braiding simulation.

Hamiltonian (Eq. 2 of the paper):
  H_EM(t) = i Ed γ_a γ_b + i E1 γ_1 γ_2
           + i|t2(t)| γ_a γ_2 - i|t1(t)| γ_b γ_1
           - i|t3(t)| γ_a γ_3

We represent the 5 Majoranas using 3 fermionic modes (dim=8 Hilbert space):
  c1: γ_1 = c1 + c1†,  γ_2 = i(c1 - c1†)
  c2: γ_3 = c2 + c2†
  c3: γ_a = c3 + c3†,  γ_b = i(c3 - c3†)

Braiding protocol (3 steps, repeated for 2 swaps = 6τ total):
  Step 1 (0→τ):   G1 OFF → t2 rises 0→tc, Ed falls E0→0
  Step 2 (τ→2τ):   G1 ON → t3 rises 0→tc, t2 falls tc→0
  Step 3 (2τ→3τ):   G2 ON → t3 falls tc→0, Ed rises 0→E0
  Steps 4-6 repeat for the second swap.

Initial state: |ψ₁₋(0)> = ground of iγ₁γ₂ = |0>₁ ⊗ |0>₂ ⊗ |0>₃
Target  state: |ψ₁₊(0)> = excited of iγ₁γ₂ = |1>₁ ⊗ |0>₂ ⊗ |0>₃

Figure 1(d) plots |<ψ₁₋(6τ)|ψ₁₊(0)>| as a color map
  x-axis: braiding time τ
  y-axis: coupling t1 (max value)
  fixed: E1 = 0.01 meV, tc = E0 = 0.3 meV
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import time

pi = np.pi

# ============================================================
# Physical parameters (all in meV)
# ============================================================
tc = 0.3       # max coupling t2, t3
E0 = 0.3       # initial QD on-site energy Ed
E1_FIXED = 0.01  # fixed E1 for Fig 1(d) (caption says 0.01 meV)

# ============================================================
# 1. Construct fermionic operators in 8-dim Fock space
#    using Jordan-Wigner transformation to ensure proper
#    fermionic anti-commutation relations.
# ============================================================
# Basis: |n₁, n₂, n₃>, n_i ∈ {0,1}
# Index: n₁·4 + n₂·2 + n₃·1

# Pauli matrices
sigma_plus  = np.array([[0., 0.], [1., 0.]])   # σ⁺ = |1><0|
sigma_minus = np.array([[0., 1.], [0., 0.]])   # σ⁻ = |0><1|
sigma_z     = np.array([[1., 0.], [0., -1.]])  # σ_z
I2          = np.eye(2)

# Jordan-Wigner fermionic operators
# c₁ = σ⁻ ⊗ I ⊗ I,         c₁† = σ⁺ ⊗ I ⊗ I
# c₂ = σ_z ⊗ σ⁻ ⊗ I,       c₂† = σ_z ⊗ σ⁺ ⊗ I
# c₃ = σ_z ⊗ σ_z ⊗ σ⁻,     c₃† = σ_z ⊗ σ_z ⊗ σ⁺
c1    = np.kron(np.kron(sigma_minus, I2), I2)
c1_dag = np.kron(np.kron(sigma_plus,  I2), I2)
c2    = np.kron(np.kron(sigma_z, sigma_minus), I2)
c2_dag = np.kron(np.kron(sigma_z, sigma_plus),  I2)
c3    = np.kron(np.kron(sigma_z, sigma_z), sigma_minus)
c3_dag = np.kron(np.kron(sigma_z, sigma_z), sigma_plus)

# Verify fermionic anti-commutation
ops = [("c₁", c1, c1_dag), ("c₂", c2, c2_dag), ("c₃", c3, c3_dag)]
all_ok = True
for ni, ci, cid in ops:
    for nj, cj, cjd in ops:
        # {c_i, c_j†} = δ_ij
        ac = ci @ cjd + cjd @ ci
        expected = np.eye(8) if ni == nj else np.zeros((8, 8))
        if not np.allclose(ac, expected, atol=1e-10):
            print(f"  FAIL: {{{ni}, {nj}†}} ≠ δ")
            all_ok = False
        # {c_i, c_j} = 0
        ac2 = ci @ cj + cj @ ci
        if not np.allclose(ac2, np.zeros((8, 8)), atol=1e-10):
            print(f"  FAIL: {{{ni}, {nj}}} ≠ 0")
            all_ok = False
if all_ok:
    print("Fermionic anti-commutation relations verified ✓")

# Majorana operators: γ = c + c†  or  γ = i(c - c†)
gamma_1 = c1 + c1_dag              # γ₁
gamma_2 = 1j * (c1 - c1_dag)       # γ₂
gamma_3 = c2 + c2_dag              # γ₃
gamma_a = c3 + c3_dag              # γ_a
gamma_b = 1j * (c3 - c3_dag)       # γ_b

# Verify Majorana Hermiticity
for name, op in [("γ₁", gamma_1), ("γ₂", gamma_2), ("γ₃", gamma_3),
                  ("γ_a", gamma_a), ("γ_b", gamma_b)]:
    assert np.allclose(op, op.conj().T), f"{name} not Hermitian!"

# Verify Majorana anti-commutation: {γ_i, γ_j} = 2 δ_{ij}
all_gammas = [gamma_1, gamma_2, gamma_3, gamma_a, gamma_b]
all_names  = ["γ₁", "γ₂", "γ₃", "γ_a", "γ_b"]
all_ok = True
for i, (gi, ni) in enumerate(zip(all_gammas, all_names)):
    for j, (gj, nj) in enumerate(zip(all_gammas, all_names)):
        ac = gi @ gj + gj @ gi
        expected = 2.0 * np.eye(8) if i == j else np.zeros((8, 8))
        if not np.allclose(ac, expected, atol=1e-10):
            print(f"  FAIL: {{{ni}, {nj}}} ≠ expected "
                  f"(diag = {np.round(np.diag(ac)[:4], 6)})")
            all_ok = False
if all_ok:
    print("Majorana anti-commutation relations verified ✓")


def build_hamiltonian(Ed, E1, t1, t2, t3):
    """
    Build H_EM(t) = i Ed γ_a γ_b + i E1 γ_1 γ_2
                   + i t2 γ_a γ_2 - i t1 γ_b γ_1 - i t3 γ_a γ_3

    All parameters are real scalars (coupling magnitudes).
    """
    H = (1j * Ed * (gamma_a @ gamma_b)
         + 1j * E1 * (gamma_1 @ gamma_2)
         + 1j * t2 * (gamma_a @ gamma_2)
         - 1j * t1 * (gamma_b @ gamma_1)
         - 1j * t3 * (gamma_a @ gamma_3))
    # Ensure exact Hermiticity (remove numerical imaginary parts)
    return (H + H.conj().T) * 0.5


# ============================================================
# 2. Braiding protocol: time-dependent couplings
# ============================================================

def pulse_rise(t, tau):
    """Rising pulse: 0 → 1, shape (1-cos(πt/τ))/2"""
    return 0.5 * (1.0 - np.cos(pi * t / tau))

def pulse_fall(t, tau):
    """Falling pulse: 1 → 0, shape (1+cos(πt/τ))/2"""
    return 0.5 * (1.0 + np.cos(pi * t / tau))


def couplings_at_time(t, tau, t1_max):
    """
    Return (Ed, t1, t2, t3) at time t for the 6-step braiding protocol.

    Protocol (each step duration = τ):
      Step 0 (G1 OFF): t2 0→tc, Ed E0→0
      Step 1 (G1 ON):  t3 0→tc, t2 tc→0
      Step 2 (G2 ON):  t3 tc→0, Ed 0→E0
      Steps 3-5: repeat for second swap.

    t1 follows the G1 gate:
      - When G1 is OFF (step 0,3): t1 rises as pulse_rise
      - When G1 is ON  (step 1,4): t1 falls as pulse_fall
      - Otherwise (step 2,5): t1 = 0
    """
    period = 6.0 * tau
    t_mod = t % period
    step = int(t_mod // tau)  # 0, 1, 2, 3, 4, 5
    s = t_mod - step * tau    # local time within step ∈ [0, τ)

    if step == 0 or step == 3:
        # G1 OFF: t2 rises, Ed falls
        Ed_val  = E0 * pulse_fall(s, tau)
        t2_val  = tc * pulse_rise(s, tau)
        t3_val  = 0.0
        t1_val  = t1_max * pulse_rise(s, tau)
    elif step == 1 or step == 4:
        # G1 ON: t3 rises, t2 falls
        Ed_val  = 0.0
        t2_val  = tc * pulse_fall(s, tau)
        t3_val  = tc * pulse_rise(s, tau)
        t1_val  = t1_max * pulse_fall(s, tau)
    else:  # step == 2 or step == 5
        # G2 ON: t3 falls, Ed rises
        Ed_val  = E0 * pulse_rise(s, tau)
        t2_val  = 0.0
        t3_val  = tc * pulse_fall(s, tau)
        t1_val  = 0.0

    return Ed_val, t1_val, t2_val, t3_val


# ============================================================
# 3. Time evolution — RK4 integrator (fast, no matrix exp per step)
# ============================================================

def _schrodinger_rhs(t, psi, tau, t1_max):
    """Right-hand side of i dψ/dt = H(t) ψ → dψ/dt = -i H(t) ψ"""
    Ed, t1, t2, t3 = couplings_at_time(t, tau, t1_max)
    H = build_hamiltonian(Ed, E1_FIXED, t1, t2, t3)
    return -1j * (H @ psi)


def evolve_braiding_rk4(psi0, tau, t1_max, n_steps_per_tau=200):
    """
    Evolve using classical RK4 on the Schrödinger ODE.
    Much faster than expm at each step for 8×8 matrices.
    """
    total_time = 6.0 * tau
    dt = tau / n_steps_per_tau
    n_total = 6 * n_steps_per_tau

    psi = psi0.copy()
    for i in range(n_total):
        t = i * dt
        # RK4 step
        k1 = dt * _schrodinger_rhs(t, psi, tau, t1_max)
        k2 = dt * _schrodinger_rhs(t + 0.5*dt, psi + 0.5*k1, tau, t1_max)
        k3 = dt * _schrodinger_rhs(t + 0.5*dt, psi + 0.5*k2, tau, t1_max)
        k4 = dt * _schrodinger_rhs(t + dt, psi + k3, tau, t1_max)
        psi = psi + (k1 + 2*k2 + 2*k3 + k4) / 6.0

        # Renormalize (RK4 is not strictly unitary)
        psi = psi / np.sqrt(np.real(psi.conj() @ psi))

    return psi


# ============================================================
# 4. Initial and target states
# ============================================================
# The initial state |ψ₁₋(0)> is the ground state of
#   H(0) = i E0 γ_a γ_b + i E1 γ₁ γ₂  (γ₃ is decoupled).
#
#   i γ₁ γ₂ = 1 - 2 c₁†c₁:  eigenvalue -1 → |1>₁ (ground, E1>0)
#   i γ_a γ_b = 2 c₃†c₃ - 1: eigenvalue -1 → |0>₃ (ground, E0>0)
#   γ₃ decoupled: |0>₂ (ground of iγ₃γ₄ with E₂=0)
#
# So initial = |1, 0, 0> = index 4 (binary: n₁=1, n₂=0, n₃=0).
#
# After two swaps of γ₂ and γ₃, the c₁ fermion flips:
#   |ψ₁₋> → |ψ₁₊> = |0>₁ (eigenvalue +1 of iγ₁γ₂).
#
# The fidelity = |<ψ₁₊(0)|ψ₁₋(6τ)>|² is computed by projecting
# |Ψ(6τ)> onto the c₁=|0> subspace (tracing out c₂, c₃):
#   fidelity = Σ_{n₂,n₃} |<0, n₂, n₃| Ψ(6τ)>|²

DIM = 8
# Initial state |Ψ(0)> = |1, 0, 0> = index 4
psi_initial = np.zeros(DIM, dtype=complex)
psi_initial[4] = 1.0  # |100>

# Projector onto c₁ = |0> (indices with n₁=0: 0, 1, 2, 3)
P_c1_zero = np.zeros((DIM, DIM))
for i in [0, 1, 2, 3]:
    P_c1_zero[i, i] = 1.0


def compute_fidelity(psi_final):
    """
    Compute |<ψ₁₊(0)|ψ₁₋(6τ)>|² = probability that c₁ is in |0> state.

    This is done by projecting onto the c₁=|0> subspace and tracing
    out c₂, c₃: ⟨Ψ| (|0><0|₁ ⊗ I₂₃) |Ψ⟩
    """
    return np.real(psi_final.conj() @ P_c1_zero @ psi_final)


# ============================================================
# 5. Quick sanity check: perfect MZM braiding (E1=0, t1=0)
# ============================================================
def sanity_check():
    """Verify that with E1=t1=0, braiding gives perfect NOT gate."""
    tau_test = 500.0   # meV⁻¹
    global E1_FIXED
    E1_saved = E1_FIXED
    E1_FIXED = 0.0     # temporarily set to 0

    psi_final = evolve_braiding_rk4(psi_initial, tau_test, t1_max=0.0,
                                     n_steps_per_tau=300)
    fidelity = compute_fidelity(psi_final)

    E1_FIXED = E1_saved  # restore
    print(f"Sanity check (E1=0, t1=0, τ={tau_test:.0f} meV⁻¹): "
          f"fidelity = {fidelity:.6f}")
    return fidelity


print("Running sanity check for perfect MZM braiding...")
fid_sanity = sanity_check()
print(f"Fidelity (should be ~1): {fid_sanity:.6f}")


# ============================================================
# 6. Parameter sweep for Fig 1(d)
# ============================================================

# τ range: paper uses units of "τ (100/meV)", so τ=1 → 100 meV⁻¹
# We sweep actual τ and label in units of 100/meV.
tau_vals_actual = np.linspace(60.0, 4000.0, 40)   # meV⁻¹
tau_vals_display = tau_vals_actual / 100.0          # units of 100/meV

# t1 range (meV)
t1_vals = np.linspace(0.0, 0.03, 30)

# Use coarser resolution for initial testing, then refine
# For speed: coarser grid first
USE_COARSE = True
if USE_COARSE:
    tau_vals_actual = np.linspace(80.0, 4000.0, 25)
    tau_vals_display = tau_vals_actual / 100.0
    t1_vals = np.linspace(0.0, 0.03, 20)

# Store results
fidelity_grid = np.zeros((len(t1_vals), len(tau_vals_actual)))

print(f"\nSweeping: {len(t1_vals)} t1 values × {len(tau_vals_actual)} τ values...")
print(f"Using RK4 integrator")
t_start = time.time()

for j, t1_max in enumerate(t1_vals):
    for i, tau in enumerate(tau_vals_actual):
        # Use fewer steps for large τ (where dynamics is slower)
        # and more steps for small τ (where fast oscillations matter)
        n_steps = max(100, int(tau / 1.0))  # dt ≈ 1 meV⁻¹
        n_steps = min(n_steps, 2000)  # cap at 2000 steps per τ
        psi_final = evolve_braiding_rk4(psi_initial, tau, t1_max,
                                         n_steps_per_tau=n_steps)
        fidelity_grid[j, i] = compute_fidelity(psi_final)

    # Progress
    if (j + 1) % 5 == 0:
        elapsed = time.time() - t_start
        print(f"  t1 = {t1_max:.4f} meV ({j+1}/{len(t1_vals)}), "
              f"elapsed: {elapsed:.1f}s")

elapsed_total = time.time() - t_start
print(f"Total sweep time: {elapsed_total:.1f}s")


# ============================================================
# 7. Plot — reproduce Fig 1(d) style
# ============================================================

# Custom colormap: dark blue → white → dark red (diverging)
# but paper likely uses a sequential colormap from 0→1
colors_list = [
    (0.00, '#0b1a3d'),   # deep navy
    (0.15, '#1a3a7a'),
    (0.30, '#2d5ec4'),
    (0.45, '#5b9cf5'),
    (0.60, '#8cc4ff'),
    (0.70, '#b8dcff'),
    (0.80, '#d6edff'),
    (0.90, '#f0f8ff'),
    (0.95, '#fff8e0'),
    (1.00, '#ffe080'),
]
cmap = LinearSegmentedColormap.from_list('fig1d', colors_list)

fig, ax = plt.subplots(figsize=(8, 5.5))

# Use pcolormesh for proper grid alignment
X, Y = np.meshgrid(tau_vals_display, t1_vals)
im = ax.pcolormesh(X, Y, fidelity_grid, cmap=cmap, shading='auto',
                    vmin=0.0, vmax=1.0, rasterized=True)

# Colorbar
cbar = plt.colorbar(im, ax=ax, label=r'$|\langle\psi_{1-}(6\tau)|\psi_{1+}(0)\rangle|$',
                     pad=0.02)
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

# Labels
ax.set_xlabel(r'$\tau$ $(100/\mathrm{meV})$', fontsize=13)
ax.set_ylabel(r'$t_1$ (meV)', fontsize=13)
ax.set_title('Fig 1(d) reproduction: Braiding fidelity vs $\\tau$ and $t_1$\n'
             f'$E_1 = {E1_FIXED}$ meV, $t_c = E_0 = {tc}$ meV',
             fontsize=12)

ax.set_xlim(tau_vals_display[0], tau_vals_display[-1])
ax.set_ylim(t1_vals[0], t1_vals[-1])

plt.tight_layout()
plt.savefig('fig1d_reproduction.png', dpi=200)
plt.savefig('fig1d_reproduction.pdf')
print("\nFigure saved: fig1d_reproduction.png, fig1d_reproduction.pdf")
