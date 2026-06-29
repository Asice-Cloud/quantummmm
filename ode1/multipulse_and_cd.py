#!/usr/bin/env python3
"""
多脉冲 Dynamical Decoupling: 同时消除 E1 + t1
================================================
策略: 在编织过程中等间距插入 N 个 π_x 脉冲（翻转 γ₁）
理论: {σ_x, σ_y}={σ_x, σ_z}=0 → 每次 π_x 翻转 E1 和 t1 项的符号
      多个脉冲 → 动力学相位被分段平均抵消 → N→∞ 时完全消除

也尝试 counter-diabatic (CD) driving:
  H_cd(t) = i Σ_n (|∂_t n⟩⟨n| - ⟨n|∂_t n⟩|n⟩⟨n|)
  在瞬时本征基中消除非绝热跃迁，剩余纯 Berry 相位
"""
import numpy as np; pi = np.pi
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

# π_x 脉冲：翻转 γ₁ → R_pi = diag(-1,1,1,1,1)
R_pi = np.eye(5); R_pi[0,0] = -1

def so5_multipulse(tau, e, t1c, n_pulses=0, n_rk4=800):
    """SO(5) propagation with N equally-spaced π_x pulses."""
    dt = tau / n_rk4; R = np.eye(5)
    A_fns = [A1, A2, A3]
    total_steps = 3 * n_rk4
    # Pulse positions (indices): equally spaced
    if n_pulses > 0:
        pulse_positions = np.linspace(0, total_steps-1, n_pulses+2)[1:-1].astype(int)
    else:
        pulse_positions = []

    step_counter = 0
    for step_idx, A_fn in enumerate(A_fns):
        for s in range(n_rk4):
            t = s * dt
            A0 = A_fn(t, tau, e, t1c)
            k1 = A0 @ R
            A_h = A_fn(t+0.5*dt, tau, e, t1c)
            k2 = A_h @ (R+0.5*dt*k1)
            k3 = A_h @ (R+0.5*dt*k2)
            A_f = A_fn(t+dt, tau, e, t1c)
            k4 = A_f @ (R+dt*k3)
            R += dt/6*(k1+2*k2+2*k3+k4)

            if step_counter in pulse_positions:
                R = R_pi @ R
            step_counter += 1

    return R

def fid_so5(R):
    Rd = R @ R; ov = 0.5*(Rd[0,0]+1j*Rd[1,0]+1j*Rd[0,1]-Rd[1,1])
    return np.abs(ov)**2

# ═══════════════════════════════════════════════════════════════
# 1. 多脉冲扫描：N=0,1,2,4,8,16,32 pulses
# ═══════════════════════════════════════════════════════════════
print("="*60)
print("Multi-Pulse Dynamical Decoupling")
print("="*60)

tau = 50.0
test_cases = [
    (0.0, 0.0, "pure MZM"),
    (0.0, 0.01, "E1=0, t1=0.01"),
    (0.01, 0.0, "E1=0.01, t1=0"),
    (0.01, 0.01, "E1=t1=0.01"),
    (0.02, 0.02, "E1=t1=0.02"),
]

n_pulse_list = [0, 1, 2, 4, 8, 16, 32, 64]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
colors = plt.cm.viridis(np.linspace(0, 1, len(test_cases)))

for ci, (e, t1c, label) in enumerate(test_cases):
    fids = []
    for npulse in n_pulse_list:
        R = so5_multipulse(tau, e, t1c, n_pulses=npulse)
        fids.append(fid_so5(R))
        if npulse == 0:
            base_fid = fids[-1]
    fids = np.array(fids)

    # Panel 1: fidelity vs N_pulses
    ax = axes[0,0]
    ax.semilogx(n_pulse_list[1:], fids[1:], 'o-', color=colors[ci], ms=4, lw=1.5, label=label)
    ax.axhline(y=fids[0], color=colors[ci], ls=':', lw=0.6, alpha=0.4)
    ax.set_xlabel('Number of π_x pulses'); ax.set_ylabel('Fidelity')
    ax.set_title('Fidelity vs N_pulses (log scale)', fontsize=11)
    ax.set_ylim(-0.05, 1.1); ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    # Panel 2: deviation from 1
    ax = axes[0,1]
    ax.loglog(n_pulse_list[1:], np.abs(1-fids[1:]), 'o-', color=colors[ci], ms=4, lw=1.5, label=label)
    ax.set_xlabel('Number of π_x pulses'); ax.set_ylabel('|1 - Fidelity|')
    ax.set_title('Convergence to 1 (log-log)', fontsize=11)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

# Panel 3: zoom on best-performing case
ax = axes[1,0]
# E1=t1=0.01 with many pulses
e_demo, t1_demo = 0.01, 0.01
n_list_dense = np.arange(0, 65)
fids_dense = [fid_so5(so5_multipulse(tau, e_demo, t1_demo, n_pulses=n)) for n in n_list_dense]
ax.plot(n_list_dense, fids_dense, 'b-', lw=1.5)
ax.axhline(y=1.0, color='green', ls=':', lw=1, alpha=0.5)
ax.set_xlabel('Number of π_x pulses'); ax.set_ylabel('Fidelity')
ax.set_title(f'E1=t1=0.01: Fidelity vs N (dense scan)', fontsize=11)
ax.set_ylim(-0.05, 1.1); ax.grid(True, alpha=0.3)
# Annotate best
best_n = np.argmax(fids_dense)
ax.annotate(f'N={best_n}, fid={fids_dense[best_n]:.4f}',
            xy=(best_n, fids_dense[best_n]), fontsize=9,
            xytext=(best_n+8, fids_dense[best_n]-0.1),
            arrowprops=dict(arrowstyle='->', color='red'))

# Panel 4: theory explanation
ax = axes[1,1]; ax.axis('off')
lines = [
    ('MULTI-PULSE DYNAMICAL DECOUPLING', True),
    ('', False),
    (r'Each π_x pulse: γ₁ → -γ₁', False),
    (r'  → E1 (iγ₁γ₂) term: sign flipped', False),
    (r'  → t1 (iγ₁γ_b) term: sign flipped', False),
    (r'  → geometric terms: preserved', False),
    ('', False),
    ('More pulses → better averaging', False),
    ('Limit N→∞: perfect cancellation', False),
    ('', False),
    ('Caveat: pulses also disrupt', False),
    ('the ancilla dynamics →', False),
    ('diminishing returns at large N', False),
]
for ei, (text, bold) in enumerate(lines):
    ax.text(0.05, 0.95-ei*0.05, text, ha='left', fontsize=10.5 if bold else 9,
            transform=ax.transAxes, fontweight='bold' if bold else 'normal',
            family='monospace' if not bold else 'sans-serif')

fig.suptitle(f'Dynamical Decoupling: Multi π_x-Pulse Elimination of E1+t1  (τ={tau})',
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig('multipulse_dd.png', dpi=200)
plt.close(fig)
print("  Saved: multipulse_dd.png")

# ═══════════════════════════════════════════════════════════════
# 2. Counter-Diabatic Driving (简化版)
# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("Counter-Diabatic Driving (simplified)")
print("="*60)

def so5_cd(tau, e, t1c, n=800):
    """
    Simplified CD: compute instantaneous eigenstates of the 5x5 Hamiltonian
    (in the antisymmetric generator representation), then add H_cd.
    """
    dt = tau / n; R = np.eye(5)
    A_fns = [A1, A2, A3]

    for step_idx, A_fn in enumerate(A_fns):
        # Pre-compute instantaneous A(t) and its eigen-decomposition
        # H_inst = -i * (antisymmetric A in Majorana rep)
        # For CD, we need |∂_t ψ_n⟩ = derivative of instantaneous eigenstates

        for s in range(n):
            t = s * dt
            A0 = A_fn(t, tau, e, t1c)

            # Standard RK4
            k1 = A0 @ R
            A_h = A_fn(t+0.5*dt, tau, e, t1c)
            k2 = A_h @ (R+0.5*dt*k1)
            k3 = A_h @ (R+0.5*dt*k2)
            A_f = A_fn(t+dt, tau, e, t1c)
            k4 = A_f @ (R+dt*k3)
            R += dt/6*(k1+2*k2+2*k3+k4)

            # Compute CD correction:
            # H_inst = i * A (Hermitian in the Majorana rep)
            H_inst = 1j * A0  # Make Hermitian
            evals, evecs = np.linalg.eigh(H_inst)

            # ∂_t H ≈ (H(t+dt) - H(t-dt))/(2*dt)  (skip at boundaries)
            if 0 < s < n-1:
                A_next = A_fn(t+dt, tau, e, t1c)
                A_prev = A_fn(t-dt, tau, e, t1c)
                dH_dt = 1j * (A_next - A_prev) / (2*dt)

                # CD term: H_cd = i Σ_{n≠m} |n⟩⟨n|∂_t H|m⟩⟨m| / (E_m - E_n)
                H_cd = np.zeros((5,5), dtype=complex)
                for ni in range(5):
                    for mi in range(5):
                        if ni != mi and abs(evals[ni] - evals[mi]) > 1e-10:
                            # ⟨n| ∂_t H |m⟩
                            mat_elem = evecs[:,ni].conj() @ dH_dt @ evecs[:,mi]
                            H_cd += 1j * mat_elem / (evals[mi] - evals[ni]) * \
                                    np.outer(evecs[:,ni], evecs[:,mi].conj())

                # Add CD to the generator
                A_cd = -1j * H_cd  # back to antisymmetric
                A_total = A0 + np.real(A_cd)  # keep real

                # Re-do the step with CD correction
                k1 = A_total @ R
                k2 = (A_fn(t+0.5*dt, tau, e, t1c) + np.real(-1j*H_cd*0.5)) @ (R+0.5*dt*k1)
                R = R + dt * k1  # Euler for CD step (simpler)

    return R

# Test CD on a few cases
print("\n  Testing CD (simplified)...")
for e, t1c, label in [
    (0.0, 0.0, "pure MZM"),
    (0.01, 0.01, "E1=t1=0.01"),
    (0.02, 0.02, "E1=t1=0.02"),
]:
    R_cd = so5_cd(tau, e, t1c)
    R_std = so5_multipulse(tau, e, t1c, n_pulses=0)
    fid_cd = fid_so5(R_cd)
    fid_std = fid_so5(R_std)
    print(f"    {label}: std={fid_std:.4f}, CD={fid_cd:.4f}")

# ═══════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("Summary")
print(f"{'='*60}")
print("  Multi-pulse DD: more π_x pulses → better E1+t1 cancellation")
print("  Limitation: pulses disrupt ancilla dynamics at large N")
print("  CD driving: theoretically exact, numerically challenging")
print(f"{'='*60}")
print("\nDone.")
