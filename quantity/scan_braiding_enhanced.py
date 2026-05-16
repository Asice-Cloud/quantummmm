"""
Enhanced braiding fidelity scanner with multiple initial states and E1=0 baseline.
Tests both computational basis and adiabatic ground state initialization.
"""
import numpy as np
from scipy.linalg import eigh, expm
import matplotlib.pyplot as plt
from importlib import util

# Load the braiding module
spec = util.spec_from_file_location('mod', 'quantity/reproduce_effective_braiding_majorana.py')
mod = util.module_from_spec(spec)
spec.loader.exec_module(mod)


def run_braiding_adiabatic(params, tau=10.0, steps_per_tau=100, repeat=2, tc=0.5):
    """
    Modified version: Initialize in adiabatic ground state of H(t=0).
    Measure overlap of final state with initial ground state.
    This is closer to the paper's experimental preparation.
    """
    n_modes = 3
    c_ops = mod.fermionic_operators(n_modes)
    gammas = []
    for c in c_ops:
        cd = c.conj().T
        gammas.append(c + cd)
        gammas.append(-1j * (c - cd))
    dim = 2 ** n_modes

    # Initialize in ground state of H(t=0)
    t1_0, t2_0, t3_0, Ed0 = mod.time_profiles(0.0, tau, tc)
    H0 = mod.build_H(gammas, params, t1_0, t2_0, t3_0, Ed0)
    evals, evecs = eigh(H0)
    psi_initial = evecs[:, 0]  # ground state at t=0
    psi = psi_initial.copy()

    total_time = repeat * 3 * tau
    steps = int(repeat * 3 * tau * steps_per_tau)
    dt = total_time / steps

    # Time evolution
    t = 0.0
    for k in range(steps):
        t += dt
        t1v, t2v, t3v, Edv = mod.time_profiles(t, tau, tc)
        Ht = mod.build_H(gammas, params, t1v, t2v, t3v, Edv)
        U = expm(-1j * Ht * dt)
        psi = U @ psi

    # Measure fidelity as overlap with initial ground state
    fidelity = abs(np.vdot(psi_initial, psi))
    return fidelity


def run_parallel_braiding(params, tau=10.0, steps_per_tau=100, repeat=2, tc=0.5):
    """
    Run both computational basis and adiabatic initializations, return both.
    """
    f_comp, _, _, _ = mod.run_braiding(params, tau, steps_per_tau, repeat, tc)
    f_adiab = run_braiding_adiabatic(params, tau, steps_per_tau, repeat, tc)
    return f_comp, f_adiab


# Main scan
print("=" * 80)
print("Enhanced Braiding Fidelity Analysis")
print("=" * 80)

taus = np.linspace(5, 150, 20)
E1_vals = [0.0, 1e-4, 1e-3, 1e-2]  # Including E1=0 baseline
tc = 0.03
steps_per_tau = 80
repeat = 2

results = {
    'comp': {E1: [] for E1 in E1_vals},
    'adiab': {E1: [] for E1 in E1_vals}
}

print(f"\nRunning scans: τ ∈ [5, 150] ({len(taus)} points)")
print(f"              E1 ∈ {E1_vals}")
print(f"              tc = {tc}, steps_per_tau = {steps_per_tau}, repeat = {repeat}")
print()

for E1 in E1_vals:
    print(f"Processing E1 = {E1:.0e}")
    for i, tau in enumerate(taus):
        f_c, f_a = run_parallel_braiding({'E1': E1}, tau=float(tau), 
                                         steps_per_tau=steps_per_tau, 
                                         repeat=repeat, tc=tc)
        results['comp'][E1].append(f_c)
        results['adiab'][E1].append(f_a)
        status = f"  τ={tau:6.1f}: comp={f_c:.4f}, adiab={f_a:.4f}"
        if i % 5 == 0:
            print(status)

# Save results
print("\nSaving results...")
for E1 in E1_vals:
    np.savetxt(f'quantity/braiding_enhanced_E1_{E1:.0e}_comp.txt',
               np.column_stack((taus, results['comp'][E1])))
    np.savetxt(f'quantity/braiding_enhanced_E1_{E1:.0e}_adiab.txt',
               np.column_stack((taus, results['adiab'][E1])))

# Plotting
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Computational basis initialization
ax = axes[0]
for E1 in E1_vals:
    label = f'E1={E1:.0e}' if E1 > 0 else 'E1=0 (ideal MZM)'
    ax.plot(taus, results['comp'][E1], '-o', label=label, markersize=4)
ax.set_xlabel('Braiding time τ')
ax.set_ylabel('Fidelity')
ax.set_title('Computational Basis Initialization |0⟩')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_ylim([0, 1.05])

# Plot 2: Adiabatic ground state initialization
ax = axes[1]
for E1 in E1_vals:
    label = f'E1={E1:.0e}' if E1 > 0 else 'E1=0 (ideal MZM)'
    ax.plot(taus, results['adiab'][E1], '-o', label=label, markersize=4)
ax.set_xlabel('Braiding time τ')
ax.set_ylabel('Fidelity')
ax.set_title('Adiabatic Ground State Initialization')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_ylim([0, 1.05])

plt.tight_layout()
plt.savefig('quantity/braiding_enhanced_comparison.png', dpi=200)
print("Saved: quantity/braiding_enhanced_comparison.png")

# Summary statistics
print("\n" + "=" * 80)
print("Summary Statistics")
print("=" * 80)
for E1 in E1_vals:
    f_c = np.array(results['comp'][E1])
    f_a = np.array(results['adiab'][E1])
    print(f"\nE1 = {E1:.0e}:")
    print(f"  Computational basis: mean={f_c.mean():.4f}, std={f_c.std():.4f}, "
          f"min={f_c.min():.4f}, max={f_c.max():.4f}")
    print(f"  Adiabatic state:     mean={f_a.mean():.4f}, std={f_a.std():.4f}, "
          f"min={f_a.min():.4f}, max={f_a.max():.4f}")

print("\n" + "=" * 80)
print("Interpretation:")
print("=" * 80)
print("""
Computational basis (|0⟩):
  - Shows oscillatory fidelity evolution
  - E1=0 (ideal): Baseline for pure phase accumulation
  - E1>0 (with ABS): How real impurity perturbs ideal braiding

Adiabatic ground state:
  - Tests preparation robustness
  - Should preserve fidelity better if initialized properly
  - Measures decoherence from time-dependent driving
""")
