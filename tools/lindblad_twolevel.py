#!/usr/bin/env python3
"""Simple Lindblad / Kraus-style decoherence utilities for a two-level (Pauli) simulator.

Provides a time-stepping routine that alternates unitary evolution and simple
Kraus channels (amplitude damping + dephasing) to mimic T1/T2 effects.
"""
import numpy as np
from scipy.linalg import expm


sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
identity2 = np.eye(2, dtype=complex)


def kraus_dephasing(rho, p):
    """Apply dephasing channel with probability p using Kraus ops K0=sqrt(1-p) I, K1=sqrt(p) σ_z."""
    if p <= 0:
        return rho
    k0 = np.sqrt(1.0 - p) * identity2
    k1 = np.sqrt(p) * sigma_z
    return k0 @ rho @ k0.conj().T + k1 @ rho @ k1.conj().T


def kraus_amplitude_damping(rho, p):
    """Amplitude damping Kraus: relax |1>->|0> with prob p."""
    if p <= 0:
        return rho
    k0 = np.array([[1.0, 0.0], [0.0, np.sqrt(1.0 - p)]], dtype=complex)
    k1 = np.array([[0.0, np.sqrt(p)], [0.0, 0.0]], dtype=complex)
    return k0 @ rho @ k0.conj().T + k1 @ rho @ k1.conj().T


def bloch_from_rho(rho):
    bx = np.real(np.trace(rho @ sigma_x))
    by = np.real(np.trace(rho @ sigma_y))
    bz = np.real(np.trace(rho @ sigma_z))
    return np.array([bx, by, bz], dtype=float)


def run_lindblad_time_series(H_list, dt, rho0=None, gamma_deph=0.0, gamma_relax=0.0):
    """Evolve a density matrix list under unitary steps followed by Kraus channels.

    H_list : list/array of 2x2 Hamiltonians (one per timestep)
    dt     : timestep (scalar)
    rho0   : initial density matrix (defaults to ground state of H_list[0])
    gamma_deph, gamma_relax : rates (1/time) for dephasing and relaxation
    """
    nsteps = len(H_list)
    # initialize rho0
    if rho0 is None:
        # ground state projector of first Hamiltonian
        evals, evecs = np.linalg.eigh(H_list[0])
        psi0 = evecs[:, 0]
        rho = np.outer(psi0, psi0.conj())
    else:
        rho = np.array(rho0, dtype=complex)

    bloch = np.zeros((nsteps, 3), dtype=float)
    for i, H in enumerate(H_list):
        # unitary step
        U = expm(-1j * H * dt)
        rho = U @ rho @ U.conj().T

        # Kraus steps: convert rates -> per-step probabilities
        p_deph = 1.0 - np.exp(-gamma_deph * dt)
        p_relax = 1.0 - np.exp(-gamma_relax * dt)

        # apply amplitude damping then dephasing
        rho = kraus_amplitude_damping(rho, p_relax)
        rho = kraus_dephasing(rho, p_deph)

        bloch[i, :] = bloch_from_rho(rho)

    return bloch


def demo_static(H, total_T=10.0, nsteps=1000, gamma_deph=0.0, gamma_relax=0.0):
    dt = total_T / nsteps
    H_list = [H] * nsteps
    bloch = run_lindblad_time_series(H_list, dt, gamma_deph=gamma_deph, gamma_relax=gamma_relax)
    return np.linspace(0, total_T, nsteps), bloch


if __name__ == '__main__':
    # small self-test: precessing σ_x Hamiltonian with dephasing
    H = 2.0 * np.pi * 0.1 * sigma_x
    t, bloch = demo_static(H, total_T=10.0, nsteps=400, gamma_deph=0.2, gamma_relax=0.05)
    print('Demo run, final Bloch:', bloch[-1])
