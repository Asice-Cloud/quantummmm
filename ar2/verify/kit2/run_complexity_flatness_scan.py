import numpy as np
import os
import sys
import matplotlib.pyplot as plt

"""Complexity / flatness layer test (experiment 3, toy version).

We quantify how hard it is for a fixed shallow LQC+permutation ansatz to
reproduce the Berry holonomy when the geometric path is gradually
"curved" away from the exactly-solvable, highly symmetric case.

Setup (building on run_dehn_twist_like_path.py and run_lqc_permutation_fit.py):

- Start from the 4-Majorana toy model and the geometric loop H(phi)
  used to define U_Berry in run_dehn_twist_like_path.py.
- Introduce a perturbation parameter eps ≥ 0 that adds a non-commuting
  Majorana bilinear to H(phi), mimicking an increase in Berry/parameter
  curvature and deviation from the YBE-like flat structure.
- For each eps in a small grid, we
    1) compute the Berry holonomy U_Berry(eps) for the perturbed path;
    2) take the 2D ground-space basis at phi=0 as logical subspace;
    3) search, within the same shallow ansatz

         U_full = (U1a ⊗ U2a) · SWAP · (U1b ⊗ U2b),

       for the best-fit logical action U_LQC(eps) on the ground
       subspace, measured by

         F(eps) = |Tr(U_Berry(eps)^† U_LQC(eps))| / 2.

The function F(eps) thus plays the role of an operational
"flatness/complexity" indicator: at fixed ansatz depth, if curvature is
small (eps ≈ 0), we expect F(eps) to be close to 1; as eps grows and the
path becomes more distorted, the best achievable F(eps) with the same
ansatz should decrease, signaling that more LQC layers / corrections are
needed to realize the same Berry holonomy.
"""

# Allow importing sibling modules when run from the repo root.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if THIS_DIR not in sys.path:
    sys.path.append(THIS_DIR)

import run_dehn_twist_like_path as dehn  # type: ignore


def su2_from_euler(alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Return an SU(2) matrix using Z-Y-Z Euler angles.

    U = exp(-i alpha σ_z / 2) exp(-i beta σ_y / 2) exp(-i gamma σ_z / 2).
    """
    # exp(-i alpha σ_z / 2)
    U1 = np.array([[np.exp(-1j * alpha / 2.0), 0.0],
                   [0.0, np.exp(1j * alpha / 2.0)]], dtype=complex)
    # exp(-i beta σ_y / 2)
    cy, sy = np.cos(beta / 2.0), np.sin(beta / 2.0)
    U2 = np.array([[cy, -sy],
                   [sy, cy]], dtype=complex)
    # exp(-i gamma σ_z / 2)
    U3 = np.array([[np.exp(-1j * gamma / 2.0), 0.0],
                   [0.0, np.exp(1j * gamma / 2.0)]], dtype=complex)

    return U1 @ U2 @ U3


def random_su2(rng: np.random.Generator) -> np.ndarray:
    """Sample a (roughly) random SU(2) via random Euler angles."""
    alpha = rng.uniform(0.0, 2.0 * np.pi)
    beta = rng.uniform(0.0, np.pi)
    gamma = rng.uniform(0.0, 2.0 * np.pi)
    return su2_from_euler(alpha, beta, gamma)


def build_H_phi_eps(phi: float, eps: float) -> np.ndarray:
    """Perturbed geometric path H(φ; eps) in the 4-Majorana model.

    We start from the original loop

        H0(φ) = (i/2)[ cos φ * γ₂γ₃ + sin φ * γ₁γ₄ ],

    and add a small non-commuting bilinear term on link (1-2):

        H_pert = eps * (i/2) γ₁γ₂.

    For eps = 0 we recover the original, more symmetric loop used in
    run_dehn_twist_like_path.py; as eps grows, the spectrum and ground
    subspace twist more nontrivially along the path, mimicking increased
    Berry/parameter curvature away from the "flat" case.
    """
    g1, g2, g3, g4 = dehn.gamma_matrices()
    H0 = 0.5j * (np.cos(phi) * (g2 @ g3) + np.sin(phi) * (g1 @ g4))
    H_pert = eps * 0.5j * (g1 @ g2)
    return H0 + H_pert


def compute_holonomy_eps(eps: float, num_steps: int = 200) -> np.ndarray:
    """Discrete Berry holonomy for the perturbed path at strength eps."""
    phis = np.linspace(0.0, 2.0 * np.pi, num_steps + 1)

    # Ground-space basis at the first point
    _, V_prev = dehn.ground_subspace(build_H_phi_eps(phis[0], eps))

    U_holo = np.eye(2, dtype=complex)

    for k in range(num_steps):
        phi_next = phis[k + 1]
        _, V_next = dehn.ground_subspace(build_H_phi_eps(phi_next, eps))

        W = V_prev.conj().T @ V_next  # 2x2 overlap
        U_step = dehn.project_to_unitary(W)
        U_holo = U_step @ U_holo

        V_prev = V_next

    return U_holo


def best_lqc_swap_fit_for_eps(eps: float,
                              num_samples: int = 3000,
                              seed: int | None = 0):
    """For given eps, find best shallow LQC+SWAP approximation to U_Berry(eps).

    Returns (U_Berry, best_F, best_U_logical).
    """
    U_berry = compute_holonomy_eps(eps, num_steps=400)

    # Ground-space basis at phi = 0 for this eps.
    _, V0 = dehn.ground_subspace(build_H_phi_eps(0.0, eps))

    SWAP = np.array([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
    ], dtype=complex)

    rng = np.random.default_rng(seed)

    best_F = -1.0
    best_U_logical: np.ndarray | None = None

    for _ in range(num_samples):
        U1a = random_su2(rng)
        U2a = random_su2(rng)
        U1b = random_su2(rng)
        U2b = random_su2(rng)

        U_a = np.kron(U1a, U2a)
        U_b = np.kron(U1b, U2b)

        U_full = U_a @ SWAP @ U_b
        U_logical = V0.conj().T @ U_full @ V0

        overlap = np.trace(U_berry.conj().T @ U_logical) / 2.0
        F = np.abs(overlap)

        if F > best_F:
            best_F = F
            best_U_logical = U_logical

    return U_berry, best_F, best_U_logical


def main():
    eps_list = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

    print("eps  |  F_max   |  1-F_max  |  Tr(U_Berry)        |  Tr(U_LQC) (best)")
    print(" --+   +   -+       +       ")

    F_list: list[float] = []

    for eps in eps_list:
        U_berry, Fmax, U_lqc = best_lqc_swap_fit_for_eps(eps)

        if U_lqc is None:
            print(f"{eps:3.1f} |  (no candidate found)")
            F_list.append(np.nan)
            continue

        tr_berry = np.trace(U_berry)
        tr_lqc = np.trace(U_lqc)

        print(f"{eps:3.1f} | {Fmax:7.4f} | {1.0-Fmax:8.4f} | "
              f"{tr_berry.real:7.4f}+{tr_berry.imag:7.4f}i | "
              f"{tr_lqc.real:7.4f}+{tr_lqc.imag:7.4f}i")

        F_list.append(Fmax)

    # Plot F_max(eps) as a simple curve and save to PNG in this directory.
    fig, ax = plt.subplots(figsize=(4.0, 3.0))
    ax.plot(eps_list, F_list, "o-", lw=1.5)
    ax.set_xlabel("epsilon (perturbation strength)")
    ax.set_ylabel("F_max(epsilon)")
    ax.set_title("Best LQC+SWAP fidelity vs perturbation strength")
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, ls=":", alpha=0.5)
    fig.tight_layout()

    out_path = os.path.join(THIS_DIR, "complexity_flatness_F_vs_eps.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    print(f"\nSaved F_max(eps) curve to {out_path}")


if __name__ == "__main__":
    main()
