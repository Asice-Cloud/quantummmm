import numpy as np

"""Small-system numerical test: Berry holonomy vs. LQC+permutation circuit.

This script implements the toy experiment sketched in the notes:

1. Take the 4-Majorana model and geometric loop H(phi) used in
   verify/run_dehn_twist_like_path.py, and compute the Berry holonomy
   U_Berry acting on the 2D ground-state subspace.
2. View the 4-Majorana system as two physical qubits (4D Hilbert space),
   and consider simple constant-depth circuits of the form

      U_full = (U1a ⊗ U2a) · SWAP · (U1b ⊗ U2b),

   where each Ujk is a single-qubit SU(2) rotation and SWAP is the
   two-qubit swap gate. This has depth O(1) and matches the
   "LQC+permutation" structure at the smallest nontrivial size.
3. Project U_full to the ground-state subspace using the same basis V0
   that defines U_Berry, obtaining a 2x2 logical matrix U_logical.
4. Compare U_logical with U_Berry up to an overall phase via the
   fidelity-like quantity

      F = |Tr(U_Berry^† U_logical)| / 2,

   where 2 is the Hilbert-space dimension. F=1 means they agree up to a
   global phase; 1-F can be viewed as an error measure.

We perform a simple random search over Euler angles for U1a, U2a, U1b,
U2b to find a low-error constant-depth realization and report the best
F and the corresponding invariants (eigenvalues, trace) for comparison.
"""

import os
import sys


# Allow importing run_dehn_twist_like_path when this script is executed
# from the repository root. We append the verify/ directory to sys.path
# explicitly instead of relying on it being a package.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if THIS_DIR not in sys.path:
    sys.path.append(THIS_DIR)

import run_dehn_twist_like_path as dehn


def su2_from_euler(alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Return an SU(2) matrix using Z-Y-Z Euler angles.

    U = exp(-i alpha σ_z / 2) exp(-i beta σ_y / 2) exp(-i gamma σ_z / 2).
    """
    cz1, sz1 = np.cos(alpha / 2.0), np.sin(alpha / 2.0)
    cy, sy = np.cos(beta / 2.0), np.sin(beta / 2.0)
    cz2, sz2 = np.cos(gamma / 2.0), np.sin(gamma / 2.0)

    # exp(-i alpha σ_z / 2)
    U1 = np.array([[np.exp(-1j * alpha / 2.0), 0.0],
                   [0.0, np.exp(1j * alpha / 2.0)]], dtype=complex)
    # exp(-i beta σ_y / 2)
    U2 = np.array([[cy, -sy],
                   [sy, cy]], dtype=complex)
    # exp(-i gamma σ_z / 2)
    U3 = np.array([[np.exp(-1j * gamma / 2.0), 0.0],
                   [0.0, np.exp(1j * gamma / 2.0)]], dtype=complex)

    return U1 @ U2 @ U3


def random_su2(rng: np.random.Generator) -> np.ndarray:
    """Sample a (roughly) uniform random SU(2) using random Euler angles."""
    alpha = rng.uniform(0.0, 2.0 * np.pi)
    beta = rng.uniform(0.0, np.pi)
    gamma = rng.uniform(0.0, 2.0 * np.pi)
    return su2_from_euler(alpha, beta, gamma)


def best_lqc_swap_fit(num_samples: int = 5000, seed: int | None = 0):
    # 1. Berry holonomy on the ground subspace for the geometric loop H(phi).
    U_berry = dehn.compute_holonomy(num_steps=400)

    # 2. Ground-space basis at phi = 0: columns of V0 span the logical subspace.
    _, V0 = dehn.ground_subspace(dehn.build_H_phi(0.0))

    # 3. Two-qubit SWAP gate acting on the 4D Hilbert space.
    SWAP = np.array([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
    ], dtype=complex)

    rng = np.random.default_rng(seed)

    best_F = -1.0
    best_params: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None = None
    best_U_logical: np.ndarray | None = None

    for _ in range(num_samples):
        # Sample four single-qubit SU(2) gates.
        U1a = random_su2(rng)
        U2a = random_su2(rng)
        U1b = random_su2(rng)
        U2b = random_su2(rng)

        U_a = np.kron(U1a, U2a)
        U_b = np.kron(U1b, U2b)

        # Constant-depth LQC+permutation circuit.
        U_full = U_a @ SWAP @ U_b

        # Logical action on the 2D ground-space subspace.
        U_logical = V0.conj().T @ U_full @ V0

        # Remove overall phase: fidelity-like overlap.
        overlap = np.trace(U_berry.conj().T @ U_logical) / 2.0
        F = np.abs(overlap)

        if F > best_F:
            best_F = F
            best_params = (U1a, U2a, U1b, U2b)
            best_U_logical = U_logical

    return U_berry, best_F, best_params, best_U_logical


def main():
    U_berry, best_F, best_params, best_U_logical = best_lqc_swap_fit()

    print("U_Berry (2x2) =")
    print(np.round(U_berry, 6))
    print("Eigenvalues U_Berry:", np.round(np.sort_complex(np.linalg.eigvals(U_berry)), 6))
    print("Trace U_Berry:", np.round(np.trace(U_berry), 6))

    if best_U_logical is None or best_params is None:
        print("No candidate circuit found (this should not happen).")
        return

    print("\nBest LQC+SWAP logical U (2x2) =")
    print(np.round(best_U_logical, 6))
    print("Eigenvalues U_LQC:", np.round(np.sort_complex(np.linalg.eigvals(best_U_logical)), 6))
    print("Trace U_LQC:", np.round(np.trace(best_U_logical), 6))

    print("\nBest fidelity F = |Tr(U_Berry^† U_LQC)|/2 =", best_F)
    print("1 - F (error measure) =", 1.0 - best_F)


if __name__ == "__main__":
    main()
