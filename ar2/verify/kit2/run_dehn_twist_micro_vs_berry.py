import numpy as np

"""Compare geometric Berry holonomy with a microscopic Dehn-twist candidate
built from half-twist gates in the 4-Majorana toy model.

We interpret:
  - Half twist: U_half = exp((pi/4) * gamma_i gamma_j), here for (i,j) = (2,3).
  - Dehn twist: U_Dehn ~ F^{-1} R^2 F, which in this toy becomes
        U_micro = F^{-1} U_half^2 F
    where F is the change of basis from the microscopic Majorana basis to
    the instantaneous ground-subspace basis used in the Berry holonomy
    computation.

We then compare U_micro to the Berry holonomy U_Berry obtained from the
closed loop H(phi) in coupling space defined in run_dehn_twist_like_path.
"""

from run_dehn_twist_like_path import (
    gamma_matrices,
    build_H_phi,
    ground_subspace,
    compute_holonomy,
)


def half_twist_matrix(i: int, j: int) -> np.ndarray:
    """Return the 4x4 half-twist unitary exp((pi/4) * gamma_i gamma_j).

    For Majorana operators gamma_i with {gamma_i, gamma_j} = 2 delta_ij and
    gamma_i^2 = 1, the bilinear K = gamma_i gamma_j is anti-Hermitian and
    squares to -I. Thus exp(theta K) can be written in closed form as

        exp(theta K) = cos(theta) * I + sin(theta) * K.
    """
    gammas = gamma_matrices()
    Gi = gammas[i - 1]
    Gj = gammas[j - 1]
    K = Gi @ Gj  # 4x4

    # Sanity check: K^2 should be -I (up to numerical precision).
    K2 = K @ K
    I4 = np.eye(4, dtype=complex)
    err = np.linalg.norm(K2 + I4)
    if err > 1e-10:
        print(f"Warning: (gamma_{i} gamma_{j})^2 deviates from -I, ||K^2 + I|| = {err:.3e}")

    theta = np.pi / 4.0
    U_half = np.cos(theta) * I4 + np.sin(theta) * K
    return U_half


def project_to_ground(U: np.ndarray, H_ref: np.ndarray) -> np.ndarray:
    """Project a 4x4 unitary U onto the 2D ground subspace of H_ref.

    This implements U_micro = F^{-1} U F with F given by the ground-subspace
    eigenvectors V_gs of H_ref, i.e.

        U_micro = V_gs^\dagger U V_gs  (2x2).
    """
    _, V_gs = ground_subspace(H_ref)
    return V_gs.conj().T @ U @ V_gs


def su2_fidelity(U_target: np.ndarray, U_candidate: np.ndarray):
    """Compute phase-agnostic overlap between two 2x2 unitaries.

    We use F = |Tr(U_target^\dagger U_candidate)| / 2, which equals 1 when the
    two unitaries agree up to an overall U(1) phase.
    """
    M = U_target.conj().T @ U_candidate
    tr = np.trace(M)
    F = np.abs(tr) / 2.0
    return F, tr


def normalize_to_su2(U: np.ndarray) -> np.ndarray:
    """Strip global phase so that det(U) = 1.

    This identifies the corresponding element in SU(2) up to an overall
    U(1) phase, which is the natural equivalence class in TQFT.
    """
    det = np.linalg.det(U)
    # np.sqrt chooses a branch of the square root; either branch differs
    # only by an overall phase, which is irrelevant for our comparison.
    sqrt_det = np.sqrt(det)
    return U / sqrt_det


def dehn_twist_ising_matrix() -> np.ndarray:
    """Abstract Dehn twist for two Ising anyons σ×σ=1+ψ in fusion basis.

    We take R-symbols
        R^{σσ}_1 = exp(-iπ/8),   R^{σσ}_ψ = exp(3iπ/8),
    so that (R^{σσ})^2 has eigenvalues
        e^{-iπ/4}, e^{3iπ/4}.

    Up to an overall phase this is diag(1, -1), i.e. a π-rotation about
    some axis; after det-normalization it lies in SU(2) and is conjugate
    to the SU(2) representation of a Dehn twist.
    """
    lam_1 = np.exp(-0.25j * np.pi)
    lam_psi = np.exp(0.75j * np.pi)
    U_R2 = np.diag([lam_1, lam_psi])  # (R^{σσ})^2 in fusion basis
    return normalize_to_su2(U_R2)


def main():
    # 1) Berry holonomy for the geometric loop H(phi) from run_dehn_twist_like_path
    U_Berry = compute_holonomy(num_steps=400)
    evals_Berry, _ = np.linalg.eig(U_Berry)

    print("=== Geometric Berry holonomy from H(phi) loop ===")
    print("U_Berry (2x2) =")
    print(np.round(U_Berry, 6))
    print("Eigenvalues of U_Berry:")
    print(np.round(np.sort_complex(evals_Berry), 6))

    # Normalize Berry holonomy to SU(2) by stripping global phase
    U_Berry_su2 = normalize_to_su2(U_Berry)
    print("\nU_Berry normalized to SU(2) (det=1):")
    print(np.round(U_Berry_su2, 6))

    # 1b) Abstract Ising Dehn twist from R-symbols
    U_Dehn_ising = dehn_twist_ising_matrix()
    evals_Dehn, _ = np.linalg.eig(U_Dehn_ising)
    print("\n=== Abstract Ising Dehn twist (two σ anyons) ===")
    print("U_Dehn_ising (2x2, SU(2)-normalized) =")
    print(np.round(U_Dehn_ising, 6))
    print("Eigenvalues of U_Dehn_ising:")
    print(np.round(np.sort_complex(evals_Dehn), 6))

    # Check SU(2) conjugacy: find V with V^† U_Berry_su2 V ≈ U_Dehn_ising
    wB, vB = np.linalg.eig(U_Berry_su2)
    wD, vD = np.linalg.eig(U_Dehn_ising)

    # Order eigenvectors consistently by sorting eigenvalues
    idxB = np.argsort(np.angle(wB))
    idxD = np.argsort(np.angle(wD))
    vB_ord = vB[:, idxB]
    vD_ord = vD[:, idxD]

    V = vB_ord @ np.linalg.inv(vD_ord)
    conj_diff = V.conj().T @ U_Berry_su2 @ V - U_Dehn_ising
    print("\n=== SU(2) conjugacy check ===")
    print("|| V^† U_Berry_su2 V - U_Dehn_ising ||_F =", np.linalg.norm(conj_diff))

    # 2) Reference ground subspace at phi=0
    H0 = build_H_phi(0.0)

    # 3) Scan all distinct Majorana pairs (i, j) and compare
    pairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    print("\n=== Scan over half-twist pairs (i, j) ===")
    for (i, j) in pairs:
        U_half_ij = half_twist_matrix(i, j)
        U_half2_ij = U_half_ij @ U_half_ij

        U_half_proj = project_to_ground(U_half_ij, H0)
        U_half2_proj = project_to_ground(U_half2_ij, H0)

        F_half, tr_half = su2_fidelity(U_Berry, U_half_proj)
        F_half2, tr_half2 = su2_fidelity(U_Berry, U_half2_proj)

        print(f"\n  Pair (i, j) = ({i}, {j})  ")
        print("U_half_proj (2x2) =")
        print(np.round(U_half_proj, 6))
        print("U_half2_proj (2x2) =")
        print(np.round(U_half2_proj, 6))
        print("Tr(U_Berry^dagger U_half_proj) =", tr_half)
        print("F_half  =", F_half)
        print("Tr(U_Berry^dagger U_half2_proj) =", tr_half2)
        print("F_half2 =", F_half2)


if __name__ == "__main__":
    main()
