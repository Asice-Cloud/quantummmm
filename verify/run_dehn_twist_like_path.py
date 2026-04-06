import numpy as np

"""Toy "Dehn-twist-like" geometric path experiment in a 4-Majorana model.

We work with the same 4x4 Majorana representation as in
run_instantaneous_braid_crosscheck.py and run_interacting_braid_check.py.

Goal:
  - Construct a closed loop H(φ) in the space of quadratic Majorana
    Hamiltonians, with φ ∈ [0, 2π], such that the spectrum is gapped and
    there is a two-dimensional ground-state subspace for all φ.
  - Discretize this loop and compute the non-Abelian Berry holonomy of the
    ground-state subspace around the loop, using parallel transport based
    on overlaps between neighboring eigenbases.
  - Interpret the resulting 2x2 holonomy matrix as a simple, geometry-based
    "Dehn-twist-like" operation on the encoded qubit.

This is not a literal spatial Dehn twist on a genus>0 surface, but a
closed geometric path in coupling space that mimics the idea of taking
zero modes around a nontrivial loop and reading off the holonomy.
"""


def gamma_matrices():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)

    g1 = np.kron(sx, I2)
    g2 = np.kron(sy, I2)
    g3 = np.kron(sz, sx)
    g4 = np.kron(sz, sy)
    return [g1, g2, g3, g4]


def build_H_phi(phi: float):
    """Quadratic Hamiltonian H(φ) defining a closed loop in coupling space.

    We take a unit-amplitude combination of two commuting Majorana bilinears:
      H(φ) = (i/2)[ cos φ * γ₂γ₃ + sin φ * γ₁γ₄ ].

    For φ=0 this reduces to the braid generator Hamiltonian (i/2)γ₂γ₃.
    For φ=π/2 it is (i/2)γ₁γ₄. As φ goes from 0 to 2π, (cos φ, sin φ)
    winds once around the unit circle and H(φ) returns to itself.

    At all φ the spectrum consists of two degenerate negative and two
    degenerate positive eigenvalues, so the ground space is always
    two-dimensional and gapped from excitations, suitable for a simple
    non-Abelian Berry holonomy computation.
    """
    g1, g2, g3, g4 = gamma_matrices()
    H = 0.5j * (np.cos(phi) * (g2 @ g3) + np.sin(phi) * (g1 @ g4))
    return H


def build_H_handle(theta1: float, theta2: float):
    """"Handle-like" path Hamiltonian H(θ1, θ2) on a coupling torus.

    We view the four Majoranas γ₁,…,γ₄ as sitting on a ring with
    couplings on the four links (1-2), (2-3), (3-4), (4-1). We then
    build two independent "cycles" of couplings:

      - A "meridian" cycle on links (2-3) and (4-1):
        H_M(θ1) = (i/2)[ cos θ1 * γ₂γ₃ + sin θ1 * γ₄γ₁ ],

      - A "longitude" cycle on links (1-2) and (3-4):
        H_L(θ2) = (i/2)[ cos θ2 * γ₁γ₂ + sin θ2 * γ₃γ₄ ].

    Each pair in H_M and H_L uses bilinears that commute within the
    pair and correspond to moving a "strong bond" around the ring along
    two independent directions. The full Hamiltonian is their sum

      H(θ1, θ2) = H_M(θ1) + H_L(θ2).

    The parameter space (θ1, θ2) is a 2-torus. A path that winds
    nontrivially in both θ1 and θ2 directions mimics going once around a
    noncontractible cycle on a genus-1 handle, and the Berry holonomy in
    the ground subspace along such a path provides a simple
    "Dehn-twist-like" geometric operation.
    """
    g1, g2, g3, g4 = gamma_matrices()

    # Basic bilinears (i/2) γ_a γ_b on each link of the ring
    B12 = 0.5j * (g1 @ g2)
    B23 = 0.5j * (g2 @ g3)
    B34 = 0.5j * (g3 @ g4)
    B41 = 0.5j * (g4 @ g1)

    H_M = np.cos(theta1) * B23 + np.sin(theta1) * B41
    H_L = np.cos(theta2) * B12 + np.sin(theta2) * B34

    # Optionally one could introduce different amplitudes for H_M, H_L
    # to optimize the gap; here we keep them equal and rely on the small
    # size of the model to keep a well-defined two-dimensional ground
    # subspace, which we always select numerically via ground_subspace.
    H = H_M + H_L
    return H


def ground_subspace(H: np.ndarray, dim_gs: int = 2):
    """Return eigenvalues and eigenvectors of the lowest dim_gs states.

    H is 4x4 Hermitian. We sort eigenvalues in ascending order and
    return the first dim_gs eigenvalues and the corresponding eigenvectors
    as columns of a 4 x dim_gs matrix.
    """
    w, v = np.linalg.eigh(H)
    idx = np.argsort(w)
    w_gs = w[idx[:dim_gs]]
    V_gs = v[:, idx[:dim_gs]]
    return w_gs, V_gs


def project_to_unitary(W: np.ndarray) -> np.ndarray:
    """Given an overlap matrix W, project it to the closest unitary via
    polar decomposition using SVD: W = U S V^† -> U_polar = U V^†.
    """
    U, s, Vh = np.linalg.svd(W)
    return U @ Vh


def compute_holonomy(num_steps: int = 200):
    """Discretize φ ∈ [0, 2π] and compute Berry holonomy on ground subspace.

    We use parallel transport: for neighboring φ_k, φ_{k+1}, with ground
    subspace bases V_k and V_{k+1} (4 x 2), define overlap
      W_k = V_k^† V_{k+1} (2 x 2),
    project W_k to the closest unitary U_k, and accumulate
      U_holo = ∏_k U_k.

    The resulting 2x2 U_holo approximates the non-Abelian Berry holonomy
    around the closed loop in coupling space.
    """
    phis = np.linspace(0.0, 2.0 * np.pi, num_steps + 1)

    # Ground-space basis at the first point
    _, V_prev = ground_subspace(build_H_phi(phis[0]))

    U_holo = np.eye(2, dtype=complex)

    for k in range(num_steps):
        phi_next = phis[k + 1]
        _, V_next = ground_subspace(build_H_phi(phi_next))

        # Overlap between consecutive ground-space bases
        W = V_prev.conj().T @ V_next  # 2x2
        U_step = project_to_unitary(W)
        U_holo = U_step @ U_holo

        V_prev = V_next

    return U_holo


def compute_holonomy_handle(num_steps: int = 400):
    """Compute Berry holonomy along a "handle-like" loop on (θ1, θ2).

    We take a closed path on the parameter torus (θ1, θ2) defined by a
    single angular parameter λ ∈ [0, 2π]:

      θ1(λ) = λ,   θ2(λ) = λ.

    As λ goes from 0 to 2π, the point (θ1, θ2) winds once along the
    diagonal of the 2-torus, which is homotopically nontrivial and
    represents a loop that threads both fundamental cycles of the torus
    simultaneously. In the language of the toy ring, this corresponds to
    dragging the pattern of strong bonds once around the ring in a way
    that intertwines the "meridian" and "longitude" cycles, analogously
    to going once around a handle.

    The Berry holonomy of the two-dimensional ground subspace along this
    loop is computed with the same parallel-transport procedure used for
    compute_holonomy.
    """
    lambdas = np.linspace(0.0, 2.0 * np.pi, num_steps + 1)

    # Ground-space basis at the first point
    theta1_0 = lambdas[0]
    theta2_0 = lambdas[0]
    _, V_prev = ground_subspace(build_H_handle(theta1_0, theta2_0))

    U_holo = np.eye(2, dtype=complex)

    for k in range(num_steps):
      lam_next = lambdas[k + 1]
      theta1_next = lam_next
      theta2_next = lam_next

      _, V_next = ground_subspace(build_H_handle(theta1_next, theta2_next))

      W = V_prev.conj().T @ V_next  # 2x2 overlap
      U_step = project_to_unitary(W)
      U_holo = U_step @ U_holo

      V_prev = V_next

    return U_holo


def main():
  # Original single-parameter loop in coupling space
  U_holo = compute_holonomy(num_steps=400)
  evals, _ = np.linalg.eig(U_holo)

  print("=== Single-parameter loop H(phi) ===")
  print("Holonomy U_holo (2x2) =")
  print(np.round(U_holo, 6))
  print("Eigenvalues of U_holo:")
  print(np.round(np.sort_complex(evals), 6))

  # New handle-like loop on the (θ1, θ2) coupling torus
  U_handle = compute_holonomy_handle(num_steps=400)
  evals_handle, _ = np.linalg.eig(U_handle)

  print("\n=== Handle-like loop H(theta1, theta2) ===")
  print("Holonomy U_handle (2x2) =")
  print(np.round(U_handle, 6))
  print("Eigenvalues of U_handle:")
  print(np.round(np.sort_complex(np.sort_complex(evals_handle)), 6))


if __name__ == "__main__":
    main()
