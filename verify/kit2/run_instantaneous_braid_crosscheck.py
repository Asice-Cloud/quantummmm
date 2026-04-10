import numpy as np

"""Numerical cross-check between our R(a,b,c,d)->Kitaev mapping
and the "instantaneous braid" generator U = exp((pi/4) gamma_2 gamma_3)
used in Zhu–Lavasani–Barkeshli (arXiv:1806.06078).

We:
1. Build a 4-Majorana Hamiltonian H from (a,b,c,d) via the R_to_Kitaev
   mapping (restricted to a two-site Kitaev chain).
2. Time evolve for tau = pi/2 and compare
      U_time = exp(-i H tau)
   to the ideal braid generator
      U_ideal = exp((pi/4) gamma_2 gamma_3).

For the special choice (a,b,c,d) = (0,1,0,0), the notes show that
  t = 1, Delta = 1, mu = 0
  A_23 = (t+Delta)/2 = 1, A_14 = 0, A_12 = A_34 = 0,
so that
  H = (i/2) * A_23 * gamma_2 gamma_3 = (i/2) gamma_2 gamma_3,
which implies
  exp(-i H tau) = exp((tau A_23 / 2) gamma_2 gamma_3).
Setting tau = pi/2 yields exactly U = exp((pi/4) gamma_2 gamma_3),
matching the braid generator used in the paper.
"""


def gamma_matrices():
    """Return a concrete 4x4 matrix representation of 4 Majorana operators.

    We use a standard Clifford representation for 2 complex fermions:
      gamma1 = sigma_x (first qubit) tensor I (second qubit)
      gamma2 = sigma_y (first qubit) tensor I (second qubit)
      gamma3 = sigma_z (first qubit) tensor sigma_x (second qubit)
      gamma4 = sigma_z (first qubit) tensor sigma_y (second qubit)
    which satisfies {gamma_a, gamma_b} = 2 delta_{ab}.
    """
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)

    g1 = np.kron(sx, I2)
    g2 = np.kron(sy, I2)
    g3 = np.kron(sz, sx)
    g4 = np.kron(sz, sy)
    return [g1, g2, g3, g4]


def build_majorana_H_from_R(a=0.0, b=0.0, c=0.0, d=0.0, mu0=0.0):
    """Build 4x4 Majorana Hamiltonian H from (a,b,c,d) using R_to_Kitaev mapping.

    For a two-site Kitaev chain we use the standard relations (see R_to_Kitaev.md):
      t = b + c
      Delta = b - c
      mu = 4 d + mu0

    The quadratic Majorana Hamiltonian in matrix form is
      H = (i/2) * sum_{a<b} A_ab gamma_a gamma_b,
    where the nonzero A_ab for the two-site Kitaev chain are
      A_23 = (t + Delta) / 2
      A_14 = (Delta - t) / 2
      A_12 = A_34 = -mu / 2.
    """
    g1, g2, g3, g4 = gamma_matrices()

    t = b + c
    Delta = b - c
    mu = 4.0 * d + mu0

    A_23 = (t + Delta) / 2.0
    A_14 = (Delta - t) / 2.0
    A_12 = -mu / 2.0
    A_34 = -mu / 2.0

    H = (1j / 2.0) * (
        A_12 * (g1 @ g2)
        + A_23 * (g2 @ g3)
        + A_34 * (g3 @ g4)
        + A_14 * (g1 @ g4)
    )
    return H, {
        "t": t,
        "Delta": Delta,
        "mu": mu,
        "A_12": A_12,
        "A_23": A_23,
        "A_34": A_34,
        "A_14": A_14,
    }


def expm_hermitian(H, tau):
    """Compute U = exp(-i H tau) for Hermitian 4x4 H using eigendecomposition."""
    w, v = np.linalg.eigh(H)
    phases = np.exp(-1j * w * tau)
    U = (v * phases) @ v.conj().T
    return U


def exp_gamma_pair(gamma_a, gamma_b, theta):
    """Compute U = exp(theta * gamma_a gamma_b).

    For Majoranas with {gamma_a, gamma_b} = 2 delta_ab,
    the product gamma_a gamma_b is anti-Hermitian, so exp(theta * gamma_a gamma_b)
    is unitary for real theta.
    """
    M = gamma_a @ gamma_b
    # Diagonalize M and exponentiate
    w, v = np.linalg.eig(M)
    phases = np.exp(theta * w)
    U = (v * phases) @ v.conj().T
    return U


def run_crosscheck_example():
    g1, g2, g3, g4 = gamma_matrices()

    # Example from notes: a=0, d=0, b=1, c=0, mu0=0
    a = 0.0
    b = 1.0
    c = 0.0
    d = 0.0
    mu0 = 0.0

    H, params = build_majorana_H_from_R(a, b, c, d, mu0)

    # Time evolution for tau = pi/2
    tau = np.pi / 2.0
    U_time = expm_hermitian(H, tau)

    # Ideal braid generator U = exp((pi/4) gamma2 gamma3)
    theta = np.pi / 4.0
    U_ideal = exp_gamma_pair(g2, g3, theta)

    # Compare
    diff = U_time - U_ideal
    max_abs_diff = np.max(np.abs(diff))

    print("R parameters: a=%.3f, b=%.3f, c=%.3f, d=%.3f, mu0=%.3f" % (a, b, c, d, mu0))
    print("Derived parameters:", params)
    print("Max |U_time - U_ideal| =", max_abs_diff)

    # Also print eigenvalues to verify spectrum agreement
    evals_time = np.linalg.eigvals(U_time)
    evals_ideal = np.linalg.eigvals(U_ideal)
    print("eig(U_time):", np.round(np.sort_complex(evals_time), 6))
    print("eig(U_ideal):", np.round(np.sort_complex(evals_ideal), 6))


def run_param_variants():
    """Scan a few parameter choices, including non-zero d, and
    compare U_time with the ideal braid generator.

    Cases:
      1) baseline: d = 0, mu0 = 0 (exact match expected)
      2) d != 0 but mu adjusted so that mu = 0 (quadratic H unchanged,
         still exact in this approximation -> tests "mu renormalization")
      3) d != 0 with mu0 = 0 so mu != 0 (chemical potential present,
         should spoil the pure braid and give a non-zero mismatch).
    """
    _, g2, g3, _ = gamma_matrices()

    theta = np.pi / 4.0
    U_ideal = exp_gamma_pair(g2, g3, theta)
    tau = np.pi / 2.0

    cases = [
        {"name": "baseline d=0, mu0=0", "a": 0.0, "b": 1.0, "c": 0.0, "d": 0.0, "mu0": 0.0},
        # choose mu0 so that mu = 4d + mu0 = 0
        {"name": "d!=0, mu adjusted so mu=0", "a": 0.0, "b": 1.0, "c": 0.0, "d": 0.2, "mu0": -0.8},
        {"name": "d!=0, mu0=0 (mu!=0)", "a": 0.0, "b": 1.0, "c": 0.0, "d": 0.2, "mu0": 0.0},
    ]

    for case in cases:
        a = case["a"]
        b = case["b"]
        c = case["c"]
        d = case["d"]
        mu0 = case["mu0"]

        H, params = build_majorana_H_from_R(a, b, c, d, mu0)
        U_time = expm_hermitian(H, tau)

        diff = U_time - U_ideal
        max_abs_diff = np.max(np.abs(diff))

        print("\nCase:", case["name"])
        print("  R params: a=%.3f, b=%.3f, c=%.3f, d=%.3f, mu0=%.3f" % (a, b, c, d, mu0))
        print("  Derived:", params)
        print("  Max |U_time - U_ideal| =", max_abs_diff)


if __name__ == "__main__":
    run_crosscheck_example()
    run_param_variants()
