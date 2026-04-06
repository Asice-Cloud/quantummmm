import numpy as np
import matplotlib.pyplot as plt
import os

"""Many-body check of the effect of a true interaction term (n1 n2)
 on the instantaneous braid generator.

We work with 2 complex fermions (4-dimensional Fock space), represented
via 4 Majorana operators gamma_a (a=1..4) as in run_instantaneous_braid_crosscheck.py.

- Free (quadratic) Hamiltonian H_free corresponds to a 2-site Kitaev chain
  with parameters (t, Delta, mu) chosen so that
      U_free(tau=pi/2) = exp((pi/4) gamma_2 gamma_3)
  holds exactly (this is the d=0 example from the cross-check script).

- We then add a genuine density-density interaction term
      H_int = U * n1 n2,
  with n_j = c_j^\dagger c_j, where c_j are constructed from gamma's.

- For several values of U we compute
      U_full = exp(-i (H_free + H_int) tau)
  and compare it to the ideal braid generator
      U_ideal = exp((pi/4) gamma_2 gamma_3).

This isolates the effect of a quartic interaction term on the braid,
separate from the chemical-potential (mu) channel that was already
captured by the BdG mapping.
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


def fermion_operators_from_gamma():
    """Construct complex fermion operators c1,c2 and densities n1,n2
    from the Majoranas.

    Standard relations:
      c1 = (gamma1 + i gamma2)/2
      c2 = (gamma3 + i gamma4)/2
      n_j = c_j^\dagger c_j.
    """
    g1, g2, g3, g4 = gamma_matrices()
    c1 = (g1 + 1j * g2) / 2.0
    c2 = (g3 + 1j * g4) / 2.0
    c1_dag = c1.conj().T
    c2_dag = c2.conj().T
    n1 = c1_dag @ c1
    n2 = c2_dag @ c2
    return (c1, c1_dag, n1), (c2, c2_dag, n2), (g1, g2, g3, g4)


def build_H_free():
    """Construct the free (quadratic) Hamiltonian in the 4-dim Fock space.

    We choose parameters matching the cross-check example:
      t = 1, Delta = 1, mu = 0.

    In Majorana form for 2 sites this reduces to
      H_free = (i/2) * gamma_2 gamma_3.
    """
    g1, g2, g3, g4 = gamma_matrices()
    H_free = 0.5j * (g2 @ g3)
    return H_free, (g1, g2, g3, g4)


def expm_hermitian(H, tau):
    w, v = np.linalg.eigh(H)
    phases = np.exp(-1j * w * tau)
    U = (v * phases) @ v.conj().T
    return U


def exp_gamma_pair(gamma_a, gamma_b, theta):
    M = gamma_a @ gamma_b
    w, v = np.linalg.eig(M)
    phases = np.exp(theta * w)
    U = (v * phases) @ v.conj().T
    return U


def run_interacting_scan():
    # Build objects
    (c1, c1_dag, n1), (c2, c2_dag, n2), (g1, g2, g3, g4) = fermion_operators_from_gamma()

    H_free, (g1, g2, g3, g4) = build_H_free()

    # Ideal braid unitary
    theta = np.pi / 4.0
    U_ideal = exp_gamma_pair(g2, g3, theta)

    tau = np.pi / 2.0

    print("H_free check (U_free vs ideal):")
    U_free = expm_hermitian(H_free, tau)
    diff_free = U_free - U_ideal
    print("  Max |U_free - U_ideal| =", np.max(np.abs(diff_free)))

    # Now add interaction U * n1 n2
    U_list = [0.0, 0.1, 0.2, 0.5, 1.0]
    errors = []
    for U_int in U_list:
      H_int = U_int * (n1 @ n2)
      H_full = H_free + H_int
      U_full = expm_hermitian(H_full, tau)
      diff = U_full - U_ideal
      max_abs_diff = np.max(np.abs(diff))
      errors.append(max_abs_diff)
      print("U_int = %.3f -> Max |U_full - U_ideal| = %.6f" % (U_int, max_abs_diff))

    # Plot braid deviation vs U_int
    out_dir = os.path.dirname(__file__)
    plt.figure()
    plt.plot(U_list, errors, marker='o')
    plt.xlabel(r'$U$ (interaction strength)')
    plt.ylabel(r'$\max\,\|U_{\mathrm{full}}-U_{\mathrm{ideal}}\|$')
    plt.title('Braid deviation vs interaction U (mu=0)')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'braid_deviation_vs_U.png'), dpi=150)
    plt.close()


def run_mu_U_scan():
    """Scan a small grid in (mu, U_int) and report braid deviation.

    Model:
      H = H_free + H_mu + H_int,
      H_free = (i/2) gamma2 gamma3,
      H_mu   = -(i mu / 4) (gamma1 gamma2 + gamma3 gamma4),
      H_int  = U_int n1 n2.
    """
    (c1, c1_dag, n1), (c2, c2_dag, n2), (g1, g2, g3, g4) = fermion_operators_from_gamma()
    H_free, (g1, g2, g3, g4) = build_H_free()

    theta = np.pi / 4.0
    U_ideal = exp_gamma_pair(g2, g3, theta)
    tau = np.pi / 2.0

    mu_vals = [-1.0, -0.5, 0.0, 0.5, 1.0]
    U_list = [0.0, 0.1, 0.2, 0.5, 1.0]

    errors = np.zeros((len(mu_vals), len(U_list)))

    print("\n(mu, U_int) scan:")
    for i_mu, mu in enumerate(mu_vals):
      for i_U, U_int in enumerate(U_list):
        H_mu = -(1j * mu / 4.0) * (g1 @ g2 + g3 @ g4)
        H_int = U_int * (n1 @ n2)
        H_full = H_free + H_mu + H_int
        U_full = expm_hermitian(H_full, tau)
        diff = U_full - U_ideal
        max_abs_diff = np.max(np.abs(diff))
        errors[i_mu, i_U] = max_abs_diff
        print("  mu = %+4.1f, U_int = %.3f -> Max |U - U_ideal| = %.6f" % (mu, U_int, max_abs_diff))

    # 2D heatmap of error vs (mu, U)
    out_dir = os.path.dirname(__file__)
    plt.figure()
    im = plt.imshow(
      errors,
      origin='lower',
      aspect='auto',
      extent=[min(U_list), max(U_list), min(mu_vals), max(mu_vals)],
      cmap='viridis'
    )
    plt.colorbar(im, label=r'$\max\,\|U_{\mathrm{full}}-U_{\mathrm{ideal}}\|$')
    plt.xlabel(r'$U$')
    plt.ylabel(r'$\mu$')
    plt.title('Braid deviation in $(\mu,U)$ plane')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'braid_deviation_mu_U_heatmap.png'), dpi=150)
    plt.close()


if __name__ == "__main__":
  run_interacting_scan()
  run_mu_U_scan()
