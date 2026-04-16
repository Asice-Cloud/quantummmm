import numpy as np
import matplotlib.pyplot as plt

"""Scan parameter(s) in the 2D p+ip vortex Berry experiment and
compute the fidelity F_Dehn(mu) between the Berry holonomy and the
Ising Dehn twist (R^{σσ})^2 on the zero-mode subspace.

This script reuses the building blocks from run_pip_vortex_berry.py
and produces both a printed table and a PNG plot of F_Dehn vs mu.
"""

from run_pip_vortex_berry import (
    compute_zero_mode_holonomy,
    normalize_to_su2,
    ising_R_and_dehn_su2,
    su2_fidelity,
)


def scan_mu_F_Dehn(Lx: int = 20, Ly: int = 20,
                   mu_values=None,
                   t: float = 1.0,
                   delta0: float = 0.5,
                   dim_sub: int = 2,
                   n_steps_segment: int = 12):
    """Scan chemical potential mu and compute F_Dehn(mu).

    mu_values: iterable of mu's to scan. If None, use a default list.
    """
    if mu_values is None:
        # A simple set of representative values in the topological regime
        mu_values = np.linspace(-3.0, -0.5, 6)

    # Ising Dehn twist in SU(2)-normalized form
    _, Dehn_ising_su2 = ising_R_and_dehn_su2()

    results = []

    print("=== F_Dehn(mu) scan for 2D p+ip vortex loop ===")
    print("Lx = {}, Ly = {}, t = {}, delta0 = {}, dim_sub = {}, n_steps_segment = {}".format(
        Lx, Ly, t, delta0, dim_sub, n_steps_segment
    ))
    print("mu\tF_Dehn")

    for mu in mu_values:
        try:
            U_holo = compute_zero_mode_holonomy(
                Lx=Lx,
                Ly=Ly,
                mu=float(mu),
                t=t,
                delta0=delta0,
                dim_sub=dim_sub,
                n_steps_segment=n_steps_segment,
            )
        except Exception as e:
            print(f"mu={mu:.3f}: FAILED ({e})")
            results.append((mu, np.nan))
            continue

        U_holo_su2 = normalize_to_su2(U_holo)
        F_Dehn, _ = su2_fidelity(Dehn_ising_su2, U_holo_su2)
        results.append((mu, F_Dehn))
        print(f"{mu:.3f}\t{F_Dehn:.6f}")

    return np.array(results)


def plot_F_Dehn_vs_mu(results: np.ndarray, out_path: str = "pip_vortex_F_Dehn_vs_mu.png"):
    """Make a simple plot of F_Dehn vs mu and save it as a PNG."""
    mus = results[:, 0]
    Fvals = results[:, 1]

    plt.figure(figsize=(5, 3.5))
    plt.plot(mus, Fvals, marker="o")
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$F_{\mathrm{Dehn}}(\mu)$")
    plt.ylim(0.0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved F_Dehn vs mu plot to {out_path}")


def main():
    results = scan_mu_F_Dehn()
    # Filter out NaNs before plotting
    results_clean = results[~np.isnan(results[:, 1])]
    if results_clean.size > 0:
        plot_F_Dehn_vs_mu(results_clean)


if __name__ == "__main__":
    main()
