#!/usr/bin/env python3
"""
Subspace LDOS-style proof for the eight-vertex model.

This script compares:
1. Subspace-projected spectral density of the local 4x4 H4(u, delta)
2. Full-chain zero-energy LDOS as a control

Goal:
- Show that delta changes the local subspace spectral weight pattern
- Make clear that full-chain LDOS can remain unchanged even when subspace geometry changes

Output:
results/workflow/subspace_ldos_proof/
"""

import argparse
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from config import OUTPUT_DIR, FIG_DPI, FIG_FORMAT, VERBOSE
from step1_eight_vertex import eight_vertex_H4, process_u_delta_point
from step3_full_chain_bdg import analyze_full_chain

OUT_DIR = OUTPUT_DIR / "subspace_ldos_proof"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def spectral_density(weights, energies_grid, eigvals, eigvecs, eta=1e-2):
    """Lorentz-broadened spectral density for a set of basis weights."""
    density = np.zeros_like(energies_grid, dtype=float)
    norm = eta / np.pi
    for n, en in enumerate(eigvals):
        lor = norm / ((energies_grid - en) ** 2 + eta ** 2)
        density += weights[n] * lor
    return density


def local_subspace_spectral(u, delta, energies_grid, eta=1e-2):
    """Project the 4x4 local spectrum onto the {|01>, |10>} subspace."""
    h4 = eight_vertex_H4(u, delta)
    eigvals, eigvecs = np.linalg.eigh(h4)

    # Computational basis ordering from kron: |00>, |01>, |10>, |11>
    subspace_indices = [1, 2]

    weights_sub = []
    weights_01 = []
    weights_10 = []
    for n in range(len(eigvals)):
        vec = eigvecs[:, n]
        w01 = float(np.abs(vec[1]) ** 2)
        w10 = float(np.abs(vec[2]) ** 2)
        weights_01.append(w01)
        weights_10.append(w10)
        weights_sub.append(w01 + w10)

    weights_sub = np.array(weights_sub)
    weights_01 = np.array(weights_01)
    weights_10 = np.array(weights_10)

    rho_sub = spectral_density(weights_sub, energies_grid, eigvals, eigvecs, eta=eta)
    rho_01 = spectral_density(weights_01, energies_grid, eigvals, eigvecs, eta=eta)
    rho_10 = spectral_density(weights_10, energies_grid, eigvals, eigvecs, eta=eta)

    return {
        "eigvals": eigvals,
        "eigvecs": eigvecs,
        "rho_sub": rho_sub,
        "rho_01": rho_01,
        "rho_10": rho_10,
        "weights_sub": weights_sub,
        "weights_01": weights_01,
        "weights_10": weights_10,
    }


def make_figure(u, delta_list, energies_grid, eta=1e-2):
    fig, axes = plt.subplots(3, len(delta_list), figsize=(5.5 * len(delta_list), 12), sharex=True)
    if len(delta_list) == 1:
        axes = np.array([[axes[0]], [axes[1]], [axes[2]]])

    for col, delta in enumerate(delta_list):
        local = local_subspace_spectral(u, delta, energies_grid, eta=eta)
        full = analyze_full_chain(u, delta, L_list=[160])

        ax0 = axes[0, col]
        ax0.plot(energies_grid, local["rho_sub"], color="black", lw=2, label="subspace total")
        ax0.plot(energies_grid, local["rho_01"], color="tab:blue", ls="--", label="|01> channel")
        ax0.plot(energies_grid, local["rho_10"], color="tab:orange", ls=":", label="|10> channel")
        ax0.set_title(f"Local H4 spectral weight\n$u={u:.3f}$, $\\delta={delta:.3f}$")
        ax0.set_ylabel("spectral weight")
        ax0.grid(True, alpha=0.3)
        ax0.legend(fontsize=9)

        ax1 = axes[1, col]
        polarization = local["rho_01"] - local["rho_10"]
        ax1.plot(energies_grid, polarization, color="tab:purple", lw=2)
        ax1.axhline(0.0, color="gray", lw=1, ls="--")
        ax1.set_title(f"Subspace polarization $\\rho_{{01}}-\\rho_{{10}}$\n$u={u:.3f}$, $\\delta={delta:.3f}$")
        ax1.set_ylabel("polarization")
        ax1.grid(True, alpha=0.3)

        ax2 = axes[2, col]
        ldos = full["ldos"]
        ldos_es = full["ldos_energies"]
        zero_idx = int(np.argmin(np.abs(ldos_es)))
        ldos_zero = ldos[:, zero_idx]
        sites = np.arange(full["L_max"])
        ax2.bar(sites, ldos_zero, color="steelblue", alpha=0.75)
        ax2.set_title(f"Full-chain LDOS at E≈0\n$u={u:.3f}$, $\\delta={delta:.3f}$")
        ax2.set_xlabel("site index")
        ax2.set_ylabel("LDOS(E≈0)")
        ax2.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    out_file = OUT_DIR / f"subspace_ldos_proof.{FIG_FORMAT}"
    plt.savefig(out_file, dpi=FIG_DPI)
    plt.close(fig)
    return out_file


def write_summary(u, delta_list, out_file):
    lines = []
    lines.append("Subspace LDOS-style proof for the eight-vertex model")
    lines.append(f"u = {u:.6f}")
    lines.append(f"delta list = {', '.join(f'{d:.3f}' for d in delta_list)}")
    lines.append("")
    for delta in delta_list:
        pv = process_u_delta_point(u, delta)
        d = pv["d"]
        lines.append(f"delta={delta:.3f}: d=({d[0]:.6f}, {d[1]:.6f}, {d[2]:.6f}), |d|={pv['|d|']:.6f}")
        lines.append(f"           t={pv['t']:.6f}, Delta={pv['Delta']:.6f}, mu={pv['mu']:.6f}")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("- The upper panels compare subspace-projected spectral weight of H4(u, delta).")
    lines.append("- The middle panels show the subspace polarization rho_01 - rho_10.")
    lines.append("- The lower panels show the full-chain LDOS at E≈0 as a control.")
    lines.append("- If the subspace claim is the target, the local polarization is the most discriminating diagnostic.")
    lines.append("- Full-chain LDOS can remain unchanged because mu ≡ 0 in this model.")
    out_file.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Subspace LDOS proof for the eight-vertex model")
    parser.add_argument("--u", type=float, default=np.pi / 4, help="Path parameter u")
    parser.add_argument("--delta", type=str, default="0.0,0.1", help="Comma-separated delta values")
    parser.add_argument("--eta", type=float, default=1e-2, help="Lorentzian broadening for local spectral density")
    parser.add_argument("--emin", type=float, default=-2.5, help="Energy min for local spectrum")
    parser.add_argument("--emax", type=float, default=2.5, help="Energy max for local spectrum")
    parser.add_argument("--ne", type=int, default=1001, help="Energy grid size")
    args = parser.parse_args()

    delta_list = [float(x) for x in args.delta.split(",") if x.strip()]
    energies_grid = np.linspace(args.emin, args.emax, args.ne)

    if VERBOSE:
        print("[Subspace LDOS] Building comparison figure...")
        print(f"[Subspace LDOS] u = {args.u:.6f}")
        print(f"[Subspace LDOS] delta = {delta_list}")

    out_fig = make_figure(args.u, delta_list, energies_grid, eta=args.eta)
    summary_file = OUT_DIR / "subspace_ldos_summary.txt"
    write_summary(args.u, delta_list, summary_file)

    print(f"[Subspace LDOS] Saved figure: {out_fig}")
    print(f"[Subspace LDOS] Saved summary: {summary_file}")


if __name__ == "__main__":
    main()
