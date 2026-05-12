#!/usr/bin/env python3
"""Verify the path-derived braid gate used in ver5.md.

This script:
- builds the exact two-spin Hamiltonian H_4(t) from the designed g_i(t) path
- integrates each of the three path segments by a midpoint product formula
- composes the six-step total gate R_{4,F}
- projects the result to the logical subspace span{|01>, |10>}
- compares the logical gate to the canonical braid targets exp(-i pi/4 sigma_a)
- checks braid-relation and constant-YBE residuals for the two-body gate
- saves diagnostics to results/path_braid/ by default and emits a compact PNG summary

The path conventions match the formulas in ver5.md:
  g^(1)(u) = (u, 0, 1-u, 1)
  g^(2)(u) = (1-u, u, 0, 1)
  g^(3)(u) = (0, 1-u, u, 1)

The Hamiltonian is written in the Pauli basis as
  H_4(t) = t_c [ -g1 g2 ZI - g1 g3 YX + g1 g4 YY - g2 g3 XX - g2 g4 XY - g3 g4 IZ ].

If you want the convention with an explicit +h.c. factor, multiply the whole
Hamiltonian by 2 via --hc-factor 2.0.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from scipy.linalg import expm, norm


# Pauli matrices.
I2 = np.array([[1, 0], [0, 1]], dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


LOGICAL_IDXS = (1, 2)  # |01>, |10>


def kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.kron(a, b)


def h4_from_path(g1: float, g2: float, g3: float, g4: float, tc: float = 1.0, hc_factor: float = 1.0) -> np.ndarray:
    """Return the 4x4 Hamiltonian from the designed path."""
    h4 = (
        -(g1 * g2) * kron(SZ, I2)
        -(g1 * g3) * kron(SY, SX)
        +(g1 * g4) * kron(SY, SY)
        -(g2 * g3) * kron(SX, SX)
        -(g2 * g4) * kron(SX, SY)
        -(g3 * g4) * kron(I2, SZ)
    )
    return hc_factor * tc * h4


def path_segment(k: int, u: float) -> Tuple[float, float, float, float]:
    """Return the segment-k control point at local parameter u in [0, 1]."""
    if k == 1:
        return u, 0.0, 1.0 - u, 1.0
    if k == 2:
        return 1.0 - u, u, 0.0, 1.0
    if k == 3:
        return 0.0, 1.0 - u, u, 1.0
    raise ValueError(f"segment index must be 1, 2, or 3; got {k}")


def integrate_segment(k: int, steps: int, T: float, tc: float, hc_factor: float) -> np.ndarray:
    """Midpoint product approximation to Texp[-i T ∫_0^1 H_k(u) du]."""
    dt = T / steps
    U = np.eye(4, dtype=complex)
    # Multiply in reverse time order so later times appear on the left.
    for u in ((np.arange(steps) + 0.5) / steps)[::-1]:
        g1, g2, g3, g4 = path_segment(k, float(u))
        H = h4_from_path(g1, g2, g3, g4, tc=tc, hc_factor=hc_factor)
        U = expm(-1j * dt * H) @ U
    return U


def compose_cycle(steps: int, T: float, tc: float, hc_factor: float) -> Dict[str, np.ndarray]:
    """Return the three segment gates and the six-step total gate."""
    R1 = integrate_segment(1, steps=steps, T=T, tc=tc, hc_factor=hc_factor)
    R2 = integrate_segment(2, steps=steps, T=T, tc=tc, hc_factor=hc_factor)
    R3 = integrate_segment(3, steps=steps, T=T, tc=tc, hc_factor=hc_factor)
    R_cycle = R3 @ R2 @ R1
    R_total = R_cycle @ R_cycle
    return {"R1": R1, "R2": R2, "R3": R3, "R_cycle": R_cycle, "R_total": R_total}


def logical_project(U4: np.ndarray) -> np.ndarray:
    """Project a 4x4 gate to the logical subspace span{|01>, |10>}."""
    idx = np.array(LOGICAL_IDXS)
    return U4[np.ix_(idx, idx)]


def phase_aligned_residual(U: np.ndarray, V: np.ndarray) -> float:
    """Minimize || e^{i phi} U - V ||_F over the global phase phi."""
    overlap = np.vdot(U, V)
    if abs(overlap) < 1e-15:
        return float(norm(U - V, ord="fro"))
    phase = np.exp(1j * np.angle(overlap))
    return float(norm(phase * U - V, ord="fro"))


def phase_aligned_gate(U: np.ndarray, V: np.ndarray) -> np.ndarray:
    """Return e^{i phi} U with phi chosen to align U with V."""
    overlap = np.vdot(U, V)
    if abs(overlap) < 1e-15:
        return U.copy()
    phase = np.exp(1j * np.angle(overlap))
    return phase * U


def unitary_error(U: np.ndarray) -> float:
    return float(norm(U.conj().T @ U - np.eye(U.shape[0]), ord="fro"))


def logical_braid_target(axis: str) -> np.ndarray:
    axis = axis.lower()
    if axis == "x":
        sigma = SX
    elif axis == "y":
        sigma = SY
    elif axis == "z":
        sigma = SZ
    else:
        raise ValueError(f"axis must be one of x/y/z, got {axis}")
    return expm(-1j * (np.pi / 4.0) * sigma)


def best_logical_axis(RL: np.ndarray) -> Tuple[str, Dict[str, float]]:
    residuals = {axis: phase_aligned_residual(RL, logical_braid_target(axis)) for axis in ("x", "y", "z")}
    best_axis = min(residuals, key=residuals.get)
    return best_axis, residuals


def swap_23() -> np.ndarray:
    """Swap qubits 2 and 3 in the 3-qubit basis |q1 q2 q3>."""
    S = np.zeros((8, 8), dtype=complex)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                in_idx = (a << 2) | (b << 1) | c
                out_idx = (a << 2) | (c << 1) | b
                S[out_idx, in_idx] = 1.0
    return S


def embed_two_body_gate(R4: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Embed a 4x4 two-body gate on qubits (1,2), (1,3), (2,3)."""
    R12 = np.kron(R4, I2)
    R23 = np.kron(I2, R4)
    S23 = swap_23()
    R13 = S23 @ R12 @ S23
    return R12, R13, R23


def braid_relation_residual(R4: np.ndarray) -> float:
    R12, _, R23 = embed_two_body_gate(R4)
    lhs = R12 @ R23 @ R12
    rhs = R23 @ R12 @ R23
    return float(norm(lhs - rhs, ord="fro"))


def ybe_residual(R4: np.ndarray) -> float:
    R12, R13, R23 = embed_two_body_gate(R4)
    lhs = R12 @ R13 @ R23
    rhs = R23 @ R13 @ R12
    return float(norm(lhs - rhs, ord="fro"))


def default_stem(steps: int, T: float, tc: float, axis: str, hc_factor: float) -> str:
    def slug(value: float) -> str:
        return format(value, "g").replace("-", "m").replace(".", "p")

    return f"path_braid_steps{steps}_T{slug(T)}_tc{slug(tc)}_axis{axis}_hc{slug(hc_factor)}"


def save_summary_figure(fig_path: Path, R_logical: np.ndarray, target: np.ndarray, axis_residuals: Dict[str, float], braid_res: float, ybe_res: float, best_axis: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    aligned = phase_aligned_gate(R_logical, target)
    diff = np.abs(aligned - target)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), dpi=200)

    labels = ["x", "y", "z", "braid", "YBE"]
    values = [axis_residuals["x"], axis_residuals["y"], axis_residuals["z"], braid_res, ybe_res]
    colors = ["#4e79a7", "#59a14f", "#f28e2b", "#e15759", "#76b7b2"]
    axes[0].bar(labels, values, color=colors)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Frobenius residual")
    axes[0].set_title(f"Residuals (best axis: {best_axis})")
    axes[0].tick_params(axis="x", rotation=20)
    for idx, value in enumerate(values):
        axes[0].text(idx, value * 1.08, f"{value:.2e}", ha="center", va="bottom", fontsize=7)

    im = axes[1].imshow(diff, cmap="magma")
    axes[1].set_title(r"$|e^{i\phi}R_L-U_{\mathrm{braid}}|$")
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels(["01", "10"])
    axes[1].set_yticks([0, 1])
    axes[1].set_yticklabels(["01", "10"])
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)

    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--steps", type=int, default=400, help="midpoint steps per segment")
    p.add_argument("--T", type=float, default=1.0, help="physical duration of each segment")
    p.add_argument("--tc", type=float, default=1.0, help="base coupling t_c")
    p.add_argument("--hc-factor", type=float, default=1.0, help="overall factor for +h.c. convention")
    p.add_argument("--axis", choices=("x", "y", "z"), default="x", help="logical braid axis to compare against")
    p.add_argument("--outdir", type=str, default="results/path_braid", help="directory for saved diagnostics")
    p.add_argument("--stem", type=str, default="", help="optional file stem override")
    p.add_argument("--save", type=str, default="", help="optional .npz output path override")
    p.add_argument("--fig", type=str, default="", help="optional .png figure path override")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    print("Building path-derived gates from ver5.md path")
    gates = compose_cycle(steps=args.steps, T=args.T, tc=args.tc, hc_factor=args.hc_factor)
    R1 = gates["R1"]
    R2 = gates["R2"]
    R3 = gates["R3"]
    R_cycle = gates["R_cycle"]
    R_total = gates["R_total"]

    print(f"unitarity ||R1^†R1-I||_F = {unitary_error(R1):.3e}")
    print(f"unitarity ||R2^†R2-I||_F = {unitary_error(R2):.3e}")
    print(f"unitarity ||R3^†R3-I||_F = {unitary_error(R3):.3e}")
    print(f"unitarity ||R_cycle^†R_cycle-I||_F = {unitary_error(R_cycle):.3e}")
    print(f"unitarity ||R_total^†R_total-I||_F = {unitary_error(R_total):.3e}")

    R_logical = logical_project(R_total)
    print(f"logical gate ||R_L^†R_L-I||_F = {unitary_error(R_logical):.3e}")

    target = logical_braid_target(args.axis)
    target_residual = phase_aligned_residual(R_logical, target)
    best_axis, axis_residuals = best_logical_axis(R_logical)

    print(f"logical target axis = {args.axis}")
    print(f"phase-aligned residual to exp(-i pi/4 sigma_{args.axis}) = {target_residual:.3e}")
    print("residual scan over x/y/z ->", {k: f"{v:.3e}" for k, v in axis_residuals.items()})
    print(f"best matching logical axis = {best_axis}")

    braid_res = braid_relation_residual(R_total)
    ybe_res = ybe_residual(R_total)
    print(f"braid relation residual ||R12 R23 R12 - R23 R12 R23||_F = {braid_res:.3e}")
    print(f"constant YBE residual ||R12 R13 R23 - R23 R13 R12||_F = {ybe_res:.3e}")

    stem = args.stem.strip() or default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor)
    outdir = Path(args.outdir)
    npz_path = Path(args.save) if args.save.strip() else outdir / f"{stem}.npz"
    fig_path = Path(args.fig) if args.fig.strip() else outdir / f"{stem}.png"

    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        npz_path,
        R1=R1,
        R2=R2,
        R3=R3,
        R_cycle=R_cycle,
        R_total=R_total,
        R_logical=R_logical,
        target=target,
        best_axis=np.array(best_axis),
        braid_residual=braid_res,
        ybe_residual=ybe_res,
        target_residual=target_residual,
        axis_residuals=np.array([axis_residuals["x"], axis_residuals["y"], axis_residuals["z"]]),
    )
    print(f"saved diagnostics to {npz_path}")

    try:
        fig_path.parent.mkdir(parents=True, exist_ok=True)
        save_summary_figure(fig_path, R_logical, target, axis_residuals, braid_res, ybe_res, best_axis)
        print(f"saved summary figure to {fig_path}")
    except Exception as exc:
        print(f"warning: figure generation failed: {exc}")


if __name__ == "__main__":
    main()