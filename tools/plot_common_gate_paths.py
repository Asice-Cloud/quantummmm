#!/usr/bin/env python3
"""Plot x-y-x control paths for common single-qubit gates.

This script uses the fitted Euler-angle parameters from `try_common_gates.py`
and visualizes the corresponding three-segment control path:

    segment 1: H_x = J sigma_x
    segment 2: H_y = J sigma_y
    segment 3: H_x = J sigma_x

For each gate it saves a PNG showing:
  - the piecewise-constant control axis label u(t)
  - the active Hamiltonian coefficients d_x(t), d_y(t)

The output directory is `results/common_gate_paths/`.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.append(str(SCRIPT_DIR))

from try_common_gates import default_gates, fit_gate  # noqa: E402


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def step_path(fractions: np.ndarray, values: np.ndarray):
    """Return stepwise sample points suitable for a piecewise-constant plot."""
    xs = [0.0]
    ys = [float(values[0])]
    cum = 0.0
    for frac, val in zip(fractions, values):
        cum_next = cum + float(frac)
        xs.extend([cum_next, cum_next])
        ys.extend([float(val), float(val)])
        cum = cum_next
    if xs[-1] < 1.0:
        xs.append(1.0)
        ys.append(float(values[-1]))
    return np.asarray(xs), np.asarray(ys)


def plot_gate(gate_name: str, out_dir: Path, J: float) -> Path:
    result = fit_gate(gate_name, J)

    durations = np.array([result.T1, result.T2, result.T3], dtype=float)
    total = float(durations.sum())

    if total <= 0.0:
        fractions = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        fractions = durations / total

    # Control labels and Hamiltonian coefficients for the x-y-x path.
    # The control parameter u(t) is 0 on x-steps and pi/2 on the y-step.
    u_values = np.array([0.0, np.pi / 2.0, 0.0], dtype=float)
    dx_values = np.array([J, 0.0, J], dtype=float)
    dy_values = np.array([0.0, J, 0.0], dtype=float)

    xu, yu = step_path(fractions, u_values)
    xd, yd = step_path(fractions, dx_values)
    _, ydy = step_path(fractions, dy_values)

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 6.5), sharex=True, constrained_layout=True)

    axes[0].plot(xu, yu, color="tab:blue", lw=2.2)
    axes[0].set_ylabel(r"$u(t)$")
    axes[0].set_yticks([0.0, np.pi / 2.0])
    axes[0].set_yticklabels([r"$0$", r"$\pi/2$"])
    axes[0].grid(True, alpha=0.25)
    axes[0].set_title(
        f"{gate_name} gate: x-y-x path | "
        f"(T1,T2,T3)=({result.T1:.4g}, {result.T2:.4g}, {result.T3:.4g}), "
        f"err={result.error:.2e}"
    )

    axes[1].plot(xd, yd, color="tab:red", lw=2.2, label=r"$d_x(t)$")
    axes[1].plot(xd, ydy, color="tab:green", lw=2.2, label=r"$d_y(t)$")
    axes[1].set_xlabel(r"normalized time $t/T$")
    axes[1].set_ylabel(r"Hamiltonian coefficients")
    axes[1].set_yticks([0.0, J])
    axes[1].set_yticklabels([r"$0$", rf"${J:.2f}$"])
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(loc="upper right")

    # Visualize segment boundaries.
    for ax in axes:
        ax.axvline(fractions[0], color="0.65", ls="--", lw=1)
        ax.axvline(fractions[0] + fractions[1], color="0.65", ls="--", lw=1)

    out_path = out_dir / f"{gate_name.lower()}_path.png"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot x-y-x control paths for common single-qubit gates.")
    parser.add_argument("--gate", type=str, default="ALL", help="Gate name or ALL.")
    parser.add_argument("--J", type=float, default=2.3, help="Control strength.")
    parser.add_argument(
        "--outdir",
        type=str,
        default=str(REPO_ROOT / "results" / "common_gate_paths"),
        help="Output directory for figures.",
    )
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    ensure_dir(out_dir)

    names = default_gates() if args.gate.upper() == "ALL" else [args.gate.upper()]

    print("Saving figures to", out_dir)
    for gate_name in names:
        out_path = plot_gate(gate_name, out_dir, args.J)
        print(f"{gate_name}: {out_path}")


if __name__ == "__main__":
    main()