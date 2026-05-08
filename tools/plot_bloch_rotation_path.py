#!/usr/bin/env python3
"""Visualize the constructive x-y-x Bloch-rotation path on the Bloch sphere.

The script shows two linked views:
  1. the Bloch-sphere trajectory of an initial state under the piecewise
     control path
  2. the corresponding control schedule u(t) = 0 -> pi/2 -> 0

By default it can either:
  - use explicit Euler angles (alpha, beta, gamma), or
  - fit a common gate via `tools/try_common_gates.py` and visualize that gate.

Output files are written to `results/bloch_rotation_paths/`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.append(str(SCRIPT_DIR))

from try_common_gates import fit_gate  # noqa: E402
from verify_su2_x_y_x import I2, SX, SY, SZ, path_gate, rotation  # noqa: E402


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def bloch_vector(psi: np.ndarray) -> np.ndarray:
    """Return the Bloch vector of a normalized qubit state."""
    psi = psi / np.linalg.norm(psi)
    x = np.vdot(psi, SX @ psi).real
    y = np.vdot(psi, SY @ psi).real
    z = np.vdot(psi, SZ @ psi).real
    return np.array([x, y, z], dtype=float)


def state_after_sequence(alpha: float, beta: float, gamma: float, J: float, n: int = 200):
    """Sample the Bloch vector trajectory under the x-y-x path.

    Returns arrays (times, traj, controls), where traj has shape (N, 3).
    """
    T1 = alpha / (2.0 * J)
    T2 = beta / (2.0 * J)
    T3 = gamma / (2.0 * J)

    # Initial qubit state |0>.
    psi0 = np.array([1.0 + 0.0j, 0.0 + 0.0j], dtype=complex)

    def rot(axis: str, angle: float) -> np.ndarray:
        return rotation(axis, angle)

    # Sample times in each segment.
    n1 = max(2, n // 3)
    n2 = max(2, n // 3)
    n3 = max(2, n - n1 - n2)

    t1 = np.linspace(0.0, T1, n1, endpoint=False)
    t2 = np.linspace(0.0, T2, n2, endpoint=False)
    t3 = np.linspace(0.0, T3, n3)

    times = []
    traj = []
    controls = []

    # Segment 1: x-rotation by alpha.
    for t in t1:
        psi = rot("x", 2.0 * J * t) @ psi0
        times.append(t)
        traj.append(bloch_vector(psi))
        controls.append(0.0)

    U1 = rot("x", alpha)

    # Segment 2: y-rotation by beta.
    for t in t2:
        psi = rot("y", 2.0 * J * t) @ (U1 @ psi0)
        times.append(T1 + t)
        traj.append(bloch_vector(psi))
        controls.append(np.pi / 2.0)

    U2 = rot("y", beta) @ U1

    # Segment 3: x-rotation by gamma.
    for t in t3:
        psi = rot("x", 2.0 * J * t) @ (U2 @ psi0)
        times.append(T1 + T2 + t)
        traj.append(bloch_vector(psi))
        controls.append(0.0)

    times = np.asarray(times, dtype=float)
    traj = np.asarray(traj, dtype=float)
    controls = np.asarray(controls, dtype=float)
    return times, traj, controls, (T1, T2, T3)


def draw_bloch_sphere(ax) -> None:
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, color="lightgray", alpha=0.12, linewidth=0)

    # Great circles.
    t = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(t), np.sin(t), 0 * t, color="0.55", lw=1.0, alpha=0.7)
    ax.plot(np.cos(t), 0 * t, np.sin(t), color="0.55", lw=1.0, alpha=0.7)
    ax.plot(0 * t, np.cos(t), np.sin(t), color="0.55", lw=1.0, alpha=0.7)

    # Axes.
    ax.quiver(0, 0, 0, 1.05, 0, 0, color="tab:red", linewidth=1.2, arrow_length_ratio=0.08)
    ax.quiver(0, 0, 0, 0, 1.05, 0, color="tab:green", linewidth=1.2, arrow_length_ratio=0.08)
    ax.quiver(0, 0, 0, 0, 0, 1.05, color="tab:blue", linewidth=1.2, arrow_length_ratio=0.08)

    ax.text(1.12, 0, 0, r"$x$", color="tab:red")
    ax.text(0, 1.12, 0, r"$y$", color="tab:green")
    ax.text(0, 0, 1.12, r"$z$", color="tab:blue")

    ax.set_xlim([-1.1, 1.1])
    ax.set_ylim([-1.1, 1.1])
    ax.set_zlim([-1.1, 1.1])
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel(r"$\langle\sigma_x\rangle$")
    ax.set_ylabel(r"$\langle\sigma_y\rangle$")
    ax.set_zlabel(r"$\langle\sigma_z\rangle$")
    ax.view_init(elev=20, azim=35)


def plot_path(alpha: float, beta: float, gamma: float, J: float, title: str, out_path: Path) -> Path:
    times, traj, controls, (T1, T2, T3) = state_after_sequence(alpha, beta, gamma, J)
    target = path_gate(alpha, beta, gamma, J)[0]
    psi_final = target @ np.array([1.0 + 0.0j, 0.0 + 0.0j], dtype=complex)
    bloch_final = bloch_vector(psi_final)

    fig = plt.figure(figsize=(10.5, 6.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.0, 1.0])

    ax3d = fig.add_subplot(gs[:, 0], projection="3d")
    draw_bloch_sphere(ax3d)
    ax3d.plot(traj[:, 0], traj[:, 1], traj[:, 2], color="tab:purple", lw=2.5)
    ax3d.scatter([traj[0, 0]], [traj[0, 1]], [traj[0, 2]], color="black", s=35, label="start")
    ax3d.scatter([bloch_final[0]], [bloch_final[1]], [bloch_final[2]], color="tab:orange", s=40, label="final")
    ax3d.set_title(title)
    ax3d.legend(loc="upper left")

    ax_u = fig.add_subplot(gs[0, 1])
    ax_u.plot(times, controls, color="tab:blue", lw=2)
    ax_u.axvline(T1, color="0.6", ls="--", lw=1)
    ax_u.axvline(T1 + T2, color="0.6", ls="--", lw=1)
    ax_u.set_ylabel(r"$u(t)$")
    ax_u.set_yticks([0.0, np.pi / 2.0])
    ax_u.set_yticklabels([r"$0$", r"$\pi/2$"])
    ax_u.set_xlim([times.min(), times.max()])
    ax_u.grid(True, alpha=0.25)

    ax_b = fig.add_subplot(gs[1, 1])
    ax_b.plot(times, traj[:, 0], label=r"$\langle\sigma_x\rangle$", lw=1.8)
    ax_b.plot(times, traj[:, 1], label=r"$\langle\sigma_y\rangle$", lw=1.8)
    ax_b.plot(times, traj[:, 2], label=r"$\langle\sigma_z\rangle$", lw=1.8)
    ax_b.set_xlabel("time")
    ax_b.set_ylabel("Bloch components")
    ax_b.grid(True, alpha=0.25)
    ax_b.legend(loc="best", fontsize=9)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize constructive x-y-x Bloch-rotation paths.")
    parser.add_argument("--gate", type=str, default="H", help="Common gate name or CUSTOM.")
    parser.add_argument("--alpha", type=float, default=1.1, help="Custom alpha angle.")
    parser.add_argument("--beta", type=float, default=0.7, help="Custom beta angle.")
    parser.add_argument("--gamma", type=float, default=0.4, help="Custom gamma angle.")
    parser.add_argument("--J", type=float, default=2.3, help="Control strength.")
    parser.add_argument(
        "--outdir",
        type=str,
        default=str(REPO_ROOT / "results" / "bloch_rotation_paths"),
        help="Directory to store the figure.",
    )
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    ensure_dir(out_dir)

    if args.gate.upper() == "CUSTOM":
        alpha, beta, gamma = args.alpha, args.beta, args.gamma
        title = rf"Custom x-y-x Bloch rotation: $({alpha:.3g}, {beta:.3g}, {gamma:.3g})$"
        out_path = out_dir / "custom_bloch_path.png"
    else:
        gate = args.gate.upper()
        result = fit_gate(gate, args.J)
        alpha, beta, gamma = result.alpha, result.beta, result.gamma
        title = (
            rf"{gate} gate via x-y-x path: "
            rf"$({alpha:.3g}, {beta:.3g}, {gamma:.3g})$, "
            rf"err={result.error:.1e}"
        )
        out_path = out_dir / f"{gate.lower()}_bloch_path.png"

    saved = plot_path(alpha, beta, gamma, args.J, title, out_path)
    print("Saved:", saved)


if __name__ == "__main__":
    main()