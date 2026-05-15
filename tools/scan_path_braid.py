#!/usr/bin/env python3
"""Grid scan over (steps, T) to minimize braid residual for the designed path.

Saves a .npz and a PNG heatmap under results/path_braid/.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple

import numpy as np
from scipy.linalg import expm, norm


# Pauli
I2 = np.array([[1, 0], [0, 1]], dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


LOGICAL_IDXS = (1, 2)


def kron(a, b):
    return np.kron(a, b)


def h4_from_path(g1, g2, g3, g4, tc=1.0, hc_factor=1.0):
    return (
        hc_factor
        * tc
        * (
            -(g1 * g2) * kron(SZ, I2)
            -(g1 * g3) * kron(SY, SX)
            +(g1 * g4) * kron(SY, SY)
            -(g2 * g3) * kron(SX, SX)
            -(g2 * g4) * kron(SX, SY)
            -(g3 * g4) * kron(I2, SZ)
        )
    )


def path_segment(k: int, u: float) -> Tuple[float, float, float, float]:
    if k == 1:
        return u, 0.0, 1.0 - u, 1.0
    if k == 2:
        return 1.0 - u, u, 0.0, 1.0
    if k == 3:
        return 0.0, 1.0 - u, u, 1.0
    raise ValueError("bad segment")


def integrate_segment(k: int, steps: int, T: float, tc: float, hc_factor: float):
    dt = T / steps
    U = np.eye(4, dtype=complex)
    for u in ((np.arange(steps) + 0.5) / steps)[::-1]:
        g1, g2, g3, g4 = path_segment(k, float(u))
        H = h4_from_path(g1, g2, g3, g4, tc=tc, hc_factor=hc_factor)
        U = expm(-1j * dt * H) @ U
    return U


def compose_cycle(steps: int, T: float, tc: float, hc_factor: float):
    R1 = integrate_segment(1, steps, T, tc, hc_factor)
    R2 = integrate_segment(2, steps, T, tc, hc_factor)
    R3 = integrate_segment(3, steps, T, tc, hc_factor)
    R_cycle = R3 @ R2 @ R1
    R_total = R_cycle @ R_cycle
    return R_total


def logical_project(U4):
    idx = np.array(LOGICAL_IDXS)
    return U4[np.ix_(idx, idx)]


def phase_aligned_residual(U, V):
    overlap = np.vdot(U, V)
    if abs(overlap) < 1e-16:
        return float(norm(U - V, ord="fro"))
    phase = np.exp(1j * np.angle(overlap))
    return float(norm(phase * U - V, ord="fro"))


def logical_braid_target(axis: str):
    if axis == "x":
        sigma = SX
    elif axis == "y":
        sigma = SY
    elif axis == "z":
        sigma = SZ
    else:
        raise ValueError(axis)
    return expm(-1j * (np.pi / 4.0) * sigma)


def embed_two_body_gate(R4):
    R12 = np.kron(R4, I2)
    R23 = np.kron(I2, R4)
    S = np.zeros((8, 8), dtype=complex)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                in_idx = (a << 2) | (b << 1) | c
                out_idx = (a << 2) | (c << 1) | b
                S[out_idx, in_idx] = 1.0
    R13 = S @ R12 @ S
    return R12, R13, R23


def braid_relation_residual(R4):
    R12, _, R23 = embed_two_body_gate(R4)
    return float(norm(R12 @ R23 @ R12 - R23 @ R12 @ R23, ord="fro"))


def ybe_residual(R4):
    R12, R13, R23 = embed_two_body_gate(R4)
    return float(norm(R12 @ R13 @ R23 - R23 @ R13 @ R12, ord="fro"))


def parse_list(text: str) -> List[float]:
    return [float(x) for x in text.split(",") if x.strip()]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--steps-list", type=str, default="40,80,120")
    p.add_argument("--T-list", type=str, default="0.5,1.0,2.0")
    p.add_argument("--tc", type=float, default=1.0)
    p.add_argument("--hc-factor", type=float, default=1.0)
    p.add_argument("--axis", choices=("x", "y", "z"), default="x")
    p.add_argument("--outdir", type=str, default="results/path_braid")
    return p.parse_args()


def main():
    args = parse_args()
    steps_list = [int(x) for x in parse_list(args.steps_list)]
    T_list = parse_list(args.T_list)
    tc = args.tc
    hc_factor = args.hc_factor
    axis = args.axis

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    residuals = np.zeros((len(steps_list), len(T_list)), dtype=float)
    braid_res = np.zeros_like(residuals)
    ybe_res = np.zeros_like(residuals)
    best_axes = np.empty(residuals.shape, dtype="U1")

    target = logical_braid_target(axis)

    best_val = 1e9
    best_meta = None

    for i, steps in enumerate(steps_list):
        for j, T in enumerate(T_list):
            R_total = compose_cycle(steps=steps, T=T, tc=tc, hc_factor=hc_factor)
            R_log = logical_project(R_total)
            res = phase_aligned_residual(R_log, target)
            residuals[i, j] = res
            braid_res[i, j] = braid_relation_residual(R_total)
            ybe_res[i, j] = ybe_residual(R_total)
            axis_vals = {a: phase_aligned_residual(R_log, logical_braid_target(a)) for a in ("x", "y", "z")}
            best_axes[i, j] = min(axis_vals, key=axis_vals.get)

            if res < best_val:
                best_val = res
                best_meta = dict(steps=steps, T=T, tc=tc, hc_factor=hc_factor, axis=axis, R_total=R_total, R_log=R_log)
            print(f"steps={steps}, T={T}: res={res:.3e}, braid={braid_res[i,j]:.3e}, ybe={ybe_res[i,j]:.3e}")

    stem = f"scan_steps_{'_'.join(map(str,steps_list))}_T_{'_'.join([s.replace('.', 'p') for s in map(str,T_list)])}_tc{tc}"
    npz = outdir / f"{stem}.npz"
    np.savez(npz, steps_list=steps_list, T_list=T_list, residuals=residuals, braid_res=braid_res, ybe_res=ybe_res, best_axes=best_axes)
    print(f"saved scan results to {npz}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(4.5, 3.5), dpi=150)
        im = ax.imshow(residuals, origin="lower", cmap="magma", aspect="auto")
        ax.set_xticks(range(len(T_list)))
        ax.set_xticklabels([str(t) for t in T_list])
        ax.set_yticks(range(len(steps_list)))
        ax.set_yticklabels([str(s) for s in steps_list])
        ax.set_xlabel("T")
        ax.set_ylabel("steps")
        ax.set_title("phase-aligned residual to target axis={}".format(axis))
        fig.colorbar(im, ax=ax, fraction=0.045)
        fig_path = outdir / f"{stem}.png"
        fig.tight_layout()
        fig.savefig(fig_path)
        print(f"saved heatmap to {fig_path}")
    except Exception as e:
        print("plotting failed:", e)

    if best_meta is not None:
        print("best:", best_meta["steps"], best_meta["T"], "res=", best_val)


if __name__ == "__main__":
    main()
