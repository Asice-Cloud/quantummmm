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
- optionally adds a Zeeman-like logical bias term to test whether a diagonal
    detuning improves the braid/target tradeoff
- optionally adds a segment-dependent phase compensation term to mimic
    dynamical-phase cancellation across the three path pieces
- optionally adds a small XX/YY cross-coupling term to test whether a more
    isotropic exchange-like perturbation improves the braid/target tradeoff

The path conventions match the formulas in ver5.md:
  g^(1)(u) = (u, 0, 1-u, 1)
  g^(2)(u) = (1-u, u, 0, 1)
  g^(3)(u) = (0, 1-u, u, 1)

The Hamiltonian is written in the Pauli basis as
  H_4(t) = t_c [ -g1 g2 ZI - g1 g3 YX + g1 g4 YY - g2 g3 XX - g2 g4 XY - g3 g4 IZ ].

If you want the convention with an explicit +h.c. factor, multiply the whole
Hamiltonian by 2 via --hc-factor 2.0. A Zeeman-like bias is modeled as
    H_Z = zeeman * (ZI - IZ),
which acts as a logical detuning between |01> and |10>.
The paper-motivated dynamic Zeeman compensation is modeled as
    H_Z^drive(k,u) = zeeman_drive_k * cos(pi u) * (ZI - IZ),
with independent amplitudes per path segment.
The dynamic Zeeman profile is selectable via --zeeman-drive-shape, with the
default cosine envelope matching the paper-style Vx(t)=Vx0+Vx1 cos(pi u).
The phase compensation term is modeled as a time-dependent ramp
    H_phi^(k)(u) = phase_comp * p_k * (2u - 1) * (ZI - IZ),
with p_k = (+1, -1, +1) for the three path segments.
The cross-coupling term can be assigned independently per segment as
    H_x^(k) = cross_comp_k * (XX + YY).
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
PHASE_PROFILE = (1.0, -1.0, 1.0)


def kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.kron(a, b)


def zeeman_drive_profile(u: float, shape: str) -> float:
    shape = shape.lower()
    if shape == "cos":
        return math.cos(math.pi * u)
    if shape == "sin":
        return math.sin(math.pi * u)
    if shape == "const":
        return 1.0
    if shape == "ramp":
        return 2.0 * u - 1.0
    raise ValueError(f"zeeman-drive shape must be one of cos/sin/const/ramp, got {shape}")


def h4_from_path(
    g1: float,
    g2: float,
    g3: float,
    g4: float,
    tc: float = 1.0,
    hc_factor: float = 1.0,
    zeeman: float = 0.0,
    zeeman_drive: float = 0.0,
    phase_comp: float = 0.0,
    phase_ramp: float = 0.0,
    cross_comp: float = 0.0,
    zeeman_u: float = 0.0,
    zeeman_drive_shape: str = "cos",
) -> np.ndarray:
    """Return the 4x4 Hamiltonian from the designed path."""
    h4 = (
        -(g1 * g2) * kron(SZ, I2)
        -(g1 * g3) * kron(SY, SX)
        +(g1 * g4) * kron(SY, SY)
        -(g2 * g3) * kron(SX, SX)
        -(g2 * g4) * kron(SX, SY)
        -(g3 * g4) * kron(I2, SZ)
    )
    zeeman_term = (zeeman + zeeman_drive * zeeman_drive_profile(zeeman_u, zeeman_drive_shape)) * (kron(SZ, I2) - kron(I2, SZ))
    phase_term = phase_comp * phase_ramp * (kron(SZ, I2) - kron(I2, SZ))
    cross_term = cross_comp * (kron(SX, SX) + kron(SY, SY))
    return hc_factor * tc * h4 + zeeman_term + phase_term + cross_term


def path_segment(k: int, u: float) -> Tuple[float, float, float, float]:
    """Return the segment-k control point at local parameter u in [0, 1]."""
    if k == 1:
        return u, 0.0, 1.0 - u, 1.0
    if k == 2:
        return 1.0 - u, u, 0.0, 1.0
    if k == 3:
        return 0.0, 1.0 - u, u, 1.0
    raise ValueError(f"segment index must be 1, 2, or 3; got {k}")


def integrate_segment(
    k: int,
    steps: int,
    T: float,
    tc: float,
    hc_factor: float,
    zeeman: float,
    zeeman_drive: float,
    zeeman_drive_shape: str,
    phase_comp: float,
    cross_comp: float,
) -> np.ndarray:
    """Midpoint product approximation to Texp[-i T ∫_0^1 H_k(u) du]."""
    dt = T / steps
    U = np.eye(4, dtype=complex)
    # Multiply in reverse time order so later times appear on the left.
    for u in ((np.arange(steps) + 0.5) / steps)[::-1]:
        g1, g2, g3, g4 = path_segment(k, float(u))
        H = h4_from_path(
            g1,
            g2,
            g3,
            g4,
            tc=tc,
            hc_factor=hc_factor,
            zeeman=zeeman,
            zeeman_drive=zeeman_drive,
            zeeman_drive_shape=zeeman_drive_shape,
            phase_comp=phase_comp,
            phase_ramp=PHASE_PROFILE[k - 1] * (2.0 * float(u) - 1.0),
            cross_comp=cross_comp,
            zeeman_u=float(u),
        )
        U = expm(-1j * dt * H) @ U
    return U


def compose_cycle(
    steps: int,
    T: float,
    tc: float,
    hc_factor: float,
    zeeman: float,
    zeeman_drive_1: float,
    zeeman_drive_2: float,
    zeeman_drive_3: float,
    zeeman_drive_shape: str,
    phase_comp: float,
    cross_comp_1: float,
    cross_comp_2: float,
    cross_comp_3: float,
) -> Dict[str, np.ndarray]:
    """Return the three segment gates and the six-step total gate."""
    R1 = integrate_segment(1, steps=steps, T=T, tc=tc, hc_factor=hc_factor, zeeman=zeeman, zeeman_drive=zeeman_drive_1, zeeman_drive_shape=zeeman_drive_shape, phase_comp=phase_comp, cross_comp=cross_comp_1)
    R2 = integrate_segment(2, steps=steps, T=T, tc=tc, hc_factor=hc_factor, zeeman=zeeman, zeeman_drive=zeeman_drive_2, zeeman_drive_shape=zeeman_drive_shape, phase_comp=phase_comp, cross_comp=cross_comp_2)
    R3 = integrate_segment(3, steps=steps, T=T, tc=tc, hc_factor=hc_factor, zeeman=zeeman, zeeman_drive=zeeman_drive_3, zeeman_drive_shape=zeeman_drive_shape, phase_comp=phase_comp, cross_comp=cross_comp_3)
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


def analyze_total_gate(R_total: np.ndarray, axis: str) -> Dict[str, object]:
    """Return logical projection diagnostics for a total 4x4 gate."""
    R_logical = logical_project(R_total)
    target = logical_braid_target(axis)
    target_residual = phase_aligned_residual(R_logical, target)
    best_axis, axis_residuals = best_logical_axis(R_logical)
    braid_res = braid_relation_residual(R_total)
    ybe_res = ybe_residual(R_total)
    score = target_residual + braid_res + ybe_res
    return {
        "R_logical": R_logical,
        "target": target,
        "target_residual": target_residual,
        "best_axis": best_axis,
        "axis_residuals": axis_residuals,
        "braid_res": braid_res,
        "ybe_res": ybe_res,
        "score": score,
    }


def scan_T_values(steps: int, tc: float, hc_factor: float, zeeman: float, zeeman_drive_1: float, zeeman_drive_2: float, zeeman_drive_3: float, zeeman_drive_shape: str, phase_comp: float, cross_comp_1: float, cross_comp_2: float, cross_comp_3: float, axis: str, T_values: np.ndarray) -> Dict[str, object]:
    """Scan a grid of T values and return the best candidate by total score."""
    records = []
    for T in T_values:
        gates = compose_cycle(
            steps=steps,
            T=float(T),
            tc=tc,
            hc_factor=hc_factor,
            zeeman=zeeman,
            zeeman_drive_1=zeeman_drive_1,
            zeeman_drive_2=zeeman_drive_2,
            zeeman_drive_3=zeeman_drive_3,
            zeeman_drive_shape=zeeman_drive_shape,
            phase_comp=phase_comp,
            cross_comp_1=cross_comp_1,
            cross_comp_2=cross_comp_2,
            cross_comp_3=cross_comp_3,
        )
        diag = analyze_total_gate(gates["R_total"], axis)
        diag.update({"T": float(T)})
        records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {
        "records": records,
        "best": best,
        "T_values": np.asarray(T_values, dtype=float),
    }


def scan_zeeman_values(steps: int, T: float, tc: float, hc_factor: float, axis: str, zeeman_values: np.ndarray, zeeman_drive_shape: str = "cos") -> Dict[str, object]:
    """Scan a grid of Zeeman values and return the best candidate by total score."""
    records = []
    for zeeman in zeeman_values:
        gates = compose_cycle(steps=steps, T=float(T), tc=tc, hc_factor=hc_factor, zeeman=float(zeeman), zeeman_drive_1=0.0, zeeman_drive_2=0.0, zeeman_drive_3=0.0, zeeman_drive_shape=zeeman_drive_shape, phase_comp=0.0, cross_comp_1=0.0, cross_comp_2=0.0, cross_comp_3=0.0)
        diag = analyze_total_gate(gates["R_total"], axis)
        diag.update({"zeeman": float(zeeman)})
        records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {
        "records": records,
        "best": best,
        "zeeman_values": np.asarray(zeeman_values, dtype=float),
    }


def scan_phase_values(steps: int, T: float, tc: float, hc_factor: float, axis: str, phase_values: np.ndarray, zeeman_drive_shape: str = "cos") -> Dict[str, object]:
    """Scan a grid of phase-compensation values and return the best candidate by total score."""
    records = []
    for phase_comp in phase_values:
        gates = compose_cycle(steps=steps, T=float(T), tc=tc, hc_factor=hc_factor, zeeman=0.0, zeeman_drive_1=0.0, zeeman_drive_2=0.0, zeeman_drive_3=0.0, zeeman_drive_shape=zeeman_drive_shape, phase_comp=float(phase_comp), cross_comp_1=0.0, cross_comp_2=0.0, cross_comp_3=0.0)
        diag = analyze_total_gate(gates["R_total"], axis)
        diag.update({"phase_comp": float(phase_comp)})
        records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {
        "records": records,
        "best": best,
        "phase_values": np.asarray(phase_values, dtype=float),
    }


def scan_cross_values(steps: int, T: float, tc: float, hc_factor: float, axis: str, cross_values: np.ndarray, zeeman_drive_shape: str = "cos") -> Dict[str, object]:
    """Scan a grid of XX/YY cross-coupling values and return the best candidate by total score."""
    records = []
    for cross_comp in cross_values:
        gates = compose_cycle(
            steps=steps,
            T=float(T),
            tc=tc,
            hc_factor=hc_factor,
            zeeman=0.0,
            zeeman_drive_1=0.0,
            zeeman_drive_2=0.0,
            zeeman_drive_3=0.0,
            zeeman_drive_shape=zeeman_drive_shape,
            phase_comp=0.0,
            cross_comp_1=float(cross_comp),
            cross_comp_2=float(cross_comp),
            cross_comp_3=float(cross_comp),
        )
        diag = analyze_total_gate(gates["R_total"], axis)
        diag.update({"cross_comp": float(cross_comp)})
        records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {
        "records": records,
        "best": best,
        "cross_values": np.asarray(cross_values, dtype=float),
    }


def scan_cross_segment_values(steps: int, T: float, tc: float, hc_factor: float, axis: str, cross_values: np.ndarray, segment: int, zeeman_drive_shape: str = "cos") -> Dict[str, object]:
    """Scan one segment's XX/YY cross-coupling while keeping the other two at zero."""
    records = []
    for cross_comp in cross_values:
        cross_profile = [0.0, 0.0, 0.0]
        cross_profile[segment - 1] = float(cross_comp)
        gates = compose_cycle(
            steps=steps,
            T=float(T),
            tc=tc,
            hc_factor=hc_factor,
            zeeman=0.0,
            zeeman_drive_1=0.0,
            zeeman_drive_2=0.0,
            zeeman_drive_3=0.0,
            zeeman_drive_shape=zeeman_drive_shape,
            phase_comp=0.0,
            cross_comp_1=cross_profile[0],
            cross_comp_2=cross_profile[1],
            cross_comp_3=cross_profile[2],
        )
        diag = analyze_total_gate(gates["R_total"], axis)
        diag.update({"cross_comp": float(cross_comp), "segment": int(segment)})
        records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {"records": records, "best": best, "cross_values": np.asarray(cross_values, dtype=float), "segment": int(segment)}


def scan_phase_cross_values(steps: int, T: float, tc: float, hc_factor: float, axis: str, phase_values: np.ndarray, cross_values: np.ndarray, zeeman_drive_shape: str = "cos") -> Dict[str, object]:
    """Scan a 2D grid of phase compensation and cross coupling values."""
    records = []
    for phase_comp in phase_values:
        for cross_comp in cross_values:
            gates = compose_cycle(
                steps=steps,
                T=float(T),
                tc=tc,
                hc_factor=hc_factor,
                zeeman=0.0,
                zeeman_drive_1=0.0,
                zeeman_drive_2=0.0,
                zeeman_drive_3=0.0,
                zeeman_drive_shape=zeeman_drive_shape,
                phase_comp=float(phase_comp),
                cross_comp_1=float(cross_comp),
                cross_comp_2=float(cross_comp),
                cross_comp_3=float(cross_comp),
            )
            diag = analyze_total_gate(gates["R_total"], axis)
            diag.update({"phase_comp": float(phase_comp), "cross_comp": float(cross_comp)})
            records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {
        "records": records,
        "best": best,
        "phase_values": np.asarray(phase_values, dtype=float),
        "cross_values": np.asarray(cross_values, dtype=float),
    }


def scan_phase_zeeman_values(
    steps: int,
    T: float,
    tc: float,
    hc_factor: float,
    axis: str,
    phase_values: np.ndarray,
    drive_values: np.ndarray,
    zeeman_drive_shape: str = "cos",
) -> Dict[str, object]:
    """Scan a 2D grid of phase compensation and (global) dynamic Zeeman-drive amplitude."""
    records = []
    for phase_comp in phase_values:
        for drive_amp in drive_values:
            gates = compose_cycle(
                steps=steps,
                T=float(T),
                tc=tc,
                hc_factor=hc_factor,
                zeeman=0.0,
                zeeman_drive_1=float(drive_amp),
                zeeman_drive_2=float(drive_amp),
                zeeman_drive_3=float(drive_amp),
                zeeman_drive_shape=zeeman_drive_shape,
                phase_comp=float(phase_comp),
                cross_comp_1=0.0,
                cross_comp_2=0.0,
                cross_comp_3=0.0,
            )
            diag = analyze_total_gate(gates["R_total"], axis)
            diag.update({"phase_comp": float(phase_comp), "zeeman_drive": float(drive_amp)})
            records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {
        "records": records,
        "best": best,
        "phase_values": np.asarray(phase_values, dtype=float),
        "drive_values": np.asarray(drive_values, dtype=float),
    }


def scan_phase_zeeman_segment_values(
    steps: int,
    T: float,
    tc: float,
    hc_factor: float,
    axis: str,
    phase_values: np.ndarray,
    drive_values: np.ndarray,
    segment: int,
    zeeman_drive_shape: str = "cos",
) -> Dict[str, object]:
    """Scan phase × dynamic Zeeman-drive amplitude for one specified segment."""
    records = []
    for phase_comp in phase_values:
        for drive_amp in drive_values:
            drives = [0.0, 0.0, 0.0]
            drives[segment - 1] = float(drive_amp)
            gates = compose_cycle(
                steps=steps,
                T=float(T),
                tc=tc,
                hc_factor=hc_factor,
                zeeman=0.0,
                zeeman_drive_1=drives[0],
                zeeman_drive_2=drives[1],
                zeeman_drive_3=drives[2],
                zeeman_drive_shape=zeeman_drive_shape,
                phase_comp=float(phase_comp),
                cross_comp_1=0.0,
                cross_comp_2=0.0,
                cross_comp_3=0.0,
            )
            diag = analyze_total_gate(gates["R_total"], axis)
            diag.update({"phase_comp": float(phase_comp), "zeeman_drive": float(drive_amp), "segment": int(segment)})
            records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {
        "records": records,
        "best": best,
        "phase_values": np.asarray(phase_values, dtype=float),
        "drive_values": np.asarray(drive_values, dtype=float),
        "segment": int(segment),
    }


def scan_phase_cross_segment_values(
    steps: int,
    T: float,
    tc: float,
    hc_factor: float,
    axis: str,
    phase_values: np.ndarray,
    cross_values: np.ndarray,
    segment: int,
    zeeman_drive_shape: str = "cos",
) -> Dict[str, object]:
    """Scan a 2D phase-cross grid for one segment while keeping the other two off."""
    records = []
    for phase_comp in phase_values:
        for cross_comp in cross_values:
            cross_profile = [0.0, 0.0, 0.0]
            cross_profile[segment - 1] = float(cross_comp)
            gates = compose_cycle(
                steps=steps,
                T=float(T),
                tc=tc,
                hc_factor=hc_factor,
                zeeman=0.0,
                zeeman_drive_1=0.0,
                zeeman_drive_2=0.0,
                zeeman_drive_3=0.0,
                zeeman_drive_shape=zeeman_drive_shape,
                phase_comp=float(phase_comp),
                cross_comp_1=cross_profile[0],
                cross_comp_2=cross_profile[1],
                cross_comp_3=cross_profile[2],
            )
            diag = analyze_total_gate(gates["R_total"], axis)
            diag.update({"phase_comp": float(phase_comp), "cross_comp": float(cross_comp), "segment": int(segment)})
            records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {
        "records": records,
        "best": best,
        "phase_values": np.asarray(phase_values, dtype=float),
        "cross_values": np.asarray(cross_values, dtype=float),
        "segment": int(segment),
    }


def scan_zeeman_drive_segment_values(steps: int, T: float, tc: float, hc_factor: float, axis: str, drive_values: np.ndarray, segment: int, zeeman_drive_shape: str = "cos") -> Dict[str, object]:
    """Scan one segment's dynamic Zeeman drive amplitude with other segments kept static."""
    records = []
    for drive_amp in drive_values:
        drives = [0.0, 0.0, 0.0]
        drives[segment - 1] = float(drive_amp)
        gates = compose_cycle(
            steps=steps,
            T=float(T),
            tc=tc,
            hc_factor=hc_factor,
            zeeman=0.0,
            zeeman_drive_1=drives[0],
            zeeman_drive_2=drives[1],
            zeeman_drive_3=drives[2],
            zeeman_drive_shape=zeeman_drive_shape,
            phase_comp=0.0,
            cross_comp_1=0.0,
            cross_comp_2=0.0,
            cross_comp_3=0.0,
        )
        diag = analyze_total_gate(gates["R_total"], axis)
        diag.update({"zeeman_drive": float(drive_amp), "segment": int(segment)})
        records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {"records": records, "best": best, "drive_values": np.asarray(drive_values, dtype=float), "segment": int(segment)}


def scan_zeeman_drive_values(steps: int, T: float, tc: float, hc_factor: float, axis: str, drive_values: np.ndarray, zeeman_drive_shape: str = "cos") -> Dict[str, object]:
    """Scan a global dynamic Zeeman drive with the same amplitude on all segments."""
    records = []
    for drive_amp in drive_values:
        gates = compose_cycle(
            steps=steps,
            T=float(T),
            tc=tc,
            hc_factor=hc_factor,
            zeeman=0.0,
            zeeman_drive_1=float(drive_amp),
            zeeman_drive_2=float(drive_amp),
            zeeman_drive_3=float(drive_amp),
            zeeman_drive_shape=zeeman_drive_shape,
            phase_comp=0.0,
            cross_comp_1=0.0,
            cross_comp_2=0.0,
            cross_comp_3=0.0,
        )
        diag = analyze_total_gate(gates["R_total"], axis)
        diag.update({"zeeman_drive": float(drive_amp)})
        records.append(diag)

    best_idx = min(range(len(records)), key=lambda i: records[i]["score"])
    best = records[best_idx]
    return {"records": records, "best": best, "drive_values": np.asarray(drive_values, dtype=float)}


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


def default_stem(steps: int, T: float, tc: float, axis: str, hc_factor: float, zeeman: float, cross_comp: float) -> str:
    def slug(value: float) -> str:
        return format(value, "g").replace("-", "m").replace(".", "p")

    return f"path_braid_steps{steps}_T{slug(T)}_tc{slug(tc)}_axis{axis}_hc{slug(hc_factor)}_ze{slug(zeeman)}_cr{slug(cross_comp)}"


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
    p.add_argument("--zeeman", type=float, default=0.0, help="Zeeman-like bias coefficient for ZI-IZ")
    p.add_argument("--zeeman-drive-1", type=float, default=0.0, help="segment-1 dynamic Zeeman drive amplitude")
    p.add_argument("--zeeman-drive-2", type=float, default=0.0, help="segment-2 dynamic Zeeman drive amplitude")
    p.add_argument("--zeeman-drive-3", type=float, default=0.0, help="segment-3 dynamic Zeeman drive amplitude")
    p.add_argument("--zeeman-drive-shape", choices=("cos", "sin", "const", "ramp"), default="cos", help="shape selector for dynamic Zeeman drive envelope")
    p.add_argument("--phase-comp", type=float, default=0.0, help="segment-dependent phase compensation coefficient for ZI-IZ")
    p.add_argument("--cross-comp", type=float, default=0.0, help="XX/YY cross-coupling coefficient")
    p.add_argument("--cross-comp-1", type=float, default=0.0, help="segment-1 XX/YY cross-coupling coefficient")
    p.add_argument("--cross-comp-2", type=float, default=0.0, help="segment-2 XX/YY cross-coupling coefficient")
    p.add_argument("--cross-comp-3", type=float, default=0.0, help="segment-3 XX/YY cross-coupling coefficient")
    p.add_argument("--axis", choices=("x", "y", "z"), default="x", help="logical braid axis to compare against")
    p.add_argument("--scan-T", action="store_true", help="scan a T grid and report the best candidate")
    p.add_argument("--T-min", type=float, default=0.25, help="minimum T value when scanning")
    p.add_argument("--T-max", type=float, default=4.0, help="maximum T value when scanning")
    p.add_argument("--T-num", type=int, default=25, help="number of T samples when scanning")
    p.add_argument("--scan-zeeman", action="store_true", help="scan a Zeeman grid and report the best candidate")
    p.add_argument("--zeeman-min", type=float, default=-0.5, help="minimum Zeeman value when scanning")
    p.add_argument("--zeeman-max", type=float, default=0.5, help="maximum Zeeman value when scanning")
    p.add_argument("--zeeman-num", type=int, default=25, help="number of Zeeman samples when scanning")
    p.add_argument("--scan-zeeman-drive", action="store_true", help="scan a global dynamic Zeeman-drive grid")
    p.add_argument("--scan-zeeman-drive-segments", action="store_true", help="scan each segment's dynamic Zeeman-drive amplitude independently")
    p.add_argument("--zeeman-drive-min", type=float, default=-0.5, help="minimum dynamic Zeeman-drive value when scanning")
    p.add_argument("--zeeman-drive-max", type=float, default=0.5, help="maximum dynamic Zeeman-drive value when scanning")
    p.add_argument("--zeeman-drive-num", type=int, default=25, help="number of dynamic Zeeman-drive samples when scanning")
    p.add_argument("--scan-phase", action="store_true", help="scan a phase-compensation grid and report the best candidate")
    p.add_argument("--phase-min", type=float, default=-0.5, help="minimum phase compensation value when scanning")
    p.add_argument("--phase-max", type=float, default=0.5, help="maximum phase compensation value when scanning")
    p.add_argument("--phase-num", type=int, default=25, help="number of phase-compensation samples when scanning")
    p.add_argument("--scan-cross", action="store_true", help="scan an XX/YY cross-coupling grid and report the best candidate")
    p.add_argument("--scan-cross-segments", action="store_true", help="scan each segment's XX/YY cross-coupling independently")
    p.add_argument("--scan-phase-cross", action="store_true", help="scan a 2D phase-compensation and cross-coupling grid")
    p.add_argument("--scan-phase-cross-segments", action="store_true", help="scan a 2D phase-cross grid independently for each segment")
    p.add_argument("--cross-min", type=float, default=-0.5, help="minimum cross-coupling value when scanning")
    p.add_argument("--cross-max", type=float, default=0.5, help="maximum cross-coupling value when scanning")
    p.add_argument("--cross-num", type=int, default=25, help="number of cross-coupling samples when scanning")
    p.add_argument("--outdir", type=str, default="results/path_braid", help="directory for saved diagnostics")
    p.add_argument("--stem", type=str, default="", help="optional file stem override")
    p.add_argument("--save", type=str, default="", help="optional .npz output path override")
    p.add_argument("--fig", type=str, default="", help="optional .png figure path override")
    p.add_argument("--scan-phase-vx", action="store_true", help="scan a 2D grid of phase-compensation and dynamic Zeeman-drive amplitude")
    p.add_argument("--scan-phase-vx-seg2", action="store_true", help="scan phase × dynamic Zeeman-drive amplitude only on segment-2")
    p.add_argument("--scan-phase-vx-all", action="store_true", help="scan phase × dynamic Zeeman-drive with all segments using the same amplitude")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.scan_T:
        T_values = np.linspace(args.T_min, args.T_max, args.T_num)
        scan = scan_T_values(
            steps=args.steps,
            tc=args.tc,
            hc_factor=args.hc_factor,
            zeeman=args.zeeman,
            zeeman_drive_1=args.zeeman_drive_1,
            zeeman_drive_2=args.zeeman_drive_2,
            zeeman_drive_3=args.zeeman_drive_3,
            zeeman_drive_shape=args.zeeman_drive_shape,
            phase_comp=args.phase_comp,
            cross_comp_1=args.cross_comp_1,
            cross_comp_2=args.cross_comp_2,
            cross_comp_3=args.cross_comp_3,
            axis=args.axis,
            T_values=T_values,
        )
        best = scan["best"]
        print("T-scan enabled")
        print("candidate Ts ->", np.array2string(scan["T_values"], precision=4, separator=", "))
        print(f"best T = {best['T']:.6g}")
        print(f"best score = {best['score']:.3e}")
        print(f"best target residual = {best['target_residual']:.3e}")
        print(f"best braid residual = {best['braid_res']:.3e}")
        print(f"best YBE residual = {best['ybe_res']:.3e}")
        print(f"best axis = {best['best_axis']}")

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_Tscan.npz"
        np.savez(
            scan_path,
            T_values=scan["T_values"],
            best_T=np.array(best["T"]),
            best_score=np.array(best["score"]),
            best_target_residual=np.array(best["target_residual"]),
            best_braid_residual=np.array(best["braid_res"]),
            best_ybe_residual=np.array(best["ybe_res"]),
            best_axis=np.array(best["best_axis"]),
        )
        print(f"saved T-scan diagnostics to {scan_path}")
        return

    if args.scan_cross:
        cross_values = np.linspace(args.cross_min, args.cross_max, args.cross_num)
        scan = scan_cross_values(steps=args.steps, T=args.T, tc=args.tc, hc_factor=args.hc_factor, axis=args.axis, cross_values=cross_values, zeeman_drive_shape=args.zeeman_drive_shape)
        best = scan["best"]
        print("XX/YY cross-coupling scan enabled")
        print("candidate cross-comp values ->", np.array2string(scan["cross_values"], precision=4, separator=", "))
        print(f"best cross_comp = {best['cross_comp']:.6g}")
        print(f"best score = {best['score']:.3e}")
        print(f"best target residual = {best['target_residual']:.3e}")
        print(f"best braid residual = {best['braid_res']:.3e}")
        print(f"best YBE residual = {best['ybe_res']:.3e}")
        print(f"best axis = {best['best_axis']}")

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_Cscan.npz"
        np.savez(
            scan_path,
            cross_values=scan["cross_values"],
            best_cross_comp=np.array(best["cross_comp"]),
            best_score=np.array(best["score"]),
            best_target_residual=np.array(best["target_residual"]),
            best_braid_residual=np.array(best["braid_res"]),
            best_ybe_residual=np.array(best["ybe_res"]),
            best_axis=np.array(best["best_axis"]),
        )
        print(f"saved cross-coupling diagnostics to {scan_path}")
        return

    if args.scan_cross_segments:
        cross_values = np.linspace(args.cross_min, args.cross_max, args.cross_num)
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for segment in (1, 2, 3):
            scan = scan_cross_segment_values(steps=args.steps, T=args.T, tc=args.tc, hc_factor=args.hc_factor, axis=args.axis, cross_values=cross_values, segment=segment, zeeman_drive_shape=args.zeeman_drive_shape)
            best = scan["best"]
            print(f"segment-{segment} cross-coupling scan enabled")
            print("candidate cross-comp values ->", np.array2string(scan["cross_values"], precision=4, separator=", "))
            print(f"best cross_comp[{segment}] = {best['cross_comp']:.6g}")
            print(f"best score = {best['score']:.3e}")
            print(f"best target residual = {best['target_residual']:.3e}")
            print(f"best braid residual = {best['braid_res']:.3e}")
            print(f"best YBE residual = {best['ybe_res']:.3e}")
            print(f"best axis = {best['best_axis']}")

            scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_Cscan_seg{segment}.npz"
            np.savez(
                scan_path,
                cross_values=scan["cross_values"],
                best_cross_comp=np.array(best["cross_comp"]),
                best_score=np.array(best["score"]),
                best_target_residual=np.array(best["target_residual"]),
                best_braid_residual=np.array(best["braid_res"]),
                best_ybe_residual=np.array(best["ybe_res"]),
                best_axis=np.array(best["best_axis"]),
                segment=np.array(segment),
            )
            print(f"saved segment-{segment} cross-coupling diagnostics to {scan_path}")
        return

    if args.scan_phase_cross_segments:
        phase_values = np.linspace(args.phase_min, args.phase_max, args.phase_num)
        cross_values = np.linspace(args.cross_min, args.cross_max, args.cross_num)
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for segment in (1, 2, 3):
            scan = scan_phase_cross_segment_values(
                steps=args.steps,
                T=args.T,
                tc=args.tc,
                hc_factor=args.hc_factor,
                axis=args.axis,
                phase_values=phase_values,
                cross_values=cross_values,
                segment=segment,
                zeeman_drive_shape=args.zeeman_drive_shape,
            )
            best = scan["best"]
            print(f"segment-{segment} phase × cross joint scan enabled")
            print("candidate phase values ->", np.array2string(scan["phase_values"], precision=4, separator=", "))
            print("candidate cross values ->", np.array2string(scan["cross_values"], precision=4, separator=", "))
            print(f"best phase_comp = {best['phase_comp']:.6g}")
            print(f"best cross_comp[{segment}] = {best['cross_comp']:.6g}")
            print(f"best score = {best['score']:.3e}")
            print(f"best target residual = {best['target_residual']:.3e}")
            print(f"best braid residual = {best['braid_res']:.3e}")
            print(f"best YBE residual = {best['ybe_res']:.3e}")
            print(f"best axis = {best['best_axis']}")

            scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_PXscan_seg{segment}.npz"
            np.savez(
                scan_path,
                phase_values=scan["phase_values"],
                cross_values=scan["cross_values"],
                best_phase_comp=np.array(best["phase_comp"]),
                best_cross_comp=np.array(best["cross_comp"]),
                best_score=np.array(best["score"]),
                best_target_residual=np.array(best["target_residual"]),
                best_braid_residual=np.array(best["braid_res"]),
                best_ybe_residual=np.array(best["ybe_res"]),
                best_axis=np.array(best["best_axis"]),
                segment=np.array(segment),
            )
            print(f"saved segment-{segment} phase-cross diagnostics to {scan_path}")
        return

    if args.scan_zeeman_drive_segments:
        drive_values = np.linspace(args.zeeman_drive_min, args.zeeman_drive_max, args.zeeman_drive_num)
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for segment in (1, 2, 3):
            scan = scan_zeeman_drive_segment_values(steps=args.steps, T=args.T, tc=args.tc, hc_factor=args.hc_factor, axis=args.axis, drive_values=drive_values, segment=segment, zeeman_drive_shape=args.zeeman_drive_shape)
            best = scan["best"]
            print(f"segment-{segment} dynamic Zeeman scan enabled")
            print("candidate drive values ->", np.array2string(scan["drive_values"], precision=4, separator=", "))
            print(f"best zeeman_drive[{segment}] = {best['zeeman_drive']:.6g}")
            print(f"best score = {best['score']:.3e}")
            print(f"best target residual = {best['target_residual']:.3e}")
            print(f"best braid residual = {best['braid_res']:.3e}")
            print(f"best YBE residual = {best['ybe_res']:.3e}")
            print(f"best axis = {best['best_axis']}")

            scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_Vscan_seg{segment}.npz"
            np.savez(
                scan_path,
                drive_values=scan["drive_values"],
                best_zeeman_drive=np.array(best["zeeman_drive"]),
                best_score=np.array(best["score"]),
                best_target_residual=np.array(best["target_residual"]),
                best_braid_residual=np.array(best["braid_res"]),
                best_ybe_residual=np.array(best["ybe_res"]),
                best_axis=np.array(best["best_axis"]),
                segment=np.array(segment),
            )
            print(f"saved segment-{segment} dynamic Zeeman diagnostics to {scan_path}")
        return

    if args.scan_zeeman_drive:
        drive_values = np.linspace(args.zeeman_drive_min, args.zeeman_drive_max, args.zeeman_drive_num)
        scan = scan_zeeman_drive_values(steps=args.steps, T=args.T, tc=args.tc, hc_factor=args.hc_factor, axis=args.axis, drive_values=drive_values, zeeman_drive_shape=args.zeeman_drive_shape)
        best = scan["best"]
        print("dynamic Zeeman-drive scan enabled")
        print("candidate drive values ->", np.array2string(scan["drive_values"], precision=4, separator=", "))
        print(f"best zeeman_drive = {best['zeeman_drive']:.6g}")
        print(f"best score = {best['score']:.3e}")
        print(f"best target residual = {best['target_residual']:.3e}")
        print(f"best braid residual = {best['braid_res']:.3e}")
        print(f"best YBE residual = {best['ybe_res']:.3e}")
        print(f"best axis = {best['best_axis']}")

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_Vscan.npz"
        np.savez(
            scan_path,
            drive_values=scan["drive_values"],
            best_zeeman_drive=np.array(best["zeeman_drive"]),
            best_score=np.array(best["score"]),
            best_target_residual=np.array(best["target_residual"]),
            best_braid_residual=np.array(best["braid_res"]),
            best_ybe_residual=np.array(best["ybe_res"]),
            best_axis=np.array(best["best_axis"]),
        )
        print(f"saved dynamic Zeeman diagnostics to {scan_path}")
        return

    if args.scan_phase_cross:
        phase_values = np.linspace(args.phase_min, args.phase_max, args.phase_num)
        cross_values = np.linspace(args.cross_min, args.cross_max, args.cross_num)
        scan = scan_phase_cross_values(steps=args.steps, T=args.T, tc=args.tc, hc_factor=args.hc_factor, axis=args.axis, phase_values=phase_values, cross_values=cross_values, zeeman_drive_shape=args.zeeman_drive_shape)
        best = scan["best"]
        print("phase × cross joint scan enabled")
        print("candidate phase values ->", np.array2string(scan["phase_values"], precision=4, separator=", "))
        print("candidate cross values ->", np.array2string(scan["cross_values"], precision=4, separator=", "))
        print(f"best phase_comp = {best['phase_comp']:.6g}")
        print(f"best cross_comp = {best['cross_comp']:.6g}")
        print(f"best score = {best['score']:.3e}")
        print(f"best target residual = {best['target_residual']:.3e}")
        print(f"best braid residual = {best['braid_res']:.3e}")
        print(f"best YBE residual = {best['ybe_res']:.3e}")
        print(f"best axis = {best['best_axis']}")

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_PXscan.npz"
        np.savez(
            scan_path,
            phase_values=scan["phase_values"],
            cross_values=scan["cross_values"],
            best_phase_comp=np.array(best["phase_comp"]),
            best_cross_comp=np.array(best["cross_comp"]),
            best_score=np.array(best["score"]),
            best_target_residual=np.array(best["target_residual"]),
            best_braid_residual=np.array(best["braid_res"]),
            best_ybe_residual=np.array(best["ybe_res"]),
            best_axis=np.array(best["best_axis"]),
        )
        print(f"saved phase-cross diagnostics to {scan_path}")
        return

    if args.scan_zeeman:
        zeeman_values = np.linspace(args.zeeman_min, args.zeeman_max, args.zeeman_num)
        scan = scan_zeeman_values(steps=args.steps, T=args.T, tc=args.tc, hc_factor=args.hc_factor, axis=args.axis, zeeman_values=zeeman_values, zeeman_drive_shape=args.zeeman_drive_shape)
        best = scan["best"]
        print("Zeeman-scan enabled")
        print("candidate Zeeman values ->", np.array2string(scan["zeeman_values"], precision=4, separator=", "))
        print(f"best zeeman = {best['zeeman']:.6g}")
        print(f"best score = {best['score']:.3e}")
        print(f"best target residual = {best['target_residual']:.3e}")
        print(f"best braid residual = {best['braid_res']:.3e}")
        print(f"best YBE residual = {best['ybe_res']:.3e}")
        print(f"best axis = {best['best_axis']}")

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_Zscan.npz"
        np.savez(
            scan_path,
            zeeman_values=scan["zeeman_values"],
            best_zeeman=np.array(best["zeeman"]),
            best_score=np.array(best["score"]),
            best_target_residual=np.array(best["target_residual"]),
            best_braid_residual=np.array(best["braid_res"]),
            best_ybe_residual=np.array(best["ybe_res"]),
            best_axis=np.array(best["best_axis"]),
        )
        print(f"saved Zeeman-scan diagnostics to {scan_path}")
        return

    if args.scan_phase:
        phase_values = np.linspace(args.phase_min, args.phase_max, args.phase_num)
        scan = scan_phase_values(steps=args.steps, T=args.T, tc=args.tc, hc_factor=args.hc_factor, axis=args.axis, phase_values=phase_values, zeeman_drive_shape=args.zeeman_drive_shape)
        best = scan["best"]
        print("phase-compensation scan enabled")
        print("candidate phase-comp values ->", np.array2string(scan["phase_values"], precision=4, separator=", "))
        print(f"best phase_comp = {best['phase_comp']:.6g}")
        print(f"best score = {best['score']:.3e}")
        print(f"best target residual = {best['target_residual']:.3e}")
        print(f"best braid residual = {best['braid_res']:.3e}")
        print(f"best YBE residual = {best['ybe_res']:.3e}")
        print(f"best axis = {best['best_axis']}")

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_Pscan.npz"
        np.savez(
            scan_path,
            phase_values=scan["phase_values"],
            best_phase_comp=np.array(best["phase_comp"]),
            best_score=np.array(best["score"]),
            best_target_residual=np.array(best["target_residual"]),
            best_braid_residual=np.array(best["braid_res"]),
            best_ybe_residual=np.array(best["ybe_res"]),
            best_axis=np.array(best["best_axis"]),
        )
        print(f"saved phase-compensation diagnostics to {scan_path}")
        return

    if args.scan_phase_vx:
        phase_values = np.linspace(args.phase_min, args.phase_max, args.phase_num)
        drive_values = np.linspace(args.zeeman_drive_min, args.zeeman_drive_max, args.zeeman_drive_num)
        scan = scan_phase_zeeman_values(
            steps=args.steps,
            T=args.T,
            tc=args.tc,
            hc_factor=args.hc_factor,
            axis=args.axis,
            phase_values=phase_values,
            drive_values=drive_values,
            zeeman_drive_shape=args.zeeman_drive_shape,
        )
        best = scan["best"]
        print("phase × dynamic-Zeeman joint scan enabled")
        print("candidate phase values ->", np.array2string(scan["phase_values"], precision=4, separator=", "))
        print("candidate drive values ->", np.array2string(scan["drive_values"], precision=4, separator=", "))
        print(f"best phase_comp = {best['phase_comp']:.6g}")
        print(f"best zeeman_drive = {best['zeeman_drive']:.6g}")
        print(f"best score = {best['score']:.3e}")
        print(f"best target residual = {best['target_residual']:.3e}")
        print(f"best braid residual = {best['braid_res']:.3e}")
        print(f"best YBE residual = {best['ybe_res']:.3e}")
        print(f"best axis = {best['best_axis']}")

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_PVscan.npz"
        np.savez(
            scan_path,
            phase_values=scan["phase_values"],
            drive_values=scan["drive_values"],
            best_phase_comp=np.array(best["phase_comp"]),
            best_zeeman_drive=np.array(best["zeeman_drive"]),
            best_score=np.array(best["score"]),
            best_target_residual=np.array(best["target_residual"]),
            best_braid_residual=np.array(best["braid_res"]),
            best_ybe_residual=np.array(best["ybe_res"]),
            best_axis=np.array(best["best_axis"]),
        )
        print(f"saved phase-Vx diagnostics to {scan_path}")
        return

    if args.scan_phase_vx_seg2:
        phase_values = np.linspace(args.phase_min, args.phase_max, args.phase_num)
        drive_values = np.linspace(args.zeeman_drive_min, args.zeeman_drive_max, args.zeeman_drive_num)
        scan = scan_phase_zeeman_segment_values(
            steps=args.steps,
            T=args.T,
            tc=args.tc,
            hc_factor=args.hc_factor,
            axis=args.axis,
            phase_values=phase_values,
            drive_values=drive_values,
            segment=2,
            zeeman_drive_shape=args.zeeman_drive_shape,
        )
        best = scan["best"]
        print("segment-2 phase × dynamic-Zeeman joint scan enabled")
        print("candidate phase values ->", np.array2string(scan["phase_values"], precision=4, separator=", "))
        print("candidate drive values ->", np.array2string(scan["drive_values"], precision=4, separator=", "))
        print(f"best phase_comp = {best['phase_comp']:.6g}")
        print(f"best zeeman_drive[2] = {best['zeeman_drive']:.6g}")
        print(f"best score = {best['score']:.3e}")
        print(f"best target residual = {best['target_residual']:.3e}")
        print(f"best braid residual = {best['braid_res']:.3e}")
        print(f"best YBE residual = {best['ybe_res']:.3e}")
        print(f"best axis = {best['best_axis']}")

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_PV_seg2_scan.npz"
        np.savez(
            scan_path,
            phase_values=scan["phase_values"],
            drive_values=scan["drive_values"],
            best_phase_comp=np.array(best["phase_comp"]),
            best_zeeman_drive=np.array(best["zeeman_drive"]),
            best_score=np.array(best["score"]),
            best_target_residual=np.array(best["target_residual"]),
            best_braid_residual=np.array(best["braid_res"]),
            best_ybe_residual=np.array(best["ybe_res"]),
            best_axis=np.array(best["best_axis"]),
            segment=np.array(2),
        )
        print(f"saved segment-2 phase-Vx diagnostics to {scan_path}")
        return

    if args.scan_phase_vx_all:
        phase_values = np.linspace(args.phase_min, args.phase_max, args.phase_num)
        drive_values = np.linspace(args.zeeman_drive_min, args.zeeman_drive_max, args.zeeman_drive_num)
        scan = scan_phase_zeeman_values(
            steps=args.steps,
            T=args.T,
            tc=args.tc,
            hc_factor=args.hc_factor,
            axis=args.axis,
            phase_values=phase_values,
            drive_values=drive_values,
            zeeman_drive_shape=args.zeeman_drive_shape,
        )
        best = scan["best"]
        print("multi-segment phase × dynamic-Zeeman joint scan enabled (all segments unified)")
        print("candidate phase values ->", np.array2string(scan["phase_values"], precision=4, separator=", "))
        print("candidate drive values ->", np.array2string(scan["drive_values"], precision=4, separator=", "))
        print(f"best phase_comp = {best['phase_comp']:.6g}")
        print(f"best zeeman_drive (all segments) = {best['zeeman_drive']:.6g}")
        print(f"best score = {best['score']:.3e}")
        print(f"best target residual = {best['target_residual']:.3e}")
        print(f"best braid residual = {best['braid_res']:.3e}")
        print(f"best YBE residual = {best['ybe_res']:.3e}")
        print(f"best axis = {best['best_axis']}")

        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        scan_path = outdir / f"{default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)}_PV_all_scan.npz"
        np.savez(
            scan_path,
            phase_values=scan["phase_values"],
            drive_values=scan["drive_values"],
            best_phase_comp=np.array(best["phase_comp"]),
            best_zeeman_drive=np.array(best["zeeman_drive"]),
            best_score=np.array(best["score"]),
            best_target_residual=np.array(best["target_residual"]),
            best_braid_residual=np.array(best["braid_res"]),
            best_ybe_residual=np.array(best["ybe_res"]),
            best_axis=np.array(best["best_axis"]),
        )
        print(f"saved multi-segment phase-Vx diagnostics to {scan_path}")
        return

    print("Building path-derived gates from ver5.md path")
    gates = compose_cycle(
        steps=args.steps,
        T=args.T,
        tc=args.tc,
        hc_factor=args.hc_factor,
        zeeman=args.zeeman,
        zeeman_drive_1=args.zeeman_drive_1,
        zeeman_drive_2=args.zeeman_drive_2,
        zeeman_drive_3=args.zeeman_drive_3,
        zeeman_drive_shape=args.zeeman_drive_shape,
        phase_comp=args.phase_comp,
        cross_comp_1=args.cross_comp_1,
        cross_comp_2=args.cross_comp_2,
        cross_comp_3=args.cross_comp_3,
    )
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

    diag = analyze_total_gate(R_total, args.axis)
    R_logical = diag["R_logical"]
    target = diag["target"]
    target_residual = diag["target_residual"]
    best_axis = diag["best_axis"]
    axis_residuals = diag["axis_residuals"]
    braid_res = diag["braid_res"]
    ybe_res = diag["ybe_res"]

    print(f"logical gate ||R_L^†R_L-I||_F = {unitary_error(R_logical):.3e}")

    print(f"logical target axis = {args.axis}")
    print(f"phase-aligned residual to exp(-i pi/4 sigma_{args.axis}) = {target_residual:.3e}")
    print("residual scan over x/y/z ->", {k: f"{v:.3e}" for k, v in axis_residuals.items()})
    print(f"best matching logical axis = {best_axis}")

    print(f"braid relation residual ||R12 R23 R12 - R23 R12 R23||_F = {braid_res:.3e}")
    print(f"constant YBE residual ||R12 R13 R23 - R23 R13 R12||_F = {ybe_res:.3e}")

    stem = args.stem.strip() or default_stem(args.steps, args.T, args.tc, args.axis, args.hc_factor, args.zeeman, args.cross_comp_1)
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