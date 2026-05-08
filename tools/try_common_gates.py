#!/usr/bin/env python3
"""Fit and verify common single-qubit gates with the x-y-x control path.

This script uses the constructive SU(2) control model from `ver2.md` and
`tools/verify_su2_x_y_x.py`.

For each target gate, it searches for angles (alpha, beta, gamma) such that

    U_path = exp(-i H_x T3) exp(-i H_y T2) exp(-i H_x T1)

with H_x = J sigma_x and H_y = J sigma_y matches the gate up to a global
phase, where T1 = alpha/(2J), T2 = beta/(2J), T3 = gamma/(2J).

The gates tested by default are: I, X, Y, Z, H, S, T.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize

from verify_su2_x_y_x import I2, SX, SY, SZ, best_global_phase, path_gate, lie_closure_rank


@dataclass(frozen=True)
class GateResult:
    name: str
    alpha: float
    beta: float
    gamma: float
    J: float
    T1: float
    T2: float
    T3: float
    phase: float
    error: float
    unitarity_error: float
    success: bool
    nit: int


def rotation(axis: str, theta: float) -> np.ndarray:
    sigma = {"x": SX, "y": SY, "z": SZ}[axis]
    return np.cos(theta / 2.0) * I2 - 1j * np.sin(theta / 2.0) * sigma


def gate_matrix(name: str) -> np.ndarray:
    """Return a canonical matrix for a common single-qubit gate."""
    name = name.upper()
    if name == "I":
        return I2
    if name == "X":
        return SX
    if name == "Y":
        return SY
    if name == "Z":
        return SZ
    if name == "H":
        return np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2.0)
    if name == "S":
        return np.array([[1, 0], [0, 1j]], dtype=complex)
    if name == "T":
        return np.array([[1, 0], [0, np.exp(1j * np.pi / 4.0)]], dtype=complex)
    raise ValueError(f"Unknown gate name: {name}")


def default_gates() -> List[str]:
    return ["I", "X", "Y", "Z", "H", "S", "T"]


def target_error(target: np.ndarray, alpha: float, beta: float, gamma: float, J: float) -> float:
    path, _ = path_gate(alpha, beta, gamma, J)
    phase = best_global_phase(path, target)
    aligned = np.exp(-1j * phase) * path
    return float(np.linalg.norm(aligned - target))


def optimize_gate(target: np.ndarray, J: float, x0: Sequence[float]) -> Tuple[np.ndarray, float, int, bool]:
    """Minimize the phase-aligned Frobenius error over (alpha, beta, gamma)."""

    def objective(x: np.ndarray) -> float:
        return target_error(target, float(x[0]), float(x[1]), float(x[2]), J)

    res = minimize(
        objective,
        np.asarray(x0, dtype=float),
        method="Nelder-Mead",
        options={"xatol": 1e-12, "fatol": 1e-12, "maxiter": 5000},
    )
    return np.asarray(res.x, dtype=float), float(res.fun), int(res.nit), bool(res.success)


def candidate_starts(name: str) -> List[Tuple[float, float, float]]:
    """A small multi-start set for the optimizer."""
    pi = np.pi
    starts = [
        (0.0, 0.0, 0.0),
        (pi, 0.0, 0.0),
        (0.0, pi, 0.0),
        (0.0, 0.0, pi),
        (pi / 2.0, pi / 2.0, pi / 2.0),
        (-pi / 2.0, pi / 2.0, pi / 2.0),
        (pi / 2.0, pi, -pi / 2.0),
        (pi, pi, pi),
        (pi / 2.0, pi / 4.0, pi / 2.0),
        (-pi / 2.0, pi / 4.0, pi / 2.0),
    ]

    if name == "H":
        starts.extend([
            (pi / 2.0, pi / 2.0, pi / 2.0),
            (pi / 2.0, pi, 0.0),
            (0.0, pi / 2.0, pi),
        ])
    elif name == "S":
        starts.extend([
            (-pi / 2.0, pi / 2.0, pi / 2.0),
            (pi / 2.0, pi / 2.0, -pi / 2.0),
        ])
    elif name == "T":
        starts.extend([
            (-pi / 2.0, pi / 4.0, pi / 2.0),
            (pi / 2.0, pi / 4.0, -pi / 2.0),
        ])
    return starts


def fit_gate(name: str, J: float) -> GateResult:
    target = gate_matrix(name)

    best_x = None
    best_fun = np.inf
    best_nit = 0
    best_success = False

    for x0 in candidate_starts(name):
        x, fun, nit, success = optimize_gate(target, J, x0)
        if fun < best_fun:
            best_x = x
            best_fun = fun
            best_nit = nit
            best_success = success

    assert best_x is not None
    # Normalize angles to [0, 2π) so the corresponding segment times are non-negative.
    alpha, beta, gamma = (float(np.mod(best_x[0], 2.0 * np.pi)),
                          float(np.mod(best_x[1], 2.0 * np.pi)),
                          float(np.mod(best_x[2], 2.0 * np.pi)))
    path, (T1, T2, T3) = path_gate(alpha, beta, gamma, J)
    phase = best_global_phase(path, target)
    aligned = np.exp(-1j * phase) * path
    error = float(np.linalg.norm(aligned - target))
    unitarity_error = float(np.linalg.norm(path.conj().T @ path - I2))

    return GateResult(
        name=name,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        J=J,
        T1=T1,
        T2=T2,
        T3=T3,
        phase=phase,
        error=error,
        unitarity_error=unitarity_error,
        success=best_success,
        nit=best_nit,
    )


def print_result(result: GateResult) -> None:
    print(f"Gate: {result.name}")
    print("  alpha,beta,gamma =", (result.alpha, result.beta, result.gamma))
    print("  T1,T2,T3 =", (result.T1, result.T2, result.T3))
    print("  best global phase phi =", result.phase)
    print("  Frobenius error after phase alignment =", result.error)
    print("  unitarity error =", result.unitarity_error)
    print("  optimizer success =", result.success, "nit =", result.nit)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit common single-qubit gates with the x-y-x control path.")
    parser.add_argument("--gate", type=str, default="ALL", help="Gate name or ALL.")
    parser.add_argument("--J", type=float, default=2.3, help="Control strength.")
    args = parser.parse_args()

    print("Lie closure rank from {iσx, iσy, i[σx,σy]} =", lie_closure_rank())
    names = default_gates() if args.gate.upper() == "ALL" else [args.gate.upper()]

    results = [fit_gate(name, args.J) for name in names]
    for result in results:
        print_result(result)


if __name__ == "__main__":
    main()