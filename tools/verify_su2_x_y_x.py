#!/usr/bin/env python3
"""Numerically verify an x-y-x control path for arbitrary SU(2) rotations.

This script implements the constructive proof used in `ver2.md`:

    U_path = exp(-i H_x T3) exp(-i H_y T2) exp(-i H_x T1)

with H_x = J sigma_x and H_y = J sigma_y.

Given Euler angles (alpha, beta, gamma), we set

    T1 = alpha / (2J),  T2 = beta / (2J),  T3 = gamma / (2J)

and verify that the resulting unitary matches the target gate

    U_target = e^{-i gamma sigma_x / 2} e^{-i beta sigma_y / 2} e^{-i alpha sigma_x / 2}

up to a global phase.

The script can also run a random sweep over many angle triples and report the
maximum / average Frobenius error after global-phase alignment.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Iterable, Tuple

import numpy as np


SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


@dataclass(frozen=True)
class VerificationResult:
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


def rotation(axis: str, theta: float) -> np.ndarray:
    """Return exp(-i theta sigma_axis / 2)."""
    sigma = {"x": SX, "y": SY, "z": SZ}[axis]
    return np.cos(theta / 2.0) * I2 - 1j * np.sin(theta / 2.0) * sigma


def target_gate(alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Target SU(2) gate in the same time-ordering convention as the path."""
    return rotation("x", gamma) @ rotation("y", beta) @ rotation("x", alpha)


def path_gate(alpha: float, beta: float, gamma: float, J: float) -> Tuple[np.ndarray, Tuple[float, float, float]]:
    """Return the x-y-x path unitary and the corresponding durations."""
    T1 = alpha / (2.0 * J)
    T2 = beta / (2.0 * J)
    T3 = gamma / (2.0 * J)

    Hx = J * SX
    Hy = J * SY

    U1 = np.cos(J * T1) * I2 - 1j * np.sin(J * T1) * SX
    U2 = np.cos(J * T2) * I2 - 1j * np.sin(J * T2) * SY
    U3 = np.cos(J * T3) * I2 - 1j * np.sin(J * T3) * SX

    # Time ordering: first Hx for T1, then Hy for T2, then Hx for T3.
    U_path = U3 @ U2 @ U1
    return U_path, (T1, T2, T3)


def best_global_phase(U: np.ndarray, V: np.ndarray) -> float:
    """Phase phi minimizing ||e^{-i phi} U - V||_F.

    For 2x2 unitaries, the optimum is phi = -arg(tr(U^† V)/2).
    """
    return float(-np.angle(np.trace(U.conj().T @ V) / 2.0))


def verify_one(alpha: float, beta: float, gamma: float, J: float) -> VerificationResult:
    target = target_gate(alpha, beta, gamma)
    path, (T1, T2, T3) = path_gate(alpha, beta, gamma, J)

    phase = best_global_phase(path, target)
    aligned = np.exp(-1j * phase) * path
    error = float(np.linalg.norm(aligned - target))
    unitarity_error = float(np.linalg.norm(path.conj().T @ path - I2))

    return VerificationResult(
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
    )


def lie_closure_rank() -> int:
    """Return the rank of {i sigma_x, i sigma_y, i [sigma_x, sigma_y]}."""
    comm = SX @ SY - SY @ SX
    basis = np.stack([1j * SX.reshape(-1), 1j * SY.reshape(-1), 1j * comm.reshape(-1)], axis=1)
    return int(np.linalg.matrix_rank(np.real_if_close(basis)))


def random_angles(rng: np.random.Generator, low: float, high: float) -> Tuple[float, float, float]:
    return tuple(float(x) for x in rng.uniform(low, high, size=3))


def sweep(n: int, J: float, seed: int, low: float, high: float) -> Iterable[VerificationResult]:
    rng = np.random.default_rng(seed)
    for _ in range(n):
        alpha, beta, gamma = random_angles(rng, low, high)
        yield verify_one(alpha, beta, gamma, J)


def print_one(result: VerificationResult) -> None:
    print("alpha,beta,gamma =", (result.alpha, result.beta, result.gamma))
    print("J =", result.J)
    print("T1,T2,T3 =", (result.T1, result.T2, result.T3))
    print("best global phase phi =", result.phase)
    print("Frobenius error after phase alignment =", result.error)
    print("unitarity error =", result.unitarity_error)


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify x-y-x control paths for SU(2) gates.")
    parser.add_argument("--alpha", type=float, default=1.1, help="First x-rotation angle.")
    parser.add_argument("--beta", type=float, default=0.7, help="Middle y-rotation angle.")
    parser.add_argument("--gamma", type=float, default=0.4, help="Last x-rotation angle.")
    parser.add_argument("--J", type=float, default=2.3, help="Control strength.")
    parser.add_argument("--sweep", type=int, default=0, help="Number of random samples to test.")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed for sweeps.")
    parser.add_argument("--low", type=float, default=-np.pi, help="Lower bound for random angles.")
    parser.add_argument("--high", type=float, default=np.pi, help="Upper bound for random angles.")
    args = parser.parse_args()

    print("Lie closure rank from {iσx, iσy, i[σx,σy]} =", lie_closure_rank())

    if args.sweep > 0:
        results = list(sweep(args.sweep, args.J, args.seed, args.low, args.high))
        errors = np.array([r.error for r in results], dtype=float)
        unitarity = np.array([r.unitarity_error for r in results], dtype=float)
        print(f"Random sweep size = {args.sweep}")
        print("max error =", float(errors.max()))
        print("mean error =", float(errors.mean()))
        print("max unitarity error =", float(unitarity.max()))
        print("mean unitarity error =", float(unitarity.mean()))
        print("first sample:")
        print_one(results[0])
        return

    result = verify_one(args.alpha, args.beta, args.gamma, args.J)
    print_one(result)

    target = target_gate(args.alpha, args.beta, args.gamma)
    path, _ = path_gate(args.alpha, args.beta, args.gamma, args.J)

    np.set_printoptions(precision=6, suppress=True)
    print("U_target =\n", target)
    print("U_path =\n", path)


if __name__ == "__main__":
    main()