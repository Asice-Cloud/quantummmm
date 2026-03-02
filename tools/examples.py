#!/usr/bin/env python3
"""Example anyon models (Ising, Fibonacci) used for tests and demos.

Contains small functions returning (fusion_rules, F_dict, R_dict) compatible
with `tools/anyon_solver.py` helpers.
"""
from __future__ import annotations
import cmath
from typing import Dict, Tuple, List, Any

Anyon = str
FusionRules = Dict[Tuple[Anyon, Anyon], List[Anyon]]


def ising_model() -> Tuple[FusionRules, Dict[str, Any], Dict[str, Any]]:
    # labels: 1, sigma, psi
    rules = {
        ("1", "1"): ["1"],
        ("1", "sigma"): ["sigma"],
        ("sigma", "1"): ["sigma"],
        ("1", "psi"): ["psi"],
        ("psi", "1"): ["psi"],
        ("sigma", "sigma"): ["1", "psi"],
        ("sigma", "psi"): ["sigma"],
        ("psi", "sigma"): ["sigma"],
        ("psi", "psi"): ["1"],
    }
    # F^{sigma}_{sigma sigma sigma} in basis [1,psi]
    F = {
        "F^sigma_sigma_sigma_sigma": [[1 / 2**0.5, 1 / 2**0.5], [1 / 2**0.5, -1 / 2**0.5]]
    }
    # R^e_{sigma sigma}
    R = {"R^1_sigma_sigma": cmath.exp(-1j * cmath.pi / 8), "R^psi_sigma_sigma": cmath.exp(3j * cmath.pi / 8)}
    return rules, F, R


def fibonacci_model() -> Tuple[FusionRules, Dict[str, Any], Dict[str, Any]]:
    # labels: 1, tau
    phi = (1 + 5 ** 0.5) / 2
    rules = {
        ("1", "1"): ["1"],
        ("1", "tau"): ["tau"],
        ("tau", "1"): ["tau"],
        ("tau", "tau"): ["1", "tau"],
    }
    F = {
        # F^tau_{tau tau tau} in basis [1,tau]
        "F^tau_tau_tau_tau": [[1 / phi, 1 / phi**0.5], [1 / phi**0.5, -1 / phi]]
    }
    R = {"R^1_tau_tau": cmath.exp(-4j * cmath.pi / 5), "R^tau_tau_tau": cmath.exp(3j * cmath.pi / 5)}
    return rules, F, R


__all__ = ["ising_model", "fibonacci_model"]
