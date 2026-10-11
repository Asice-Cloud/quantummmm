#!/usr/bin/env python3
"""Check the paper's overlap, parameter expansion, and adiabatic formulas.

Run with .venv/bin/python verify_so5_fidelity.py. This performs a small set
of verification calculations, not a parameter-map scan. Energies are meV,
tau is inverse meV, and hbar=1. The paper gates t1 with t2 and repeats the
three-step protocol twice. Dimensionless perturbation amplitudes are
x=2*tau*E1 and y=2*tau*abs(t1).
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.special import ellipe, ellipk


def generator(i, j, dimension=5):
    out = np.zeros((dimension, dimension))
    out[i, j], out[j, i] = 1.0, -1.0
    return out


LE = generator(0, 1)
LT = generator(0, 4)
CM = np.array([1, -1j, 0, 0, 0]) / np.sqrt(2)
CP = CM.conj()
C = float(ellipk(0.5) / np.pi)
A = float(ellipe(0.5) / np.pi)
B = 0.5 + A
RIDGE = (1 + 2 * C) / B


def controls(s):
    segment = min(int(s), 5)
    local = s - segment
    rise = (1 - np.cos(np.pi * local)) / 2
    fall = 1 - rise
    if segment % 3 == 0:
        return 0.3 * rise, 0.0, 0.3 * fall, rise
    if segment % 3 == 1:
        return 0.3 * fall, 0.3 * rise, 0.0, fall
    return 0.0, 0.3 * fall, 0.3 * rise, 0.0


def k0(s, tau):
    t2, t3, ed, gate = controls(s)
    mat = 2 * tau * (
        t2 * generator(3, 1)
        - t3 * generator(3, 2)
        + ed * generator(3, 4)
    )
    return mat, gate


def integrate_segments(rhs, initial, count=6):
    y = initial.ravel()
    for step in range(count):
        sol = solve_ivp(rhs, (step, step + 1), y, method="DOP853",
                        rtol=2e-12, atol=2e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
    return y.reshape(initial.shape)


def endpoint(tau, x, y):
    def rhs(s, flat):
        mat, gate = k0(s, tau)
        return ((mat + x * LE + y * gate * LT) @ flat.reshape(5, 5)).ravel()
    return integrate_segments(rhs, np.eye(5))


def amplitude(r):
    return np.vdot(CP, r @ CM)


def fidelity_entries(r):
    return float(((r[0, 0] - r[1, 1]) ** 2
                  + (r[0, 1] + r[1, 0]) ** 2) / 4)


def tangent_expansion(tau):
    # Taylor coefficients, not second derivatives: mixed term is x*y*Cxy.
    initial = np.zeros((6, 5, 5))
    initial[0] = np.eye(5)

    def rhs(s, flat):
        r = flat.reshape(6, 5, 5)
        mat, gate = k0(s, tau)
        out = mat @ r
        out[1] += LE @ r[0]
        out[2] += gate * LT @ r[0]
        out[3] += LE @ r[1]
        out[4] += LE @ r[2] + gate * LT @ r[1]
        out[5] += gate * LT @ r[2]
        return out.ravel()

    r = integrate_segments(rhs, initial)
    a0, ax, ay, axx, axy, ayy = map(amplitude, r)
    hessian = np.array([
        [abs(ax) ** 2 + 2 * (a0.conjugate() * axx).real,
         (ax.conjugate() * ay + a0.conjugate() * axy).real],
        [(ax.conjugate() * ay + a0.conjugate() * axy).real,
         abs(ay) ** 2 + 2 * (a0.conjugate() * ayy).real],
    ])
    omega_x = r[0].T @ r[1]
    omega_y = r[0].T @ r[2]
    linear = [2 * (a0.conjugate() * ax).real,
              2 * (a0.conjugate() * ay).real]
    return float(abs(a0) ** 2), hessian, omega_x[0, 1:], omega_y[0, 1:], linear


def alpha(v):
    rise = (1 - np.cos(np.pi * v)) / 2
    return np.arctan2(rise, 1 - rise)


def slow_endpoint(x, y):
    l12 = generator(0, 1, 3)
    l13 = generator(0, 2, 3)
    j = np.array([[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]])
    s1 = expm((C * x - A * y) * l12)
    s3 = expm(x * l13)

    def rhs(v, flat):
        fall = (1 + np.cos(np.pi * v)) / 2
        k = -y * fall * l12 + x * np.sin(alpha(v)) * l13
        return (k @ flat.reshape(3, 3)).ravel()

    s2 = integrate_segments(rhs, np.eye(3), count=1)
    cycle = j @ s3 @ s2 @ s1
    return cycle @ cycle


def main():
    results = {"conventions": {"tc_meV": 0.3, "E0_meV": 0.3,
                                "hbar": 1, "x": "2*tau*E1",
                                "y": "2*tau*abs(t1)"},
               "constants": {"C": C, "A": A, "B": B, "ridge_t1_over_E1": RIDGE}}
    random = np.random.default_rng(71).normal(size=(5, 5))
    r = expm(random - random.T)
    f = fidelity_entries(r)
    f_overlap = float(abs(amplitude(r)) ** 2)
    returning = abs(np.vdot(CM, r @ CM)) ** 2
    leakage = np.linalg.norm((r @ CM)[2:]) ** 2
    assert abs(f - f_overlap) < 1e-13
    assert abs(f + returning + leakage - 1) < 1e-13
    results["exact_overlap_identity_residual"] = abs(f - f_overlap)
    results["channel_conservation_residual"] = float(abs(f + returning + leakage - 1))

    asymptotic_hessian = -0.5 * np.outer([1 + 2 * C, -B], [1 + 2 * C, -B])
    coefficients = []
    for tau in [100., 300., 1000.]:
        f0, h, ax, ay, linear = tangent_expansion(tau)
        assert max(abs(v) for v in linear) < 1e-12
        coefficients.append({"tau": tau, "F0": f0, "quadratic_matrix_xy": h.tolist(),
                             "normalized_response_E": ax.tolist(),
                             "normalized_response_t1": ay.tolist(),
                             "linear_terms": linear})
    results["finite_time_coefficients"] = coefficients
    results["adiabatic_quadratic_matrix_xy"] = asymptotic_hessian.tolist()
    assert np.max(np.abs(np.array(coefficients[-1]["quadratic_matrix_xy"])
                         - asymptotic_hessian)) < 3e-4

    tau = 100.
    f0, h, _, _, _ = tangent_expansion(tau)
    expansion_checks = []
    for scale in [1., .5, .25]:
        xy = np.array([.08, .032]) * scale
        actual = fidelity_entries(endpoint(tau, *xy))
        predicted = f0 + xy @ h @ xy
        expansion_checks.append({"x": float(xy[0]), "y": float(xy[1]),
                                 "F_exact": actual, "F_quadratic": float(predicted),
                                 "absolute_residual": abs(actual - predicted)})
    residuals = [row["absolute_residual"] for row in expansion_checks]
    ratios = [residuals[i] / residuals[i + 1] for i in range(2)]
    assert all(14 < ratio < 18 for ratio in ratios)
    results["finite_time_expansion_checks"] = expansion_checks
    results["residual_ratios_on_halving"] = ratios

    pure_checks = []
    for x, y in [(.8, 0.), (0., .8)]:
        slow = slow_endpoint(x, y)
        f_slow = (1 - slow[2, 2]) / 2
        expected = np.cos(((1 + 2 * C) * x - B * y) / 2) ** 2
        assert abs(f_slow - expected) < 1e-11
        full_f = np.sqrt(fidelity_entries(endpoint(1000., x, y)))
        assert abs(full_f - expected) < 2e-4
        pure_checks.append({"x": x, "y": y, "f_slow": float(f_slow),
                            "f_analytic": float(expected), "f_full_tau1000": float(full_f)})
    results["pure_parameter_checks"] = pure_checks

    ordered_integral = (2 * quad(lambda v: (v / 2 + np.sin(np.pi * v) / (2 * np.pi))
                               * np.sin(alpha(v)), 0, 1, epsabs=1e-13)[0] - C / 2)
    area_cross = A * (1 + C) + .5 + ordered_integral
    area_on_ridge = (1 + 2 * C) ** 2 - (area_cross + C * B) * RIDGE
    quartic_coefficient = .5 * area_on_ridge ** 2
    ridge_checks = []
    for x in [.1, .05, .025]:
        slow = slow_endpoint(x, RIDGE * x)
        infidelity = 1 - ((1 - slow[2, 2]) / 2) ** 2
        ridge_checks.append({"x": x, "one_minus_F": float(infidelity),
                             "one_minus_F_over_x4": float(infidelity / x ** 4)})
    assert abs(ridge_checks[-1]["one_minus_F_over_x4"] - quartic_coefficient) < 3e-5
    results["ridge_quartic_coefficient"] = float(quartic_coefficient)
    results["ridge_slow_model_checks"] = ridge_checks

    output = Path(__file__).with_name("so5_fidelity_validation.json")
    output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "constants": results["constants"],
                      "residual_ratios_on_halving": ratios,
                      "pure_parameter_checks": pure_checks,
                      "ridge_quartic_coefficient": quartic_coefficient}, indent=2))


if __name__ == "__main__":
    main()
