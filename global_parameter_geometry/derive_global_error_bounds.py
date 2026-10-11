#!/usr/bin/env python3
"""Evaluate analytic parameter-interpolation bounds; no time evolution is run.

The bounds assume exact time integration at the Gauss nodes, K=2A, and the
existing six-stage protocol. They separate parameter truncation from time
integration error. See global_peak_bounds.md for the derivation.
"""

import argparse
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial import legendre as leg


mp.mp.dps = 70


def endpoint_derivative(degree, derivative):
    """max |P_degree^(derivative)| on [-1,1]."""
    if derivative > degree:
        return mp.mpf(0)
    return mp.factorial(degree + derivative) / (
        mp.power(2, derivative)
        * mp.factorial(derivative)
        * mp.factorial(degree - derivative)
    )


def interpolation_norm(degree, derivative):
    return mp.sqrt(sum(
        (2 * k + 1) * endpoint_derivative(k, derivative) ** 2
        for k in range(degree + 1)
    ))


def remainder_bound(nodes, derivative, phase_width):
    """One-dimensional Gauss interpolation derivative remainder bound."""
    leading = mp.binomial(2 * nodes, nodes) / mp.power(2, nodes)
    return sum(
        math.comb(derivative, j)
        * endpoint_derivative(nodes, derivative - j)
        / leading
        * mp.factorial(j)
        * phase_width ** (nodes + j)
        / mp.factorial(nodes + j)
        for j in range(derivative + 1)
    )


def tensor_bound(P, Q, r, s, chi_E, chi_eta):
    error_x = remainder_bound(P + 1, r, chi_E)
    error_y = remainder_bound(Q + 1, s, chi_eta)
    return min(
        chi_eta ** s * error_x + interpolation_norm(P, r) * error_y,
        chi_E ** r * error_y + interpolation_norm(Q, s) * error_x,
    )


def fidelity_bounds(errors, chi_E, chi_eta):
    e0, ex, ey = errors['00'], errors['10'], errors['01']

    def second(error, ei, ej, ci, cj):
        return (mp.mpf('1.5') * error + e0 * (ci * cj + error)
                + cj * ei + ci * ej + ei * ej)

    return {
        'F': mp.mpf('1.5') * e0 + e0 ** 2 / 2,
        'Fx': mp.mpf('1.5') * ex + e0 * (chi_E + ex),
        'Fy': mp.mpf('1.5') * ey + e0 * (chi_eta + ey),
        'Fxx': second(errors['20'], ex, ex, chi_E, chi_E),
        'Fxy': second(errors['11'], ex, ey, chi_E, chi_eta),
        'Fyy': second(errors['02'], ey, ey, chi_eta, chi_eta),
    }


def gauss_basis(degree):
    nodes, weights = leg.leggauss(degree + 1)
    weights = weights / 2
    basis = leg.legvander(nodes, degree).T * np.sqrt(
        2 * np.arange(degree + 1) + 1
    )[:, None]
    U = basis * np.sqrt(weights)[None, :]
    indices = np.arange(1, degree + 1)
    off = indices / np.sqrt((2 * indices - 1) * (2 * indices + 1))
    J = np.diag(off, 1) + np.diag(off, -1)
    return nodes, weights, U, J


def adjoint_matrix(K):
    pairs = np.triu_indices(5, 1)
    result = np.zeros((10, 10))
    for col, (i, j) in enumerate(zip(*pairs)):
        B = np.zeros((5, 5))
        B[i, j], B[j, i] = 1., -1.
        result[:, col] = (K @ B - B @ K)[pairs]
    return result


def algebra_checks():
    """Check the reduction and the bound formula, independently of a solve."""
    def L(i, j):
        result = np.zeros((5, 5))
        result[i, j], result[j, i] = 1., -1.
        return result

    P, Q = 2, 3
    xn, wx, Ux, Jx = gauss_basis(P)
    yn, wy, Uy, Jy = gauss_basis(Q)
    Mc = adjoint_matrix(2 * (-.17 * L(1, 3) + .29 * L(2, 3) + .13 * L(3, 4)))
    ME = adjoint_matrix(2 * L(0, 1))
    Meta = adjoint_matrix(2 * .37 * L(0, 4))
    center = Mc + .01 * ME + .025 * Meta
    tensor = (np.kron(np.eye(P + 1), np.kron(np.eye(Q + 1), center))
              + .005 * np.kron(Jx, np.kron(np.eye(Q + 1), ME))
              + .015 * np.kron(np.eye(P + 1), np.kron(Jy, Meta)))
    U = np.kron(Ux, np.kron(Uy, np.eye(10)))
    expected = np.zeros_like(tensor)
    for i, x in enumerate(xn):
        for j, y in enumerate(yn):
            index = 10 * (i * (Q + 1) + j)
            expected[index:index + 10, index:index + 10] = (
                center + .005 * x * ME + .015 * y * Meta
            )
    coefficient_initial = np.zeros((P + 1, Q + 1, 10))
    coefficient_initial[0, 0, 0] = 1
    initial_transformed = (U.T @ coefficient_initial.ravel()).reshape(
        coefficient_initial.shape
    )
    expected_initial = np.zeros_like(initial_transformed)
    expected_initial[:, :, 0] = np.sqrt(wx[:, None] * wy[None, :])
    checks = {
        'gauss_orthogonality': float(np.max(np.abs(U.T @ U - np.eye(len(U))))),
        'jacobi_diagonalization': float(max(
            np.max(np.abs(Jx @ Ux - Ux * xn)),
            np.max(np.abs(Jy @ Uy - Uy * yn)),
        )),
        'full_generator_similarity': float(np.max(np.abs(U.T @ tensor @ U - expected))),
        'initial_value_transformation': float(np.max(np.abs(initial_transformed - expected_initial))),
    }
    assert max(checks.values()) < 3e-14, checks

    # The test field has ||d^k b/dx^k||=a^k, matching the hypothesis.
    # These samples check implementation of the formula; the proof is in the note.
    sampled_ratios = []
    xs = np.linspace(-1, 1, 301)
    phase = 6.
    for degree in (4, 8, 12):
        nodes, weights, _, _ = gauss_basis(degree)
        nodal = np.column_stack((np.cos(phase * nodes), np.sin(phase * nodes)))
        coefficients = (leg.legvander(nodes, degree).T @ (weights[:, None] * nodal)
                        * (2 * np.arange(degree + 1) + 1)[:, None])
        for derivative in (0, 1, 2):
            approx = leg.legval(xs, leg.legder(coefficients, derivative, axis=0)).T
            exact = phase ** derivative * np.column_stack((
                np.cos(phase * xs + derivative * np.pi / 2),
                np.sin(phase * xs + derivative * np.pi / 2),
            ))
            bound = float(remainder_bound(degree + 1, derivative, mp.mpf(phase)))
            ratio = float(np.max(np.linalg.norm(exact - approx, axis=1)) / bound)
            assert ratio <= 1 + 1e-10, (degree, derivative, ratio)
            sampled_ratios.append(ratio)
    checks['largest_sampled_remainder_to_bound_ratio'] = max(sampled_ratios)
    checks['scope'] = 'Algebra and interpolation checks; no six-stage evolution or peak computation.'
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tau', default='100')
    parser.add_argument('--delta-E', dest='delta_E', default='0.005')
    parser.add_argument('--delta-eta', dest='delta_eta', default='0.015')
    parser.add_argument('--degrees', type=int, nargs='+', default=[8, 12, 16, 20, 24])
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_name('bounds_at_tau100.json'))
    args = parser.parse_args()
    tau, delta_E, delta_eta = map(mp.mpf, (args.tau, args.delta_E, args.delta_eta))
    if min(tau, delta_E, delta_eta) <= 0 or min(args.degrees) < 0:
        parser.error('tau and parameter half-widths must be positive; degrees must be nonnegative')

    chi_E, chi_eta = 12 * tau * delta_E, 4 * tau * delta_eta
    derivatives = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
    result = {
        'scope': 'A priori parameter-truncation bounds under exact time integration',
        'tau': str(tau), 'delta_E': str(delta_E), 'delta_eta': str(delta_eta),
        'chi_E_end': str(chi_E), 'chi_eta_end': str(chi_eta),
        'derivative_coordinates': 'x=(E-Ec)/delta_E; y=(eta-etac)/delta_eta',
        'precision_decimal_digits': mp.mp.dps,
        'numeric_note': '70-digit evaluation; stored decimal strings are rounded, not interval arithmetic.',
        'algebra_checks': algebra_checks(),
        'rows': [],
    }
    print('degree    state bound       F bound           Fx bound          Fxx bound         Fxy bound')
    for degree in args.degrees:
        errors = {f'{r}{s}': tensor_bound(degree, degree, r, s, chi_E, chi_eta)
                  for r, s in derivatives}
        fidelity = fidelity_bounds(errors, chi_E, chi_eta)
        row = {
            'P': degree, 'Q': degree,
            'state_derivative_bounds': {key: mp.nstr(value, 30) for key, value in errors.items()},
            'fidelity_derivative_bounds': {key: mp.nstr(value, 30) for key, value in fidelity.items()},
            'nodal_time_error_amplification_bounds': {
                f'{r}{s}': mp.nstr(interpolation_norm(degree, r) * interpolation_norm(degree, s), 30)
                for r, s in derivatives
            },
        }
        result['rows'].append(row)
        values = [errors['00'], *[fidelity[key] for key in ('F', 'Fx', 'Fxx', 'Fxy')]]
        print(f'{degree:6d}    ' + '    '.join(f'{float(value):.8e}' for value in values))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    print(f'Saved {args.output}')


if __name__ == '__main__':
    main()
