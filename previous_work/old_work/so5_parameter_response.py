#!/usr/bin/env python3
"""Exact first and second parameter derivatives in the SO(5) adjoint.

The physical propagation is not truncated in E1 or abs(t1). Derivative
fields are propagated at the actual parameters over two complete cycles.
Energies are meV, tau is inverse meV, and hbar=1.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from so5_full_fidelity import controls, parameter_gradients, skew


def commutator(a, b):
    return a @ b-b @ a


def parameter_response(tau, e1, eta):
    # Fields: B, B_E, B_eta, B_EE, B_Eeta, B_etaeta.
    fields = np.zeros((6, 5, 5))
    fields[0, 0, 1], fields[0, 1, 0] = 1., -1.
    ke = 2*tau*skew(1., 0., 0., 0., 0.)
    kt = 2*tau*skew(0., 1., 0., 0., 0.)

    for segment in range(6):
        def rhs(v, flat):
            b = flat.reshape(6, 5, 5)
            t2, t3, ed, gate = controls(segment % 3, v)
            k = 2*tau*skew(e1, eta*gate, t2, t3, ed)
            keta = gate*kt
            out = commutator(k, b)
            out[1] += commutator(ke, b[0])
            out[2] += commutator(keta, b[0])
            out[3] += 2*commutator(ke, b[1])
            out[4] += commutator(ke, b[2])+commutator(keta, b[1])
            out[5] += 2*commutator(keta, b[2])
            return out.ravel()

        sol = solve_ivp(rhs, (0., 1.), fields.ravel(), method="DOP853",
                        rtol=2e-12, atol=2e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        fields = sol.y[:, -1].reshape(6, 5, 5)

    b = fields[0]
    fidelity = (np.sum(b[:2]**2)-2*b[0, 1])/4
    gradient = np.array([.5*(np.sum(b[:2]*p[:2])-p[0, 1])
                         for p in fields[1:3]])
    hessian = np.empty((2, 2))
    for i, j, index in ((0, 0, 3), (0, 1, 4), (1, 1, 5)):
        value = .5*(np.sum(fields[i+1, :2]*fields[j+1, :2])
                     + np.sum(b[:2]*fields[index, :2])-fields[index, 0, 1])
        hessian[i, j] = hessian[j, i] = value
    assert np.max(np.abs(fields+fields.swapaxes(-1, -2))) < 1e-8
    assert abs(np.linalg.norm(b, "fro")**2-2) < 1e-9
    assert np.max(np.abs(b @ b @ b+b)) < 1e-9
    return float(fidelity), gradient, hessian


def main():
    samples = []
    for tau, e1, eta in ((100., .01, .01), (100., .01, .02344645627)):
        fidelity, gradient, hessian = parameter_response(tau, e1, eta)
        ref_f, ref_g, _ = parameter_gradients(tau, e1, eta)
        assert abs(fidelity-ref_f) < 1e-10
        assert np.max(np.abs(gradient-ref_g[:2])) < 2e-7
        samples.append({"tau": tau, "E1": e1, "t1": eta,
                        "F": fidelity, "gradient_E1_t1": gradient.tolist(),
                        "hessian_E1_t1": hessian.tolist(),
                        "hessian_eigenvalues": np.linalg.eigvalsh(hessian).tolist(),
                        "independent_gradient_residual": float(np.max(np.abs(gradient-ref_g[:2])))})

    tau, e1, eta = 100., .01, .01
    _, _, h = parameter_response(tau, e1, eta)
    finite = np.empty((2, 2))
    for index in range(2):
        epsilon = 5e-8
        plus, minus = [e1, eta], [e1, eta]
        plus[index] += epsilon
        minus[index] -= epsilon
        gp = parameter_gradients(tau, *plus)[1][:2]
        gm = parameter_gradients(tau, *minus)[1][:2]
        finite[:, index] = (gp-gm)/(2*epsilon)
    discrepancy = float(np.max(np.abs(finite-h)))
    assert discrepancy < 1e-3

    # Follow a local peak near the former quadratic cancellation line.
    def slope(value):
        return parameter_gradients(100., .01, value)[1][1]
    root = brentq(slope, .020, .02344645627, xtol=1e-13)
    fidelity, gradient, hessian = parameter_response(100., .01, root)
    assert hessian[1, 1] < 0
    assert abs(gradient[1]) < 1e-6
    ridge_slope = -hessian[0, 1]/hessian[1, 1]
    epsilon_e = 1e-7
    shifted_roots = [brentq(lambda t1: parameter_gradients(100., .01+shift, t1)[1][1],
                           .020, .024, xtol=1e-13)
                     for shift in (-epsilon_e, epsilon_e)]
    finite_ridge_slope = (shifted_roots[1]-shifted_roots[0])/(2*epsilon_e)
    assert abs(ridge_slope-finite_ridge_slope) < 1e-5

    report = {"conventions": {"hbar": 1., "tc_meV": .3, "E0_meV": .3,
                               "gradient_units": "meV^-1", "hessian_units": "meV^-2"},
              "samples": samples,
              "hessian_check": {"exact": h.tolist(), "finite_difference": finite.tolist(),
                                  "max_absolute_residual": discrepancy},
              "local_ridge": {"tau": 100., "E1": .01, "t1": root,
                               "t1_over_E1": root/.01, "F": fidelity,
                               "gradient": gradient.tolist(), "hessian": hessian.tolist(),
                               "dt1_dE1": float(ridge_slope),
                               "independent_dt1_dE1": finite_ridge_slope}}
    output = Path(__file__).with_name("so5_parameter_response.json")
    output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
