"""Cross-check the scalar continuation formula through total degree four.

The symbolic route starts with C,S integral moments and the beta readout,
independently of the Pauli-product derivation in perturbation.md. Numerical
checks integrate the defining C,S matrix equation and test the sixth-order
Taylor remainder. This checks the projected model, not its physical error
relative to the full five-Majorana model. No normalization rescaling is used.
"""

import json
import math
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.special import ellipe, ellipk, ellipkinc


DATA = Path(__file__).resolve().parents[1] / "data"


def symbolic_check():
    alpha, nu, t, j, k, J, T, U = sp.symbols(
        "alpha nu t j k J_p T_p U_p", real=True)
    A, B = sp.Rational(1, 2) + j, 1 + 2*k
    m1 = (k-J)/2
    m2 = T/2 + 3*k/16 - J/4
    d1 = k**2/8 - U/2
    C = 1 - (t*alpha)**2*k**2/2 - 2*sp.I*t**3*alpha**2*nu*d1
    S = t*alpha*k + 2*sp.I*t**2*alpha*nu*m1 \
        - 2*t**3*alpha*nu**2*m2 - t**3*alpha**3*k**3/6

    def trunc(expr, degree):
        expanded = sp.expand(expr)
        return sp.Add(*(expanded.coeff(t, n)*t**n for n in range(degree+1)))

    def exp_poly(argument, degree):
        return sum(argument**n/sp.factorial(n) for n in range(degree+1))

    # z = i*sqrt(2)*beta. Terms of degree four in z cannot contribute to
    # degree four of f: the readout has zero first derivative at z=1,delta=0.
    z = trunc(exp_poly(sp.I*t*nu/2, 3)*sp.conjugate(C)
              - sp.I*exp_poly(-sp.I*t*nu/2, 3)*S, 3)
    delta = t*((1+k)*alpha-j*nu)
    rotated = trunc(exp_poly(-sp.I*delta, 4)*z, 4)
    real_rotated = sp.expand((rotated+sp.conjugate(rotated))/2)
    norm_z = trunc(z*sp.conjugate(z), 4)
    direct_f = trunc((2-norm_z)*trunc(real_rotated**2, 4), 4)

    dz = z-1
    log_z = trunc(dz-dz**2/2+dz**3/3, 3)
    phase = sp.expand((log_z-sp.conjugate(log_z))/(2*sp.I)-delta)
    visibility = trunc(1-(norm_z-1)**2, 4)
    v = B*alpha-A*nu
    a4 = U-k**2/4-(1+k)*J
    b4 = k/8+j*J-T
    expected_phase = -t*v-t**3*alpha*nu*(J*v+a4*alpha+b4*nu)
    expected_visibility = 1-4*t**4*J**2*alpha**2*nu**2
    p4 = v**4/3-2*alpha*nu*J*v**2-4*alpha**2*nu**2*J**2 \
        - 2*alpha*nu*v*(a4*alpha+b4*nu)
    expected_f = 1-t**2*v**2+t**4*p4

    identities = {
        "phase_through_degree_three": sp.expand(phase-expected_phase),
        "visibility_through_degree_four": sp.expand(visibility-expected_visibility),
        "f_through_degree_four": sp.expand(direct_f-expected_f),
        "F_through_degree_four": sp.expand(
            trunc(direct_f**2, 4)-(1-2*t**2*v**2+t**4*(v**4+2*p4))),
    }
    for name, residual in identities.items():
        if residual != 0:
            raise AssertionError(f"Symbolic identity failed: {name}: {residual}")

    coefficient_rows = []
    for degree in (2, 4):
        direct = sp.Poly(sp.expand(direct_f.coeff(t, degree)), alpha, nu)
        target = sp.Poly(sp.expand(expected_f.coeff(t, degree)), alpha, nu)
        for power_alpha in range(degree, -1, -1):
            monomial = alpha**power_alpha*nu**(degree-power_alpha)
            coefficient = direct.coeff_monomial(monomial)
            difference = sp.expand(coefficient-target.coeff_monomial(monomial))
            assert difference == 0
            coefficient_rows.append({
                "alpha_power": power_alpha, "nu_power": degree-power_alpha,
                "coefficient": str(sp.factor(coefficient)), "symbolic_difference": "0",
            })
    return ({"identities": {key: str(value) for key, value in identities.items()},
             "coefficients": coefficient_rows,
             "odd_degrees_one_and_three_zero": all(
                 sp.expand(direct_f.coeff(t, n)) == 0 for n in (1, 3))},
            (j, k, J, T, U))


def pulse(s):
    a = (1-np.cos(np.pi*s))/2
    b = 1-a
    return b, a/np.hypot(a, b), s/2+np.sin(np.pi*s)/(2*np.pi)


def constants():
    k, j = float(ellipk(.5)/np.pi), float(ellipe(.5)/np.pi)

    def G(s):
        return (ellipkinc(np.pi*s, .5)
                - np.sqrt(2)*np.arcsin(np.sin(np.pi*s)/np.sqrt(2)))/(2*np.pi)

    def integral(fun):
        return float(quad(fun, 0., 1., epsabs=2e-13, epsrel=2e-13)[0])

    J = integral(lambda s: pulse(s)[1]*(1-2*pulse(s)[2]))
    T = integral(lambda s: .5*pulse(s)[1]*(2*pulse(s)[2]-.5)**2)
    U = integral(lambda s: .5*pulse(s)[0]*(2*G(s)-k)**2)
    m1 = integral(lambda s: pulse(s)[1]*pulse(s)[2])
    m2 = integral(lambda s: pulse(s)[1]*pulse(s)[2]**2)
    d1 = integral(lambda s: pulse(s)[0]*G(s)*(k-G(s)))
    residuals = {
        "m1_minus_reduction": m1-(k-J)/2,
        "m2_minus_reduction": m2-(T/2+3*k/16-J/4),
        "D1_minus_reduction": d1-(k*k/8-U/2),
    }
    assert max(abs(value) for value in residuals.values()) < 2e-12
    return {"j": j, "k": k, "A": .5+j, "B": 1+2*k,
            "J_p": J, "T_p": T, "U_p": U,
            "a4": U-k*k/4-(1+k)*J, "b4": k/8+j*J-T,
            "moment_identity_residuals": residuals}


def scalar_reference(alpha, nu, c):
    """Integrate the normalized-function IVP, then use the exact beta readout."""
    def rhs(s, flat):
        _, g, P = pulse(s)
        kernel = g*np.exp(2j*nu*P)
        matrix = alpha*np.array([[0., -kernel], [kernel.conjugate(), 0.]])
        return (matrix@flat.reshape(2, 2)).ravel()

    sol = solve_ivp(rhs, (0., 1.), np.eye(2, dtype=complex).ravel(),
                    method="DOP853", rtol=3e-14, atol=1e-15)
    if not sol.success:
        raise RuntimeError(sol.message)
    W = sol.y[:, -1].reshape(2, 2)
    C, S = W[0, 0], -W[0, 1]
    beta = -(np.exp(-.5j*nu)*S+1j*np.exp(.5j*nu)*C.conjugate())/np.sqrt(2)
    delta = alpha*(1+c["k"])-nu*c["j"]
    return float(4*(1-abs(beta)**2)*np.imag(np.exp(-1j*delta)*beta)**2)


def numerical_check(c):
    rstar = c["B"]/c["A"]
    rays = [("alpha_axis", (1., 0.)), ("nu_axis", (0., 1.)),
            ("r=0.1", (1., .1)), ("r=1", (1., 1.)),
            ("r=r_star", (1/rstar, 1.)), ("r=10", (.1, 1.))]
    reports = []
    max_parity_error = 0.
    for name, direction in rays:
        rows = []
        for h in (.4, .2, .1, .05):
            alpha, nu = h*np.array(direction)
            v = c["B"]*alpha-c["A"]*nu
            f4 = (1-v*v+v**4/3-2*alpha*nu*c["J_p"]*v*v
                  - 4*alpha**2*nu**2*c["J_p"]**2
                  - 2*alpha*nu*v*(c["a4"]*alpha+c["b4"]*nu))
            exact = scalar_reference(alpha, nu, c)
            parity_error = abs(exact-scalar_reference(-alpha, -nu, c))
            max_parity_error = max(max_parity_error, parity_error)
            residual = exact-f4
            Q = c["B"]*abs(alpha)+c["A"]*abs(nu)
            # Conservative rigorous bound for the Taylor remainder in the
            # projected model; not a bound for the ODE numerical error.
            bound = math.exp(4*Q)*(4*Q)**6/math.factorial(6)
            assert abs(residual) <= bound+2e-13
            row = {"h": h, "alpha": float(alpha), "nu": float(nu),
                   "f_reference": exact, "f_fourth_order": float(f4),
                   "signed_remainder": float(residual),
                   "remainder_divided_by_h6": float(residual/h**6),
                   "conservative_analytic_remainder_bound": bound}
            if rows:
                row["observed_order_on_halving"] = float(np.log2(
                    abs(rows[-1]["signed_remainder"]/residual)))
            rows.append(row)
        assert 5.8 < rows[-1]["observed_order_on_halving"] < 6.2
        reports.append({"ray": name, "direction_alpha_nu": list(direction),
                        "samples": rows})
    assert max_parity_error < 2e-12
    return reports, max_parity_error


def main():
    symbolic, symbols = symbolic_check()
    c = constants()
    substitutions = dict(zip(symbols, (c["j"], c["k"], c["J_p"], c["T_p"], c["U_p"])))
    local_symbols = {str(symbol): symbol for symbol in symbols}
    for row in symbolic["coefficients"]:
        coefficient = sp.sympify(row["coefficient"], locals=local_symbols)
        row["numeric_coefficient"] = float(coefficient.subs(substitutions))
    numerical, parity_error = numerical_check(c)
    report = {
        "model": "First-order projected model; original K=2A convention; no factor-2 trial",
        "comparison": "C,S moment expansion plus beta readout vs independent Pauli-product perturbation polynomial",
        "symbolic": symbolic, "constants": c,
        "numerical_reference": "DOP853 integration of normalized C,S IVP, rtol=3e-14, atol=1e-15",
        "numerical_samples": numerical, "max_simultaneous_sign_reversal_error": parity_error,
        "scope": "Algebraic consistency and local Taylor remainder only; no claim of full-model projection error certification",
    }
    DATA.mkdir(parents=True, exist_ok=True)
    output = DATA/"continuation_fourth_order_check.json"
    output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"output": str(output), "symbolic_identities": symbolic["identities"],
                      "constants": c, "coefficients": symbolic["coefficients"],
                      "observed_final_orders": {row["ray"]: row["samples"][-1]["observed_order_on_halving"]
                                                for row in numerical}}, indent=2))


if __name__ == "__main__":
    main()
