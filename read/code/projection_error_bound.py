"""Evaluate an explicit bound between the scalar formula and the five-mode model.

Only pulse derivatives, small matrices, and one-dimensional quadratures are
used. No dynamics are propagated. The integral inequality is rigorous;
floating-point quadrature values below are not interval-certified bounds.
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad_vec

from full_model import segment_coefficients


DATA = Path(__file__).resolve().parents[1] / "data"


def make_pulses():
    s, x, y, ell = sp.symbols("s x y ell", real=True)
    a = (1-sp.cos(sp.pi*s))/2
    b = 1-a
    q = sp.sqrt(a*a+b*b)
    speed = sp.diff(a, s)/(2*ell*q**3)
    # Coordinates are (gamma_1, f_1, f_2, gamma_a, gamma_c).
    rows = [
        ([2*ell*(x*b-y*a*a)/q, 0], [a*(x+y*b)/q**2, -speed, 0]),
        ([-2*ell*y*b, 2*ell*x*a/q], [x*b/q**2, 0, speed]),
        ([0, 2*ell*x], [0, speed, 0]),
    ]
    functions = []
    for entries, u in rows:
        expressions = (entries + [sp.diff(z, s) for z in entries]
                       + u + [sp.diff(z, s) for z in u]
                       + [sp.diff(z, s, 2) for z in u]
                       + [2*ell*q, 2*ell*sp.diff(q, s)])
        functions.append(sp.lambdify((s, x, y, ell), expressions, "numpy", cse=True))
    return functions


PULSES = make_pulses()


def skew_slow(entries):
    p, q = entries
    return np.array([[0., p, q], [-p, 0., 0.], [-q, 0., 0.]])


def blocks(stage, s, x, y, ell):
    values = np.asarray(PULSES[stage-1](s, x, y, ell), dtype=float)
    A, Ap = skew_slow(values[:2]), skew_slow(values[2:4])
    u, up, upp = values[4:7], values[7:10], values[10:13]
    omega0, omega0p = values[13:]
    sigma = np.sqrt(1+u@u)
    sigmap = u@up/sigma
    factor = 1/(sigma*(sigma+1))
    factorp = -factor*sigmap*(1/sigma+1/(sigma+1))
    M = np.eye(3)-factor*np.outer(u, u)
    Mp = -factorp*np.outer(u, u)-factor*(np.outer(up, u)+np.outer(u, up))
    wedge = np.outer(u, up)-np.outer(up, u)
    wedgep = np.outer(u, upp)-np.outer(upp, u)
    conn = factor*wedge
    connp = factorp*wedge+factor*wedgep
    H = M@A@M-conn
    Hp = Mp@A@M+M@Ap@M+M@A@Mp-connp
    drive = A@u-up
    drivep = Ap@u+A@up-upp
    c = M@drive/sigma
    cp = (Mp@drive+M@drivep)/sigma-c*sigmap/sigma
    omega = omega0*sigma
    omegap = omega0p*sigma+omega0*sigmap
    if omega <= np.linalg.norm(H, 2):
        raise ValueError("The sufficient spectral separation condition fails.")
    denominator = omega**2*np.eye(3)+H@H
    denominatorp = 2*omega*omegap*np.eye(3)+Hp@H+H@Hp
    h = np.linalg.solve(denominator, c)
    hp = np.linalg.solve(denominator, cp-denominatorp@h)
    X = np.column_stack((H@h, omega*h))
    Xp = np.column_stack((Hp@h+H@hp, omegap*h+omega*hp))
    Q = np.eye(5)
    Q[:3, :3], Q[:3, 3], Q[3, :3], Q[3, 3] = M, u/sigma, -u/sigma, 1/sigma
    Qp = np.zeros((5, 5))
    Qp[:3, :3] = Mp
    Qp[:3, 3] = up/sigma-u*sigmap/sigma**2
    Qp[3, :3] = -Qp[:3, 3]
    Qp[3, 3] = -sigmap/sigma**2
    return dict(A=A, H=H, c=c, omega=omega, X=X, Xp=Xp, Q=Q, Qp=Qp,
                u=u, omega0=omega0)


def uniform_gap_lower_bound(epsilon, ell):
    """Analytic sufficient gap, uniform over all three segments."""
    U0 = np.hypot(np.sqrt(2)*epsilon, np.pi/(np.sqrt(2)*ell))
    U1 = np.hypot(2*np.sqrt(2)*np.pi*epsilon, 5*np.pi**2/(2*np.sqrt(2)*ell))
    return (np.sqrt(2)-2*epsilon)*ell-U0*U1/2


def evaluate_bound(x, y, ell):
    if not np.isfinite([x, y, ell]).all() or ell <= 0:
        raise ValueError("Finite parameters and L > 0 are required.")
    epsilon = np.hypot(x, y)
    gap = uniform_gap_lower_bound(epsilon, ell)
    if gap <= 0:
        raise ValueError("Uniform gap certificate fails; use a sharper gap check.")
    reports = []
    total = np.zeros(4)
    quadrature_error = 0.
    for stage in (1, 2, 3):
        def integrand(s):
            z = blocks(stage, s, x, y, ell)
            return np.array([np.linalg.norm(z["H"]-z["A"], 2),
                             np.linalg.norm(z["Xp"], 2),
                             np.linalg.norm(z["c"])*np.linalg.norm(z["X"], 2)])
        integrals, err = quad_vec(integrand, 0., 1., epsabs=1e-10, epsrel=2e-9)
        endpoints = sum(np.linalg.norm(blocks(stage, s, x, y, ell)["X"], 2)
                        for s in (0., 1.))
        parts = np.r_[integrals, endpoints]
        total += 2*parts  # Two cycles; internal endpoint terms are retained.
        # quad_vec estimates the Euclidean error of three integrals; convert
        # that estimate to their sum and include the two cycles.
        quadrature_error += 2*np.sqrt(3)*err
        reports.append(dict(stage=stage, diagonal_integral=float(parts[0]),
                            S_derivative_integral=float(parts[1]),
                            V_S_integral=float(parts[2]), endpoint_norms=float(parts[3])))
    raw = float(total.sum())
    return dict(e_over_C=x, eta_over_C=y, L=ell, epsilon=float(epsilon),
                uniform_gap_lower_bound=float(gap),
                two_cycle_components=dict(zip(
                    ("diagonal", "S_derivative", "V_S", "endpoints"), map(float, total))),
                state_error_bound_integral_value=raw, f_error_bound=min(1., raw),
                F_error_bound=min(1., 2*raw),
                quadrature_error_estimate=float(quadrature_error), stages=reports)


def physical_frame(stage, s):
    """Independent physical moving frame and its derivative."""
    a, ap = (1-np.cos(np.pi*s))/2, np.pi*np.sin(np.pi*s)/2
    b, bp = 1-a, -ap
    q = np.hypot(a, b)
    qp = (a*ap+b*bp)/q
    aq, bq = a/q, b/q
    aqp, bqp = ap/q-a*qp/q**2, bp/q-b*qp/q**2
    if stage == 1:
        vectors = ([bq, aq, 0], [-aq, bq, 0], [0, 0, 1])
        derivatives = ([bqp, aqp, 0], [-aqp, bqp, 0], [0, 0, 0])
    elif stage == 2:
        vectors = ([0, bq, -aq], [-1, 0, 0], [0, aq, bq])
        derivatives = ([0, bqp, -aqp], [0, 0, 0], [0, aqp, bqp])
    else:
        vectors = ([aq, 0, -bq], [-bq, 0, -aq], [0, 1, 0])
        derivatives = ([aqp, 0, -bqp], [-bqp, 0, -aqp], [0, 0, 0])
    T, Tp = np.zeros((5, 5)), np.zeros((5, 5))
    T[0, 0], T[3, 3] = 1., 1.
    for col, index in ((4, 0), (1, 1), (2, 2)):
        T[[4, 1, 2], col] = vectors[index]
        Tp[[4, 1, 2], col] = derivatives[index]
    return T, Tp


def structural_checks():
    maxima = dict(frame_identity=0., orthogonality=0., sylvester=0., S_derivative=0.,
                  endpoint_matching=0.)
    for x, y, ell in ((0., 0., 15.), (1/30, 1/30, 30.), (1/30, 1/3, 120.)):
        for stage in (1, 2):
            Tend, _ = physical_frame(stage, 1.)
            Tstart, _ = physical_frame(stage+1, 0.)
            jump = (Tend@blocks(stage, 1., x, y, ell)["Q"]
                    - Tstart@blocks(stage+1, 0., x, y, ell)["Q"])
            maxima["endpoint_matching"] = max(maxima["endpoint_matching"],
                                               np.linalg.norm(jump, 2))
        for stage, s in ((1, 0.), (3, 1.)):
            maxima["endpoint_matching"] = max(maxima["endpoint_matching"],
                np.linalg.norm(blocks(stage, s, x, y, ell)["Q"]-np.eye(5), 2))
        for stage in (1, 2, 3):
            for s in (.13, .47, .82):
                z = blocks(stage, s, x, y, ell)
                T, Tp = physical_frame(stage, s)
                k0, k1 = segment_coefficients(stage, y, tc=1., E0=1., E=x)
                moving = T.T@(ell*(k0+k1*np.cos(np.pi*s)))@T-T.T@Tp
                expected = np.zeros((5, 5))
                expected[:3, :3] = z["H"]
                expected[:3, 3], expected[3, :3] = z["c"], -z["c"]
                expected[3, 4], expected[4, 3] = z["omega"], -z["omega"]
                transformed = z["Q"].T@moving@z["Q"]-z["Q"].T@z["Qp"]
                maxima["frame_identity"] = max(maxima["frame_identity"],
                    np.linalg.norm(transformed-expected, 2))
                maxima["orthogonality"] = max(maxima["orthogonality"],
                    np.linalg.norm(z["Q"].T@z["Q"]-np.eye(5), 2))
                D = np.array([[0., z["omega"]], [-z["omega"], 0.]])
                B = np.column_stack((z["c"], np.zeros(3)))
                maxima["sylvester"] = max(maxima["sylvester"],
                    np.linalg.norm(z["H"]@z["X"]-z["X"]@D-B, 2))
                step = 2e-6
                derivative = (blocks(stage, s+step, x, y, ell)["X"]
                              - blocks(stage, s-step, x, y, ell)["X"])/(2*step)
                maxima["S_derivative"] = max(maxima["S_derivative"],
                    np.linalg.norm(derivative-z["Xp"], 2))
    assert maxima["frame_identity"] < 2e-10
    assert maxima["orthogonality"] < 2e-14
    assert maxima["sylvester"] < 2e-12
    assert maxima["S_derivative"] < 2e-8
    assert maxima["endpoint_matching"] < 2e-14
    return {key: float(value) for key, value in maxima.items()}


def main():
    checks = structural_checks()
    records = []
    cache = DATA/"solution_vs_scalar_formula.npz"
    cached = None
    if cache.exists():
        with np.load(cache) as saved:
            cached = {key: saved[key].copy() for key in saved.files}
    for ratio in (.1, 1., 2.3446, 10.):
        for tau in (50., 100., 400., 1200.):
            report = evaluate_bound(1/30, ratio/30, .3*tau)
            report.update(ratio=ratio, tau=tau, lambda_value=.01*tau)
            if cached is not None:
                index = int(np.argmin(abs(cached["lambda"]-.01*tau)))
                assert abs(cached["lambda"][index]-.01*tau) < 1e-12
                f_scalar = float(cached[f"f_scalar_r_{ratio}"][index])
                f_full = float(cached[f"f_full_r_{ratio}"][index])
                difference = abs(f_full-f_scalar)
                assert difference <= report["f_error_bound"]+1e-8
                report["existing_array_comparison"] = dict(
                    f_scalar=f_scalar, f_full=f_full, absolute_difference=difference,
                    source=str(cache.relative_to(DATA.parent)))
            records.append(report)
            print(f"r={ratio:g}, tau={tau:g}, bound={report['f_error_bound']:.6g}, "
                  f"raw={report['state_error_bound_integral_value']:.6g}", flush=True)
    result = dict(method="Exact fast-direction rotation and Sylvester integration by parts",
                  normalization="Original K=2A, no rescaling", structural_checks=checks,
                  parameters=dict(E1=.01, C=.3, hbar=1., cycles=2),
                  uniform_domain=dict(epsilon_max=float(np.sqrt(101)/30), L_min=15.,
                      gap_over_L_lower_bound=float(
                          uniform_gap_lower_bound(np.sqrt(101)/30, 15.)/15.)),
                  points=records,
                  caveat="Integral inequality is rigorous; decimal quadratures are not interval-certified.")
    output = DATA/"projection_error_bound.json"
    output.write_text(json.dumps(result, indent=2)+"\n")
    print("Saved:", output, flush=True)


if __name__ == "__main__":
    main()
