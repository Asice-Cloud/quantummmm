"""Check the unified function, its parameter response, and spectral channels.

Results are stored in read/data and read/figures/fidelity_formula.
The original gate convention is K=2A.
"""

import json
import os
from pathlib import Path
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/quantumsss_special_functions_mpl")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar

if __package__:
    from .full_model import (L, LE, LH, LP, LR, LD, cosine_connection,
                                        direct_protocol, readout, segment_coefficients)
    from .unified_matrix_response import (INDICES, LABELS, ZETA, angle_response,
                                         conditioning_scales, jet_product, peak_response,
                                         protocol_response, so5_channels)
else:
    from full_model import (L, LE, LH, LP, LR, LD, cosine_connection,
                                       direct_protocol, readout, segment_coefficients)
    from unified_matrix_response import (INDICES, LABELS, ZETA, angle_response,
                                        conditioning_scales, jet_product, peak_response,
                                        protocol_response, so5_channels)


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = HERE.parent / "figures/fidelity_formula"


def independent_derivatives(tau, eta, tc=.3, E0=.3, E=.01):
    """Six continuous stages, explicit differentiated gate equations.

    Does not use the z recurrence, coefficient jets, jet_product, or W^2
    composition.  Coordinate scales are undone before returning.
    """
    st, se, sh = conditioning_scales(tau)
    state = np.zeros((10, 5, 5))
    state[0] = np.eye(5)
    for kind in (1, 2, 3, 1, 2, 3):
        def rhs(u, flat):
            a, c = (1-np.cos(np.pi*u))/2, (1+np.cos(np.pi*u))/2
            if kind == 1:
                p, r, ed, g = tc*a, 0., E0*c, a
            elif kind == 2:
                p, r, ed, g = tc*c, tc*a, 0., c
            else:
                p, r, ed, g = 0., tc*c, E0*a, 0.
            k = 2*(E*LE+eta*g*LH-p*LP+r*LR+ed*LD)
            ke, kh = 2*LE, 2*g*LH
            v = flat.reshape(10, 5, 5)
            result = tau*k@v
            result[1] += st*k@v[0]
            result[2] += se*tau*ke@v[0]
            result[3] += sh*tau*kh@v[0]
            result[4] += 2*st*k@v[1]
            result[5] += st*k@v[2]+se*tau*ke@v[1]+st*se*ke@v[0]
            result[6] += st*k@v[3]+sh*tau*kh@v[1]+st*sh*kh@v[0]
            result[7] += 2*se*tau*ke@v[2]
            result[8] += se*tau*ke@v[3]+sh*tau*kh@v[2]
            result[9] += 2*sh*tau*kh@v[3]
            return result.ravel()
        sol = solve_ivp(rhs, (0., 1.), state.ravel(), method="DOP853",
                        rtol=2e-12, atol=2e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        state = sol.y[:, -1].reshape(10, 5, 5)
    factors = np.prod(np.array([st, se, sh])[None, :]**np.array(INDICES), axis=1)
    return state/factors[:, None, None]


def connection(p, q):
    return cosine_connection(p, q, [1.], order=32, local_phase=1.)[0][0]


def universal_identities():
    x, y = segment_coefficients(2, .023446, E=.013, tc=.27)
    p, q = 5*x, 5*y
    base = connection(p, q)
    change = expm(.37*L(0, 2)-.23*L(1, 4))
    checks = {
        "Q_zero_exponential": float(abs(connection(p, 0*q)-expm(p)).max()),
        "P_zero_endpoint_identity": float(abs(connection(0*p, q)-np.eye(5)).max()),
        "reverse_P_inverse": float(abs(connection(-p, q)@base-np.eye(5)).max()),
        "reverse_P_transpose": float(abs(connection(-p, q)-base.T).max()),
        "constant_basis_covariance": float(abs(
            connection(change@p@change.T, change@q@change.T)-change@base@change.T).max()),
    }
    basis = [L(i, j) for i in range(5) for j in range(i+1, 5)]
    pairs = [(i, j) for i in range(5) for j in range(i+1, 5)]
    def adjoint(a):
        return np.array([[(a@b-b@a)[i, j] for b in basis] for i, j in pairs])
    represented = connection(adjoint(p), adjoint(q))
    expected = np.array([[(base@b@base.T)[i, j] for b in basis] for i, j in pairs])
    checks["adjoint_representation_same_function"] = float(abs(represented-expected).max())
    assert max(checks.values()) < 2e-10, checks
    return checks


def verify_response():
    records = []
    cases = [
        {"tau": 0., "eta": .01}, {"tau": 8., "eta": .1},
        {"tau": 236., "eta": .01}, {"tau": 600., "eta": .023446},
        {"tau": 160., "eta": .02, "E": .013, "tc": .27, "E0": .41},
        {"tau": 80., "eta": 0.}, {"tau": 80., "eta": .01, "E": 0.},
    ]
    for case in cases:
        result = protocol_response(**case)
        ode = independent_derivatives(**case)
        errors = {}
        for j, label in enumerate(LABELS):
            absolute = float(abs(result["R"][j]-ode[j]).max())
            errors[label] = {"max_abs": absolute,
                             "scaled_error": absolute/max(1., float(abs(ode[j]).max()))}
        maximum = max(v["scaled_error"] for v in errors.values())
        assert maximum < 3e-8, (case, errors)
        row = {"parameters": case, "R_derivatives_against_six_stage_ODE": errors,
               "imaginary_scaled_R": result["imaginary_scaled_R"],
               "panels": result["panels"]}
        b = result["B"]
        one_minus_d = -b[:, 0, 1].copy()
        one_minus_d[0] += 1.
        through_b = jet_product(one_minus_d[:, None, None],
                                one_minus_d[:, None, None])[:, 0, 0]
        for i, j in ((2, 3), (2, 4), (3, 4)):
            component = b[:, i, j, None, None]
            through_b -= jet_product(component, component)[:, 0, 0]
        through_b /= 4
        row["F_and_derivatives_through_B_scaled_error"] = float(np.max(
            abs(through_b-result["F"])/np.maximum(1., abs(result["F"]))))
        assert row["F_and_derivatives_through_B_scaled_error"] < 2e-8, row
        if case["tau"]:
            refined = protocol_response(**case, order=34, local_phase=1.2)
            row["refinement_max_scaled_derivative_error"] = max(
                float(abs(result["R"][j]-refined["R"][j]).max())
                /max(1., float(abs(refined["R"][j]).max())) for j in range(10))
            assert row["refinement_max_scaled_derivative_error"] < 2e-9, row
            ch = so5_channels(result["W"][0])
            row["channel_amplitude_error"] = ch["amplitude_error"]
            row["channel_F_error"] = float(abs(ch["F_reconstructed"]-result["F"][0]))
            row["projector_idempotency_error"] = float(max(
                abs(p@p-p).max() for p in ch["projectors"]))
            assert max(row["channel_amplitude_error"], row["channel_F_error"],
                       row["projector_idempotency_error"]) < 2e-7, row
        records.append(row)
        print("response check", json.dumps(case), maximum, flush=True)
    # A physical finite difference separately checks the coefficient/jet wiring.
    tau, E, eta = 90., .013, .023446
    result = protocol_response(tau, eta, E=E)
    fd_errors = {}
    for name, index, step in (("tau", 1, 2e-4), ("E", 2, 2e-7), ("eta", 3, 2e-7)):
        args = {"tau": tau, "E": E, "eta": eta}
        plus, minus = args.copy(), args.copy()
        plus[name] += step
        minus[name] -= step
        derivative = (direct_protocol(**plus)-direct_protocol(**minus))/(2*step)
        error = float(abs(derivative-result["R"][index]).max())/max(1., float(abs(derivative).max()))
        fd_errors[name] = error
        assert error < 2e-7, (name, error)
    return records, fd_errors


def verify_angle_response():
    tau, E, eta = 160., .013, .023446
    base = protocol_response(tau, eta, E=E)
    errors = {}
    for name, index, step in (("tau", 1, 1e-4), ("E", 2, 1e-7), ("eta", 3, 1e-7)):
        def angles(delta):
            args = {"tau": tau, "E": E, "eta": eta}
            args[name] += delta
            stages = [cosine_connection(*segment_coefficients(kind, args["eta"], E=args["E"]),
                                        [args["tau"]])[0][0] for kind in (1, 2, 3)]
            return so5_channels(stages[2]@stages[1]@stages[0])["angles"]
        fd = (angles(step)-angles(-step))/(2*step)
        exact = angle_response(base["W"][0], base["W"][index])
        error = float(abs(fd-exact).max())/max(1., float(abs(exact).max()))
        errors[name] = error
        assert error < 2e-7, errors
    return errors


def peak_prediction():
    E, eta = .01, .01
    root = brentq(lambda tau: protocol_response(tau, eta, E=E)["F"][1],
                  236., 236.4, xtol=2e-9)
    result = protocol_response(root, eta, E=E)
    slopes = peak_response(result)
    assert abs(slopes["F_tau"]) < 1e-10 and slopes["F_tau_tau"] < 0
    trials = []
    for name in ("E", "eta"):
        for delta in (1e-5, 5e-6):
            args = {"E": E, "eta": eta}
            args[name] += delta
            predicted_tau = root+slopes["tau_"+name]*delta
            predicted_height = slopes["F"]+slopes["height_"+name]*delta
            # Independent full ODE, not the jet or its derivative root finder.
            def objective(tau):
                return -float(readout(direct_protocol(tau, **args))[1])
            optimum = minimize_scalar(objective, bounds=(root-1.2, root+1.2),
                                      method="bounded", options={"xatol": 2e-8})
            if not optimum.success:
                raise RuntimeError(optimum.message)
            actual_tau, actual_height = float(optimum.x), float(-optimum.fun)
            trials.append({"parameter": name, "delta": delta,
                           "predicted_tau": predicted_tau, "actual_tau": actual_tau,
                           "tau_residual": actual_tau-predicted_tau,
                           "predicted_F_height": predicted_height,
                           "actual_F_height": actual_height,
                           "F_height_residual": actual_height-predicted_height})
    # A first-order implicit prediction should have quadratic local error.
    for name in ("E", "eta"):
        pair = [row for row in trials if row["parameter"] == name]
        for key in ("tau_residual", "F_height_residual"):
            ratio = abs(pair[0][key]/pair[1][key])
            assert 2.5 < ratio < 5.5, (name, key, ratio)
    return {"base": {"tau": root, "E": E, "eta": eta, **slopes},
            "independent_recomputed_peaks": trials,
            "prediction_type": "local first-order implicit-function formula; no global ridge claim"}


def channel_figure():
    taus = np.linspace(40., 600., 281)
    figure, axes = plt.subplots(2, 2, figsize=(11.4, 6.8), constrained_layout=True)
    records, arrays = [], {"tau": taus}
    for col, ratio in enumerate((1., 10.)):
        stages = [cosine_connection(*segment_coefficients(kind, ratio*.01), taus,
                                   order=28, local_phase=2.)[0] for kind in (1, 2, 3)]
        cycles = stages[2]@stages[1]@stages[0]
        f, F, _ = readout(cycles@cycles)
        channel_values = np.full((len(taus), 2), np.nan+1j*np.nan)
        reconstruction_errors = []
        for n, w in enumerate(cycles):
            try:
                channel = so5_channels(w)
            except ValueError:
                continue
            if channel["amplitude_error"] > 2e-7:
                raise RuntimeError((ratio, taus[n], channel))
            channel_values[n] = channel["channels"]
            reconstruction_errors.append(channel["amplitude_error"])
        diagonal = abs(channel_values)**2
        interference = 2*np.real(channel_values[:, 0]*channel_values[:, 1].conj())
        bound = np.minimum(1., abs(channel_values).sum(axis=1))
        axes[0, col].plot(taus, f, color="#1556a4", label=r"Exact $f=\sqrt{F}$")
        axes[0, col].plot(taus, bound, "--", color="#c67b14", linewidth=1.2,
                          label=r"Channel bound $\min(1,|a_1|+|a_2|)$")
        axes[0, col].set(title=rf"$|t_1|/E_1={ratio:g}$", ylabel=r"$f$", ylim=(-.03, 1.04))
        axes[1, col].plot(taus, F, color="black", linewidth=1.6, label="Exact F")
        # The individual plane labels can exchange at spectral crossings.
        # Their summed intensity and interference are label independent.
        axes[1, col].plot(taus, diagonal.sum(axis=1), "--", color="#30934a",
                          label=r"$D=|a_1|^2+|a_2|^2$")
        axes[1, col].plot(taus, interference, color="#d16036", linewidth=1.,
                          label=r"$2\,\mathrm{Re}(a_1 a_2^*)$")
        axes[1, col].set(xlabel=r"$\tau\ (\mathrm{meV}^{-1})$", ylabel="F and channel terms")
        for row in range(2):
            axes[row, col].grid(alpha=.22)
            axes[row, col].set_xlim(40, 600)
        records.append({"ratio": ratio, "regular_chart_points": len(reconstruction_errors),
                        "total_points": len(taus),
                        "channel_amplitude_error_max": max(reconstruction_errors),
                        "triangle_bound_minus_f_max": float(np.nanmax(bound-f)),
                        "bound_is": "pointwise triangle bound; not an established slow envelope"})
        arrays[f"f_ratio_{ratio}"] = f
        arrays[f"F_ratio_{ratio}"] = F
        arrays[f"channels_ratio_{ratio}"] = channel_values
        arrays[f"bound_ratio_{ratio}"] = bound
    axes[0, 0].legend(fontsize=8, loc="lower left")
    axes[1, 0].legend(fontsize=8, loc="best")
    figure.suptitle(r"Unified matrix function: exact SO(5) channel decomposition"
                   "\n"+r"$E_1=0.01,\ t_c=E_0=0.3\ \mathrm{meV}$; original gates, $K=2A$, two cycles")
    figure.savefig(OUT/"unified_channels.png", dpi=170)
    plt.close(figure)
    np.savez_compressed(DATA/"unified_channels.npz", **arrays)
    return records


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    report = {"physical_protocol": "original cosine gates; K=2A; six continuous stages",
              "jet_variables": ["tau", "E1", "eta=abs(t1)"],
              "jet_order": 2, "base_evolution_truncated_in_parameters": False}
    report["universal_function_identities"] = universal_identities()
    report["parameter_response_checks"], report["physical_finite_difference_checks"] = verify_response()
    report["angle_response_finite_difference_checks"] = verify_angle_response()
    report["local_peak_prediction"] = peak_prediction()
    report["channel_figure"] = channel_figure()
    report["elapsed_seconds"] = time.monotonic()-start
    (DATA/"unified_function_verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"identities": report["universal_function_identities"],
                      "peak_prediction": report["local_peak_prediction"],
                      "channel_figure": report["channel_figure"],
                      "elapsed_seconds": report["elapsed_seconds"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
