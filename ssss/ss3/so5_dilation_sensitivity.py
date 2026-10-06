"""Partial parameter derivatives of the subspace-dilation maps.

Question: for the PRB111 three-stage protocol of model.md / n7.md, may we take
partial derivatives with respect to E1 and t1, and how much do higher-order
parameter terms contribute to the "space expansion"?

Space expansion is measured by the principal angles of the physical 3-frame
W = R(T)[I3;0] = [M; N] inside R^5 (n5.md, n7.md):

    cos(theta_i) = singular values of M,     N^T N = I_3 - M^T M,
    volume factor J = det M = prod cos(theta_i),
    expansion cost E_leak = (3 - sum cos(theta_i))/6
                          = (||S-I_3||_F^2 + ||N||_F^2)/12,  S = (M^T M)^(1/2).

Method.  With x = E1/e_scale and y = t1/t1_scale the propagator obeys

    dR/ds = tau [o(s) + x oe(s) + y og(s)] R,

and o, oe, og are independent of (x, y).  Hence R is analytic in (x, y), and the
bivariate propagator coefficients c[m, n] (Taylor coefficients, factorials
already divided out) at (x, y) = (0, 0) follow from a single linear hierarchy.
The observables are then computed by exact polynomial algebra: the matrix square
root is expanded as a series about M^T M = I_3, and det(I + X) as exp tr log.
Every mixed partial derivative is therefore available in closed form.

Nothing is fitted.  The coefficients are cross-checked against parameter finite
differences of the full nonlinear propagator, and the truncated series are
compared with the exact evolution over the physical parameter range.

The protocol is gated, not adiabatic, so at (E1, t1) = (0, 0) the frame is not
yet the ideal braid B when tau is small.  The series therefore keeps its
constant term, which is reported separately as ideal_protocol_baseline, and the
maps show the E1, t1-induced change with that baseline subtracted.

Units: energies in meV, time in meV^-1, SO(5) normalization factor 2 retained.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from verify_manifold_expansion import B, Model


# --------------------------------------------------------------------------- #
# Propagator and its parameter coefficients
# --------------------------------------------------------------------------- #
def integrate(model, rhs, y0, n_out=None):
    """Run the six protocol segments; optionally sample the trajectory."""
    y = np.asarray(y0, dtype=float).ravel()
    traj, times = [], []
    for segment in range(6):
        t_eval = None if not n_out else np.linspace(0., 1., n_out, endpoint=segment == 5)
        sol = solve_ivp(lambda s, v: rhs(segment, s, v), (0., 1.), y,
                        method="DOP853", rtol=model.rtol, atol=model.atol, t_eval=t_eval)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
        if n_out:
            traj.append(sol.y)
            times.append(segment + sol.t)
    if not n_out:
        return y, None, None
    return y, np.concatenate(times), np.concatenate(traj, axis=1)


def propagator_coefficients(model):
    """Bivariate Taylor coefficients c[m, n] of R(x, y) about (0, 0)."""
    p = model.series
    initial = np.zeros((p.count, 5, 5))
    initial[p.index[(0, 0)]] = np.eye(5)
    with_x = [k for k, (m, _) in enumerate(p.keys) if m]
    from_x = [p.index[(m-1, n)] for m, n in p.keys if m]
    with_y = [k for k, (_, n) in enumerate(p.keys) if n]
    from_y = [p.index[(m, n-1)] for m, n in p.keys if n]

    def rhs(segment, s, vector):
        r = vector.reshape(p.count, 5, 5)
        o, oe, og = model.generators(segment, s)
        result = o @ r
        result[with_x] += oe @ r[from_x]
        result[with_y] += og @ r[from_y]
        return result.ravel()
    return model.integrate(rhs, initial).reshape(p.count, 5, 5)


# --------------------------------------------------------------------------- #
# Polynomial algebra (Taylor coefficients in x, y)
# --------------------------------------------------------------------------- #
def scalar_mul(p, a, b):
    """Product of two scalar Taylor series in (x, y)."""
    out = np.zeros_like(a)
    for k, terms in enumerate(p.pairs):
        out[k] = sum(a[i]*b[j] for i, j in terms)
    return out


def determinant_3x3(p, m):
    """det of a 3x3 matrix Taylor series: cubic polynomial in the entries."""
    e = [[m[:, i, j] for j in range(3)] for i in range(3)]
    mul = lambda a, b: scalar_mul(p, a, b)
    minor = lambda a, b, c, d: mul(a, b) - mul(c, d)
    return (mul(e[0][0], minor(e[1][1], e[2][2], e[1][2], e[2][1]))
            - mul(e[0][1], minor(e[1][0], e[2][2], e[1][2], e[2][0]))
            + mul(e[0][2], minor(e[1][0], e[2][1], e[1][1], e[2][0])))


def expansion_series(coefficients, model):
    """Parameter series of the dilation observables, exact to the model order.

    S = (M^T M)^(1/2) is obtained by the Sylvester recursion of Series.sqrt_spd,
    which makes no smallness assumption on M^T M - I_3, and det M is a cubic
    polynomial in the entries.  The constant term is retained: the protocol at
    (E1, t1) = (0, 0) is adiabatic only for large tau, so E_leak has a baseline
    (see ideal_protocol_baseline in the report).
    """
    p = model.series
    m_block = coefficients[:, :3, :3]
    n_block = coefficients[:, 3:, :3]
    gram = p.matmul(m_block.transpose(0, 2, 1), m_block)
    stretch = p.sqrt_spd(gram)
    trace_s = np.trace(stretch, axis1=1, axis2=2)
    leak = -trace_s/6.
    leak[p.index[(0, 0)]] += .5                     # (3 - tr S)/6 at the base point
    leakage = np.trace(p.matmul(n_block.transpose(0, 2, 1), n_block),
                       axis1=1, axis2=2)
    return {"E_leak": leak, "det_M": determinant_3x3(p, m_block),
            "leak": leakage, "trace_S": trace_s}


def series_derivative(p, a, axis):
    """Partial derivative of a scalar Taylor series w.r.t. x or y."""
    out = np.zeros_like(a)
    for k, (m, n) in enumerate(p.keys):
        if (axis == "x" and m) or (axis == "y" and n):
            key = (m-1, n) if axis == "x" else (m, n-1)
            out[p.index[key]] = (m if axis == "x" else n)*a[k]
    return out


def physical_powers(p, a, e_scale, g_scale):
    return {f"{m},{n}": float(a[k]/e_scale**m/g_scale**n)
            for k, (m, n) in enumerate(p.keys)}


# --------------------------------------------------------------------------- #
# Observables of the exact propagator
# --------------------------------------------------------------------------- #
def dilation(R):
    """Principal angles, volume factor and expansion costs of the frame R[:, :3]."""
    w = R[:, :3]
    m_block, n_block = w[:3], w[3:]
    cos_theta = np.clip(np.linalg.svd(m_block, compute_uv=False), 0., 1.)
    theta = np.arccos(cos_theta)
    sin_theta = np.linalg.svd(n_block, compute_uv=False)
    return {
        "cos_theta": cos_theta,
        "theta": theta,
        "sin_theta_from_N": sin_theta,
        "volume": float(np.linalg.det(m_block)),
        "E_leak": float((3.-cos_theta.sum())/6.),
        "leak": float(np.sum(n_block*n_block)),
        "grassmann_distance": float(np.sqrt(np.sum(theta**2))),
        "orthogonality_defect": float(np.linalg.norm(R.T @ R - np.eye(5))),
    }


def frame(R):
    w = R[:, :3]
    return w, w[:3], w[3:]# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #
def fd_derivatives(model, h=1.e-2):
    """Richardson-extrapolated parameter derivatives of the exact propagator."""
    def first(key, step):
        if key == (1, 0):
            return (model.exact(step, 0.) - model.exact(-step, 0.))/(2*step)
        return (model.exact(0., step) - model.exact(0., -step))/(2*step)

    def pure(key, step):
        centre = model.exact(0., 0.)
        if key == (2, 0):
            return (model.exact(step, 0.) - 2*centre + model.exact(-step, 0.))/step**2
        return (model.exact(0., step) - 2*centre + model.exact(0., -step))/step**2

    def mixed(step):
        return (model.exact(step, step) - model.exact(step, -step)
                - model.exact(-step, step) + model.exact(-step, -step))/(4*step**2)

    coarse, fine = h, h/2
    return {
        (1, 0): (4*first((1, 0), fine) - first((1, 0), coarse))/3,
        (0, 1): (4*first((0, 1), fine) - first((0, 1), coarse))/3,
        (2, 0): (4*pure((2, 0), fine) - pure((2, 0), coarse))/3,
        (0, 2): (4*pure((0, 2), fine) - pure((0, 2), coarse))/3,
        (1, 1): (4*mixed(fine) - mixed(coarse))/3,
    }


def propagator_checks(model, coefficients, h=1.e-2):
    """Series coefficients against finite differences of the nonlinear flow."""
    p = model.series
    reference = fd_derivatives(model, h)
    order = {(1, 0): 1., (0, 1): 1., (2, 0): 2., (0, 2): 2., (1, 1): 1.}
    out = {}
    for key, derivative in reference.items():
        coefficient = coefficients[p.index[key]]*order[key]
        scale = max(float(np.max(np.abs(derivative))), 1e-300)
        out[f"d{key[0]}x_d{key[1]}y"] = {
            "finite_difference_max_abs": scale,
            "coefficient_max_abs": float(np.max(np.abs(coefficient))),
            "max_abs_difference": float(np.max(np.abs(coefficient - derivative))),
            "max_rel_difference": float(np.max(np.abs(coefficient - derivative))/scale),
        }
    return out


def series_truncation(model, coefficients, values, points, order):
    """Truncated series against the exact evolution at selected parameters."""
    p = model.series
    rows = []
    for x, y in points:
        exact_r = model.exact(x, y)
        exact = dilation(exact_r)
        row = {"x": float(x), "y": float(y), "E_leak_exact": exact["E_leak"],
               "leak_exact": exact["leak"], "volume_exact": exact["volume"]}
        for degree in range(2, order + 1, 2):
            r_truncated = p.evaluate(coefficients, x, y, degree)
            row[f"E_leak_degree{degree}"] = p.evaluate(values["E_leak"], x, y, degree)
            row[f"E_leak_degree{degree}_error"] = abs(
                row[f"E_leak_degree{degree}"] - exact["E_leak"])
            row[f"E_leak_truncated_propagator_degree{degree}"] = \
                dilation(r_truncated)["E_leak"]
            row[f"propagator_degree{degree}_error"] = float(
                np.max(np.abs(r_truncated - exact_r)))
        row["E_leak_series_full"] = p.evaluate(values["E_leak"], x, y, order)
        row["volume_series"] = p.evaluate(values["det_M"], x, y, order)
        row["leak_series"] = p.evaluate(values["leak"], x, y, order)
        rows.append(row)
    return rows


def observable_identities(model, points):
    """Internal consistency of the exact dilation observables."""
    worst = 0.
    rows = []
    for x, y in points:
        d = dilation(model.exact(x, y))
        # N has two columns, so only two principal sines are non-zero.
        padded = np.concatenate((d["sin_theta_from_N"],
                                 np.zeros(3 - d["sin_theta_from_N"].size)))
        residual = max(abs(d["volume"] - float(np.prod(d["cos_theta"]))),
                       abs(d["leak"] - float(np.sum(np.sin(d["theta"])**2))),
                       float(np.max(np.abs(np.sort(padded) - np.sort(np.sin(d["theta"]))))))
        worst = max(worst, residual)
        rows.append({"x": float(x), "y": float(y), "residual": residual,
                     "orthogonality_defect": d["orthogonality_defect"]})
    return worst, rows


def gradient_checks(model, values, points, order, h=1.e-2):
    """Parameter partial derivatives of E_leak: series against differences.

    The analytic identity behind the first-order result is d tr S / dx =
    tr(U^T dM/dx) with U the polar factor of M; the series result is checked
    against Richardson-extrapolated differences of the exact observable.
    """
    p = model.series
    rows = []
    for x, y in points:
        exact = lambda a, b: dilation(model.exact(a, b))["E_leak"]

        def difference(step, axis, mixed=False):
            if mixed:
                return (exact(x+step, y+step) - exact(x+step, y-step)
                        - exact(x-step, y+step) + exact(x-step, y-step))/(4*step**2)
            if axis == "x":
                return (exact(x+step, y) - exact(x-step, y))/(2*step)
            return (exact(x, y+step) - exact(x, y-step))/(2*step)

        richardson = lambda axis: (4*difference(h/2, axis) - difference(h, axis))/3
        mixed = lambda: (4*difference(h/2, "x", True) - difference(h, "x", True))/3
        d_x = series_derivative(p, values["E_leak"], "x")
        d_y = series_derivative(p, values["E_leak"], "y")
        d_xy = series_derivative(p, series_derivative(p, values["E_leak"], "x"), "y")
        rows.append({
            "x": float(x), "y": float(y),
            "dE_dx_series": float(p.evaluate(d_x, x, y, order)),
            "dE_dx_finite_difference": richardson("x"),
            "dE_dy_series": float(p.evaluate(d_y, x, y, order)),
            "dE_dy_finite_difference": richardson("y"),
            "d2E_dxdy_series": float(p.evaluate(d_xy, x, y, order)),
            "d2E_dxdy_finite_difference": mixed(),
        })
    return rows


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tau", type=float, default=1.)
    parser.add_argument("--e-scale", type=float, default=.02)
    parser.add_argument("--t1-scale", type=float, default=.03)
    parser.add_argument("--order", type=int, default=6)
    parser.add_argument("--grid", type=int, default=21)
    parser.add_argument("--map-range", type=float, default=1.,
                        help="half width of the (x, y) map; the repository scans "
                             "t1/E1 up to 10, i.e. y up to 10*t1_scale/E1")
    parser.add_argument("--physical", type=float, nargs=2, default=[.01, .01],
                        metavar=("E1", "T1"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).parent/"dilation_sensitivity")
    args = parser.parse_args()
    if min(args.tau, args.e_scale, args.t1_scale) <= 0:
        parser.error("tau and parameter scales must be positive")

    model = Model(args.tau, args.e_scale, args.t1_scale, args.order, 2e-13, 2e-15)
    p = model.series
    coefficients = propagator_coefficients(model)
    values = expansion_series(coefficients, model)

    ideal = dilation(model.exact(0., 0.))
    ideal_frame = model.exact(0., 0.)[:3, :3]
    baseline = {
        "E_leak_at_ideal_protocol": ideal["E_leak"],
        "leak_at_ideal_protocol": ideal["leak"],
        "max_abs_M_minus_B": float(np.max(np.abs(ideal_frame - B))),
        "orthogonality_defect": ideal["orthogonality_defect"],
    }

    checks = {
        "propagator_coefficients_vs_finite_differences": propagator_checks(model, coefficients),
        "observable_identities": {},
        "parity_of_expansion_series": {},
    }
    worst_identity, identity_rows = observable_identities(
        model, [(0., 0.), (args.e_scale, args.t1_scale), (-args.e_scale, 0.)])
    checks["observable_identities"] = {
        "worst_residual": worst_identity, "rows": identity_rows,
        "scope": "det M = prod cos theta; ||N||_F^2 = sum sin^2 theta; N singular values.",
    }
    odd = [k for k, (m, n) in enumerate(p.keys) if (m + n) % 2]
    checks["parity_of_expansion_series"] = {
        "max_abs_odd_degree_coefficient": float(np.max(np.abs(values["E_leak"][odd]))),
        "note": "(E1, t1) -> (-E1, -t1) is a Majorana relabelling; E_leak must be even.",
    }

    test_points = [(args.e_scale, args.t1_scale), (0.1, 0.1), (-0.15, 0.2),
                   (0.5, 1/3), (0.5, 1/30), (0.05, 1/3),
                   (0.5, 10/3), (0.5, -10/3), (1., 10/3), (2., 2.)]
    checks["parameter_gradients_vs_finite_differences"] = gradient_checks(
        model, values, test_points, args.order)
    truncation = series_truncation(model, coefficients, values, test_points, args.order)

    # Parameter map: exact against the truncated series.
    n = args.grid
    xs = np.linspace(-args.map_range, args.map_range, n)
    ys = np.linspace(-args.map_range, args.map_range, n)
    exact_e = np.empty((n, n))
    degree2 = np.empty((n, n))
    degree4 = np.empty((n, n))
    degree6 = np.empty((n, n))
    leak = np.empty((n, n))
    volume = np.empty((n, n))
    theta = np.empty((3, n, n))
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            d = dilation(model.exact(x, y))
            exact_e[i, j] = d["E_leak"]
            leak[i, j] = d["leak"]
            volume[i, j] = d["volume"]
            theta[:, i, j] = d["theta"]
            degree2[i, j] = p.evaluate(values["E_leak"], x, y, 2)
            degree4[i, j] = p.evaluate(values["E_leak"], x, y, 4)
            degree6[i, j] = p.evaluate(values["E_leak"], x, y, 6)
        if j % 5 == 0:
            print(f"map column {j+1}/{n}", flush=True)

    # Time-resolved expansion at the physical point.
    px, py = args.physical[0]/args.e_scale, args.physical[1]/args.t1_scale
    n_out = 40

    def exact_rhs(segment, s, vector):
        o, oe, og = model.generators(segment, s)
        return ((o + px*oe + py*og) @ vector.reshape(5, 5)).ravel()

    _, times, states = integrate(model, exact_rhs, np.eye(5), n_out=n_out)

    def coefficient_rhs(segment, s, vector):
        r = vector.reshape(p.count, 5, 5)
        o, oe, og = model.generators(segment, s)
        result = o @ r
        result[[k for k, (m, _) in enumerate(p.keys) if m]] += \
            oe @ r[[p.index[(m-1, nn)] for m, nn in p.keys if m]]
        result[[k for k, (_, nn) in enumerate(p.keys) if nn]] += \
            og @ r[[p.index[(m, nn-1)] for m, nn in p.keys if nn]]
        return result.ravel()

    _, _, c_traj = integrate(model, coefficient_rhs,
                             np.concatenate((np.eye(5).ravel(),
                                             np.zeros((p.count-1)*25))), n_out=n_out)

    timeline = {"u": times, "E_leak_exact": [], "E_leak_degree2": [],
                "E_leak_degree4": [], "E_leak_degree6": [],
                "dE_leak_dx": [], "dE_leak_dy": []}
    for k in range(times.size):
        c = c_traj[:, k].reshape(p.count, 5, 5)
        local = expansion_series(c, model)
        timeline["E_leak_exact"].append(dilation(states[:, k].reshape(5, 5))["E_leak"])
        for degree in (2, 4, 6):
            timeline[f"E_leak_degree{degree}"].append(
                p.evaluate(local["E_leak"], px, py, degree))
        timeline["dE_leak_dx"].append(p.evaluate(
            series_derivative(p, local["E_leak"], "x"), px, py, args.order))
        timeline["dE_leak_dy"].append(p.evaluate(
            series_derivative(p, local["E_leak"], "y"), px, py, args.order))
    timeline = {key: np.asarray(value) for key, value in timeline.items()}

    origin_gradient = {
        axis: float(p.evaluate(series_derivative(p, values["E_leak"], axis), 0., 0.))
        for axis in ("x", "y")
    }
    physical_row = truncation[3]
    result = {
        "parameters": {
            "tau_meV_inverse": args.tau, "total_time_meV_inverse": 6*args.tau,
            "E1_scale_meV": args.e_scale, "t1_scale_meV": args.t1_scale,
            "order": args.order, "rtol": model.rtol, "atol": model.atol,
            "physical_E1_meV": args.physical[0], "physical_t1_meV": args.physical[1],
            "physical_x": px, "physical_y": py,
        },
        "ideal_protocol_baseline": baseline,
        "checks": checks,
        "origin_gradient": origin_gradient,
        "expansion_series_physical_units": {
            name: physical_powers(p, poly, args.e_scale, args.t1_scale)
            for name, poly in values.items()},
        "physical_point": physical_row,
        "truncation_study": truncation,
        "scope": ("Fixed three-stage protocol; expansion about the ideal protocol "
                  "in the smooth positive-polar branch; coefficients from time "
                  "integration, not parameter fitting."),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    np.savez(args.output_dir/"maps.npz", xs=xs, ys=ys, E_leak_exact=exact_e,
             E_leak_degree2=degree2, E_leak_degree4=degree4, E_leak_degree6=degree6,
             leak=leak, volume=volume, theta=theta, timeline_u=timeline["u"],
             timeline_E_leak_exact=timeline["E_leak_exact"],
             timeline_E_leak_degree2=timeline["E_leak_degree2"],
             timeline_E_leak_degree4=timeline["E_leak_degree4"],
             timeline_E_leak_degree6=timeline["E_leak_degree6"],
             timeline_dE_leak_dx=timeline["dE_leak_dx"],
             timeline_dE_leak_dy=timeline["dE_leak_dy"],
             E1_scale=args.e_scale, t1_scale=args.t1_scale, tau=args.tau)

    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), constrained_layout=True)
    extent = (xs[0], xs[-1], ys[0], ys[-1])
    base = exact_e[0, 0]
    im = axes[0].imshow((exact_e - base).T, origin="lower", extent=extent, cmap="magma")
    axes[0].set(title=rf"exact $\Delta E_{{\rm leak}}$ (baseline ${base:.4g}$)",
                xlabel=r"$x=E_1/0.02$", ylabel=r"$y=t_1/0.03$")
    fig.colorbar(im, ax=axes[0])
    for ax, field, label in ((axes[1], degree2, "degree 2"), (axes[2], degree4, "degree 4")):
        im = ax.imshow(np.abs(field - exact_e).T, origin="lower", extent=extent,
                       cmap="viridis")
        ax.set(title=rf"$|{label} - $ exact$|$", xlabel=r"$x$", ylabel=r"$y$")
        fig.colorbar(im, ax=ax)
    fig.suptitle(rf"subspace dilation, $\tau={args.tau:g}$ meV$^{{-1}}$, "
                 rf"series degree {args.order}")
    fig.savefig(args.output_dir/"dilation_map.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    axes[0].plot(timeline["u"], timeline["E_leak_exact"], "k", lw=1.8, label="exact")
    for degree in (2, 4, 6):
        axes[0].plot(timeline["u"], timeline[f"E_leak_degree{degree}"], "--",
                     label=f"series degree {degree}")
    axes[0].set(title=rf"$E_{{\rm leak}}(u)$ at $E_1={args.physical[0]:g}$, $t_1={args.physical[1]:g}$ meV",
                xlabel=r"$u=t/\tau$", ylabel=r"$E_{\rm leak}$")
    axes[0].grid(alpha=.25)
    axes[0].legend()
    axes[1].plot(timeline["u"], timeline["dE_leak_dx"], label=r"$\partial E_{\rm leak}/\partial x$")
    axes[1].plot(timeline["u"], timeline["dE_leak_dy"], label=r"$\partial E_{\rm leak}/\partial y$")
    axes[1].axhline(0., color="k", lw=.6)
    axes[1].set(title="expansion sensitivity along the protocol", xlabel=r"$u=t/\tau$",
                ylabel="partial derivative")
    axes[1].grid(alpha=.25)
    axes[1].legend()
    fig.savefig(args.output_dir/"dilation_timeline.png", dpi=200)
    plt.close(fig)

    print(json.dumps({"ideal_protocol_baseline": baseline,
                      "checks": {k: v for k, v in checks.items()
                                 if k != "propagator_coefficients_vs_finite_differences"},
                      "origin_gradient": origin_gradient,
                      "physical_point": physical_row}, indent=2), flush=True)
    print("\npropagator coefficients vs finite differences:", flush=True)
    for name, row in checks["propagator_coefficients_vs_finite_differences"].items():
        print(f"  {name:14s} max|FD|={row['finite_difference_max_abs']:.6e} "
              f"max|diff|={row['max_abs_difference']:.3e} "
              f"rel={row['max_rel_difference']:.3e}", flush=True)
    print("\nparameter gradients vs finite differences:", flush=True)
    for row in checks["parameter_gradients_vs_finite_differences"]:
        print(f"  x={row['x']:.4f} y={row['y']:.4f} "
              f"dE/dx={row['dE_dx_series']:+.6e} (fd {row['dE_dx_finite_difference']:+.6e}) "
              f"dE/dy={row['dE_dy_series']:+.6e} (fd {row['dE_dy_finite_difference']:+.6e}) "
              f"d2E/dxdy={row['d2E_dxdy_series']:+.6e} "
              f"(fd {row['d2E_dxdy_finite_difference']:+.6e})", flush=True)
    print("\ntruncation study:", flush=True)
    for row in truncation:
        print(f"  x={row['x']:.4f} y={row['y']:.4f} E_leak={row['E_leak_exact']:.8e} "
              f"err2={row['E_leak_degree2_error']:.3e} "
              f"err4={row['E_leak_degree4_error']:.3e} "
              f"err6={row['E_leak_degree6_error']:.3e} "
              f"volume={row['volume_exact']:+.6f} (series {row['volume_series']:+.6f})",
              flush=True)
    print("\nE_leak series (physical units, per meV powers):",
          json.dumps(result["expansion_series_physical_units"]["E_leak"], indent=2), flush=True)
    print(f"\nOutputs: {args.output_dir.resolve()}", flush=True)


if __name__ == "__main__":
    main()
