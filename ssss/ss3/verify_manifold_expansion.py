"""Fourth-order parameter expansion in a moving Grassmann chart.

Retained as a model/series utility and local coefficient cross-check.
The current separate-error analysis is in verify_error_channels.py and n6.md.
These local polynomials are not a validated global braiding map.

No fit to an E1,t1 grid is used. A triangular Riccati coefficient hierarchy
produces a bivariate polynomial, including the polar factor as a matrix series.
Independent linear SO(5) coefficient and full propagator solves validate it.

Units: tau in meV^-1, E1 and t1 in meV, hbar=1. Two three-stage cycles.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_sylvester


B = np.diag([1., -1., -1.])
W0 = np.eye(5)[:, :3]


class Series:
    """Ordinary Taylor coefficients, with factorials already divided out."""

    def __init__(self, order):
        self.order = order
        self.keys = [(d-n, n) for d in range(order+1) for n in range(d+1)]
        self.index = {key: k for k, key in enumerate(self.keys)}
        self.count = len(self.keys)
        self.pairs = [
            [(self.index[(i, j)], self.index[(m-i, n-j)])
             for i in range(m+1) for j in range(n+1)]
            for m, n in self.keys
        ]

    def matmul(self, a, b):
        out = np.zeros((self.count, a.shape[1], b.shape[2]))
        for k, terms in enumerate(self.pairs):
            for i, j in terms:
                out[k] += a[i] @ b[j]
        return out

    def scalar_product(self, a, b):
        return np.array([sum(a[i]*b[j] for i, j in terms)
                         for terms in self.pairs])

    def sqrt_spd(self, a):
        vals, vecs = np.linalg.eigh(a[0])
        if vals.min() <= 0:
            raise ValueError("Singular polar base point: change expansion/chart.")
        h = np.zeros_like(a)
        h[0] = (vecs * np.sqrt(vals)) @ vecs.T
        for k in range(1, self.count):
            rhs = a[k].copy()
            for i, j in self.pairs[k]:
                if i and j:
                    rhs -= h[i] @ h[j]
            h[k] = solve_sylvester(h[0], h[0], rhs)
        return h

    def inverse(self, a):
        out = np.zeros_like(a)
        out[0] = np.linalg.inv(a[0])
        for k in range(1, self.count):
            rhs = sum((a[i] @ out[j] for i, j in self.pairs[k] if i),
                      np.zeros_like(a[0]))
            out[k] = -out[0] @ rhs
        return out

    def evaluate(self, a, x, y, order=None):
        order = self.order if order is None else order
        return sum((a[k] * x**m * y**n for k, (m, n) in enumerate(self.keys)
                    if m+n <= order), np.zeros_like(a[0]))


class Model:
    def __init__(self, tau, e_scale, g_scale, order, rtol, atol):
        self.tau, self.e_scale, self.g_scale = tau, e_scale, g_scale
        self.series = Series(order)
        self.rtol, self.atol = rtol, atol

    def generators(self, segment, s):
        """dR/ds; segment-local s in [0,1]. Energy normalization is 2."""
        rise = (1-np.cos(np.pi*s))/2
        fall = 1-rise
        if segment % 3 == 0:
            t2, t3, ed, envelope = .3*rise, 0., .3*fall, rise
        elif segment % 3 == 1:
            t2, t3, ed, envelope = .3*fall, .3*rise, 0., fall
        else:
            t2, t3, ed, envelope = 0., .3*fall, .3*rise, 0.
        o, oe, og = np.zeros((3, 5, 5))
        o[3, 1], o[1, 3] = 2*t2, -2*t2
        o[3, 2], o[2, 3] = -2*t3, 2*t3
        o[3, 4], o[4, 3] = 2*ed, -2*ed
        oe[0, 1], oe[1, 0] = 2*self.e_scale, -2*self.e_scale
        og[0, 4], og[4, 0] = 2*self.g_scale*envelope, -2*self.g_scale*envelope
        return self.tau*o, self.tau*oe, self.tau*og

    def integrate(self, rhs, y0):
        y = y0.ravel()
        for segment in range(6):
            sol = solve_ivp(lambda s, v: rhs(segment, s, v), (0., 1.), y,
                            method="DOP853", rtol=self.rtol, atol=self.atol)
            if not sol.success:
                raise RuntimeError(sol.message)
            y = sol.y[:, -1]
        return y

    def chart_coefficients(self):
        """Integrate R0, graph Z, and top frame X as polynomial coefficients."""
        p = self.series
        z = np.zeros((p.count, 2, 3))
        x = np.zeros((p.count, 3, 3))
        x[0] = np.eye(3)
        z_end = 25 + z.size

        def rhs(segment, s, y):
            r0 = y[:25].reshape(5, 5)
            z = y[25:z_end].reshape(p.count, 2, 3)
            x = y[z_end:].reshape(p.count, 3, 3)
            o, oe, og = self.generators(segment, s)
            dz, dx = np.zeros_like(z), np.zeros_like(x)
            # V = x_parameter*A + y_parameter*C has only degree-one terms.
            for unit, source in [((1, 0), oe), ((0, 1), og)]:
                a = r0.T @ source @ r0
                a11, a12, a21, a22 = a[:3, :3], a[:3, 3:], a[3:, :3], a[3:, 3:]
                for k, (m, n) in enumerate(p.keys):
                    beta = (m-unit[0], n-unit[1])
                    if beta not in p.index:
                        continue
                    b = p.index[beta]
                    dz[k] += a22 @ z[b] - z[b] @ a11
                    dx[k] += a11 @ x[b]
                    if b == 0:
                        dz[k] += a21
                    for i, j in p.pairs[b]:
                        dz[k] -= z[i] @ a12 @ z[j]
                        dx[k] += a12 @ z[i] @ x[j]
            return np.concatenate(((o @ r0).ravel(), dz.ravel(), dx.ravel()))

        result = self.integrate(rhs, np.concatenate((np.eye(5).ravel(), z.ravel(), x.ravel())))
        r0 = result[:25].reshape(5, 5)
        z = result[25:z_end].reshape(p.count, 2, 3)
        x = result[z_end:].reshape(p.count, 3, 3)
        frame_rel = np.concatenate((x, p.matmul(z, x)), axis=1)
        frame_abs = np.einsum("ij,kjl->kil", r0, frame_rel)
        return r0, z, x, frame_rel, frame_abs

    def linear_coefficients(self):
        """Independent linear hierarchy for the absolute SO(5) frame."""
        p = self.series
        initial = np.zeros((p.count, 5, 3))
        initial[0] = W0

        def rhs(segment, s, y):
            w = y.reshape(p.count, 5, 3)
            o, oe, og = self.generators(segment, s)
            dw = np.einsum("ij,kjl->kil", o, w)
            for k, (m, n) in enumerate(p.keys):
                if m:
                    dw[k] += oe @ w[p.index[(m-1, n)]]
                if n:
                    dw[k] += og @ w[p.index[(m, n-1)]]
            return dw.ravel()
        return self.integrate(rhs, initial).reshape(p.count, 5, 3)

    def exact(self, x, y):
        """Full propagator, without SVD reorthogonalization."""
        def rhs(segment, s, vector):
            o, oe, og = self.generators(segment, s)
            return ((o+x*oe+y*og) @ vector.reshape(5, 5)).ravel()
        return self.integrate(rhs, np.eye(5)).reshape(5, 5)

    def exact_chart(self, x, y):
        """Direct nonlinear Z,U evolution: no defective full R is propagated."""
        def rhs(segment, s, vector):
            r0 = vector[:25].reshape(5, 5)
            z = vector[25:31].reshape(2, 3)
            u = vector[31:].reshape(3, 3)
            if np.linalg.norm(z) > 1e3:
                raise ValueError("Grassmann chart nearly singular; use smaller parameters or another chart.")
            o, oe, og = self.generators(segment, s)
            v = r0.T @ (x*oe+y*og) @ r0
            zdot = v[3:, :3]+v[3:, 3:]@z-z@v[:3, :3]-z@v[:3, 3:]@z
            eig, vec = np.linalg.eigh(np.eye(3)+z.T@z)
            root = (vec*np.sqrt(eig)) @ vec.T
            rootdot = solve_sylvester(root, root, zdot.T@z+z.T@zdot)
            j = (rootdot+root@(v[:3, :3]+v[:3, 3:]@z)) @ np.linalg.inv(root)
            return np.concatenate(((o@r0).ravel(), zdot.ravel(), (j@u).ravel()))
        initial = np.concatenate((np.eye(5).ravel(), np.zeros(6), np.eye(3).ravel()))
        result = self.integrate(rhs, initial)
        r0 = result[:25].reshape(5, 5)
        z, u = result[25:31].reshape(2, 3), result[31:].reshape(3, 3)
        eig, vec = np.linalg.eigh(np.eye(3)+z.T@z)
        root = (vec*np.sqrt(eig)) @ vec.T
        m = np.linalg.solve(root, u)
        frame = np.vstack((m, z@m))
        return r0, z, u, frame


def metric_series(p, frame, target):
    m = frame[:, :3, :]
    if np.linalg.det(m[0]) <= 0:
        raise ValueError("Base polar factor is not in SO(3); this local expansion is not applicable.")
    h = p.sqrt_spd(p.matmul(m, m.transpose(0, 2, 1)))
    u = p.matmul(p.inverse(h), m)
    surv = np.trace(p.matmul(m.transpose(0, 2, 1), m), axis1=1, axis2=2)/3
    rot = np.einsum("ij,kij->k", target, u)/4
    rot[0] += .25
    return {"F_surv": surv, "F_rot": rot, "F_group": p.scalar_product(surv, rot)}


def metrics(frame, target):
    m = frame[:3]
    left, values, vh = np.linalg.svd(m)
    u = left@vh
    if np.linalg.det(u) <= 0 or values.min() < 1e-10:
        raise ValueError("Singular or orientation-reversing block outside tested smooth polar branch.")
    surv = float(np.sum(m*m)/3)
    rot = float((1+np.trace(target.T@u))/4)
    return {"F_surv": surv, "F_rot": rot, "F_group": surv*rot}


def serialize_coefficients(p, values, e_scale, g_scale):
    return {
        name: {f"{m},{n}": float(a[k]/e_scale**m/g_scale**n)
               for k, (m, n) in enumerate(p.keys)}
        for name, a in values.items()
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tau", type=float, default=1.)
    parser.add_argument("--e-scale", type=float, default=.02)
    parser.add_argument("--t1-scale", type=float, default=.03)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).parent/"manifold_results_tau1")
    args = parser.parse_args()
    if min(args.tau, args.e_scale, args.t1_scale) <= 0:
        parser.error("tau and parameter scales must be positive")
    model = Model(args.tau, args.e_scale, args.t1_scale, 4, 2e-13, 2e-15)
    loose = Model(args.tau, args.e_scale, args.t1_scale, 4, 2e-11, 2e-13)
    p = model.series
    r0, z, x, relative, absolute = model.chart_coefficients()
    linear = model.linear_coefficients()
    # Nonlinear chart and linear equations are different coefficient routes.
    coefficients = {"absolute": metric_series(p, absolute, B),
                    "relative": metric_series(p, relative, np.eye(3))}
    other_coefficients = metric_series(p, linear, B)
    coeff_difference = float(np.max(np.abs(absolute-linear)))
    metric_difference = max(float(np.max(np.abs(coefficients["absolute"][k]-v)))
                            for k, v in other_coefficients.items())
    gram = p.matmul(relative.transpose(0, 2, 1), relative)
    gram[0] -= np.eye(3)
    # The graph formula for survival is independently assembled as a series.
    a = p.matmul(z.transpose(0, 2, 1), z)
    a[0] += np.eye(3)
    graph_surv = np.trace(p.inverse(a), axis1=1, axis2=2)/3
    surv_identity = float(np.max(np.abs(graph_surv-coefficients["relative"]["F_surv"])))
    # An invariant fourth-order formula, independent of expanding trace(U):
    # 1-F = ||Z||^2/3 + ||a||^2/4 - tr((Z^T Z)^2)/3
    #       - ||a||^4/48 - ||Z||^2 ||a||^2/12 + O((||Z||+||a||)^6).
    # Here hat(a)=log(U), U=sqrt(I+Z^T Z) X, all as formal series.
    u_series = p.matmul(p.sqrt_spd(a), x)
    u_minus_i = u_series.copy()
    u_minus_i[0] -= np.eye(3)
    log_u, power = np.zeros_like(u_series), u_minus_i.copy()
    for degree in range(1, p.order+1):
        log_u += (-1)**(degree+1)*power/degree
        power = p.matmul(power, u_minus_i)
    rotation_vector = np.stack((log_u[:, 2, 1], log_u[:, 0, 2], log_u[:, 1, 0]), axis=1)
    spin_squared = sum(p.scalar_product(rotation_vector[:, j], rotation_vector[:, j])
                       for j in range(3))
    ztz = a.copy()
    ztz[0] -= np.eye(3)
    z_squared = np.trace(ztz, axis1=1, axis2=2)
    z_fourth = np.trace(p.matmul(ztz, ztz), axis1=1, axis2=2)
    invariant_error = (z_squared/3+spin_squared/4-z_fourth/3
                       -p.scalar_product(spin_squared, spin_squared)/48
                       -p.scalar_product(z_squared, spin_squared)/12)
    group_error = -coefficients["relative"]["F_group"].copy()
    group_error[0] += 1
    invariant_residual = float(np.max(np.abs(invariant_error-group_error)))
    odd = [k for k, (m, n) in enumerate(p.keys) if (m+n) % 2]
    odd_residual = max(float(np.max(np.abs(v[odd])))
                       for group in coefficients.values() for v in group.values())
    r0_direct = model.exact(0., 0.)
    r = model.exact(1., 1.)
    rg0, zg, ug, wg = model.exact_chart(1., 1.)
    relative_direct = r0_direct.T@r@W0
    chart_error = float(np.linalg.norm(wg-relative_direct))
    flipped = model.exact(-1., -1.)
    sign = np.diag([-1., 1., 1., 1., 1.])
    sign_residual = float(np.linalg.norm(flipped-sign@r@sign))
    orthogonality = float(np.linalg.norm(r.T@r-np.eye(5)))
    checks = {
        "scaled_frame_coefficient_max_difference": coeff_difference,
        "scaled_fidelity_coefficient_max_difference": metric_difference,
        "frame_gram_series_residual": float(np.max(np.abs(gram))),
        "graph_survival_series_identity_residual": surv_identity,
        "invariant_relative_infidelity_degree4_residual": invariant_residual,
        "relative_rotation_log_skew_residual": float(np.max(np.abs(log_u+log_u.transpose(0, 2, 1)))),
        "odd_total_degree_fidelity_coefficient_max": odd_residual,
        "nonlinear_Z_U_vs_full_SO5_frame_norm": chart_error,
        "absolute_frame_reconstruction_norm": float(np.linalg.norm(rg0@wg-r@W0)),
        "raw_SO5_orthogonality_norm": orthogonality,
        "raw_U_orthogonality_norm": float(np.linalg.norm(ug.T@ug-np.eye(3))),
        "simultaneous_sign_reversal_norm": sign_residual,
        "relative_endpoint_Z_norm": float(np.linalg.norm(zg)),
        "absolute_base_M_min_singular_value": float(np.linalg.svd(absolute[0, :3], compute_uv=False).min()),
    }
    for label in ("scaled_frame_coefficient_max_difference",
                  "scaled_fidelity_coefficient_max_difference",
                  "frame_gram_series_residual", "graph_survival_series_identity_residual",
                  "invariant_relative_infidelity_degree4_residual",
                  "relative_rotation_log_skew_residual",
                  "nonlinear_Z_U_vs_full_SO5_frame_norm",
                  "absolute_frame_reconstruction_norm", "raw_SO5_orthogonality_norm",
                  "raw_U_orthogonality_norm", "simultaneous_sign_reversal_norm"):
        if checks[label] > 2e-8:
            raise AssertionError(f"{label}: {checks[label]}")
    scales = [1., .75, .5, .375, .25]
    scan = []
    loose_r0 = loose.exact(0., 0.)
    for s in scales:
        exact = model.exact(s, s)
        less_accurate = loose.exact(s, s)
        entry = {"scale": s, "E1_meV": s*args.e_scale, "t1_meV": s*args.t1_scale}
        for mode, frame, alt_frame, target in [
            ("absolute", exact@W0, less_accurate@W0, B),
            ("relative", r0_direct.T@exact@W0, loose_r0.T@less_accurate@W0, np.eye(3)),
        ]:
            truth, alt = metrics(frame, target), metrics(alt_frame, target)
            entry[mode] = {}
            for name, polynomial in coefficients[mode].items():
                f2 = float(p.evaluate(polynomial, s, s, 2))
                f4 = float(p.evaluate(polynomial, s, s, 4))
                entry[mode][name] = {
                    "full": truth[name], "degree2": f2, "degree4": f4,
                    "error2": abs(f2-truth[name]), "error4": abs(f4-truth[name]),
                    "tolerance_sensitivity": abs(truth[name]-alt[name]),
                }
        scan.append(entry)
    points = {(a, b): model.exact(a, b)@W0 for a, b in [(0, 0), (1, 0), (0, 1), (1, 1)]}
    mixed = {}
    for mode, target in [("absolute", B), ("relative", np.eye(3))]:
        transform = np.eye(5) if mode == "absolute" else r0_direct.T
        vals = {key: metrics(transform@w, target) for key, w in points.items()}
        mixed[mode] = {}
        for name, polynomial in coefficients[mode].items():
            delta = vals[(1, 1)][name]-vals[(1, 0)][name]-vals[(0, 1)][name]+vals[(0, 0)][name]
            second = float(polynomial[p.index[(1, 1)]])
            fourth = sum(float(polynomial[p.index[key]]) for key in [(3, 1), (2, 2), (1, 3)])
            mixed[mode][name] = {"full_four_point": delta, "degree2": second,
                                "degree4_correction": fourth, "degree4": second+fourth}
    result = {
        "parameters": {"tau_meV_inverse": args.tau, "total_time_meV_inverse": 6*args.tau,
                       "tc_meV": .3, "Ed_max_meV": .3, "normalization": 2,
                       "E1_scale_meV": args.e_scale, "t1_scale_meV": args.t1_scale,
                       "polynomial_degree": 4, "rtol": model.rtol, "atol": model.atol},
        "checks": checks,
        "coefficients_physical_meV_powers": {
            mode: serialize_coefficients(p, v, args.e_scale, args.t1_scale)
            for mode, v in coefficients.items()},
        "scaling": scan, "mixed_response": mixed,
        "scope": "Fixed protocol; local smooth chart/polar branch; coefficients use time integration, not parameter fitting.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    np.savez(args.output_dir/"coefficients.npz", powers=p.keys, R0=r0, Z=z, X=x,
             relative_frame=relative, absolute_frame=absolute,
             E1_scale=args.e_scale, t1_scale=args.t1_scale, tau=args.tau)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, mode in zip(axes, ("absolute", "relative")):
        for order in (2, 4):
            err = [row[mode]["F_group"][f"error{order}"] for row in scan]
            ax.loglog(scales, np.maximum(err, 1e-17), "o-", label=f"degree {order}")
        floor = [row[mode]["F_group"]["tolerance_sensitivity"] for row in scan]
        ax.loglog(scales, np.maximum(floor, 1e-17), "k:", label="tolerance sensitivity")
        ax.set(title=f"{mode.capitalize()} group fidelity", xlabel="parameter scale s", ylabel="absolute error")
        ax.grid(alpha=.25, which="both")
        ax.legend()
    fig.suptitle(rf"$\tau={args.tau:g}$ meV$^{{-1}}$, $E_1={args.e_scale:g}s$, $t_1={args.t1_scale:g}s$ meV")
    fig.savefig(args.output_dir/"convergence.png", dpi=200)
    plt.close(fig)
    print(json.dumps(checks, indent=2), flush=True)
    for mode in ("absolute", "relative"):
        print(f"\n{mode} F_group: s, full, error2, error4, tolerance sensitivity", flush=True)
        for row in scan:
            v = row[mode]["F_group"]
            print(f"{row['scale']:.3f} {v['full']:.12g} {v['error2']:.6e} {v['error4']:.6e} {v['tolerance_sensitivity']:.3e}")
        print("F_group physical coefficients:", result["coefficients_physical_meV_powers"][mode]["F_group"])
        print("F_group mixed response:", mixed[mode]["F_group"])
    print(f"\nOutputs: {args.output_dir.resolve()}", flush=True)


if __name__ == "__main__":
    main()
