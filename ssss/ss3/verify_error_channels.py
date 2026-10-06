"""Verify how higher Lie-generator orders affect leakage and rotation separately.

Keeps the original two group-space errors. No state overlap, parameter fitting,
or global extrapolation of a scalar fidelity polynomial is used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm, logm, polar

from verify_manifold_expansion import B, Model, Series


def bracket(a, b):
    return a @ b-b @ a


def diagonal_part(a):
    result = a.copy()
    result[:3, 3:] = 0
    result[3:, :3] = 0
    return result


def homogeneous(p, a, degree):
    result = np.zeros_like(a)
    for k, (m, n) in enumerate(p.keys):
        if m+n == degree:
            result[k] = a[k]
    return result


def inner_series(p, a, b):
    return np.trace(p.matmul(a.transpose(0, 2, 1), b), axis1=1, axis2=2)


def exponential_series(p, k):
    result = np.zeros_like(k)
    result[0] = np.eye(k.shape[1])
    term = result.copy()
    for degree in range(1, p.order+1):
        term = p.matmul(term, k)/degree
        result += term
    return result


def logarithm_series(p, s):
    delta = s.copy()
    delta[0] -= np.eye(s.shape[1])
    result, term = np.zeros_like(s), delta.copy()
    for degree in range(1, p.order+1):
        result += (-1)**(degree+1)*term/degree
        term = p.matmul(term, delta)
    return result


def frame_error_series(p, frame, target):
    m, n = frame[:, :3], frame[:, 3:]
    if np.linalg.det(m[0]) <= 0:
        raise ValueError("Outside the smooth positive-polar SO(3) branch.")
    h = p.sqrt_spd(p.matmul(m, m.transpose(0, 2, 1)))
    u = p.matmul(p.inverse(h), m)
    leak = inner_series(p, n, n)/3
    rot = -np.einsum("ij,kij->k", target, u)/4
    rot[0] += .75
    return {"leakage": leak, "rotation": rot}, u


def frame_errors(frame, target):
    m, n = frame[:3], frame[3:]
    left, singular, vh = np.linalg.svd(m)
    u = left @ vh
    if np.linalg.det(u) <= 0 or singular.min() < 1e-9:
        raise ValueError("Outside the tested smooth polar branch.")
    return {"leakage": float(np.sum(n*n)/3),
            "rotation": float((3-np.trace(target.T @ u))/4)}


def generator_error_terms(p, k):
    """Homogeneous contributions through parameter degree four, at S(0)=I."""
    a1, a2, a3 = [homogeneous(p, k[:, :3, :3], d) for d in (1, 2, 3)]
    c1, c2, c3 = [homogeneous(p, k[:, 3:, :3], d) for d in (1, 2, 3)]
    d1 = homogeneous(p, k[:, 3:, 3:], 1)
    ca, dc = p.matmul(c1, a1), p.matmul(d1, c1)
    p1 = p.matmul(c1.transpose(0, 2, 1), c1)
    a_norm = inner_series(p, a1, a1)/2  # ||vee(A1)||^2
    leak = {
        "degree2_C1_norm": inner_series(p, c1, c1)/3,
        "degree3_C1_C2": 2*inner_series(p, c1, c2)/3,
        "degree4_C2_norm": inner_series(p, c2, c2)/3,
        "degree4_C1_C3": 2*inner_series(p, c1, c3)/3,
        "degree4_internal_transport": -inner_series(p, dc-ca, dc-ca)/36,
        "degree4_sine_saturation": -np.trace(p.matmul(p1, p1), axis1=1, axis2=2)/9,
    }
    rot = {
        "degree2_A1_norm": inner_series(p, a1, a1)/8,
        "degree3_A1_A2": inner_series(p, a1, a2)/4,
        "degree4_A2_norm": inner_series(p, a2, a2)/8,
        "degree4_A1_A3": inner_series(p, a1, a3)/4,
        "degree4_rotation_saturation": -p.scalar_product(a_norm, a_norm)/48,
        "degree4_leakage_feedback": inner_series(p, ca, ca-dc)/24,
    }
    # Polar-rotation logarithm, not merely the pp block of log(S).
    a, c, d = k[:, :3, :3], k[:, 3:, :3], k[:, 3:, 3:]
    cp = p.matmul(c.transpose(0, 2, 1), c)
    polar_log = (a+(p.matmul(a, cp)+p.matmul(cp, a))/12
                 -p.matmul(p.matmul(c.transpose(0, 2, 1), d), c)/6)
    return {"leakage": leak, "rotation": rot}, polar_log


def propagator_coefficients(model):
    p = model.series
    initial = np.zeros((p.count, 5, 5))
    initial[0] = np.eye(5)
    de = np.array([j for j, (m, n) in enumerate(p.keys) if m])
    se = np.array([p.index[(m-1, n)] for m, n in p.keys if m])
    dg = np.array([j for j, (m, n) in enumerate(p.keys) if n])
    sg = np.array([p.index[(m, n-1)] for m, n in p.keys if n])

    def rhs(segment, s, vector):
        r = vector.reshape(p.count, 5, 5)
        o, oe, og = model.generators(segment, s)
        result = o @ r
        result[de] += oe @ r[se]
        result[dg] += og @ r[sg]
        return result.ravel()

    return model.integrate(rhs, initial).reshape(p.count, 5, 5)


def magnus_through_three(model):
    """Independent Magnus ODE; separately track hh, mm, hm sources of K2."""
    p = model.series
    size = p.count*25
    units = ((1, 0), (0, 1))
    initial = np.concatenate((np.eye(5).ravel(), np.zeros(4*size)))

    def rhs(segment, s, vector):
        r0 = vector[:25].reshape(5, 5)
        k, hh, mm, hm = vector[25:].reshape(4, p.count, 5, 5)
        o, oe, og = model.generators(segment, s)
        vs = (r0.T @ oe @ r0, r0.T @ og @ r0)
        dk, dhh, dmm, dhm = np.zeros((4, p.count, 5, 5))
        for unit, v in zip(units, vs):
            dk[p.index[unit]] += v
            hv = diagonal_part(v)
            mv = v-hv
            for key in p.keys:
                degree = sum(key)
                if degree not in (2, 3):
                    continue
                rest = (key[0]-unit[0], key[1]-unit[1])
                if rest not in p.index:
                    continue
                dest, prev = p.index[key], p.index[rest]
                dk[dest] -= bracket(k[prev], v)/2
                if degree == 2:
                    hk = diagonal_part(k[prev])
                    mk = k[prev]-hk
                    dhh[dest] += bracket(hv, hk)/2
                    dmm[dest] += bracket(mv, mk)/2
                    dhm[dest] += (bracket(hv, mk)+bracket(mv, hk))/2
                else:
                    for first in units:
                        second = (rest[0]-first[0], rest[1]-first[1])
                        if second in units:
                            dk[dest] += bracket(k[p.index[first]],
                                                bracket(k[p.index[second]], v))/12
        return np.concatenate(((o @ r0).ravel(), np.stack((dk, dhh, dmm, dhm)).ravel()))

    result = model.integrate(rhs, initial)
    return result[25:].reshape(4, p.count, 5, 5)


def algebra_checks():
    """Check general A,C,D, including nonzero auxiliary D at first order."""
    p = Series(4)
    rng = np.random.default_rng(42)
    k = np.zeros((p.count, 5, 5))
    for j, key in enumerate(p.keys):
        if 1 <= sum(key) <= 3:
            raw = rng.normal(size=(5, 5))
            k[j] = .15*(raw-raw.T)
    s = exponential_series(p, k)
    errors, u = frame_error_series(p, s[:, :, :3], np.eye(3))
    terms, polar_log = generator_error_terms(p, k)
    residual = max(float(np.max(abs(sum(terms[name].values())-errors[name]))) for name in errors)
    polar_residual = float(np.max(abs(logarithm_series(p, u)-polar_log)))
    # Independent finite-matrix check of the O(K^6) and O(K^5) remainders.
    raw = rng.normal(size=(5, 5))
    direction = (raw-raw.T)/2
    scaling = []
    for scale in (.2, .1, .05):
        matrix = scale*direction
        a, c, d = matrix[:3, :3], matrix[3:, :3], matrix[3:, 3:]
        cp = c.T @ c
        av = np.array((a[2, 1], a[0, 2], a[1, 0]))
        full = expm(matrix)
        eu, _ = polar(full[:3, :3], side="left")
        exact = frame_errors(full[:, :3], np.eye(3))
        leak4 = np.sum(c*c)/3-np.sum((d@c-c@a)**2)/36-np.trace(cp@cp)/9
        rot4 = av@av/4-(av@av)**2/48+np.sum((c@a)*(c@a-d@c))/24
        log3 = a+(a@cp+cp@a)/12-c.T@d@c/6
        scaling.append({"scale": scale, "leakage_remainder": abs(exact["leakage"]-leak4),
                        "rotation_remainder": abs(exact["rotation"]-rot4),
                        "polar_log_remainder": float(np.linalg.norm(logm(eu)-log3))})
    if max(residual, polar_residual) > 1e-10:
        raise AssertionError("General block formula failed.")
    return {"formal_error_coefficient_residual": residual,
            "formal_polar_log_coefficient_residual": polar_residual,
            "finite_matrix_scaling": scaling}


def physical_coefficients(p, values, es, gs):
    return {f"{m},{n}": float(values[j]/es**m/gs**n) for j, (m, n) in enumerate(p.keys)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tau", type=float, default=100.)
    parser.add_argument("--e-scale", type=float, default=.0002)
    parser.add_argument("--t1-scale", type=float, default=.0003)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).parent/"error_channels_tau100")
    args = parser.parse_args()
    if min(args.tau, args.e_scale, args.t1_scale) <= 0:
        parser.error("tau and parameter scales must be positive.")
    algebra = algebra_checks()
    model = Model(args.tau, args.e_scale, args.t1_scale, 4, 2e-13, 2e-15)
    p = model.series
    r = propagator_coefficients(model)
    r0 = r[0]
    s = np.einsum("ij,kjl->kil", r0.T, r)
    base_residual = float(np.linalg.norm(s[0]-np.eye(5)))
    s[0] = np.eye(5)  # exact base condition; record the discarded numerical drift
    k = logarithm_series(p, s)
    km, khh, kmm, khm = magnus_through_three(model)
    mask3 = np.array([sum(key) <= 3 for key in p.keys])
    relative, u = frame_error_series(p, s[:, :, :3], np.eye(3))
    absolute, _ = frame_error_series(p, r[:, :, :3], B)
    terms, polar_log = generator_error_terms(p, k)
    c1 = homogeneous(p, k[:, 3:, :3], 1)
    a1 = homogeneous(p, k[:, :3, :3], 1)
    c_norm, a_norm = inner_series(p, c1, c1), inner_series(p, a1, a1)/2
    ca_norm = inner_series(p, p.matmul(c1, a1), p.matmul(c1, a1))
    # Explicit E^2 g^2 coefficient in the model-specific rank-one formulas.
    cs = {key: k[p.index[key], 3:, :3] for key in p.keys}
    av = {key: np.array((k[p.index[key], 2, 1], k[p.index[key], 0, 2],
                         k[p.index[key], 1, 0])) for key in p.keys}
    ce, cg = cs[(1, 0)], cs[(0, 1)]
    ae, ag = av[(1, 0)], av[(0, 1)]
    uc, vc, wc = np.sum(ce*ce), np.sum(ce*cg), np.sum(cg*cg)
    ua, va, wa = ae@ae, ae@ag, ag@ag
    cross22 = uc*wa+4*vc*va+wc*ua
    leak22 = (np.sum(cs[(1, 1)]**2)+2*np.sum(cs[(2, 0)]*cs[(0, 2)])
              +2*np.sum(ce*cs[(1, 2)])+2*np.sum(cg*cs[(2, 1)]))/3
    leak22 -= cross22/36+(2*uc*wc+4*vc**2)/9
    rot22 = (av[(1, 1)]@av[(1, 1)]+2*av[(2, 0)]@av[(0, 2)]
             +2*ae@av[(1, 2)]+2*ag@av[(2, 1)])/4
    rot22 += -(2*ua*wa+4*va**2)/48+cross22/24
    odd = np.array([j for j, key in enumerate(p.keys) if sum(key) % 2])
    checks = {
        "R0_orthogonality_norm": base_residual,
        "K_antisymmetry_max_residual": float(np.max(abs(k+k.transpose(0, 2, 1)))),
        "log_exp_series_inverse_residual": float(np.max(abs(exponential_series(p, k)-s))),
        "Magnus_vs_log_Taylor_degree3_residual": float(np.max(abs(k[mask3]-km[mask3]))),
        "K2_hh_mm_hm_sum_residual": float(np.max(abs(homogeneous(p, k, 2)-khh-kmm-khm))),
        "polar_log_feedback_formula_residual": float(np.max(abs(logarithm_series(p, u)-polar_log))),
        "leakage_order4_formula_residual": float(np.max(abs(sum(terms["leakage"].values())-relative["leakage"]))),
        "rotation_order4_formula_residual": float(np.max(abs(sum(terms["rotation"].values())-relative["rotation"]))),
        "physical_D1_zero_residual": float(np.max(abs(homogeneous(p, k[:, 3:, 3:], 1)))),
        "physical_C1_columns_2_3_zero_residual": float(np.max(abs(c1[:, :, 1:]))),
        "physical_mm_second_order_A_zero_residual": float(np.max(abs(kmm[:, :3, :3]))),
        "physical_rank_one_feedback_identity_residual": float(np.max(abs(ca_norm-p.scalar_product(c_norm, a_norm)))),
        "explicit_mixed22_leakage_residual": float(abs(leak22-relative["leakage"][p.index[(2, 2)]])),
        "explicit_mixed22_rotation_residual": float(abs(rot22-relative["rotation"][p.index[(2, 2)]])),
        "odd_degree_error_residual": max(float(np.max(abs(values[odd])))
                                        for channel in (absolute, relative) for values in channel.values()),
    }
    for name, value in checks.items():
        if value > 2e-8:
            raise AssertionError(f"{name}: {value}")
    # Absolute leakage's fourth order includes the finite-time N0/N4 term.
    ns = [homogeneous(p, r[:, 3:, :3], degree) for degree in range(5)]
    absolute_leak4 = {
        "N2_norm": inner_series(p, ns[2], ns[2])/3,
        "N1_N3": 2*inner_series(p, ns[1], ns[3])/3,
        "N0_N4_baseline_interference": 2*inner_series(p, ns[0], ns[4])/3,
    }
    checks["absolute_leakage_fourth_order_residual"] = float(np.max(
        abs(sum(absolute_leak4.values())-homogeneous(p, absolute["leakage"], 4))))
    if checks["absolute_leakage_fourth_order_residual"] > 2e-8:
        raise AssertionError("Absolute baseline interference formula failed.")
    scan = []
    for scale in (1., .5, .25):
        exact = model.exact(scale, scale)
        exact_relative = r0.T @ exact
        row = {"scale": scale, "E1_meV": scale*args.e_scale, "t1_meV": scale*args.t1_scale,
               "absolute_full": frame_errors(exact[:, :3], B),
               "relative_full": frame_errors(exact_relative[:, :3], np.eye(3)),
               "generator_reconstruction": {}}
        for degree in (1, 2, 3, 4):
            generator = p.evaluate(k, scale, scale, degree)
            # Roundoff antisymmetry drift is measured above; enforce the Lie algebra.
            generator = (generator-generator.T)/2
            reconstructed = expm(generator)
            abs_errors = frame_errors((r0 @ reconstructed)[:, :3], B)
            rel_errors = frame_errors(reconstructed[:, :3], np.eye(3))
            row["generator_reconstruction"][str(degree)] = {
                "absolute": abs_errors, "relative": rel_errors,
                "absolute_difference_from_full": {name: abs_errors[name]-row["absolute_full"][name] for name in abs_errors},
                "relative_difference_from_full": {name: rel_errors[name]-row["relative_full"][name] for name in rel_errors},
                "orthogonality_norm": float(np.linalg.norm(reconstructed.T@reconstructed-np.eye(5))),
            }
        row["relative_degree4_terms"] = {
            name: {term: float(p.evaluate(value, scale, scale))
                   for term, value in pieces.items()} for name, pieces in terms.items()}
        row["relative_degree4_total_correction"] = {
            name: sum(value for term, value in pieces.items() if term.startswith("degree4"))
            for name, pieces in row["relative_degree4_terms"].items()}
        row["absolute_degree4_total_correction"] = {
            name: float(p.evaluate(homogeneous(p, value, 4), scale, scale))
            for name, value in absolute.items()}
        scan.append(row)
    second_sources = {}
    for key in ((2, 0), (1, 1), (0, 2)):
        j = p.index[key]
        second_sources[str(key)] = {
            name: {"A_Frobenius_norm": float(np.linalg.norm(value[j, :3, :3])),
                   "C_Frobenius_norm": float(np.linalg.norm(value[j, 3:, :3])),
                   "D_Frobenius_norm": float(np.linalg.norm(value[j, 3:, 3:]))}
            for name, value in (("hh", khh), ("mm", kmm), ("hm", khm), ("total", homogeneous(p, k, 2)))}
    result = {
        "parameters": {"tau_meV_inverse": args.tau, "total_time_meV_inverse": 6*args.tau,
                       "E1_scale_meV": args.e_scale, "t1_scale_meV": args.t1_scale,
                       "normalization": 2, "tc_meV": .3, "Ed_max_meV": .3},
        "checks": checks, "general_block_checks": algebra,
        "absolute_baseline_errors": {name: float(value[0]) for name, value in absolute.items()},
        "error_coefficients_physical": {
            mode: {name: physical_coefficients(p, values, args.e_scale, args.t1_scale)
                   for name, values in channel.items()} for mode, channel in
            (("absolute", absolute), ("relative", relative))},
        "relative_contribution_coefficients_physical": {
            name: {term: physical_coefficients(p, values, args.e_scale, args.t1_scale)
                   for term, values in pieces.items()} for name, pieces in terms.items()},
        "absolute_leakage_fourth_order_terms_at_scale1": {
            name: float(p.evaluate(values, 1, 1)) for name, values in absolute_leak4.items()},
        "second_order_generator_source_norms_scaled": second_sources,
        "scaling": scan,
        "scope": "Local error-channel analysis; relative and absolute errors are distinct; exponentiation preserves geometry, not a global accuracy guarantee.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    np.savez(args.output_dir/"generators.npz", powers=p.keys, R0=r0, R_coefficients=r,
             K=k, K_magnus=km, K2_hh=khh, K2_mm=kmm, K2_hm=khm,
             E1_scale=args.e_scale, t1_scale=args.t1_scale, tau=args.tau)
    print(json.dumps(checks, indent=2), flush=True)
    print("Absolute baseline errors:", result["absolute_baseline_errors"], flush=True)
    print("Absolute full errors:", scan[0]["absolute_full"], flush=True)
    print("Relative full errors:", scan[0]["relative_full"], flush=True)
    print("Relative degree-four decomposition:", json.dumps(
        scan[0]["relative_degree4_terms"], indent=2), flush=True)
    print("Absolute degree-four correction:", scan[0]["absolute_degree4_total_correction"], flush=True)
    print(f"Outputs: {args.output_dir.resolve()}", flush=True)


if __name__ == "__main__":
    main()
