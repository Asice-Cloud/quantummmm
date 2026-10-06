"""Exact frame-distance split and parameter orders for n7.md.

The total frame error is global. The pure leakage / SO(3) rotation
interpretation is restricted to nonsingular positive-determinant M.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from verify_manifold_expansion import B, Model, Series
from verify_error_channels import (exponential_series, homogeneous,
                                   inner_series, logarithm_series,
                                   propagator_coefficients)


def frame_metrics(w, target=B):
    m, n = w[..., :3, :], w[..., 3:, :]
    left, singular, vh = np.linalg.svd(m)
    u = left @ vh
    valid = (np.linalg.det(u) > 0) & (singular[..., -1] > 1e-8)
    total = (3-np.einsum("ij,...ij->...", target, m))/6
    leakage = (3-singular.sum(axis=-1))/6
    rotation = total-leakage
    # On invalid points these raw scalars only describe alignment over O(3).
    return {"total": total, "leakage_raw": leakage, "rotation_raw": rotation,
            "leakage": np.where(valid, leakage, np.nan),
            "rotation": np.where(valid, rotation, np.nan),
            "valid": valid, "sigma_min": singular[..., -1]}


def error_series(p, frame, target):
    m, n = frame[:, :3], frame[:, 3:]
    if np.linalg.det(m[0]) <= 0:
        raise ValueError("The expansion base is outside the positive SO(3) polar branch.")
    s = p.sqrt_spd(p.matmul(m.transpose(0, 2, 1), m))
    total = -np.einsum("ij,kij->k", target, m)/6
    total[0] += .5
    leakage = -np.trace(s, axis1=1, axis2=2)/6
    leakage[0] += .5
    rotation = total-leakage
    u = p.matmul(m, p.inverse(s))
    du = u.copy()
    du[0] -= target
    weighted = np.trace(p.matmul(p.matmul(du.transpose(0, 2, 1), du), s), axis1=1, axis2=2)/12
    delta = m.copy()
    delta[0] -= target
    squared = (inner_series(p, delta, delta)+inner_series(p, n, n))/12
    residuals = {"trace_vs_squared_distance": float(np.max(np.abs(total-squared))),
                 "split_vs_weighted_rotation": float(np.max(np.abs(rotation-weighted)))}
    return {"total": total, "leakage": leakage, "rotation": rotation}, residuals


def relative_terms(p, k):
    a1, a2, a3 = [homogeneous(p, k[:, :3, :3], d) for d in (1, 2, 3)]
    c1, c2, c3 = [homogeneous(p, k[:, 3:, :3], d) for d in (1, 2, 3)]
    d1 = homogeneous(p, k[:, 3:, 3:], 1)
    ca, dc = p.matmul(c1, a1), p.matmul(d1, c1)
    cp = p.matmul(c1.transpose(0, 2, 1), c1)
    anorm = inner_series(p, a1, a1)/2
    leakage = {
        "degree2_C1": inner_series(p, c1, c1)/12,
        "degree3_C1_C2": inner_series(p, c1, c2)/6,
        "degree4_C2": inner_series(p, c2, c2)/12,
        "degree4_C1_C3": inner_series(p, c1, c3)/6,
        "degree4_transport": -inner_series(p, dc-ca, dc-ca)/144,
        "degree4_leakage_saturation": -np.trace(p.matmul(cp, cp), axis1=1, axis2=2)/144,
    }
    rotation = {
        "degree2_A1": inner_series(p, a1, a1)/12,
        "degree3_A1_A2": inner_series(p, a1, a2)/6,
        "degree4_A2": inner_series(p, a2, a2)/12,
        "degree4_A1_A3": inner_series(p, a1, a3)/6,
        "degree4_rotation_saturation": -p.scalar_product(anorm, anorm)/72,
        "degree4_surviving_weight": -inner_series(p, ca, ca)/72,
        "degree4_auxiliary_rotation": -inner_series(p, ca, dc)/36,
    }
    return {"leakage": leakage, "rotation": rotation}


def general_checks():
    rng = np.random.default_rng(2706)
    p = Series(4)
    k = rng.normal(size=(p.count, 5, 5))*.15
    k = k-k.transpose(0, 2, 1)
    k[0] = 0
    frame = exponential_series(p, k)[:, :, :3]
    coefficients, algebra = error_series(p, frame, np.eye(3))
    terms = relative_terms(p, k)
    for channel in terms:
        residual = np.max(np.abs(coefficients[channel]-sum(terms[channel].values())))
        algebra[channel+"_quartic_formula"] = float(residual)
        assert residual < 1e-12
    # Separate finite-matrix test; K has nonzero A,C,D and does not commute blockwise.
    z = rng.normal(size=(5, 5))*.3
    z = z-z.T
    univariate = np.zeros_like(k)
    univariate[p.index[(1, 0)]] = z
    terms = relative_terms(p, univariate)
    remainders = []
    for scale in (.2, .1, .05):
        exact = frame_metrics(expm(scale*z)[:, :3], np.eye(3))
        item = {"scale": scale}
        for channel in terms:
            approximation = p.evaluate(sum(terms[channel].values()), scale, 0)
            item[channel+"_remainder"] = float(abs(exact[channel]-approximation))
        remainders.append(item)
    for channel in terms:
        ratio = remainders[0][channel+"_remainder"]/remainders[1][channel+"_remainder"]
        assert 50 < ratio < 80
    worst = np.zeros(4)
    for _ in range(100):
        n = rng.normal(size=(2, 3))
        n *= rng.uniform(.01, .98)/np.linalg.norm(n, 2)
        def root(x):
            v, q = np.linalg.eigh(x)
            return (q*np.sqrt(v))@q.T
        s = root(np.eye(3)-n.T@n)
        j = root(np.eye(2)-n@n.T)
        a = rng.normal(size=(3, 3))
        u = expm(a-a.T)
        leakage_step = np.block([[s, -n.T], [n, j]])
        w = np.concatenate((u@s, n))
        out = frame_metrics(w)
        direct = (np.linalg.norm(w[:3]-B)**2+np.linalg.norm(n)**2)/12
        weighted = np.linalg.norm((u-B)@root(s))**2/12
        worst = np.maximum(worst, [np.max(abs(leakage_step.T@leakage_step-np.eye(5))),
                                  abs(np.linalg.det(leakage_step)-1),
                                  abs(direct-out["total"]),
                                  abs(weighted-out["rotation"])])
    assert worst.max() < 1e-12
    algebra["random_sequential_decomposition_max"] = worst.tolist()
    return {"algebra": algebra, "sixth_order_remainders": remainders}


def physical_case(tau, e_scale, g_scale):
    model = Model(tau, e_scale, g_scale, 4, 3e-12, 3e-14)
    p = model.series
    r = propagator_coefficients(model)
    relative = r[0].T @ r
    relative[0] = np.eye(5)
    k = logarithm_series(p, relative)
    rel, relcheck = error_series(p, relative[:, :, :3], np.eye(3))
    absolute, abscheck = error_series(p, r[:, :, :3], B)
    terms = relative_terms(p, k)
    formula_checks = {ch: float(np.max(np.abs(rel[ch]-sum(parts.values()))))
                      for ch, parts in terms.items()}
    assert max(formula_checks.values()) < 1e-9
    def c(m, n):
        return k[p.index[(m, n)], 3:, :3]
    def a(m, n):
        block = k[p.index[(m, n)], :3, :3]
        return np.array([block[2, 1], block[0, 2], block[1, 0]])
    def dot(v, w):
        return float(np.sum(v*w))
    uc, vc, wc = dot(c(1, 0), c(1, 0)), dot(c(1, 0), c(0, 1)), dot(c(0, 1), c(0, 1))
    ua, va, wa = dot(a(1, 0), a(1, 0)), dot(a(1, 0), a(0, 1)), dot(a(0, 1), a(0, 1))
    q22 = uc*wa+4*vc*va+wc*ua
    dc22 = (dot(c(1, 1), c(1, 1))+2*dot(c(2, 0), c(0, 2))
            +2*dot(c(1, 0), c(1, 2))+2*dot(c(0, 1), c(2, 1)))
    da22 = (dot(a(1, 1), a(1, 1))+2*dot(a(2, 0), a(0, 2))
            +2*dot(a(1, 0), a(1, 2))+2*dot(a(0, 1), a(2, 1)))
    mixed = {"leakage": dc22/12-(q22+2*uc*wc+4*vc**2)/144,
             "rotation": da22/6-(2*ua*wa+4*va**2+q22)/72}
    for ch, value in mixed.items():
        residual = abs(value-rel[ch][p.index[(2, 2)]])
        formula_checks[ch+"_explicit_mixed22"] = float(residual)
        assert residual < 1e-10
    scales = np.array([e_scale**m*g_scale**n for m, n in p.keys])
    output = {"tau": tau, "E1": e_scale, "t1": g_scale,
              "checks": {"absolute": abscheck, "relative": relcheck,
                         "relative_formulas": formula_checks},
              "coefficients_physical": {}, "relative_term_contributions": {},
              "reconstruction": [], "local_scalar_expansion": []}
    for convention, coeffs in (("absolute", absolute), ("relative", rel)):
        output["coefficients_physical"][convention] = {
            ch: {f"{m},{n}": float(coeffs[ch][i]/scales[i]) for i, (m, n) in enumerate(p.keys)}
            for ch in coeffs}
    for ch, parts in terms.items():
        output["relative_term_contributions"][ch] = {name: float(p.evaluate(v, 1, 1))
                                                      for name, v in parts.items()}
    full = model.exact(1, 1)
    for degree in range(1, 5):
        kp = p.evaluate(k, 1, 1, degree)
        approx = frame_metrics((r[0]@expm((kp-kp.T)/2))[:, :3])
        exact = frame_metrics(full[:, :3])
        output["reconstruction"].append({"order": degree,
            "error": {ch: float(approx[ch]) for ch in ("total", "leakage", "rotation")},
            "difference_from_full": {ch: float(approx[ch]-exact[ch]) for ch in ("total", "leakage", "rotation")}})
    output["full_absolute"] = {ch: float(exact[ch]) for ch in ("total", "leakage", "rotation")}
    for scale in (1., .5, .25):
        rr = model.exact(scale, scale)
        item = {"scale": scale}
        for convention, frame, target in (("absolute", rr[:, :3], B),
                                           ("relative", (r[0].T@rr)[:, :3], np.eye(3))):
            exact = frame_metrics(frame, target)
            coeffs = absolute if convention == "absolute" else rel
            item[convention] = {ch: {
                "full": float(exact[ch]),
                "quadratic": float(p.evaluate(coeffs[ch], scale, scale, 2)),
                "quartic": float(p.evaluate(coeffs[ch], scale, scale, 4))}
                for ch in ("total", "leakage", "rotation")}
        output["local_scalar_expansion"].append(item)
    output["odd_coefficient_max"] = float(max(np.max(np.abs(v[[i for i, key in enumerate(p.keys) if sum(key)%2]]))
                                              for v in absolute.values()))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent/"additive_error_results")
    args = parser.parse_args()
    result = {"general": general_checks(), "cases": [physical_case(100, .0002, .0003),
                                                     physical_case(1, .02, .03)]}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"verification.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"general": result["general"], "cases": [{
        "tau": c["tau"], "checks": c["checks"], "full_absolute": c["full_absolute"],
        "relative_terms": c["relative_term_contributions"]} for c in result["cases"]]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
