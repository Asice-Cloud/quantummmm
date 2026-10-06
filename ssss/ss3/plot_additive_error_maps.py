"""n7: new additive frame error, with separate leakage and rotation costs.

Reuse validated generator coefficients from n6; solve full reference frames
again because the previous scalar metrics cannot determine the new metric.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import multiprocessing as mp
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_error_channel_maps import exp_skew, orthogonal_factor, panel
from verify_additive_errors import frame_metrics
from verify_error_channels import logarithm_series, propagator_coefficients
from verify_manifold_expansion import B, Model


def full_frames(x, y, energy, rtol=3e-12, atol=3e-14):
    model = Model(100*x, .01, .01, 4, rtol, atol)
    t1 = energy*10.**y
    def rhs(segment, s, vector):
        r = vector.reshape(len(y), 5, 5)
        o, oe, og = model.generators(segment, s)
        return ((o+(energy/.01)*oe+(t1/.01)[:, None, None]*og) @ r).ravel()
    raw = model.integrate(rhs, np.broadcast_to(np.eye(5), (len(y), 5, 5)).copy()).reshape(len(y), 5, 5)
    drift = np.max(abs(raw.swapaxes(-1, -2)@raw-np.eye(5)))
    return orthogonal_factor(raw)[:, :, :3], float(drift)


def reconstruct(r0, k, powers, scales, energy, y):
    degree = powers.sum(axis=1)
    weights = np.array([(energy/scales[0])**m*(energy*10.**y/scales[1])**n
                        for m, n in powers]).T
    frames, norms = [], []
    for order in range(1, 5):
        kp = np.einsum("yk,kij->yij", weights*(degree <= order), k)
        frames.append((r0@exp_skew(kp))[:, :, :3])
        kd = np.einsum("yk,kij->yij", weights*(degree == order), k)
        norms.append(np.linalg.norm(kd, axis=(-1, -2)))
    return np.array(frames), np.array(norms)


def column(job):
    i, x, y, energy, r0, k, powers, scales = job
    approx, norms = reconstruct(r0, k, powers, scales, energy, y)
    reference, drift = full_frames(x, y, energy)
    return i, np.concatenate((approx, reference[None])), norms, drift


def summarize(data):
    out = frame_metrics(data["frames"])
    common = out["valid"].all(axis=0)
    m, n = data["frames"][..., :3, :], data["frames"][..., 3:, :]
    direct = (np.sum((m-B)**2, axis=(-1, -2))+np.sum(n*n, axis=(-1, -2)))/12
    residual = float(np.max(abs(direct-out["total"])))
    assert residual < 1e-9
    assert np.min(out["total"]) >= -1e-9 and np.max(out["total"]) <= 1+1e-9
    assert np.nanmin(out["leakage"]) >= -1e-9 and np.nanmin(out["rotation"]) >= -1e-9
    assert np.nanmax(abs(out["total"]-out["leakage"]-out["rotation"])) < 1e-12
    result = {"parameters": {"E1": float(data["E1"]), "tc": .3, "Ed_max": .3,
                              "normalization": 2, "cycles": 2,
                              "x_range": [float(data["x"][0]), float(data["x"][-1])],
                              "y_range": [float(data["y"][0]), float(data["y"][-1])],
                              "grid_shape_y_x": [len(data["y"]), len(data["x"])]},
              "common_positive_polar_fraction": float(np.mean(common)),
              "distance_vs_trace_identity_max": residual,
              "full_reference_orthogonality_before_projection_max": float(data["integration_drift"].max()),
              "comparison": {}}
    for j in range(5):
        label = f"order_{j+1}" if j < 4 else "full"
        entry = {"invalid_polar_points": int((~out["valid"][j]).sum())}
        for key in ("total", "leakage", "rotation"):
            mask = np.ones_like(common) if key == "total" else common
            a, b = out[key][j][mask], out[key][4][mask]
            entry[key] = {"mean_absolute_difference": float(np.mean(abs(a-b))),
                          "max_absolute_difference": float(np.max(abs(a-b))),
                          "rmse": float(np.sqrt(np.mean((a-b)**2))),
                          "min_on_valid_points": float(np.nanmin(out[key][j])),
                          "max_on_valid_points": float(np.nanmax(out[key][j]))}
        result["comparison"][label] = entry
    return result, out


def render(data, out, directory):
    x, y = data["x"], data["y"]
    selected = [0, 3, 4]
    titles = [r"$K^{[1]}$", r"$K^{[1]}+\cdots+K^{[4]}$", "Full SO(5) evolution"]
    subtitle = (rf"$E_1={float(data['E1']):g}$ meV; $t_c=E_{{d,\max}}=0.3$ meV; "
                r"$T=6\tau$; normalization $2$")
    fidelity = 1-out["total"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), layout="constrained")
    for ax, j, title in zip(axes, selected, titles):
        im = panel(ax, x, y, fidelity[j], title)
    fig.colorbar(im, ax=axes, label=r"$F_{\rm frame}=1-\mathcal{E}_{\rm total}$", shrink=.85)
    fig.suptitle("New additive frame-distance fidelity\n"+subtitle, fontsize=12)
    fig.savefig(directory/"fidelity_comparison.png", dpi=190)
    fig.savefig(directory/"fidelity_comparison.pdf")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.8), layout="constrained")
    im = panel(ax, x, y, fidelity[3], r"$R^{(4)}=R_0\exp(K^{[1]}+\cdots+K^{[4]})$")
    fig.colorbar(im, ax=ax, label=r"$F_{\rm frame}$")
    fig.suptitle(subtitle, fontsize=10)
    fig.savefig(directory/"fidelity_order4.png", dpi=220)
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout="constrained")
    for row, (key, limit, label) in enumerate([
            ("leakage", 1/3, r"$\mathcal{E}_{\rm leakage}$"),
            ("rotation", 1, r"$\mathcal{E}_{\rm rotation}$")]):
        for col, (j, title) in enumerate(zip(selected, titles)):
            im = panel(axes[row, col], x, y, out[key][j], title, 0, limit)
        fig.colorbar(im, ax=axes[row], label=label, shrink=.85)
    fig.suptitle("Additive error components (brighter = larger error)\n"+subtitle+
                 "\nGray: pure leakage / SO(3) rotation split unavailable", fontsize=12)
    fig.savefig(directory/"error_components.png", dpi=190)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), layout="constrained")
    for ax, j in zip(axes, [0, 2, 3]):
        im = panel(ax, x, y, fidelity[j]-fidelity[4], f"Order {j+1} minus full", -1, 1, "RdBu_r")
    fig.colorbar(im, ax=axes, label=r"$\Delta F_{\rm frame}$", shrink=.85)
    fig.suptitle("Approximation error of the new fidelity\n"+subtitle, fontsize=12)
    fig.savefig(directory/"approximation_error.png", dpi=190)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), layout="constrained")
    for ax, key in zip(axes, ["leakage", "rotation", "total"]):
        diff = out[key][3]-out[key][0]
        limit = max(float(np.nanmax(abs(diff))), 1e-12)
        im = panel(ax, x, y, diff, key+": order 4 minus order 1", -limit, limit, "RdBu_r")
        fig.colorbar(im, ax=ax, label="Change in error", shrink=.8)
    fig.suptitle("Effect of generator orders 2, 3, 4 (separate color scales)", fontsize=12)
    fig.savefig(directory/"high_order_change.png", dpi=190)
    plt.close(fig)


def validate(data):
    ix = sorted(set([1, len(data["x"])//2, len(data["x"])-1]))
    iy = sorted(set([0, len(data["y"])//2, len(data["y"])-1]))
    out = []
    for i in ix:
        y = data["y"][iy]
        w, drift = full_frames(data["x"][i], y, float(data["E1"]), 1e-13, 1e-15)
        tight = frame_metrics(w)
        old = frame_metrics(data["frames"][4, iy, i])
        out.append({"x": float(data["x"][i]), "y": y.tolist(),
                    "max_full_metric_change": {key: float(np.nanmax(abs(tight[key]-old[key])))
                                                for key in ("total", "leakage", "rotation")}})
    # Recompute coefficients at the long-time endpoint; do not merely reuse K.
    x = float(data["x"][-1])
    model = Model(100*x, .01, .01, 4, 1e-13, 1e-15)
    coeff = propagator_coefficients(model)
    s = coeff[0].T@coeff
    s[0] = np.eye(5)
    k = logarithm_series(model.series, s)
    k = (k-k.swapaxes(-1, -2))/2
    ww, _ = reconstruct(orthogonal_factor(coeff[0]), k, data["powers"], data["coefficient_scales_meV"],
                         float(data["E1"]), data["y"][iy])
    tight, old = frame_metrics(ww[3]), frame_metrics(data["frames"][3, iy, -1])
    return {"tighter_full_reference": out,
            "tighter_fourth_order_at_largest_tau": {key: float(np.nanmax(abs(tight[key]-old[key])))
                                                    for key in ("total", "leakage", "rotation")}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).parent/"error_channel_maps/maps.npz")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent/"additive_error_maps")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--render-only", action="store_true")
    parser.add_argument("--validate", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache = args.output_dir/"maps.npz"
    if args.render_only:
        with np.load(cache) as saved:
            data = dict(saved)
    else:
        with np.load(args.source) as saved:
            data = {k: saved[k] for k in ("x", "y", "E1", "R0", "K", "powers", "coefficient_scales_meV")}
        nx, ny = len(data["x"]), len(data["y"])
        data["source_cache"] = str(args.source)
        data["frames"] = np.empty((5, ny, nx, 5, 3))
        data["homogeneous_K_norm"] = np.empty((4, ny, nx))
        data["integration_drift"] = np.empty(nx)
        jobs = [(i, x, data["y"], float(data["E1"]), data["R0"][i], data["K"][i],
                 data["powers"], data["coefficient_scales_meV"]) for i, x in enumerate(data["x"])]
        start = time.monotonic()
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn")) as pool:
            for count, future in enumerate(as_completed([pool.submit(column, j) for j in jobs]), 1):
                i, frames, norms, drift = future.result()
                data["frames"][:, :, i] = frames
                data["homogeneous_K_norm"][:, :, i] = norms
                data["integration_drift"][i] = drift
                if count%10 == 0 or count == nx:
                    print(f"{count}/{nx} columns; {time.monotonic()-start:.1f} seconds", flush=True)
        np.savez_compressed(cache, **data)
    report, metrics = summarize(data)
    np.savez_compressed(args.output_dir/"metrics.npz", x=data["x"], y=data["y"], **metrics)
    (args.output_dir/"summary.json").write_text(json.dumps(report, indent=2)+"\n")
    if args.validate:
        checks = validate(data)
        (args.output_dir/"validation.json").write_text(json.dumps(checks, indent=2)+"\n")
        print(json.dumps(checks, indent=2), flush=True)
    render(data, metrics, args.output_dir)
    print(json.dumps(report, indent=2), flush=True)
    print(f"Wrote {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
