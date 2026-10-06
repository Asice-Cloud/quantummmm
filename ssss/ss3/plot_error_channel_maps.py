"""Fig.1(d)-style maps of the two original group-space error channels.

Truncate log(R0.T R) in parameters, retain its matrix exponential, then
evaluate the ORIGINAL absolute leakage and polar-rotation metrics. Compare
with full SO(5) evolution. No state overlap or scalar-polynomial clipping.
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
from scipy.linalg import expm

from verify_manifold_expansion import B, Model
from verify_error_channels import propagator_coefficients, logarithm_series


def orthogonal_factor(r):
    u, _, vh = np.linalg.svd(r)
    result = u @ vh
    if np.any(np.linalg.det(result) < 0):
        raise ValueError("Unexpected improper SO(5) evolution.")
    return result


def exp_skew(k):
    """Batch exponential through Hermitian iK; retains SO(5) at large angles."""
    values, vectors = np.linalg.eigh(1j*k)
    result = (vectors*np.exp(-1j*values)[..., None, :]) @ vectors.conj().swapaxes(-1, -2)
    if np.max(np.abs(result.imag)) > 1e-7:
        raise ValueError("Exponential lost real structure.")
    return result.real


def metrics(r):
    m, n = r[..., :3, :3], r[..., 3:, :3]
    left, singular, vh = np.linalg.svd(m)
    u = left @ vh
    valid = (np.linalg.det(u) > 0) & (singular[..., -1] > 1e-8)
    fs = 1-np.sum(n*n, axis=(-1, -2))/3
    fr = (1+np.einsum("ij,...ij->...", B, u))/4
    fr = np.where(valid, fr, np.nan)
    return fs, fr, valid, singular[..., -1]


def calculate_column(job):
    index, x, ys, energy, rtol, atol = job
    model = Model(100*x, .01, .01, 4, rtol, atol)
    coefficients = propagator_coefficients(model)
    r0raw = coefficients[0]
    s = r0raw.T @ coefficients
    s[0] = np.eye(5)
    kraw = logarithm_series(model.series, s)
    skew_defect = float(np.max(np.abs(kraw+kraw.swapaxes(-1, -2))))
    k = (kraw-kraw.swapaxes(-1, -2))/2
    r0 = orthogonal_factor(r0raw)
    eratio = energy/.01
    gratios = energy*10.**ys/.01
    weights = np.array([eratio**m*gratios**n for m, n in model.series.keys]).T
    degree = np.sum(model.series.keys, axis=1)
    reconstructed = []
    knorm = []
    for order in range(1, 5):
        kp = np.einsum("yk,kij->yij", weights*(degree <= order), k)
        reconstructed.append(r0 @ exp_skew(kp))
        knorm.append(np.linalg.norm(kp, axis=(-1, -2)))

    def rhs(segment, s, vector):
        r = vector.reshape(len(ys), 5, 5)
        o, oe, og = model.generators(segment, s)
        generators = o+eratio*oe+gratios[:, None, None]*og
        return (generators @ r).ravel()

    initial = np.broadcast_to(np.eye(5), (len(ys), 5, 5)).copy()
    raw_reference = model.integrate(rhs, initial).reshape(len(ys), 5, 5)
    drift = max(float(np.max(np.abs(raw_reference.swapaxes(-1, -2) @ raw_reference-np.eye(5)))),
                float(np.max(np.abs(r0raw.T @ r0raw-np.eye(5)))))
    reference = orthogonal_factor(raw_reference)
    all_r = np.concatenate((np.array(reconstructed), reference[None]), axis=0)
    fs, fr, valid, sigma = metrics(all_r)
    return index, r0, k, fs, fr, valid, sigma, np.array(knorm), drift, skew_defect


def panel(ax, x, y, z, title, vmin=0, vmax=1, cmap="viridis"):
    palette = plt.get_cmap(cmap).copy()
    palette.set_bad("0.75")
    image = ax.pcolormesh(x, y, np.ma.masked_invalid(z), shading="auto",
                         cmap=palette, vmin=vmin, vmax=vmax, rasterized=True)
    if vmin == 0 and vmax == 1 and np.nanmax(z)-np.nanmin(z) > .05:
        ax.contour(x, y, np.ma.masked_invalid(z), levels=[.2, .4, .6, .8],
                   colors="white", linewidths=.35, alpha=.55)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(r"$\tau/(100\,\mathrm{meV}^{-1})$")
    ax.set_ylabel(r"$\log_{10}(t_1/E_1)$")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[0], y[-1])
    return image


def render(data, directory):
    x, y = data["x"], data["y"]
    fs, fr = data["F_surv"], data["F_rot"]
    cols, titles = [0, 3, 4], [r"$K^{[1]}$", r"$K^{[1]}+\cdots+K^{[4]}$", "Full SO(5) evolution"]
    subtitle = (rf"$E_1={float(data['E1']):g}$ meV; $t_c=E_{{d,\max}}=0.3$ meV; "
                r"$T=6\tau$; SO(5) factor $2$")
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout="constrained")
    for row, (array, label) in enumerate(((fs, r"$F_{\rm surv}=1-\varepsilon_{\rm leak}$"),
                                         (fr, r"$F_{\rm rot}=1-\varepsilon_{\rm rot}$"))):
        for col, (j, title) in enumerate(zip(cols, titles)):
            im = panel(axes[row, col], x, y, array[j], title)
        fig.colorbar(im, ax=axes[row], label=label, shrink=.85)
    fig.suptitle("Original group-space fidelity channels\n"+subtitle+
                 "\nGray: positive SO(3) polar branch unavailable", fontsize=12)
    fig.savefig(directory/"error_channel_comparison.png", dpi=190)
    fig.savefig(directory/"error_channel_comparison.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), layout="constrained")
    for ax, j, title in zip(axes, cols, titles):
        im = panel(ax, x, y, fs[j]*fr[j], title)
    fig.colorbar(im, ax=axes, label=r"$F_{\rm group}=F_{\rm surv}F_{\rm rot}$", shrink=.8)
    fig.suptitle("Same original product metric (not a state overlap)\n"+subtitle, fontsize=12)
    fig.savefig(directory/"group_fidelity_comparison.png", dpi=190)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.8), layout="constrained")
    im = panel(ax, x, y, fs[3]*fr[3], r"$R^{(4)}=R_0\exp(K^{[1]}+\cdots+K^{[4]})$")
    fig.colorbar(im, ax=ax, label=r"$F_{\rm group}=F_{\rm surv}F_{\rm rot}$")
    fig.suptitle(subtitle, fontsize=10)
    fig.savefig(directory/"group_fidelity_order4.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout="constrained")
    for row, (array, label) in enumerate(((fs, r"$F_{\rm surv}$"), (fr, r"$F_{\rm rot}$"))):
        for col, (j, title) in enumerate(zip([0, 2, 3], ["Order 1", "Order 3", "Order 4"])):
            delta = array[j]-array[4]
            im = panel(axes[row, col], x, y, delta, title+" minus full evolution", -1, 1, "RdBu_r")
        fig.colorbar(im, ax=axes[row], label=r"$\Delta$"+label, shrink=.85)
    fig.suptitle("Approximation error, with the fidelity definitions held fixed\n"+subtitle, fontsize=12)
    fig.savefig(directory/"approximation_error.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), layout="constrained")
    for ax, array, label in zip(axes, [fs, fr], [r"$F_{\rm surv}$", r"$F_{\rm rot}$"]):
        im = panel(ax, x, y, array[3]-array[0], "Order 4 minus order 1: "+label, -1, 1, "RdBu_r")
    fig.colorbar(im, ax=axes, label="Change after adding generator orders 2, 3, 4", shrink=.8)
    fig.suptitle("High-order correction (a change is not necessarily an improvement)", fontsize=12)
    fig.savefig(directory/"high_order_change.png", dpi=190)
    plt.close(fig)


def summarize(data):
    result = {"parameters": {"E1_meV": float(data["E1"]), "tc_meV": .3,
                              "Ed_max_meV": .3, "cycles": 2, "normalization": 2,
                              "grid_shape_y_x": [len(data["y"]), len(data["x"])]},
              "integration_orthogonality_max": float(data["integration_drift"].max()),
              "K_coefficient_skew_defect_max": float(data["skew_defect"].max()),
              "comparison": {}}
    for j in range(5):
        label = f"order_{j+1}" if j < 4 else "full_evolution"
        entry = {"invalid_rotation_fraction": float(np.mean(~data["valid_rotation"][j]))}
        for key in ("F_surv", "F_rot"):
            a = data[key][j]
            delta = a-data[key][4]
            entry[key] = {"min": float(np.nanmin(a)), "max": float(np.nanmax(a)),
                          "mean_absolute_error": float(np.nanmean(np.abs(delta))),
                          "max_absolute_error": float(np.nanmax(np.abs(delta))),
                          "rmse": float(np.sqrt(np.nanmean(delta**2)))}
        if j < 4:
            ok = (np.abs(data["F_surv"][j]-data["F_surv"][4]) < .02) & (
                np.abs(data["F_rot"][j]-data["F_rot"][4]) < .02)
            entry["both_channel_errors_below_0.02_fraction_of_entire_grid"] = float(np.mean(ok))
        result["comparison"][label] = entry
    common = np.all(data["valid_rotation"], axis=0)
    result["rotation_comparison_on_same_valid_points"] = {
        "fraction_of_grid": float(np.mean(common)),
        "mean_absolute_error_by_order": [float(np.mean(np.abs(
            data["F_rot"][j][common]-data["F_rot"][4][common]))) for j in range(4)]}
    return result


def validate(data):
    """Tighter independent solves at nine points, plus an older RK4 map."""
    result = {"tighter_integration": []}
    ix = sorted(set([1, len(data["x"])//2, len(data["x"])-1]))
    iy = sorted(set([0, len(data["y"])//2, len(data["y"])-1]))
    for i in ix:
        column = calculate_column((i, data["x"][i], data["y"][iy],
                                   float(data["E1"]), 1e-13, 1e-15))
        _, _, _, fs, fr, _, _, _, _, _ = column
        result["tighter_integration"].append({
            "x": float(data["x"][i]), "y": data["y"][iy].tolist(),
            "F_surv_max_difference_by_order_then_full":
                np.nanmax(np.abs(fs-data["F_surv"][:, iy, i]), axis=1).tolist(),
            "F_rot_max_difference_by_order_then_full":
                np.nanmax(np.abs(fr-data["F_rot"][:, iy, i]), axis=1).tolist()})
    old_path = Path(__file__).parent/"so5_group_fidelity_map.npz"
    with np.load(old_path) as old:
        if (np.array_equal(data["x"], old["x"]) and np.array_equal(data["y"], old["y"])
                and float(data["E1"]) == float(old["E1"])):
            result["older_RK4_map_max_difference_on_valid_polar_branch"] = {
                key: float(np.nanmax(np.abs(data[key][4]-old[oldkey])))
                for key, oldkey in [("F_surv", "F_surv"), ("F_rot", "F_rot_cond")]}
    weights = np.array([(float(data["E1"])/.01)**m *
                        (float(data["E1"])*10**data["y"][-1]/.01)**n
                        for m, n in data["powers"]])
    k = np.einsum("k,kij->ij", weights, data["K"][-1])
    result["largest_parameter_endpoint_eigen_expm_vs_scipy_expm_max"] = float(
        np.max(np.abs(exp_skew(k)-expm(k))))
    result["largest_parameter_endpoint_K_norm"] = float(np.linalg.norm(k))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=121)
    parser.add_argument("--ny", type=int, default=101)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--energy", type=float, default=.01)
    parser.add_argument("--xmax", type=float, default=12.)
    parser.add_argument("--ymin", type=float, default=-1.)
    parser.add_argument("--ymax", type=float, default=1.)
    parser.add_argument("--rtol", type=float, default=3e-12)
    parser.add_argument("--atol", type=float, default=3e-14)
    parser.add_argument("--render-only", action="store_true")
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent/"error_channel_maps")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache = args.output_dir/"maps.npz"
    if args.render_only:
        with np.load(cache) as saved:
            data = dict(saved)
    else:
        x = np.linspace(0, args.xmax, args.nx)
        y = np.linspace(args.ymin, args.ymax, args.ny)
        data = {"x": x, "y": y, "E1": args.energy, "rtol": args.rtol, "atol": args.atol,
                "powers": np.array(Model(0, .01, .01, 4, 1e-10, 1e-12).series.keys),
                "coefficient_scales_meV": np.array([.01, .01]),
                "R0": np.empty((args.nx, 5, 5)), "K": np.empty((args.nx, 15, 5, 5)),
                "F_surv": np.empty((5, args.ny, args.nx)), "F_rot": np.empty((5, args.ny, args.nx)),
                "valid_rotation": np.empty((5, args.ny, args.nx), dtype=bool),
                "sigma_min_M": np.empty((5, args.ny, args.nx)),
                "K_norm": np.empty((4, args.ny, args.nx)),
                "integration_drift": np.empty(args.nx), "skew_defect": np.empty(args.nx)}
        jobs = [(i, v, y, args.energy, args.rtol, args.atol) for i, v in enumerate(x)]
        start = time.monotonic()
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn")) as pool:
            pending = [pool.submit(calculate_column, job) for job in jobs]
            for count, future in enumerate(as_completed(pending), 1):
                i, r0, k, fs, fr, valid, sigma, knorm, drift, skew = future.result()
                data["R0"][i], data["K"][i] = r0, k
                for key, value in (("F_surv", fs), ("F_rot", fr), ("valid_rotation", valid),
                                   ("sigma_min_M", sigma), ("K_norm", knorm)):
                    data[key][..., i] = value
                data["integration_drift"][i], data["skew_defect"][i] = drift, skew
                if count % 10 == 0 or count == len(jobs):
                    print(f"{count}/{len(jobs)} columns; {time.monotonic()-start:.1f} seconds", flush=True)
        np.savez_compressed(cache, **data)
    report = summarize(data)
    (args.output_dir/"summary.json").write_text(json.dumps(report, indent=2)+"\n")
    if args.validate:
        checks = validate(data)
        (args.output_dir/"validation.json").write_text(json.dumps(checks, indent=2)+"\n")
        print(json.dumps(checks, indent=2), flush=True)
    render(data, args.output_dir)
    print(json.dumps(report, indent=2), flush=True)
    print(f"Wrote maps and figures to {args.output_dir}", flush=True)


if __name__ == "__main__":
    main()
