"""Overlay the scalar fidelity formula on curve samples from solution.png.

This compares the user's four existing cuts with the same scalar formula as
the heatmap. The PDF color digitization used in the earlier plot is unrelated
to these source samples. All time axes and physical parameters stay fixed.
"""

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/quantumsss_fidelity_mpl")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

from plot_fidelity_formula import formula, J, K

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "read/figures/fidelity_formula"
DATA = ROOT / "read/data"
from image_samples import image_samples
from full_model import ensure_full_model_cuts


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lam = np.linspace(0., 12., 601)
    ratios = (.1, 1., 2.3446, 10.)
    report = {"image_source": "solution.png", "fit_parameters": 0,
              "formula": "Scalar continuation formula, uniform ordered integrals through I9",
              "parameters": {"E1_meV": .01, "tc_meV": .3, "E0_meV": .3,
                             "K": "2A", "cycles": 2},
              "image_values_are_digitized": True, "curves": []}
    arrays = {"lambda": lam}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 7.8), layout="constrained")
    with np.load(ensure_full_model_cuts(DATA)) as full:
        assert np.allclose(full["tau"], 100*lam)
        for index, (ratio, ax) in enumerate(zip(ratios, axes.flat)):
            data = formula(ratio, lam)
            samples = image_samples(index)
            black, orange = samples["black"], samples["orange"]
            ax.scatter(black[:, 0]/100, black[:, 1], s=7, color=".25", alpha=.5,
                       linewidths=0, label="solution.png: Complete B (pixels)", zorder=1)
            ax.scatter(orange[:, 0]/100, orange[:, 1], s=9, marker="x", color="#e87514",
                       linewidths=.7, alpha=.8,
                       label="solution.png: Slow rotations (pixels)", zorder=2)
            ax.plot(lam, data["f"], color="#2264c3", linewidth=1.55,
                    label="Our scalar fidelity formula", zorder=3)
            ax.set(title=rf"$r=|t_1|_{{\max}}/E_1={ratio:g}$", xlim=(0, 12),
                   ylim=(-.025, 1.035), xlabel=r"$\tau\ (100/\mathrm{meV})$",
                   ylabel=r"$f_r(\tau)=\sqrt{F_r(\tau)}$")
            ax.grid(alpha=.18)
            def rho(s):
                a = (1-np.cos(np.pi*s))/2
                b = 1-a
                return np.hypot(ratio*b, a/np.hypot(a, b))
            action = quad(rho, 0., 1., epsabs=2e-13, epsrel=2e-13)[0]
            row = {"ratio": ratio, "action_S1": action,
                   "phase_coefficient_D": action+1-ratio*J+K, "image_comparison": {}}
            for name, points in samples.items():
                if len(points) == 0:
                    continue
                points = points[points[:, 0] >= 50]
                if len(points) == 0:
                    continue
                error = np.interp(points[:, 0]/100, lam, data["f"])-points[:, 1]
                row["image_comparison"][name] = {
                    "n": len(points), "rmse": float(np.sqrt(np.mean(error**2))),
                    "max_abs": float(np.max(abs(error))),
                    "caution": "Contains rasterization, finite line width, and interpolation errors."}
            exact = full[f"f_ratio_{ratio}"]
            use = lam >= .5
            delta = data["f"][use]-exact[use]
            row["scalar_vs_full_SO5_x_ge_0p5"] = {
                "max_abs": float(np.max(abs(delta))), "rmse": float(np.sqrt(np.mean(delta**2)))}
            report["curves"].append(row)
            arrays[f"f_scalar_r_{ratio}"] = data["f"]
            arrays[f"f_full_r_{ratio}"] = exact
            print(json.dumps(row), flush=True)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=9)
    fig.suptitle("Your solution.png cuts compared with the explicit scalar formula\n"
                 r"$E_1=0.01$, $t_c=E_0=0.3$ meV; two cycles; same axes", fontsize=13)
    fig.savefig(OUT/"solution_vs_scalar_formula.png", dpi=190)
    fig.savefig(OUT/"solution_vs_scalar_formula.pdf")
    plt.close(fig)
    np.savez_compressed(DATA/"solution_vs_scalar_formula.npz", **arrays)
    (DATA/"solution_vs_scalar_formula.json").write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    main()
