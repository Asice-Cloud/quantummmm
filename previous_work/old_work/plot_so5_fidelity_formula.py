#!/usr/bin/env python3
"""Plot the derived quadratic fidelity formula on Fig. 1(d)'s axes.

This evaluates an analytic formula only; no Hamiltonian propagation is run.
F2 = 1 - 2*(tau/hbar)^2*((1+2*C)*E1-B*abs(t1))^2.
The full plotting range extends beyond the formula's perturbative regime.
Unphysical negative values are masked, not clipped to zero.
"""
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/quantumsss-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.special import ellipe, ellipk


def main():
    root = Path(__file__).resolve().parent
    e1 = 0.01  # meV, following the caption of the paper's Fig. 1(d)
    hbar = 1.0
    c = float(ellipk(.5) / np.pi)
    b = float(.5 + ellipe(.5) / np.pi)
    a = 1 + 2 * c
    ridge_ratio = a / b
    ridge_y = np.log10(ridge_ratio)
    tau_axis = np.linspace(0, 12, 1201)  # tau in units of 100/meV
    log_ratio = np.linspace(-1, 1, 1001)
    physical_tau = 100 * tau_axis[None, :]
    t1 = e1 * 10 ** log_ratio[:, None]
    f2 = 1 - 2 * (physical_tau / hbar) ** 2 * (a * e1 - b * t1) ** 2
    plotted = np.ma.masked_less(f2, 0)

    np.savez_compressed(root / "so5_fidelity_formula_fig1d.npz",
                        tau_axis=tau_axis, tau_inverse_meV=100 * tau_axis,
                        log10_t1_over_E1=log_ratio, F_quadratic=f2,
                        E1_meV=e1, C=c, B=b,
                        cancellation_ratio=ridge_ratio)

    cmap = LinearSegmentedColormap.from_list(
        "paper_blue_yellow", ["#24218c", "#285ac6", "#22a6cf",
                              "#2bbf99", "#a8d94e", "#fff05a"])
    cmap.set_bad("#f4f4f4")
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.linewidth": .9}):
        fig, ax = plt.subplots(figsize=(8.6, 5.8))
        fig.subplots_adjust(left=.11, right=.88, bottom=.24, top=.82)
        im = ax.pcolormesh(tau_axis, log_ratio, plotted, shading="auto",
                           cmap=cmap, vmin=0, vmax=1, rasterized=True)
        ax.axhline(ridge_y, color="#972e25", linewidth=1, linestyle="--")
        ax.text(11.8, ridge_y + .07,
                rf"$|t_1|/E_1={ridge_ratio:.5f}$", ha="right", va="bottom",
                color="#972e25", fontsize=10,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
        ax.set(xlim=(0, 12), ylim=(-1, 1),
               xlabel=r"$\tau\ (100/\mathrm{meV})$",
               ylabel=r"$\lg(|t_1|/E_1)$")
        ax.set_xticks(np.arange(0, 13, 2))
        ax.set_yticks(np.linspace(-1, 1, 5))
        cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.035)
        cb.set_ticks(np.linspace(0, 1, 6))
        cb.set_label(r"$F^{(2)}$ (squared state fidelity)")
        fig.suptitle("Derived SO(5) fidelity formula on Fig. 1(d) axes",
                     y=.97, fontsize=13)
        fig.text(.5, .88,
                 r"$F^{(2)}=1-2(\tau/\hbar)^2"
                 r"[2.180340599\,E_1-0.929923300\,|t_1|]^2$",
                 ha="center", fontsize=13)
        fig.text(.11, .14,
                 r"$E_1=0.01$ meV; equal ideal control amplitudes; "
                 r"$\hbar=1$; two three-step exchanges.", fontsize=10)
        fig.text(.11, .092,
                 r"Blank region: $F^{(2)}<0$ (unphysical extrapolation). "
                 "Dashed line: quadratic cancellation.", fontsize=10)
        fig.text(.11, .048,
                 "The formula requires adiabatic control and small accumulated perturbations;\n"
                 "a displayed value between 0 and 1 does not by itself guarantee validity.",
                 fontsize=9, color="#555555")
        for extension in ["png", "pdf"]:
            fig.savefig(root / f"so5_fidelity_formula_fig1d.{extension}", dpi=220)
        plt.close(fig)

    # Check the plotted analytic expression at its defining special points.
    assert np.all(f2[:, 0] == 1)
    assert abs(1 - 2 * 1200 ** 2 * (a * e1 - b * ridge_ratio * e1) ** 2 - 1) < 1e-12
    print(f"Created so5_fidelity_formula_fig1d.png and .pdf; cancellation ratio={ridge_ratio:.9f}")
    print(f"Negative-formula fraction masked: {np.mean(f2 < 0):.2%}")


if __name__ == "__main__":
    main()
