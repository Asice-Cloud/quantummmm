#!/usr/bin/env python3
"""Phase-1 strict candidate reproduction for paper Fig.2.

This script keeps the panel semantics and control variables aligned with Fig.2:
(a) schematic (uniform), (b) spectrum+eta vs B, (c) E1/t1 vs B,
(d-f) braiding fidelity vs tau at three B values, (g-i) spectrum vs Ed.

Notes:
- Uses the repository's 6-Majorana effective Hamiltonian implementation.
- Provides a "nanowire-track" (full 6-Majorana evolution) and an
  "effective-track" (2-level reduced model) for direct panel comparison.
- This is a strict-candidate baseline, not final publication parity.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Tuple
from importlib import util as importlib_util

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from scipy.linalg import expm, eigh


def _load_majorana_module():
    this_dir = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(this_dir, "reproduce_effective_braiding_majorana.py")
    spec = importlib_util.spec_from_file_location("majorana_repro_mod", mod_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {mod_path}")
    module = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mod = _load_majorana_module()


@dataclass(frozen=True)
class Fig2Anchors:
    b: Tuple[float, float, float] = (2.25, 2.55, 2.75)
    e1: Tuple[float, float, float] = (3.368e-5, 5.500e-3, 4.885e-4)
    t1: Tuple[float, float, float] = (2.00e-4, 2.00e-3, 5.10e-3)


def build_interpolators(anchors: Fig2Anchors):
    b = np.array(anchors.b, dtype=float)
    e1 = np.array(anchors.e1, dtype=float)
    t1 = np.array(anchors.t1, dtype=float)
    return PchipInterpolator(b, e1), PchipInterpolator(b, t1)


def _gamma_ops(n_modes: int = 3):
    c_ops = mod.fermionic_operators(n_modes)
    gammas = []
    for c in c_ops:
        cd = c.conj().T
        gammas.append(c + cd)
        gammas.append(-1j * (c - cd))
    return gammas


def run_braiding_full_track(e1: float, tc: float, tau: float, steps_per_tau: int = 80, repeat: int = 2) -> float:
    fid, _, _, _ = mod.run_braiding({"E1": float(e1), "E2": 0.0}, tau=float(tau), steps_per_tau=steps_per_tau, repeat=repeat, tc=float(tc))
    return float(fid)


def run_braiding_effective_track(
    e1: float,
    tc: float,
    tau: float,
    steps_per_tau: int = 160,
    repeat: int = 2,
    alpha: float = 1.0,
) -> float:
    """Two-level reduced effective model for dashed-line comparison.

    H_eff(t) = e1 * sz + (t1+t2+t3)/2 * sx
    where (t1,t2,t3) follow the same gate protocol from the full model.
    """
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

    psi0 = np.array([1.0 + 0j, 0.0 + 0j])
    psi = psi0.copy()

    total_time = repeat * 3.0 * tau
    steps = max(10, int(repeat * 3.0 * tau * steps_per_tau))
    dt = total_time / steps

    t = 0.0
    for _ in range(steps):
        t += dt
        t1v, t2v, t3v, _ = mod.time_profiles(t, tau, tc)
        g = alpha * 0.5 * (t1v + t2v + t3v)
        h = e1 * sz + g * sx
        u = expm(-1j * h * dt)
        psi = u @ psi

    return float(abs(np.vdot(psi0, psi)))


def spectrum_vs_b(e1_fn, b_grid: np.ndarray, tc_ref: float = 0.03):
    gammas = _gamma_ops(3)
    all_e = []
    eta = []

    for b in b_grid:
        e1 = float(max(e1_fn(b), 0.0))
        # Static proxy point for spectrum panel: use finite couplings around protocol midpoint.
        h = mod.build_H(gammas, {"E1": e1, "E2": 0.0}, tc_ref, tc_ref, tc_ref, 0.0)
        ev = np.sort(np.real(np.linalg.eigvalsh(h)))
        all_e.append(ev)
        # Local-estimator proxy: ratio of low-energy crowding near zero.
        eta.append(float(np.exp(-40.0 * np.min(np.abs(ev)))))

    return np.array(all_e), np.array(eta)


def spectrum_vs_ed_for_b(e1: float, tc: float, ed_grid: np.ndarray):
    gammas = _gamma_ops(3)
    vals = []
    for ed in ed_grid:
        h = mod.build_H(gammas, {"E1": float(e1), "E2": 0.0}, float(tc), float(tc), float(tc), float(ed))
        ev = np.sort(np.real(np.linalg.eigvalsh(h)))
        vals.append(ev)
    return np.array(vals)


def calibrate_effective_alpha(
    e1: float,
    tc: float,
    tau_grid: np.ndarray,
    full_curve: np.ndarray,
    alpha_grid: np.ndarray,
) -> float:
    best_alpha = float(alpha_grid[0])
    best_rmse = np.inf
    for alpha in alpha_grid:
        eff = np.array([
            run_braiding_effective_track(
                e1=e1,
                tc=tc,
                tau=max(t, 0.2),
                steps_per_tau=160,
                repeat=2,
                alpha=float(alpha),
            )
            for t in tau_grid
        ])
        rmse = float(np.sqrt(np.mean((eff - full_curve) ** 2)))
        if rmse < best_rmse:
            best_rmse = rmse
            best_alpha = float(alpha)
    return best_alpha


def generate_fig2_strict_candidate(output_png: str, output_npz: str, output_txt: str):
    anchors = Fig2Anchors()
    e1_fn, t1_fn = build_interpolators(anchors)

    b_grid = np.linspace(0.0, 4.5, 220)
    evals_b, eta_b = spectrum_vs_b(e1_fn, b_grid, tc_ref=0.03)

    b_focus = np.linspace(anchors.b[0], anchors.b[-1], 120)
    e1_focus = np.maximum(e1_fn(b_focus), 0.0)
    t1_focus = np.maximum(t1_fn(b_focus), 0.0)

    tau_grid_long = np.linspace(0.0, 40.0, 120)
    tau_grid_short = np.linspace(0.0, 20.0, 90)
    tau_grids = [tau_grid_long, tau_grid_short, tau_grid_short]

    # Scale effective t1 anchors to protocol amplitude scale used by current simulator.
    # 5.1e-3 meV -> ~0.03 gives scale ~5.88.
    tc_scale = 5.88

    full_curves: Dict[float, np.ndarray] = {}
    eff_curves: Dict[float, np.ndarray] = {}
    spectrum_ed: Dict[float, np.ndarray] = {}

    alpha_fit: Dict[float, float] = {}

    for i, b0 in enumerate(anchors.b):
        e1 = float(anchors.e1[i])
        tc = float(anchors.t1[i] * tc_scale)
        tgrid = tau_grids[i]
        f_full = np.array([
            run_braiding_full_track(e1=e1, tc=tc, tau=max(t, 0.2), steps_per_tau=80, repeat=2)
            for t in tgrid
        ])

        alpha_grid = np.linspace(0.2, 2.5, 28)
        alpha_star = calibrate_effective_alpha(e1=e1, tc=tc, tau_grid=tgrid, full_curve=f_full, alpha_grid=alpha_grid)
        alpha_fit[b0] = alpha_star

        f_eff = np.array([
            run_braiding_effective_track(
                e1=e1,
                tc=tc,
                tau=max(t, 0.2),
                steps_per_tau=160,
                repeat=2,
                alpha=alpha_star,
            )
            for t in tgrid
        ])
        full_curves[b0] = f_full
        eff_curves[b0] = f_eff

    ed_grid = np.linspace(-1.0, 1.0, 220)
    for i, b0 in enumerate(anchors.b):
        e1 = float(anchors.e1[i])
        tc = float(anchors.t1[i] * tc_scale)
        spectrum_ed[b0] = spectrum_vs_ed_for_b(e1=e1, tc=tc, ed_grid=ed_grid)

    # Plot
    fig, axs = plt.subplots(3, 3, figsize=(13, 10))

    # (a) schematic
    ax = axs[0, 0]
    x = np.linspace(0, 2.3, 200)
    delta = np.full_like(x, 0.25)
    v = np.zeros_like(x)
    ax.plot(x, delta, color="tab:purple", lw=2, label=r"$\Delta(x)$")
    ax.plot(x, v, color="red", lw=2, label=r"$V(x)$")
    ax.set_xlim(0, 2.3)
    ax.set_ylim(0, 0.45)
    ax.set_title("(a) Uniform nanowire schematic")
    ax.set_xlabel("L (um)")
    ax.set_ylabel("E (meV)")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    # (b) spectrum + eta
    ax = axs[0, 1]
    for n in range(evals_b.shape[1]):
        ax.plot(b_grid, evals_b[:, n], color="blue", lw=0.9)
    ax.set_xlim(0, 4.5)
    ax.set_ylim(-0.3, 0.3)
    ax.set_xlabel("B (T)")
    ax.set_ylabel("E (meV)")
    ax.set_title("(b) Spectrum with local estimator")
    ax2 = ax.twinx()
    eta_norm = eta_b / max(np.max(eta_b), 1e-12)
    ax2.plot(b_grid, eta_norm, "--", color="brown", lw=1.8)
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel(r"$\eta$", color="brown")
    for b0, c in zip(anchors.b, ["black", "red", "limegreen"]):
        ax.axvline(b0, color=c, lw=1.2)

    # (c) E1 and t1 vs B
    ax = axs[0, 2]
    ax.plot(b_focus, e1_focus, "-o", color="tab:blue", markevery=20, ms=3, label=r"$E_1$")
    ax.plot(b_focus, t1_focus, "-o", color="tab:orange", markevery=20, ms=3, label=r"$t_1$")
    ax.set_xlim(anchors.b[0], anchors.b[-1])
    ymax = 1.15 * max(np.max(e1_focus), np.max(t1_focus))
    ax.set_ylim(0, ymax)
    ax.set_xlabel("B (T)")
    ax.set_ylabel("E (meV)")
    ax.set_title("(c) Extracted $E_1$ and effective $t_1$")
    ax.legend(frameon=False, fontsize=8)

    # (d-f) braiding vs tau
    panel_titles = ["(d)", "(e)", "(f)"]
    dots = ["black", "red", "limegreen"]
    for i, b0 in enumerate(anchors.b):
        ax = axs[1, i]
        tgrid = tau_grids[i]
        ax.plot(tgrid, eff_curves[b0], "--", color="magenta", lw=2, label="effective model")
        ax.plot(tgrid, full_curves[b0], "-", color="blue", lw=1.7, label="nanowire track")
        ax.scatter([tgrid[-10]], [0.92], s=160, color=dots[i], edgecolor="white", linewidth=0.5)
        ax.set_ylim(0, 1.02)
        ax.set_xlim(tgrid[0], tgrid[-1])
        ax.set_xlabel(r"$\tau$ (100/meV)")
        if i == 0:
            ax.set_ylabel(r"$|\langle \psi_i|\mathcal{G}(6\tau)|\psi_i^+\rangle|$")
        ax.set_title(f"{panel_titles[i]}  B={b0:.2f} T, alpha={alpha_fit[b0]:.2f}")
        if i == 0:
            ax.legend(frameon=False, fontsize=8, loc="upper left")

    # (g-i) spectrum vs Ed
    lower_titles = ["(g)", "(h)", "(i)"]
    for i, b0 in enumerate(anchors.b):
        ax = axs[2, i]
        vals = spectrum_ed[b0]
        # draw 4 states closest to zero to match paper-like low-energy focus
        idx = np.argsort(np.mean(np.abs(vals), axis=0))[:4]
        for k in idx:
            ax.plot(ed_grid, vals[:, k], color="blue", lw=1.2)
        lim = max(0.03, 1.05 * np.max(np.abs(vals[:, idx])))
        ax.set_ylim(-lim, lim)
        ax.set_xlim(ed_grid[0], ed_grid[-1])
        ax.set_xlabel(r"$E_d$ (meV)")
        ax.set_ylabel("E (meV)")
        ax.set_title(f"{lower_titles[i]}  Spectrum vs $E_d$, B={b0:.2f}")

    fig.suptitle("Fig.2 Strict-Candidate Reproduction (Phase-1)")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_png, dpi=220)
    plt.close(fig)

    # Save data for auditing
    np.savez(
        output_npz,
        b_grid=b_grid,
        evals_b=evals_b,
        eta_b=eta_b,
        b_focus=b_focus,
        e1_focus=e1_focus,
        t1_focus=t1_focus,
        tau_long=tau_grid_long,
        tau_short=tau_grid_short,
        full_b1=full_curves[anchors.b[0]],
        full_b2=full_curves[anchors.b[1]],
        full_b3=full_curves[anchors.b[2]],
        eff_b1=eff_curves[anchors.b[0]],
        eff_b2=eff_curves[anchors.b[1]],
        eff_b3=eff_curves[anchors.b[2]],
        ed_grid=ed_grid,
        spec_ed_b1=spectrum_ed[anchors.b[0]],
        spec_ed_b2=spectrum_ed[anchors.b[1]],
        spec_ed_b3=spectrum_ed[anchors.b[2]],
        alpha_b1=np.array([alpha_fit[anchors.b[0]]]),
        alpha_b2=np.array([alpha_fit[anchors.b[1]]]),
        alpha_b3=np.array([alpha_fit[anchors.b[2]]]),
        b_anchor=np.array(anchors.b),
        e1_anchor=np.array(anchors.e1),
        t1_anchor=np.array(anchors.t1),
        tc_scale=np.array([tc_scale]),
    )

    # Simple overlap metrics between full/effective curves
    def _curve_metrics(y_ref: np.ndarray, y_cmp: np.ndarray):
        rmse = float(np.sqrt(np.mean((y_ref - y_cmp) ** 2)))
        mae = float(np.mean(np.abs(y_ref - y_cmp)))
        corr = float(np.corrcoef(y_ref, y_cmp)[0, 1]) if y_ref.size > 2 else np.nan
        return rmse, mae, corr

    m1 = _curve_metrics(full_curves[anchors.b[0]], eff_curves[anchors.b[0]])
    m2 = _curve_metrics(full_curves[anchors.b[1]], eff_curves[anchors.b[1]])
    m3 = _curve_metrics(full_curves[anchors.b[2]], eff_curves[anchors.b[2]])

    with open(output_txt, "w", encoding="utf-8") as f:
        f.write("Fig.2 strict-candidate phase-1 metrics\n")
        f.write("This report compares full-track and effective-track fidelity curves.\n\n")
        f.write(f"B={anchors.b[0]:.2f}T: rmse={m1[0]:.6f}, mae={m1[1]:.6f}, corr={m1[2]:.6f}\n")
        f.write(f"  fitted alpha={alpha_fit[anchors.b[0]]:.6f}\n")
        f.write(f"B={anchors.b[1]:.2f}T: rmse={m2[0]:.6f}, mae={m2[1]:.6f}, corr={m2[2]:.6f}\n")
        f.write(f"  fitted alpha={alpha_fit[anchors.b[1]]:.6f}\n")
        f.write(f"B={anchors.b[2]:.2f}T: rmse={m3[0]:.6f}, mae={m3[1]:.6f}, corr={m3[2]:.6f}\n")
        f.write(f"  fitted alpha={alpha_fit[anchors.b[2]]:.6f}\n")
        f.write("\nAnchor values from paper caption used in this phase:\n")
        f.write(f"B anchors: {anchors.b}\n")
        f.write(f"E1 anchors (meV): {anchors.e1}\n")
        f.write(f"t1 anchors (meV): {anchors.t1}\n")


def main():
    os.makedirs("quantity", exist_ok=True)
    out_png = "quantity/fig2_strict_candidate_phase1.png"
    out_npz = "quantity/fig2_strict_candidate_phase1_data.npz"
    out_txt = "quantity/fig2_strict_candidate_phase1_metrics.txt"

    print("[Fig2-Strict] Generating phase-1 candidate...")
    generate_fig2_strict_candidate(out_png, out_npz, out_txt)
    print(f"[Fig2-Strict] Saved figure: {out_png}")
    print(f"[Fig2-Strict] Saved data:   {out_npz}")
    print(f"[Fig2-Strict] Saved report: {out_txt}")


if __name__ == "__main__":
    main()
