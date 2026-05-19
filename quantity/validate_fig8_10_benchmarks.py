#!/usr/bin/env python3
"""Lightweight validation for the comparable Fig.8-10 paper benchmarks.

This script evaluates only the observables that are comparable with the current
Pauli-tensor model:
- spectra vs magnetic field / disorder parameter
- zero-bias conductance maps
- sampled non-Abelian indicator values from the current model

It intentionally avoids paper-level strict reproduction claims.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from dataclasses import dataclass
from itertools import combinations

import numpy as np
import matplotlib.pyplot as plt


def _load_model_module():
    this_dir = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(this_dir, "reproduce_fig8_fig9_fig10_pauli.py")
    spec = importlib.util.spec_from_file_location("pauli_fig_mod", mod_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module from {mod_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


mod = _load_model_module()


@dataclass(frozen=True)
class Metric:
    name: str
    value: float
    threshold: float | None
    passed: bool | None
    detail: str


def rel_l2(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = float(np.linalg.norm(a) + np.linalg.norm(b) + 1e-15)
    return float(np.linalg.norm(a - b) / denom)


def mean_zero_bias(G, E_grid):
    idx = int(np.argmin(np.abs(E_grid)))
    return float(np.mean(G[idx, :]))


def conductance_summary(Evals, Evecs, E_grid, gamma=0.01, lead_indices=(0,)):
    G = mod.conductance_map_from_spectrum(Evals, Evecs, E_grid, gamma=gamma, lead_indices=lead_indices)
    zbp = mean_zero_bias(G, E_grid)
    return G, zbp


def fig8_metrics():
    B = np.linspace(0.0, 4.5, 120)
    E = np.linspace(-0.3, 0.3, 181)
    pars_uniform = {"delta": 0.10}
    pars_inhom = {"delta": 0.04}

    spec_u, vecs_u = mod.pauli_spectrum_vs_B(pars_uniform, B, return_vecs=True)
    spec_i, vecs_i = mod.pauli_spectrum_vs_B(pars_inhom, B, return_vecs=True)
    G_u, zbp_u = conductance_summary(spec_u, vecs_u, E, gamma=0.012)
    G_i, zbp_i = conductance_summary(spec_i, vecs_i, E, gamma=0.009)

    u_samples = np.array([0.5, 2.0, 4.0])
    N_u = np.array([mod.compute_N(pars_uniform, u=float(u), v=0.7)[0] for u in u_samples])
    N_i = np.array([mod.compute_N(pars_inhom, u=float(u), v=0.7)[0] for u in u_samples])

    metrics = [
        Metric("Fig.8 spectrum relative difference", rel_l2(spec_u, spec_i), 0.05, None, "uniform vs inhomogeneous spectrum contrast"),
        Metric("Fig.8 conductance relative difference", rel_l2(G_u, G_i), 0.05, None, "uniform vs inhomogeneous conductance contrast"),
        Metric("Fig.8 zero-bias conductance ratio", zbp_i / max(zbp_u, 1e-12), 1.0, None, "inhomogeneous / uniform at E≈0"),
        Metric("Fig.8 N-sample relative difference", rel_l2(N_u, N_i), 0.02, None, "sampled kernel indicator difference"),
        Metric("Fig.8 N monotonicity (uniform)", float(np.sum(np.diff(N_u) >= -1e-12) / max(len(N_u) - 1, 1)), 0.66, None, "fraction of nondecreasing intervals"),
        Metric("Fig.8 N monotonicity (inhomogeneous)", float(np.sum(np.diff(N_i) >= -1e-12) / max(len(N_i) - 1, 1)), 0.66, None, "fraction of nondecreasing intervals"),
    ]
    return {
        "B": B,
        "E": E,
        "spec_u": spec_u,
        "spec_i": spec_i,
        "vecs_u": vecs_u,
        "vecs_i": vecs_i,
        "G_u": G_u,
        "G_i": G_i,
        "zbp_u": zbp_u,
        "zbp_i": zbp_i,
        "u_samples": u_samples,
        "N_u": N_u,
        "N_i": N_i,
        "metrics": metrics,
    }


def fig9_metrics():
    B = np.linspace(0.0, 4.5, 110)
    E = np.linspace(-0.3, 0.3, 171)
    lead_choices = [
        (0,),
        (1,),
        (0, 1),
        (0, 3),
        (0, 1, 2, 3),
    ]
    scenarios = [
        {"name": "Uniform", "delta": 0.05},
        {"name": "Smooth distortion", "delta": 0.10},
        {"name": "QD-like", "delta": 0.14},
        {"name": "Disorder-like", "delta": 0.18},
        {"name": "Strong distortion", "delta": 0.22},
        {"name": "Bulk disorder", "delta": 0.28},
    ]
    lead_maps = {lead: [] for lead in lead_choices}
    lead_zbps = {lead: [] for lead in lead_choices}
    nvals = []
    spectra = []
    vecsets = []
    for sc in scenarios:
        spec, vecs = mod.pauli_spectrum_vs_B({"delta": sc["delta"]}, B, return_vecs=True)
        spectra.append(spec)
        vecsets.append(vecs)
        for lead in lead_choices:
            G, zbp = conductance_summary(spec, vecs, E, gamma=0.01, lead_indices=lead)
            lead_maps[lead].append(G)
            lead_zbps[lead].append(zbp)
        nvals.append(mod.compute_N({"delta": sc["delta"]}, u=1.0, v=0.7)[0])

    nvals = np.array(nvals)
    lead_pairwise = {}
    for lead in lead_choices:
        pairwise = []
        for i, j in combinations(range(len(scenarios)), 2):
            pairwise.append(rel_l2(lead_maps[lead][i], lead_maps[lead][j]))
        lead_pairwise[lead] = np.array(pairwise)

    # Sensitivity to the most separated scenario pair for each lead.
    delta_pairwise = {}
    for lead in lead_choices:
        delta_pairwise[lead] = rel_l2(lead_maps[lead][0], lead_maps[lead][-1])

    metrics = [
        Metric("Fig.9 spectrum spread across scenarios", float(np.mean([rel_l2(spectra[i], spectra[j]) for i, j in combinations(range(len(scenarios)), 2)])), 0.01, None, "mean pairwise spectral distance"),
        Metric("Fig.9 N spread across scenarios", float(np.std(nvals) / (np.mean(nvals) + 1e-15)), 0.05, None, "relative spread of sampled N values"),
    ]
    for lead in lead_choices:
        metrics.append(
            Metric(
                f"Fig.9 conductance spread lead={lead}",
                float(np.mean(lead_pairwise[lead])),
                0.03,
                None,
                f"mean pairwise conductance-map distance for lead={lead}",
            )
        )
        metrics.append(
            Metric(
                f"Fig.9 strong-pair diff lead={lead}",
                float(delta_pairwise[lead]),
                0.05,
                None,
                f"delta=0.05 vs 0.28 map distance for lead={lead}",
            )
        )
    spectra = np.array(spectra)
    for lead in lead_choices:
        lead_maps[lead] = np.array(lead_maps[lead])
        lead_zbps[lead] = np.array(lead_zbps[lead])
    return {
        "B": B,
        "E": E,
        "scenarios": scenarios,
        "lead_choices": lead_choices,
        "lead_maps": lead_maps,
        "lead_zbps": lead_zbps,
        "lead_pairwise": lead_pairwise,
        "delta_pairwise": delta_pairwise,
        "spectra": spectra,
        "nvals": nvals,
        "metrics": metrics,
    }


def fig10_metrics():
    B = np.linspace(0.0, 4.5, 120)
    E = np.linspace(-0.3, 0.3, 181)
    B_markers = [1.84, 1.93, 2.02]
    delta_values = [0.08, 0.14, 0.22]

    spec_d, vecs_d = mod.pauli_spectrum_vs_B({"delta": 0.22}, B, return_vecs=True)
    G_d, zbp_d = conductance_summary(spec_d, vecs_d, E, gamma=0.012)

    taus = np.array([0.6, 2.0, 4.0, 6.0])
    curves = []
    aucs = []
    for delta_val in delta_values:
        curve = np.array([mod.compute_N({"delta": delta_val}, u=float(u), v=0.7)[0] for u in taus])
        curves.append(curve)
        aucs.append(float(np.trapezoid(curve, taus)))
    curves = np.array(curves)
    aucs = np.array(aucs)

    delta_grid = np.linspace(0.02, 0.30, 80)
    spectrum_shift = []
    for Bval in B_markers:
        u_fixed = mod.map_B_to_u(Bval)
        spec_vs_delta = []
        for delta_val in delta_grid:
            Htmp = mod.effective_hamiltonian_pauli(u_fixed, {"delta": delta_val})
            spec_vs_delta.append(np.sort(np.real(np.linalg.eigvalsh(Htmp))))
        spec_vs_delta = np.array(spec_vs_delta)
        spectrum_shift.append(float(np.mean(np.abs(spec_vs_delta[-1] - spec_vs_delta[0]))))
    spectrum_shift = np.array(spectrum_shift)

    metrics = [
        Metric("Fig.10 zero-bias conductance", zbp_d, 0.05, None, "disorder conductance map zero-bias level"),
        Metric("Fig.10 N curve spread across delta", float(np.std(aucs) / (np.mean(aucs) + 1e-15)), 0.05, None, "AUC spread over delta values"),
        Metric("Fig.10 N ordering monotonicity", float(np.sum(np.diff(aucs) >= -1e-12) / max(len(aucs) - 1, 1)), 0.66, None, "AUC ordering across delta values"),
        Metric("Fig.10 spectrum shift vs delta", float(np.mean(spectrum_shift)), 0.01, None, "mean absolute spectral shift between delta extremes"),
    ]
    return {
        "B": B,
        "E": E,
        "B_markers": B_markers,
        "delta_values": delta_values,
        "G_d": G_d,
        "zbp_d": zbp_d,
        "taus": taus,
        "curves": curves,
        "aucs": aucs,
        "delta_grid": delta_grid,
        "spectrum_shift": spectrum_shift,
        "metrics": metrics,
    }


def render_fig8_fixed_figure(data8, output_png):
    best_lead = (1,)
    gamma = 0.01
    G_u_corr, zbp_u_corr = conductance_summary(data8["spec_u"], data8["vecs_u"], data8["E"], gamma=gamma, lead_indices=best_lead)
    G_i_corr, zbp_i_corr = conductance_summary(data8["spec_i"], data8["vecs_i"], data8["E"], gamma=gamma, lead_indices=best_lead)
    ratio_corr = zbp_i_corr / max(zbp_u_corr, 1e-12)

    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    ax = axs[0, 0]
    im0 = ax.imshow(G_u_corr, origin="lower", aspect="auto", extent=[data8["B"][0], data8["B"][-1], data8["E"][0], data8["E"][-1]], cmap="turbo", vmin=0, vmax=1)
    ax.set_title(f"Uniform conductance (lead={best_lead})")
    ax.set_xlabel("B")
    ax.set_ylabel("E")

    ax = axs[0, 1]
    im1 = ax.imshow(G_i_corr, origin="lower", aspect="auto", extent=[data8["B"][0], data8["B"][-1], data8["E"][0], data8["E"][-1]], cmap="turbo", vmin=0, vmax=1)
    ax.set_title(f"Inhomogeneous conductance (lead={best_lead})")
    ax.set_xlabel("B")
    ax.set_ylabel("E")

    cbar = fig.colorbar(im1, ax=[axs[0, 0], axs[0, 1]], fraction=0.03, pad=0.02)
    cbar.set_label("G (normalized)")

    ax = axs[1, 0]
    ax.bar(["uniform", "inhomogeneous"], [zbp_u_corr, zbp_i_corr], color=["tab:blue", "tab:orange"])
    ax.set_title(f"Zero-bias comparison, ratio={ratio_corr:.3f}")
    ax.set_ylabel("mean G(E≈0)")

    ax = axs[1, 1]
    ax.plot(data8["u_samples"], data8["N_u"], "-o", label="uniform")
    ax.plot(data8["u_samples"], data8["N_i"], "-o", label="inhomogeneous")
    ax.set_title("Sampled N(u)")
    ax.set_xlabel("u")
    ax.set_ylabel("N")
    ax.legend(frameon=False)

    fig.suptitle(f"Fig.8 corrected reproduction using lead={best_lead} and common gamma={gamma:.3f}", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def render_summary_figure(data8, data9, data10, output_png):
    fig, axs = plt.subplots(2, 3, figsize=(15, 9))

    ax = axs[0, 0]
    labels = ["spec diff", "cond diff", "zbp ratio", "N diff"]
    values = [
        data8["metrics"][0].value,
        data8["metrics"][1].value,
        data8["metrics"][2].value,
        data8["metrics"][3].value,
    ]
    ax.bar(labels, values, color=["tab:blue", "tab:blue", "tab:orange", "tab:green"])
    ax.set_title("Fig.8 validation metrics")
    ax.tick_params(axis='x', rotation=25)
    ax.set_ylabel("value")

    ax = axs[0, 1]
    lead_labels = [str(lead) for lead in data9["lead_choices"]]
    lead_sens = [data9["delta_pairwise"][lead] for lead in data9["lead_choices"]]
    ax.bar(lead_labels, lead_sens, color="tab:purple")
    ax.set_title("Fig.9 lead sensitivity (delta=0.05 vs 0.28)")
    ax.tick_params(axis='x', rotation=30)
    ax.set_ylabel("relative map diff")

    ax = axs[0, 2]
    ax.bar([str(d) for d in data10["delta_values"]], data10["aucs"], color="tab:red")
    ax.set_title("Fig.10 N-curve AUC vs delta")
    ax.set_xlabel("delta")
    ax.set_ylabel("AUC")

    ax = axs[1, 0]
    ax.plot(data8["u_samples"], data8["N_u"], "-o", label="uniform")
    ax.plot(data8["u_samples"], data8["N_i"], "-o", label="inhomogeneous")
    ax.set_title("Fig.8 sampled N(u)")
    ax.set_xlabel("u")
    ax.set_ylabel("N")
    ax.legend(frameon=False)

    ax = axs[1, 1]
    best_lead = max(data9["lead_choices"], key=lambda lead: data9["delta_pairwise"][lead])
    for i, sc in enumerate(data9["scenarios"]):
        ax.plot(data9["B"], data9["lead_maps"][best_lead][i][data9["lead_maps"][best_lead][i].shape[0] // 2], lw=1.5, label=sc["name"])
    ax.set_title(f"Fig.9 middle-line slices (best lead={best_lead})")
    ax.set_xlabel("B")
    ax.set_ylabel("G(E≈0,B)")
    ax.legend(frameon=False, fontsize=7)

    ax = axs[1, 2]
    for i, delta in enumerate(data10["delta_values"]):
        ax.plot(data10["taus"], data10["curves"][i], "-o", label=fr"$\delta={delta:.2f}$")
    ax.set_title("Fig.10 sampled N(τ)")
    ax.set_xlabel("τ")
    ax.set_ylabel("N")
    ax.legend(frameon=False, fontsize=7)

    fig.suptitle("Fig.8-10 benchmark validation summary", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def render_fig9_fixed_figure(data9, output_png):
    best_lead = max(data9["lead_choices"], key=lambda lead: data9["delta_pairwise"][lead])
    fig, axs = plt.subplots(2, 3, figsize=(12, 7))
    last_im = None

    for k, sc in enumerate(data9["scenarios"]):
        i, j = divmod(k, 3)
        ax = axs[i, j]
        last_im = ax.imshow(
            data9["lead_maps"][best_lead][k],
            origin="lower",
            aspect="auto",
            extent=[data9["B"][0], data9["B"][-1], data9["E"][0], data9["E"][-1]],
            cmap="turbo",
            vmin=0,
            vmax=1,
        )
        ax.set_title(f"{sc['name']} (lead={best_lead})")
        ax.set_xlabel("B (T)")
        ax.set_ylabel("E (meV)")

    cbar = fig.colorbar(last_im, ax=axs.ravel().tolist(), fraction=0.02, pad=0.02)
    cbar.set_label("G_{zb} (arb. units), normalized")
    fig.suptitle(f"Fig.9 corrected reproduction using best lead={best_lead}")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def main():
    os.makedirs("quantity", exist_ok=True)
    report_path = os.path.join("quantity", "fig8_10_validation_report.txt")
    figure_path = os.path.join("quantity", "fig8_10_validation_summary.png")
    fig8_path = os.path.join("quantity", "fig8_corrected_best_lead.png")
    fig9_path = os.path.join("quantity", "fig9_corrected_best_lead.png")

    data8 = fig8_metrics()
    data9 = fig9_metrics()
    data10 = fig10_metrics()

    metrics = data8["metrics"] + data9["metrics"] + data10["metrics"]

    report_lines = []
    report_lines.append("Fig.8-10 benchmark validation report")
    report_lines.append("=" * 80)
    report_lines.append("Scope: only comparable observables are validated; no strict paper-level equivalence is claimed.")
    report_lines.append("")

    for item in metrics:
        if item.threshold is None:
            status = "INFO"
        else:
            status = "PASS" if item.value >= item.threshold else "PARTIAL"
        report_lines.append(f"[{status}] {item.name}: value={item.value:.6e}, threshold={item.threshold} | {item.detail}")

    report_lines.append("")
    report_lines.append("Fig.8 details:")
    report_lines.append(f"  zero-bias conductance uniform    = {data8['zbp_u']:.6e}")
    report_lines.append(f"  zero-bias conductance inhomog.   = {data8['zbp_i']:.6e}")
    report_lines.append(f"  sampled N(u) uniform             = {np.array2string(data8['N_u'], precision=6)}")
    report_lines.append(f"  sampled N(u) inhomogeneous       = {np.array2string(data8['N_i'], precision=6)}")

    report_lines.append("")
    report_lines.append("Fig.9 details:")
    report_lines.append(f"  sampled N values                 = {np.array2string(data9['nvals'], precision=6)}")
    best_lead = max(data9["lead_choices"], key=lambda lead: data9["delta_pairwise"][lead])
    report_lines.append(f"  corrected reproduction lead      = {best_lead}")
    for lead in data9["lead_choices"]:
        report_lines.append(f"  lead={lead} mean pairwise map diff = {np.mean(data9['lead_pairwise'][lead]):.6e}")
        report_lines.append(f"  lead={lead} delta pairwise diff    = {data9['delta_pairwise'][lead]:.6e}")

    report_lines.append("")
    report_lines.append("Fig.10 details:")
    report_lines.append(f"  zero-bias conductance (disorder) = {data10['zbp_d']:.6e}")
    report_lines.append(f"  N curve AUCs                     = {np.array2string(data10['aucs'], precision=6)}")
    report_lines.append(f"  spectrum shift vs delta          = {np.array2string(data10['spectrum_shift'], precision=6)}")

    report = "\n".join(report_lines)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report)

    render_summary_figure(data8, data9, data10, figure_path)
    render_fig8_fixed_figure(data8, fig8_path)
    render_fig9_fixed_figure(data9, fig9_path)

    print(report)
    print(f"\nSaved report: {report_path}")
    print(f"Saved figure: {figure_path}")
    print(f"Saved Fig.8 correction figure: {fig8_path}")
    print(f"Saved Fig.9 correction figure: {fig9_path}")


if __name__ == "__main__":
    main()