#!/usr/bin/env python3
"""Finite-difference parameter sensitivity for verify_path_braid model.

Saves CSV and NPZ results to results/path_braid/.
"""

from __future__ import annotations

import json
from importlib import util
from pathlib import Path
import sys
import csv
import numpy as np

# Load verify_path_braid.py as a module by path (robust import)
HERE = Path(__file__).resolve().parent
spec = util.spec_from_file_location("verify_path_braid", str(HERE / "verify_path_braid.py"))
vpb = util.module_from_spec(spec)
spec.loader.exec_module(vpb)

# Baseline parameters (match recent runs)
BASELINE = {
    "steps": 300,
    "T": 3.175,
    "tc": 1.0,
    "hc_factor": 2.125,
    "zeeman": 0.0,
    "zeeman_drive_1": 0.0,
    "zeeman_drive_2": 0.0,
    "zeeman_drive_3": 0.0,
    "zeeman_drive_shape": "cos",
    "phase_comp": 0.0,
    "cross_comp_1": 0.0,
    "cross_comp_2": 0.0,
    "cross_comp_3": 0.0,
    "axis": "z",
}

# Parameters to test
PARAM_KEYS = [
    "T",
    "tc",
    "hc_factor",
    "phase_comp",
    "zeeman",
    "zeeman_drive_1",
    "zeeman_drive_2",
    "zeeman_drive_3",
    "cross_comp_1",
    "cross_comp_2",
    "cross_comp_3",
]

# Evaluate helper
def evaluate(params: dict) -> dict:
    gates = vpb.compose_cycle(
        steps=int(params["steps"]),
        T=float(params["T"]),
        tc=float(params["tc"]),
        hc_factor=float(params["hc_factor"]),
        zeeman=float(params["zeeman"]),
        zeeman_drive_1=float(params["zeeman_drive_1"]),
        zeeman_drive_2=float(params["zeeman_drive_2"]),
        zeeman_drive_3=float(params["zeeman_drive_3"]),
        zeeman_drive_shape=params["zeeman_drive_shape"],
        phase_comp=float(params["phase_comp"]),
        cross_comp_1=float(params["cross_comp_1"]),
        cross_comp_2=float(params["cross_comp_2"]),
        cross_comp_3=float(params["cross_comp_3"]),
    )
    diag = vpb.analyze_total_gate(gates["R_total"], params["axis"])
    # Keep only numeric diagnostics
    return {
        "score": float(diag["score"]),
        "target_residual": float(diag["target_residual"]),
        "braid_res": float(diag["braid_res"]),
        "ybe_res": float(diag["ybe_res"]),
        "best_axis": str(diag["best_axis"]),
    }


def main():
    outdir = Path("results/path_braid")
    outdir.mkdir(parents=True, exist_ok=True)

    baseline_diag = evaluate(BASELINE)
    base_score = baseline_diag["score"]

    results = []

    for key in PARAM_KEYS:
        v0 = float(BASELINE[key])
        # choose delta: 2% of scale, minimum absolute 0.02
        delta = 0.02 * max(1.0, abs(v0))
        p_plus = dict(BASELINE)
        p_minus = dict(BASELINE)
        p_plus[key] = v0 + delta
        p_minus[key] = v0 - delta

        diag_p = evaluate(p_plus)
        diag_m = evaluate(p_minus)

        score_p = diag_p["score"]
        score_m = diag_m["score"]
        sens = (score_p - score_m) / (2.0 * delta)

        results.append(
            {
                "param": key,
                "delta": float(delta),
                "sens_score": float(sens),
                "base_score": float(base_score),
                "score_plus": float(score_p),
                "score_minus": float(score_m),
                "target_plus": float(diag_p["target_residual"]),
                "braid_plus": float(diag_p["braid_res"]),
                "ybe_plus": float(diag_p["ybe_res"]),
                "target_minus": float(diag_m["target_residual"]),
                "braid_minus": float(diag_m["braid_res"]),
                "ybe_minus": float(diag_m["ybe_res"]),
            }
        )

    # Sort by absolute sensitivity
    results_sorted = sorted(results, key=lambda r: abs(r["sens_score"]), reverse=True)

    # Write CSV
    csv_path = outdir / "sensitivity_matrix_steps300_T3p175.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        header = [
            "param",
            "delta",
            "sens_score",
            "base_score",
            "score_plus",
            "score_minus",
            "target_plus",
            "braid_plus",
            "ybe_plus",
            "target_minus",
            "braid_minus",
            "ybe_minus",
        ]
        writer.writerow(header)
        for r in results_sorted:
            writer.writerow([r[h] for h in header])

    # Save NPZ (pickleable object array)
    npz_path = outdir / "sensitivity_matrix_steps300_T3p175.npz"
    np.savez(npz_path, baseline=baseline_diag, results=np.array(results_sorted, dtype=object))

    # Print concise table
    print(f"Baseline score = {base_score:.6g}")
    print("param, delta, sens_score, score_plus, score_minus")
    for r in results_sorted:
        print(f"{r['param']}, {r['delta']:.3g}, {r['sens_score']:.3e}, {r['score_plus']:.6g}, {r['score_minus']:.6g}")

    print("Saved:")
    print(f" - CSV -> {csv_path}")
    print(f" - NPZ -> {npz_path}")


if __name__ == '__main__':
    main()
