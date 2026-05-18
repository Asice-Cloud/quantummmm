#!/usr/bin/env python3
"""Cross-validate the candidate point (tc=1.0, cross_comp_2=0.00833333).
Evaluate combinations: steps in {80,300} and zeeman in {0, +0.01, -0.01}.
Save CSV and NPZ to results/path_braid/.
"""

from __future__ import annotations

from importlib import util
from pathlib import Path
import csv
import numpy as np

HERE = Path(__file__).resolve().parent
spec = util.spec_from_file_location("verify_path_braid", str(HERE / "verify_path_braid.py"))
vpb = util.module_from_spec(spec)
spec.loader.exec_module(vpb)

OUTDIR = Path("results/path_braid")
OUTDIR.mkdir(parents=True, exist_ok=True)

CANDIDATE = {
    "tc": 1.0,
    "cross_comp_2": 0.00833333,
    "T": 3.175,
    "hc_factor": 2.125,
    "axis": "z",
}

STEPS = [80, 300]
ZEEMANS = [0.0, 0.01, -0.01]

records = []

for steps in STEPS:
    for zeeman in ZEEMANS:
        gates = vpb.compose_cycle(
            steps=int(steps),
            T=float(CANDIDATE["T"]),
            tc=float(CANDIDATE["tc"]),
            hc_factor=float(CANDIDATE["hc_factor"]),
            zeeman=float(zeeman),
            zeeman_drive_1=0.0,
            zeeman_drive_2=0.0,
            zeeman_drive_3=0.0,
            zeeman_drive_shape="cos",
            phase_comp=0.0,
            cross_comp_1=0.0,
            cross_comp_2=float(CANDIDATE["cross_comp_2"]),
            cross_comp_3=0.0,
        )
        diag = vpb.analyze_total_gate(gates["R_total"], CANDIDATE["axis"])
        records.append(
            {
                "steps": int(steps),
                "zeeman": float(zeeman),
                "tc": float(CANDIDATE["tc"]),
                "cross_comp_2": float(CANDIDATE["cross_comp_2"]),
                "score": float(diag["score"]),
                "target_residual": float(diag["target_residual"]),
                "braid_res": float(diag["braid_res"]),
                "ybe_res": float(diag["ybe_res"]),
                "best_axis": str(diag["best_axis"]),
            }
        )

# Aggregate per steps
summary = []
for steps in STEPS:
    subset = [r for r in records if r["steps"] == steps]
    scores = np.array([r["score"] for r in subset])
    targets = np.array([r["target_residual"] for r in subset])
    braids = np.array([r["braid_res"] for r in subset])
    ybes = np.array([r["ybe_res"] for r in subset])
    summary.append(
        {
            "steps": int(steps),
            "mean_score": float(scores.mean()),
            "std_score": float(scores.std()),
            "mean_target": float(targets.mean()),
            "std_target": float(targets.std()),
            "mean_braid": float(braids.mean()),
            "std_braid": float(braids.std()),
            "mean_ybe": float(ybes.mean()),
            "std_ybe": float(ybes.std()),
        }
    )

# Save CSV
csv_path = OUTDIR / f"crossval_candidate_tc{CANDIDATE['tc']}_cross2{CANDIDATE['cross_comp_2']}_steps{STEPS[0]}_{STEPS[1]}.csv"
with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["steps", "zeeman", "tc", "cross_comp_2", "score", "target_residual", "braid_res", "ybe_res", "best_axis"])
    for r in records:
        writer.writerow([r["steps"], r["zeeman"], r["tc"], r["cross_comp_2"], r["score"], r["target_residual"], r["braid_res"], r["ybe_res"], r["best_axis"]])

# Save NPZ
npz_path = OUTDIR / f"crossval_candidate_tc{CANDIDATE['tc']}_cross2{CANDIDATE['cross_comp_2']}_steps{STEPS[0]}_{STEPS[1]}.npz"
np.savez(npz_path, records=np.array(records, dtype=object), summary=np.array(summary, dtype=object))

# Print concise summary
print(f"Performed {len(records)} evaluations")
for s in summary:
    print(f"steps={s['steps']}: mean_score={s['mean_score']:.6g} ± {s['std_score']:.3g}, mean_target={s['mean_target']:.6g} ± {s['std_target']:.3g}, mean_braid={s['mean_braid']:.6g} ± {s['std_braid']:.3g}, mean_ybe={s['mean_ybe']:.6g} ± {s['std_ybe']:.3g}")

print(f"Saved CSV -> {csv_path}")
print(f"Saved NPZ -> {npz_path}")

if __name__ == '__main__':
    pass
