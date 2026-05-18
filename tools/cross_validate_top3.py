#!/usr/bin/env python3
"""Cross-validate top-N candidates from a scan CSV.

Defaults to reading results/path_braid/scan_tc_cross2_steps300_T3.175.csv
and running steps in {80,300} and zeeman in {0, +0.01, -0.01}.
Saves per-candidate CSV/NPZ and an aggregated CSV/NPZ.
"""
from __future__ import annotations

import argparse
import csv
from importlib import util
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
spec = util.spec_from_file_location("verify_path_braid", str(HERE / "verify_path_braid.py"))
vpb = util.module_from_spec(spec)
spec.loader.exec_module(vpb)


def slug(value: float) -> str:
    return format(value, "g").replace("-", "m").replace(".", "p")


def read_scan_csv(path: Path):
    rows = []
    with open(path, "r") as f:
        r = csv.DictReader(f)
        for row in r:
            row_parsed = {k: v for k, v in row.items()}
            row_parsed["score"] = float(row_parsed["score"])
            row_parsed["tc"] = float(row_parsed["tc"])
            row_parsed["cross_comp_2"] = float(row_parsed["cross_comp_2"])
            rows.append(row_parsed)
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-csv", type=str, default="results/path_braid/scan_tc_cross2_steps300_T3.175.csv")
    p.add_argument("--top-n", type=int, default=3)
    p.add_argument("--steps", type=str, default="80,300")
    p.add_argument("--zeemans", type=str, default="0,0.01,-0.01")
    p.add_argument("--T", type=float, default=3.175)
    p.add_argument("--hc-factor", type=float, default=2.125)
    p.add_argument("--axis", type=str, default="z")
    p.add_argument("--outdir", type=str, default="results/path_braid")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    steps_list = [int(x) for x in args.steps.split(",") if x.strip()]
    zeeman_list = [float(x) for x in args.zeemans.split(",") if x.strip()]

    scan_rows = read_scan_csv(Path(args.input_csv))
    scan_sorted = sorted(scan_rows, key=lambda r: r["score"])[: args.top_n]

    agg_records = []
    per_candidate_summaries = []

    for idx, cand in enumerate(scan_sorted, start=1):
        tc = float(cand["tc"])
        cross2 = float(cand["cross_comp_2"])
        records = []
        for steps in steps_list:
            for zeeman in zeeman_list:
                gates = vpb.compose_cycle(
                    steps=int(steps),
                    T=float(args.T),
                    tc=float(tc),
                    hc_factor=float(args.hc_factor),
                    zeeman=float(zeeman),
                    zeeman_drive_1=0.0,
                    zeeman_drive_2=0.0,
                    zeeman_drive_3=0.0,
                    zeeman_drive_shape="cos",
                    phase_comp=0.0,
                    cross_comp_1=0.0,
                    cross_comp_2=float(cross2),
                    cross_comp_3=0.0,
                )
                diag = vpb.analyze_total_gate(gates["R_total"], args.axis)
                rec = {
                    "candidate_idx": idx,
                    "tc": float(tc),
                    "cross_comp_2": float(cross2),
                    "steps": int(steps),
                    "zeeman": float(zeeman),
                    "score": float(diag["score"]),
                    "target_residual": float(diag["target_residual"]),
                    "braid_res": float(diag["braid_res"]),
                    "ybe_res": float(diag["ybe_res"]),
                    "best_axis": str(diag["best_axis"]),
                }
                records.append(rec)
                agg_records.append(rec)

        # per-candidate summary
        scores = np.array([r["score"] for r in records])
        targets = np.array([r["target_residual"] for r in records])
        braids = np.array([r["braid_res"] for r in records])
        ybes = np.array([r["ybe_res"] for r in records])
        summary = {
            "candidate_idx": idx,
            "tc": float(tc),
            "cross_comp_2": float(cross2),
            "mean_score": float(scores.mean()),
            "std_score": float(scores.std()),
            "mean_target": float(targets.mean()),
            "std_target": float(targets.std()),
            "mean_braid": float(braids.mean()),
            "std_braid": float(braids.std()),
            "mean_ybe": float(ybes.mean()),
            "std_ybe": float(ybes.std()),
        }
        per_candidate_summaries.append(summary)

        # save per-candidate CSV/NPZ
        slug_tc = slug(tc)
        slug_cross = slug(cross2)
        csv_path = outdir / f"crossval_top{idx}_tc{slug_tc}_cross2{slug_cross}.csv"
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["candidate_idx", "tc", "cross_comp_2", "steps", "zeeman", "score", "target_residual", "braid_res", "ybe_res", "best_axis"])
            for r in records:
                w.writerow([r["candidate_idx"], r["tc"], r["cross_comp_2"], r["steps"], r["zeeman"], r["score"], r["target_residual"], r["braid_res"], r["ybe_res"], r["best_axis"]])

        npz_path = outdir / f"crossval_top{idx}_tc{slug_tc}_cross2{slug_cross}.npz"
        np.savez(npz_path, records=np.array(records, dtype=object), summary=summary)

        print(f"Candidate {idx}: tc={tc}, cross_comp_2={cross2} -> mean_score={summary['mean_score']:.6g} ± {summary['std_score']:.3g}")

    # Save aggregated CSV/NPZ
    agg_csv = outdir / f"crossval_top{args.top_n}_aggregate.csv"
    with open(agg_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["candidate_idx", "tc", "cross_comp_2", "steps", "zeeman", "score", "target_residual", "braid_res", "ybe_res", "best_axis"])
        for r in agg_records:
            w.writerow([r["candidate_idx"], r["tc"], r["cross_comp_2"], r["steps"], r["zeeman"], r["score"], r["target_residual"], r["braid_res"], r["ybe_res"], r["best_axis"]])

    agg_npz = outdir / f"crossval_top{args.top_n}_aggregate.npz"
    np.savez(agg_npz, records=np.array(agg_records, dtype=object), summaries=np.array(per_candidate_summaries, dtype=object))

    print(f"Saved per-candidate CSVs and NPZs to {outdir}")
    print(f"Saved aggregate CSV -> {agg_csv}")
    print(f"Saved aggregate NPZ -> {agg_npz}")


if __name__ == '__main__':
    main()
