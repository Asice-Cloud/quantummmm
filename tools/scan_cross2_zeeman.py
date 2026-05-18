#!/usr/bin/env python3
"""Scan cross_comp_2 × zeeman 2D grid and save diagnostics.

Default: cross_comp_2 ∈ [-0.02,0.02], zeeman ∈ [-0.02,0.02], 25×25, steps=300.
"""
from __future__ import annotations

import argparse
from importlib import util
from pathlib import Path
import csv
import numpy as np

HERE = Path(__file__).resolve().parent
spec = util.spec_from_file_location("verify_path_braid", str(HERE / "verify_path_braid.py"))
vpb = util.module_from_spec(spec)
spec.loader.exec_module(vpb)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--steps", type=int, default=300)
    p.add_argument("--T", type=float, default=3.175)
    p.add_argument("--tc", type=float, default=1.0)
    p.add_argument("--hc-factor", type=float, default=2.125)
    p.add_argument("--cross-min", type=float, default=-0.02)
    p.add_argument("--cross-max", type=float, default=0.02)
    p.add_argument("--cross-num", type=int, default=25)
    p.add_argument("--zeeman-min", type=float, default=-0.02)
    p.add_argument("--zeeman-max", type=float, default=0.02)
    p.add_argument("--zeeman-num", type=int, default=25)
    p.add_argument("--axis", choices=("x", "y", "z"), default="z")
    p.add_argument("--outdir", type=str, default="results/path_braid")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cross_vals = np.linspace(args.cross_min, args.cross_max, args.cross_num)
    zeeman_vals = np.linspace(args.zeeman_min, args.zeeman_max, args.zeeman_num)

    records = []
    for cross in cross_vals:
        for zeeman in zeeman_vals:
            gates = vpb.compose_cycle(
                steps=int(args.steps),
                T=float(args.T),
                tc=float(args.tc),
                hc_factor=float(args.hc_factor),
                zeeman=float(zeeman),
                zeeman_drive_1=0.0,
                zeeman_drive_2=0.0,
                zeeman_drive_3=0.0,
                zeeman_drive_shape="cos",
                phase_comp=0.0,
                cross_comp_1=0.0,
                cross_comp_2=float(cross),
                cross_comp_3=0.0,
            )
            diag = vpb.analyze_total_gate(gates["R_total"], args.axis)
            records.append(
                {
                    "cross_comp_2": float(cross),
                    "zeeman": float(zeeman),
                    "score": float(diag["score"]),
                    "target_residual": float(diag["target_residual"]),
                    "braid_res": float(diag["braid_res"]),
                    "ybe_res": float(diag["ybe_res"]),
                    "best_axis": str(diag["best_axis"]),
                }
            )

    best = min(records, key=lambda r: r["score"])

    csv_path = outdir / f"scan_cross2_zeeman_steps{args.steps}_T{args.T}.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cross_comp_2", "zeeman", "score", "target_residual", "braid_res", "ybe_res", "best_axis"])
        for r in records:
            writer.writerow([r["cross_comp_2"], r["zeeman"], r["score"], r["target_residual"], r["braid_res"], r["ybe_res"], r["best_axis"]])

    npz_path = outdir / f"scan_cross2_zeeman_steps{args.steps}_T{args.T}.npz"
    np.savez(npz_path, records=np.array(records, dtype=object), best=best)

    print(f"scan complete: {len(records)} points")
    print(f"best cross_comp_2 = {best['cross_comp_2']:.6g}")
    print(f"best zeeman = {best['zeeman']:.6g}")
    print(f"best score = {best['score']:.6g}")
    print(f"target_res = {best['target_residual']:.6g}")
    print(f"braid_res = {best['braid_res']:.6g}")
    print(f"ybe_res = {best['ybe_res']:.6g}")
    print(f"saved CSV -> {csv_path}")
    print(f"saved NPZ -> {npz_path}")


if __name__ == '__main__':
    main()
