#!/usr/bin/env python3
"""Simple standalone script to map toy `delta` to BdG `VD` scale factors.

This script is intentionally minimal: it builds toy 2x2 instantaneous Eg(t)
and BdG lowest-|E|(t) for different VD scale factors, then finds for each
`delta` the VD factor that best matches (by time-averaged Eg and by a
reference time). It saves a PNG and NPZ with results.

Usage (example):
  .venv/bin/python tools/scan_delta_vd_simple.py --outdir results/scan_simple

Defaults are moderate resolution to keep runtime reasonable.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import time
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.interpolate import interp1d

# Make repo root importable
THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import tetron_path_sim as tetron
from tools import embed_kitaev
import tools.paper_params as P


def map_gates_to_links(g1, g2, g3, g4, t0, Delta0, L, mu0, VD, qd_width):
    mu = mu0 * np.ones(L)
    for i in range(qd_width):
        mu[i] = mu[i] - VD

    t_links = t0 * np.ones(L - 1)
    t_links_mod = t_links.copy()
    t_links_mod[0] = t0 * g1
    if L > 2:
        t_links_mod[1] = t0 * g3

    Delta_links = Delta0 * np.ones(L - 1)
    Delta_mod = Delta_links.copy()
    Delta_mod[0] = Delta0 * (g1 if g1 > 0 else 1e-6)
    if L > 2:
        Delta_mod[1] = Delta0 * (g3 if g3 > 0 else 1e-6)

    return mu, t_links_mod, Delta_mod


def compute_bdg_Eg_series(vd_factor: float, T_step: float, n_per_step: int, qd_width: int | None = None):
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)
    Eg_ts = np.zeros(N)
    qd_w = qd_width if qd_width is not None else P.QD_WIDTH

    base_steps = int(max(step_idx))
    for k in range(N):
        step = int(step_idx[k])
        s = float(slist[k])
        step_base = ((step - 1) % base_steps) + 1
        g1, g2, g3, g4 = tetron.gates_at(step_base, s)

        VD_here = float(P.VD) * float(vd_factor)
        mu, t_links_mod, Delta_mod = map_gates_to_links(
            g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, VD_here, qd_w
        )
        Hfull = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
        Efull = np.linalg.eigvalsh(Hfull)
        Eg_ts[k] = float(np.min(np.abs(Efull)))

    return tlist / T_step, Eg_ts


def compute_toy_Eg_series(delta: float, T_step: float, n_per_step: int):
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)
    Eg_ts = np.zeros(N)

    base_steps = int(max(step_idx))
    for k in range(N):
        step = int(step_idx[k])
        s = float(slist[k])
        theta = tetron.theta_from_time(step, s)
        Ht = tetron.H_eff_from_theta(theta, delta=float(delta))
        vals = eigh(Ht)[0]
        Eg_ts[k] = float(np.min(np.abs(vals)))

    return tlist / T_step, Eg_ts


def fit_line(x: np.ndarray, y: np.ndarray):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return 0.0, 0.0, float('nan'), float('nan')
    a, b = np.polyfit(x[mask], y[mask], 1)
    ypred = a * x + b
    rmse = float(np.sqrt(np.mean((y[mask] - ypred[mask]) ** 2)))
    try:
        corr = float(np.corrcoef(x[mask], y[mask])[0, 1])
    except Exception:
        corr = float('nan')
    return float(a), float(b), rmse, corr


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser(description='Standalone δ ↔ VD mapping scan (simple).')
    p.add_argument('--outdir', type=str, default='results/scan_simple')
    p.add_argument('--deltas', type=float, nargs='+', default=None, help='list of toy δ values')
    p.add_argument('--vd-min', type=float, default=1.0)
    p.add_argument('--vd-max', type=float, default=8.0)
    p.add_argument('--vd-n', type=int, default=41)
    p.add_argument('--T-step', type=float, default=200.0)
    p.add_argument('--n-per-step', type=int, default=80)
    p.add_argument('--ref-time', type=float, default=0.5)
    p.add_argument('--qd-width', type=int, default=None)
    args = p.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.deltas is None:
        deltas = np.linspace(0.0, 0.2, 11)
    else:
        deltas = np.array(args.deltas, dtype=float)

    vd_factors = np.linspace(float(args.vd_min), float(args.vd_max), int(args.vd_n))

    print('Scan settings:')
    print(' deltas=', deltas)
    print(' vd_factors (len)=', len(vd_factors), f'[{vd_factors[0]}..{vd_factors[-1]}]')
    print(' T_step=', args.T_step, ' n_per_step=', args.n_per_step)

    # build time grid to find reference index
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=args.T_step, n_per_step=args.n_per_step)
    t_over_T = tlist / args.T_step
    idx_ref = int(np.argmin(np.abs(t_over_T - float(args.ref_time))))

    # compute BdG series for each vd_factor (avg and ref)
    Eg_bdg_avg = np.zeros(len(vd_factors))
    Eg_bdg_ref = np.zeros(len(vd_factors))

    t0 = time.time()
    for j, vf in enumerate(vd_factors):
        t1 = time.time()
        _, Eg_ts = compute_bdg_Eg_series(vf, T_step=args.T_step, n_per_step=args.n_per_step, qd_width=args.qd_width)
        Eg_bdg_avg[j] = float(np.mean(Eg_ts))
        Eg_bdg_ref[j] = float(Eg_ts[idx_ref])
        if (j + 1) % 5 == 0 or j == len(vd_factors) - 1:
            print(f'  vd [{j+1}/{len(vd_factors)}] vf={vf:.3f}  avgEg={Eg_bdg_avg[j]:.4e}  refEg={Eg_bdg_ref[j]:.4e}  (time {time.time()-t1:.1f}s)')
    print('BdG scan done, elapsed', time.time() - t0)

    # compute toy Eg for each delta
    Eg_toy_avg = np.zeros(len(deltas))
    Eg_toy_ref = np.zeros(len(deltas))
    for i, dd in enumerate(deltas):
        _, Eg_ts = compute_toy_Eg_series(dd, T_step=args.T_step, n_per_step=args.n_per_step)
        Eg_toy_avg[i] = float(np.mean(Eg_ts))
        Eg_toy_ref[i] = float(Eg_ts[idx_ref])

    # attempt monotonic inversion (Eg -> vd) using interpolation if possible
    best_vd_avg = np.zeros(len(deltas))
    best_vd_ref = np.zeros(len(deltas))

    # check monotonicity
    monot_inc = np.all(np.diff(Eg_bdg_avg) >= 0)
    monot_dec = np.all(np.diff(Eg_bdg_avg) <= 0)
    if monot_inc or monot_dec:
        try:
            f_avg = interp1d(Eg_bdg_avg, vd_factors, kind='linear', bounds_error=False, fill_value=(vd_factors[0], vd_factors[-1]))
            best_vd_avg = f_avg(Eg_toy_avg)
        except Exception:
            best_vd_avg = vd_factors[np.argmin(np.abs(Eg_bdg_avg[None, :] - Eg_toy_avg[:, None]), axis=1)]
    else:
        best_vd_avg = vd_factors[np.argmin(np.abs(Eg_bdg_avg[None, :] - Eg_toy_avg[:, None]), axis=1)]

    # same for reference-time
    monot_inc_r = np.all(np.diff(Eg_bdg_ref) >= 0)
    monot_dec_r = np.all(np.diff(Eg_bdg_ref) <= 0)
    if monot_inc_r or monot_dec_r:
        try:
            f_ref = interp1d(Eg_bdg_ref, vd_factors, kind='linear', bounds_error=False, fill_value=(vd_factors[0], vd_factors[-1]))
            best_vd_ref = f_ref(Eg_toy_ref)
        except Exception:
            best_vd_ref = vd_factors[np.argmin(np.abs(Eg_bdg_ref[None, :] - Eg_toy_ref[:, None]), axis=1)]
    else:
        best_vd_ref = vd_factors[np.argmin(np.abs(Eg_bdg_ref[None, :] - Eg_toy_ref[:, None]), axis=1)]

    # fits
    a_avg, b_avg, rmse_avg, corr_avg = fit_line(deltas, best_vd_avg)
    a_ref, b_ref, rmse_ref, corr_ref = fit_line(deltas, best_vd_ref)

    # save results
    outpng = outdir / 'delta_vs_vd_map_simple.png'
    outnpz = outdir / 'delta_vs_vd_map_simple.npz'

    fig, ax = plt.subplots(1, 1, figsize=(6.6, 4.2), dpi=160)
    ax.scatter(deltas, best_vd_avg, color='#1f77b4', label='best vd (avg E)')
    ax.scatter(deltas, best_vd_ref, color='#d62728', label='best vd (ref time)')
    xs = np.linspace(float(deltas.min()), float(deltas.max()), 200)
    ax.plot(xs, a_avg * xs + b_avg, '--', color='#1f77b4', alpha=0.7, label=f'fit avg: vd={a_avg:.3f}·δ+{b_avg:.3f}')
    ax.plot(xs, a_ref * xs + b_ref, '--', color='#d62728', alpha=0.7, label=f'fit ref: vd={a_ref:.3f}·δ+{b_ref:.3f}')
    ax.set_xlabel('δ (toy)')
    ax.set_ylabel('VD scale factor')
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpng, bbox_inches='tight')
    plt.close(fig)

    np.savez(
        outnpz,
        deltas=deltas,
        vd_factors=vd_factors,
        Eg_bdg_avg=Eg_bdg_avg,
        Eg_bdg_ref=Eg_bdg_ref,
        Eg_toy_avg=Eg_toy_avg,
        Eg_toy_ref=Eg_toy_ref,
        best_vd_avg=best_vd_avg,
        best_vd_ref=best_vd_ref,
        fit_avg=(a_avg, b_avg, rmse_avg, corr_avg),
        fit_ref=(a_ref, b_ref, rmse_ref, corr_ref),
    )

    print('\nScan complete. Results:')
    print('  outpng=', outpng)
    print('  outnpz=', outnpz)
    print('  fit avg: a={:.6e}, b={:.6e}, rmse={:.6e}, corr={}'.format(a_avg, b_avg, rmse_avg, corr_avg))
    print('  fit ref: a={:.6e}, b={:.6e}, rmse={:.6e}, corr={}'.format(a_ref, b_ref, rmse_ref, corr_ref))


if __name__ == '__main__':
    main()
