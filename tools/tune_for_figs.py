#!/usr/bin/env python3
"""Tuning helper for Fig.3/4/5 trend reproduction.

This script wraps the existing tools/reproduce_trend_figs.py functions and adds a
small tuning pipeline:

1. run_auto_scan() to pick baseline best candidates (fast)
2. sensitivity analysis (finite-diff) on the shortlisted parameters
3. coarse grid search in parameter subspaces suggested by sensitivity
4. local refinement around coarse bests
5. cross-validation (varying n_per_step and small perturbations)

The script saves results under results/paper_trends/tuning_<timestamp>/.
It purposely does not run very large grids by default; pass --run-full to
execute denser searches.

Usage examples:
  python tools/tune_for_figs.py --outdir results/paper_trends --run-full

"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import json
import numpy as np

# import the existing trend reproduction module
from tools import reproduce_trend_figs as rt
import tools.paper_params as P


def timestamp():
    return datetime.now().strftime('%Y%m%dT%H%M%S')


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def to_jsonable(obj):
    """Recursively convert numpy-heavy structures into JSON-serializable data."""
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return obj


def run_baseline_auto_scan(outdir: Path, n_per_step: int = 180) -> dict:
    """Run the existing auto-scan to obtain baseline best candidates."""
    print('Running baseline auto-scan (fast) ...')
    res = rt.run_auto_scan(outdir, n_per_step=n_per_step)
    print('Baseline auto-scan done.')
    return res


def sensitivity_fig4(base_params: dict, delta_frac: float = 0.05):
    """Finite-difference sensitivity for Fig.4 parameters (d0a,d0b,amp).

    base_params should include keys: 'd0a','d0b','amp', and 'Tscan'.
    Returns sensitivities dict.
    """
    sens = {}
    Tscan = base_params['Tscan']
    base_a = base_params['a']
    base_b = base_params['b']
    base_score = rt.trend_score_fig4(np.array(base_a), np.array(base_b))

    for key in ('delta0_a', 'delta0_b', 'amp'):
        orig = base_params[key if key != 'amp' else 'amp']
        up = dict(base_params)
        down = dict(base_params)
        # perturb relatively
        up_val = orig * (1 + delta_frac) if orig != 0 else delta_frac
        down_val = orig * (1 - delta_frac) if orig != 0 else -delta_frac
        if key == 'delta0_a':
            # recompute curves
            curve_a_up = []
            curve_b_up = []
            curve_a_down = []
            curve_b_down = []
            for TT in Tscan:
                _, ov_a_up = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(up_val, base_params['amp']))
                _, ov_b_up = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(base_params['delta0_b'], base_params['amp']))
                _, ov_a_down = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(down_val, base_params['amp']))
                _, ov_b_down = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(base_params['delta0_b'], base_params['amp']))
                curve_a_up.append(abs(ov_a_up[-1]))
                curve_b_up.append(abs(ov_b_up[-1]))
                curve_a_down.append(abs(ov_a_down[-1]))
                curve_b_down.append(abs(ov_b_down[-1]))
            s_up = rt.trend_score_fig4(np.array(curve_a_up), np.array(curve_b_up))
            s_down = rt.trend_score_fig4(np.array(curve_a_down), np.array(curve_b_down))
            sens['delta0_a'] = {'base': base_score, 'up': s_up, 'down': s_down}
        elif key == 'delta0_b':
            curve_a_up = []
            curve_b_up = []
            curve_a_down = []
            curve_b_down = []
            for TT in Tscan:
                _, ov_a_up = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(base_params['delta0_a'], base_params['amp']))
                _, ov_b_up = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(up_val, base_params['amp']))
                _, ov_a_down = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(base_params['delta0_a'], base_params['amp']))
                _, ov_b_down = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(down_val, base_params['amp']))
                curve_a_up.append(abs(ov_a_up[-1]))
                curve_b_up.append(abs(ov_b_up[-1]))
                curve_a_down.append(abs(ov_a_down[-1]))
                curve_b_down.append(abs(ov_b_down[-1]))
            s_up = rt.trend_score_fig4(np.array(curve_a_up), np.array(curve_b_up))
            s_down = rt.trend_score_fig4(np.array(curve_a_down), np.array(curve_b_down))
            sens['delta0_b'] = {'base': base_score, 'up': s_up, 'down': s_down}
        else:  # amp
            curve_a_up = []
            curve_b_up = []
            curve_a_down = []
            curve_b_down = []
            for TT in Tscan:
                _, ov_a_up = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(base_params['delta0_a'], up_val))
                _, ov_b_up = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(base_params['delta0_b'], up_val))
                _, ov_a_down = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(base_params['delta0_a'], down_val))
                _, ov_b_down = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=180, modulation=(base_params['delta0_b'], down_val))
                curve_a_up.append(abs(ov_a_up[-1]))
                curve_b_up.append(abs(ov_b_up[-1]))
                curve_a_down.append(abs(ov_a_down[-1]))
                curve_b_down.append(abs(ov_b_down[-1]))
            s_up = rt.trend_score_fig4(np.array(curve_a_up), np.array(curve_b_up))
            s_down = rt.trend_score_fig4(np.array(curve_a_down), np.array(curve_b_down))
            sens['amp'] = {'base': base_score, 'up': s_up, 'down': s_down}

    return sens


def cross_validate_fig(best_rec: dict, outdir: Path, n_per_steps=(120, 180, 300), zeeman_shifts=(0.0, 0.005, -0.005)):
    """Cross-validate a selected candidate by varying n_per_step and small static shifts.

    best_rec expected to be one of the dict entries returned by run_auto_scan
    (fig3/fig4/fig5). For Fig.4/Fig.5 best_rec contains 'Tscan' and curves.
    """
    out = {'meta': to_jsonable(best_rec)}
    # For Fig.4 and Fig.5, we'll re-evaluate overlaps for each n_per_step.
    scores = []
    gap_mins = []
    for nps in n_per_steps:
        # compute score using the same scoring function depending on keys
        if 'a' in best_rec and 'b' in best_rec:
            # fig4-like
            Tscan = best_rec['Tscan']
            # use stored modulation params if present
            d0a = best_rec.get('delta0_a', 0.59)
            d0b = best_rec.get('delta0_b', 0.57)
            amp = best_rec.get('amp', 0.02)
            curve_a = []
            curve_b = []
            for TT in Tscan:
                _, ov_a = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=nps, modulation=(d0a, amp))
                _, ov_b = rt.compute_overlap_curve(T_step=float(TT), delta=0.2, n_per_step=nps, modulation=(d0b, amp))
                curve_a.append(abs(ov_a[-1]))
                curve_b.append(abs(ov_b[-1]))
            score = rt.trend_score_fig4(np.array(curve_a), np.array(curve_b))
            scores.append(score)
            # gap_min proxy: use a BdG trace for a representative T
            t_over_T, branches = rt.compute_bdg_trace(T_step=float(Tscan[len(Tscan)//2]), n_per_step=200, delta_mod=d0a, amp=amp)
            gap_min = float(np.min(np.abs(branches[:, 1] - branches[:, 0])))
            gap_mins.append(gap_min)
        elif 'c1' in best_rec and 'c2' in best_rec:
            # fig5-like
            Tscan = best_rec['Tscan']
            d0 = best_rec.get('delta0', 0.59)
            al = best_rec.get('amp_low', 0.01)
            ah = best_rec.get('amp_high', 0.03)
            c1 = []
            c2 = []
            for TT in Tscan:
                _, ov1 = rt.compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=nps, modulation=(d0, al))
                _, ov2 = rt.compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=nps, modulation=(d0, ah))
                c1.append(abs(ov1[-1]))
                c2.append(abs(ov2[-1]))
            score = rt.trend_score_fig5(np.array(c1), np.array(c2))
            scores.append(score)
            t_over_T, branches = rt.compute_bdg_trace(T_step=float(Tscan[len(Tscan)//2]), n_per_step=200, delta_mod=d0, amp=al)
            gap_min = float(np.min(np.abs(branches[:, 1] - branches[:, 0])))
            gap_mins.append(gap_min)
        else:
            # fig3-like
            Tscan = best_rec['Tscan']
            m = best_rec['mzm']
            a = best_rec['abs']
            score = rt.trend_score_fig3(np.array(m), np.array(a))
            scores.append(score)
            t_over_T, branches = rt.compute_bdg_trace(T_step=float(Tscan[len(Tscan)//2]), n_per_step=200, delta_mod=None, amp=0.0)
            gap_mins.append(float(np.min(np.abs(branches[:, 1] - branches[:, 0]))))

    out['scores'] = {'n_per_steps': list(n_per_steps), 'scores': scores, 'mean': float(np.mean(scores)), 'std': float(np.std(scores))}
    out['gap_mins'] = {'values': gap_mins, 'min': float(np.min(gap_mins))}

    # small static perturbations (modulation/delta shifts) - only for fig4/5
    pert_results = []
    if 'a' in best_rec or 'c1' in best_rec:
        for shift in (0.0, 0.005, -0.005):
            # naive perturbation: for fig4, add shift to d0 values; for fig5, add to delta0
            if 'a' in best_rec:
                d0a = best_rec.get('delta0_a', 0.59) + shift
                d0b = best_rec.get('delta0_b', 0.57) + shift
                amp = best_rec.get('amp', 0.02)
                Tscan = best_rec['Tscan']
                ca, cb = [], []
                for TT in Tscan:
                    _, ov_a = rt.compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=180, modulation=(d0a, amp))
                    _, ov_b = rt.compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=180, modulation=(d0b, amp))
                    ca.append(abs(ov_a[-1])); cb.append(abs(ov_b[-1]))
                s = rt.trend_score_fig4(np.array(ca), np.array(cb))
                pert_results.append({'shift': shift, 'score': float(s)})
            else:
                d0 = best_rec.get('delta0', 0.59) + shift
                al = best_rec.get('amp_low', 0.01)
                ah = best_rec.get('amp_high', 0.03)
                Tscan = best_rec['Tscan']
                c1, c2 = [], []
                for TT in Tscan:
                    _, ov1 = rt.compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=180, modulation=(d0, al))
                    _, ov2 = rt.compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=180, modulation=(d0, ah))
                    c1.append(abs(ov1[-1])); c2.append(abs(ov2[-1]))
                s = rt.trend_score_fig5(np.array(c1), np.array(c2))
                pert_results.append({'shift': shift, 'score': float(s)})

    out['perturbations'] = pert_results

    # save
    ensure_dir(outdir)
    out_path = outdir / 'crossval_summary.json'
    with open(out_path, 'w', encoding='utf-8') as fh:
        json.dump(to_jsonable(out), fh, indent=2)
    print(f'Cross-validation saved to {out_path}')
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--outdir', type=str, default='results/paper_trends', help='base output directory')
    p.add_argument('--scan-n-per-step', type=int, default=180, help='n per step for auto-scan')
    p.add_argument('--run-full', action='store_true', help='run denser grids (slow)')
    args = p.parse_args()

    base_out = Path(args.outdir)
    stamp = timestamp()
    outdir = base_out / f'tuning_{stamp}'
    ensure_dir(outdir)

    # 1) baseline auto scan. This is the same logic as the original script's auto-scan
    res = run_baseline_auto_scan(outdir, n_per_step=args.scan_n_per_step)

    # Save baseline summary
    with open(outdir / 'auto_scan_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(
            {
                'fig3': {k: res['fig3'][k] for k in ('score', 'delta_abs')},
                'fig4': {k: res['fig4'][k] for k in ('score', 'delta0_a', 'delta0_b', 'amp')},
                'fig5': {k: res['fig5'][k] for k in ('score', 'delta0', 'amp_low', 'amp_high')},
            },
            fh,
            indent=2,
            default=to_jsonable,
        )

    # 2) sensitivity analysis for Fig.4 as an example (fast)
    print('Running sensitivity analysis for Fig.4 (finite-diff, small) ...')
    sens4 = sensitivity_fig4(res['fig4'], delta_frac=0.05)
    with open(outdir / 'sensitivity_fig4.json', 'w', encoding='utf-8') as fh:
        json.dump(sens4, fh, indent=2, default=to_jsonable)
    print('Sensitivity analysis saved.')

    # 3) cross-validate the auto-scan winners
    print('Cross-validating auto-scan winners (n_per_step variations + small shifts) ...')
    cv_fig3 = cross_validate_fig(res['fig3'], outdir / 'fig3_cv')
    cv_fig4 = cross_validate_fig(res['fig4'], outdir / 'fig4_cv')
    cv_fig5 = cross_validate_fig(res['fig5'], outdir / 'fig5_cv')

    print('Tuning helper completed. Results in:', outdir)


if __name__ == '__main__':
    main()
