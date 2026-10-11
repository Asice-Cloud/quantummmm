#!/usr/bin/env python3
"""Compare complete bivector evolution with the published vector contours.

The paper's ten color bands give intervals, not raw wavefunction data. The
literal Hamiltonian independently fixes K=2A (verify_effective_hamiltonian.py);
alternative coefficient scales below are empirical diagnostic comparisons.
All principal comparisons keep tc=E0=.3, the gates, and two cycles unchanged.
"""

import copy
import hashlib
import json
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from analyze_peak_structure import (
    ATOL, MAX_STEP, RTOL, independent_wave, propagate,
)
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.optimize import least_squares
from scipy.signal import find_peaks


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/'paper_alignment'
PAPER = ROOT.parent/'arts/prb111.pdf'
SVG_NS = 'http://www.w3.org/2000/svg'
ET.register_namespace('', SVG_NS)
DPI = 600
SCALE = DPI/72
MODEL_CASES = [
    ('caption_original', .01, 2., 1., 1.),
    ('caption_global_half', .01, 1., 1., 1.),
    ('caption_E_only_half', .01, 2., .5, 1.),
    ('caption_eta_only_half', .01, 2., 1., .5),
    ('caption_weak_half', .01, 2., .5, .5),
]


def extract_contours():
    OUTPUT.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='prb111_alignment_') as tmp:
        svg_path = Path(tmp)/'page2.svg'
        png_stem = Path(tmp)/'page2'
        subprocess.run(['pdftocairo', '-svg', '-f', '2', '-l', '2',
                        str(PAPER), str(svg_path)], check=True)
        subprocess.run(['pdftocairo', '-png', '-singlefile', '-f', '2', '-l', '2',
                        '-r', str(DPI), str(PAPER), str(png_stem)], check=True)
        root = ET.parse(svg_path).getroot()
        paths, colors, coordinates = [], [], []
        excluded = {'rgb(0%, 0%, 0%)', 'rgb(100%, 0%, 0%)', 'rgb(100%, 100%, 100%)'}
        for node in root.iter():
            fill, d = node.get('fill', ''), node.get('d', '')
            if (node.tag.endswith('path') and fill.startswith('rgb(')
                    and fill not in excluded and len(d) > 1000):
                if set(re.findall(r'[A-Za-z]', d))-{'M', 'L', 'Z'}:
                    raise ValueError('Unexpected contour path commands.')
                paths.append(node)
                colors.append([float(x)*2.55 for x in re.findall(r'[\d.]+', fill)])
                coordinates.append(np.array([float(x) for x in re.findall(
                    r'-?\d+(?:\.\d+)?', d)]).reshape(-1, 2))
        if len(paths) != 10:
            raise ValueError(f'Expected 10 contour bands, found {len(paths)}.')
        all_coords = np.concatenate(coordinates)
        x0, y0 = all_coords.min(axis=0)
        x1, y1 = all_coords.max(axis=0)
        bounds = [float(x0), float(x1), float(y0), float(y1)]
        original = Image.open(png_stem.with_suffix('.png')).convert('RGB')
        raster = np.array(original)
        original.crop(tuple(int(round(x*SCALE)) for x in (x0, y0, x1, y1))).save(
            OUTPUT/'paper_heatmap.png')
        # A standalone copy of the source contour geometry, without page text.
        contour_svg = ET.Element(f'{{{SVG_NS}}}svg', {
            'viewBox': f'{x0} {y0} {x1-x0} {y1-y0}',
            'width': '900', 'height': '758',
        })
        for node in paths:
            contour_svg.append(copy.deepcopy(node))
        ET.ElementTree(contour_svg).write(OUTPUT/'paper_contours.svg', encoding='utf-8',
                                         xml_declaration=True)
    return raster, np.array(colors), bounds


def digitize(raster, colors, bounds, tau, logarithm):
    tau, logarithm = np.broadcast_arrays(tau, logarithm)
    x0, x1, y0, y1 = bounds
    x = x0+(x1-x0)*tau/1200
    y = y0+(y1-y0)*(1-logarithm)/2
    rgb = raster[np.rint(y*SCALE).astype(int), np.rint(x*SCALE).astype(int)]
    distance = np.sum((rgb[..., None, :]-colors)**2, axis=-1)
    index = np.argmin(distance, axis=-1)
    valid = np.min(distance, axis=-1) < 100
    lower, upper = index/10, (index+1)/10
    return {'midpoint': (lower+upper)/2, 'lower': lower, 'upper': upper,
            'valid': valid, 'band_index': index}


def score(prediction, paper, mask=None):
    if mask is None:
        mask = paper['valid']
    delta = prediction-paper['midpoint']
    excess = np.maximum(abs(delta)-.05, 0)
    return {'midpoint_RMSE': float(np.sqrt(np.mean(delta[mask]**2))),
            'outside_band_RMSE': float(np.sqrt(np.mean(excess[mask]**2))),
            'inside_band_fraction': float(np.mean(excess[mask] == 0)),
            'valid_point_count': int(np.count_nonzero(mask))}


def model_prediction(case, tau, logarithm, response=False):
    _, E, factor, E_scale, eta_scale = case
    return propagate(tau, E*E_scale, E*10**logarithm*eta_scale,
                     response=response, generator_factor=factor)


def select_bands(paper, tau_axis, after=100):
    """Intervals containing the high, low-ratio peaks, not exact peak roots."""
    high = paper['valid'] & (paper['band_index'] == 9) & (tau_axis >= after)
    indices = np.flatnonzero(high)
    groups = np.split(indices, np.flatnonzero(np.diff(indices)>1)+1)
    return [{'tau_bright_interval': [float(tau_axis[g[0]]), float(tau_axis[g[-1]])],
             'f_color_interval': [.9, 1.]} for g in groups if len(g) >= 3]


def source_grid_check(tau_axis, logs, paper, curves, fit_case):
    """Test the regular mesh suggested by the source polygon vertices.

    The SVG's repeated coordinates suggest 100 duration and 101 logarithm
    nodes. These three cuts fall on logarithm nodes. This is a mesh hypothesis,
    not a reconstruction of the unpublished data or contour algorithm.
    """
    source_tau = np.linspace(0., 1200., 100)
    tau, logarithm = np.meshgrid(source_tau, logs)
    interpolated, metrics = {}, {}
    for case in (*MODEL_CASES[:2], fit_case):
        name = case[0]
        node_values = model_prediction(case, tau, logarithm)['f']
        prediction = np.array([np.interp(tau_axis, source_tau, values)
                               for values in node_values])
        interpolated[name] = prediction
        mask = paper['valid'] & (tau_axis[None, :] >= 100.)
        delta = prediction-curves[name]
        metrics[name] = {
            'continuous': score(curves[name], paper, mask),
            'linear_interpolation': score(prediction, paper, mask),
            'per_cut': [{
                'log10_ratio': float(logarithm_value),
                'interpolation_bias_RMSE': float(np.sqrt(np.mean(delta[row, mask[row]]**2))),
                'interpolation_bias_max': float(np.max(abs(delta[row, mask[row]]))),
            } for row, logarithm_value in enumerate(logs)],
        }
    return {
        'hypothesized_node_count': [101, 100],
        'scope': 'Repeated SVG coordinates suggest this regular mesh; tested as a hypothesis, not identified raw source data.',
        'cut_scores_tau_at_least_100': metrics,
    }, interpolated


def plot_results(tau_axis, logs, digitized, curves, scores, fit_curve):
    fig, axes = plt.subplots(3, 1, figsize=(11., 10.), layout='constrained')
    labels = {
        'caption_original': (r'Current $K=2A$', '#D55E00', '-'),
        'caption_global_half': (r'Diagnostic rescaling $K=A$', '#0072B2', '-'),
    }
    for row, ax in enumerate(axes):
        p = {key: value[row] for key, value in digitized.items()}
        ax.fill_between(tau_axis/100, p['lower'], p['upper'], where=p['valid'],
                        color='black', alpha=.22, label='Paper color interval')
        ax.plot(tau_axis/100, np.where(p['valid'], p['midpoint'], np.nan),
                color='black', linewidth=.7, alpha=.6)
        for name, (label, color, line) in labels.items():
            ax.plot(tau_axis/100, curves[name][row], line, color=color,
                    linewidth=1.2, label=label)
        ax.plot(tau_axis/100, fit_curve[row], '--', color='#009E73', linewidth=1.2,
                label='Diagnostic calibration of weak coefficients')
        ax.set(xlim=(0, 12), ylim=(-.03, 1.03), ylabel=r'$f=\sqrt{F}$',
               title=rf'$\log_{{10}}(|t_1|/E_1)={logs[row]:g}$')
        ax.grid(alpha=.15)
    axes[-1].set_xlabel(r'$\tau\ (100/\mathrm{meV})$')
    axes[0].legend(fontsize=8, ncol=2)
    fig.suptitle(r'Published contours versus full bivector evolution; $E_1=0.01$, $t_c=E_0=0.3$ meV')
    fig.savefig(OUTPUT/'alignment_slices.png', dpi=180)
    fig.savefig(OUTPUT/'alignment_slices.pdf')
    plt.close(fig)

    names = [row['name'] for row in scores]
    fig, ax = plt.subplots(figsize=(10., 4.5), layout='constrained')
    ax.barh(names[::-1], [row['midpoint_RMSE'] for row in scores[::-1]], color='#0072B2')
    ax.set(xlabel='RMSE against the midpoint of the published color interval',
           title='Explicit parameter and matrix-normalization candidates')
    fig.savefig(OUTPUT/'candidate_errors.png', dpi=160)
    plt.close(fig)


def main():
    raster, colors, bounds = extract_contours()
    tau_axis = np.linspace(20., 1180., 25)
    log_axis = np.linspace(-.98, .98, 17)
    tau, logarithm = np.meshgrid(tau_axis, log_axis)
    paper = digitize(raster, colors, bounds, tau, logarithm)
    scores, predictions = [], {}
    for case in MODEL_CASES:
        prediction = model_prediction(case, tau, logarithm)['f']
        predictions[case[0]] = prediction
        result = {'name': case[0], 'label_E1': case[1], 'generator_factor': case[2],
                  'E_coefficient_multiplier': case[3], 'eta_coefficient_multiplier': case[4],
                  **score(prediction, paper)}
        scores.append(result)
        print(case[0], result['midpoint_RMSE'], flush=True)

    # A calibration is an empirical diagnostic, not a Hamiltonian convention.
    # Fit only half the sample locations, with the controlled design unchanged.
    i, j = np.indices(tau.shape)
    usable = paper['valid'] & (tau >= 60.)
    training = usable & ((i+j) % 2 == 0)
    validation = usable & ~training
    def residual(x):
        prediction = propagate(tau[training], x[0], x[1]*10**logarithm[training],
                               response=False, generator_factor=1.)['f']
        return prediction-paper['midpoint'][training]
    fit = least_squares(residual, [.0094, .0105], bounds=([.006, .006], [.014, .014]),
                        diff_step=1e-4, max_nfev=15, ftol=1e-7, xtol=1e-7, gtol=1e-7)
    fit_case = ('diagnostic_fit', .01, 1., float(fit.x[0]/.01), float(fit.x[1]/.01))
    fitted = model_prediction(fit_case, tau, logarithm)['f']
    calibration = {'generator_factor': 1., 'E_matrix_coefficient_meV': float(fit.x[0]),
                   'eta_matrix_coefficient_per_labeled_ratio_meV': float(fit.x[1]),
                   'label_E1_meV': .01,
                   'E_coefficient_multiplier': fit_case[3],
                   'eta_coefficient_multiplier': fit_case[4],
                   'controlled_tc': .3, 'controlled_E0': .3,
                   'training': score(fitted, paper, training),
                   'validation': score(fitted, paper, validation),
                   'optimizer_success': bool(fit.success),
                   'scope': 'Fitted numerical correspondence only; not adopted as the model or proved to be the authors convention.'}
    print('Calibration', calibration, flush=True)

    # Save a denser set of complete-state cut lines and source color intervals.
    curve_tau = np.linspace(0., 1200., 481)
    curve_logs = np.array([-.98, 0., .98])
    cut_tau, cut_log = np.meshgrid(curve_tau, curve_logs)
    digitized = digitize(raster, colors, bounds, cut_tau, cut_log)
    curves = {case[0]: model_prediction(case, cut_tau, cut_log)['f']
              for case in MODEL_CASES[:2]}
    fit_curve = model_prediction(fit_case, cut_tau, cut_log)['f']
    grid_report, grid_curves = source_grid_check(
        curve_tau, curve_logs, digitized, {**curves, 'diagnostic_fit': fit_curve}, fit_case)
    published_bands = select_bands({key: value[0] for key, value in digitized.items()}, curve_tau)
    peak_samples = {}
    for name, values in curves.items():
        indices, _ = find_peaks(values[0], prominence=.3)
        indices = indices[curve_tau[indices] >= 100]
        peak_samples[name] = [{'tau_grid_estimate': float(curve_tau[k]),
                              'f': float(values[0, k])} for k in indices]

    checks = {}
    for case in (MODEL_CASES[0], MODEL_CASES[1], fit_case):
        _, E, factor, e_scale, eta_scale = case
        full = propagate(300., E*e_scale, E*eta_scale, generator_factor=factor)
        wave = independent_wave(300., E*e_scale, E*eta_scale, generator_factor=factor)
        checks[case[0]] = {key: float(np.max(abs(full[key]-wave[key])))
                          for key in ('b', 'F', 'F_tau', 'F_tautau')}
        if max(checks[case[0]].values()) > 2e-8:
            raise AssertionError(checks[case[0]])

    report = {
        'source': {'pdf': str(PAPER.relative_to(ROOT.parent)),
                   'sha256': hashlib.sha256(PAPER.read_bytes()).hexdigest(), 'page': 2,
                   'contour_bounds_pdf_points': bounds, 'color_count': len(colors)},
        'parameter_resolution': {'E1_meV': .01,
                                 'basis': 'Figure caption; user confirmed that the body value .001 is a typo.'},
        'literal_Hamiltonian_audit': {
            'generator_factor': 2.,
            'reference': '../effective_hamiltonian_verification.json',
            'scope': 'Clifford and Fock/BdG constructions fix K=2A for the written Hamiltonian. Other scales are diagnostic, not equivalent conventions with the same coefficients.'},
        'digitization': {'axes': 'tau/(100/meV), log10(abs(t1)/E1)', 'range': [[0, 1200], [-1, 1]],
                         'assumed_equal_color_intervals': 10, 'color_half_width': .05,
                         'scope': 'Color intervals from the vector contour image, not the authors raw data.',
                         'antialias_exclusion_max_RGB_distance': 10,
                         'grid': [len(log_axis), len(tau_axis)]},
        'solver': {'method': 'DOP853', 'rtol': RTOL, 'atol': ATOL, 'max_local_step': MAX_STEP},
        'candidates': sorted(scores, key=lambda row: row['midpoint_RMSE']),
        'diagnostic_calibration': calibration,
        'source_grid_interpolation_check': grid_report,
        'low_ratio': {'log10_ratio': float(curve_logs[0]),
                      'published_high_color_intervals': published_bands,
                      'grid_selected_model_peaks': peak_samples},
        'independent_wave_checks': checks,
        'conclusion': 'The stated-parameter candidates do not fully reproduce the published contour map; calibration is diagnostic and leaves residual differences. The current model convention is unchanged.',
    }
    (OUTPUT/'alignment.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez_compressed(OUTPUT/'alignment.npz', sample_tau=tau_axis, sample_log=log_axis,
                        sample_paper=paper['midpoint'], sample_valid=paper['valid'],
                        **predictions, fitted=fitted, colors=colors,
                        curve_tau=curve_tau, curve_log=curve_logs,
                        curve_paper_lower=digitized['lower'], curve_paper_upper=digitized['upper'],
                        curve_paper_valid=digitized['valid'],
                        curve_original=curves['caption_original'],
                        curve_global_half=curves['caption_global_half'], curve_fit=fit_curve,
                        grid_curve_original=grid_curves['caption_original'],
                        grid_curve_global_half=grid_curves['caption_global_half'],
                        grid_curve_fit=grid_curves['diagnostic_fit'])
    plot_results(curve_tau, curve_logs, digitized, curves, scores, fit_curve)
    print('Saved source contours, candidate scores, independent checks, and comparison plots.', flush=True)


if __name__ == '__main__':
    main()
