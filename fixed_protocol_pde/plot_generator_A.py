#!/usr/bin/env python3
"""Compare the K=A normalization trial with the paper and original K=2A.

Original amplitudes, cosine gates, solver settings and duration axes are kept.
At fixed normalized time, b_A(tau)=b_2A(tau/2). Reuse only exact source nodes;
all missing nodes are integrated with the complete ten-component equation.
"""

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

from analyze_peak_structure import ATOL, MAX_STEP, RTOL, independent_wave, propagate, readout
from plot_state_evolution import invariants, map_axes, palette
from plot_normalization_comparison import draw_model, save
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / 'normalization_comparison'
SOURCE = ROOT / 'state_evolution' / 'endpoint_map.npz'
DATA = OUTPUT / 'endpoint_map_KA.npz'


def compute_map(base):
    tau_axis = base['tau'].copy()
    log_axis = base['log_ratio'].copy()
    even = np.arange(0, len(tau_axis), 2)
    assert np.array_equal(tau_axis[even]/2, base['tau'][:len(even)])
    b = np.empty((len(log_axis), len(tau_axis), 10))
    b[:, even] = base['b'][:, :len(even)]
    missing = np.arange(1, len(tau_axis), 2)
    started = perf_counter()
    print(f'Reused {len(even)} exact duration nodes using f_A(tau)=f_2A(tau/2).', flush=True)
    for start in range(0, len(missing), 12):
        selected = missing[start:start+12]
        tau, logs = np.meshgrid(tau_axis[selected], log_axis)
        b[:, selected] = propagate(tau, .01, .01*10**logs, tc=.3, e0=.3,
                                   response=False, generator_factor=1.)['b']
        print(f'K=A map: {min(start+12, len(missing))}/{len(missing)} new columns; '
              f'elapsed {perf_counter()-started:.1f}s', flush=True)
    data = {'tau': tau_axis, 'log_ratio': log_axis, 'b': b, **readout(b)}
    np.savez_compressed(DATA, **data)
    return data


def verify(data):
    checks = invariants(data['b'])
    assert checks['max_abs_B_cubed_plus_B'] < 2e-8, checks
    assert checks['max_abs_component_norm_squared_minus_one'] < 2e-8, checks
    assert -2e-8 <= checks['F_min'] <= checks['F_max'] <= 1+2e-8, checks
    samples = []
    # Cover both freshly integrated nodes and an exactly reused node.
    for tau, log_ratio in [(105., -1.), (300., 0.), (995., 1.)]:
        row = int(np.argmin(abs(data['log_ratio']-log_ratio)))
        col = int(np.argmin(abs(data['tau']-tau)))
        assert data['tau'][col] == tau
        ratio = float(10**data['log_ratio'][row])
        wave = independent_wave(tau, .01, .01*ratio, tc=.3, e0=.3,
                                generator_factor=1.)
        sample = {
            'tau': tau, 'eta_over_E1': ratio,
            'map_node_source': 'integrated' if col % 2 else 'exact_scaling',
            'B_max_error': float(np.max(abs(wave['b']-data['b'][row, col]))),
            'F_error': float(abs(wave['F']-data['F'][row, col])),
            'sqrt_F_error': float(abs(np.sqrt(wave['F'])-data['f'][row, col])),
        }
        assert sample['B_max_error'] < 2e-8 and sample['F_error'] < 2e-8, sample
        samples.append(sample)
    return {'invariants': checks, 'independent_wave_checks': samples,
            'reuse_identity': 'b_A(zeta;tau)=b_2A(zeta;tau/2), at fixed E,eta,tc,E0'}


def plot(data, base):
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12,
                         'mathtext.fontset': 'stix', 'axes.linewidth': .8,
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'pdf.fonttype': 42})
    cmap = palette()
    parameters = r'$E_1=0.01,\quad t_c=E_0=0.3\ \mathrm{meV}$'
    fig, ax = plt.subplots(figsize=(7.2, 6.4), layout='constrained')
    mesh = draw_model(ax, data, cmap)
    ax.set_title(parameters+'\n'+r'$|t_1|\in[0.001,0.1]\ \mathrm{meV}$', fontsize=12)
    fig.colorbar(mesh, ax=ax, pad=.025, ticks=np.linspace(0, 1, 6),
                 label=r'$\sqrt{F(6\tau)}$')
    fig.suptitle(r'$\sqrt{F}:\quad K=A$', fontsize=18)
    fig.supxlabel('Original amplitudes and six cosine stages; generator prefactor = 1.', fontsize=10)
    save(fig, 'sqrtF_KA_original_parameters')

    fig, axes = plt.subplots(1, 3, figsize=(18.2, 6.0), layout='constrained')
    paper = np.asarray(Image.open(ROOT / 'paper_alignment' / 'paper_heatmap.png').convert('RGB'))
    axes[0].imshow(paper, origin='upper', extent=(0, 12, -1, 1),
                   aspect='auto', interpolation='nearest')
    map_axes(axes[0])
    axes[0].set_box_aspect(655/778)
    axes[0].set_title('Paper Fig. 1(d)')
    draw_model(axes[1], base, cmap)
    axes[1].set_title(r'Previous normalization: $K=2A$')
    mesh = draw_model(axes[2], data, cmap)
    axes[2].set_title(r'Current trial: $K=A$')
    for ax in axes[1:]:
        ax.set_ylabel('')
    fig.colorbar(mesh, ax=axes.tolist(), pad=.014, shrink=.88,
                 ticks=np.linspace(0, 1, 6), label=r'$f=\sqrt{F}$')
    fig.suptitle(r'$f=|\langle\psi_1^+(0)|\psi_1^-(6\tau)\rangle|$'+'; '+parameters, fontsize=16)
    fig.supxlabel(r'Exact scaling: $f_A(\tau,\eta/E_1)=f_{2A}(\tau/2,\eta/E_1)$.'+
                  ' Same amplitudes, duration axes, gates and color scale.', fontsize=11)
    save(fig, 'paper_vs_K2A_vs_KA_sqrtF')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replot', action='store_true', help='Plot and verify saved K=A data.')
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    with np.load(SOURCE) as stored:
        base = {key: stored[key] for key in stored.files}
    metadata = json.loads((SOURCE.parent / 'verification.json').read_text())
    assert metadata['generator'] == 'K=2A'
    assert metadata['E1_map_meV'] == .01
    assert metadata['tc_meV'] == metadata['E0_meV'] == .3
    if args.replot:
        with np.load(DATA) as stored:
            data = {key: stored[key] for key in stored.files}
    else:
        data = compute_map(base)
    checks = verify(data)
    report = {
        'generator': 'K=A', 'generator_factor': 1.,
        'matrix_parameters_meV': {'E1': .01, 'eta_range': [.001, .1], 'tc': .3, 'E0': .3},
        'tau_range': [0., 1200.], 'end_time': '6*tau', 'grid': list(data['F'].shape),
        'output_tau_spacing': 5.,
        'solver': {'method': 'DOP853', 'rtol': RTOL, 'atol': ATOL, 'max_local_step': MAX_STEP},
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'scope': 'Normalization trial with original amplitudes; does not change the K=2A derivation from literal H.',
        'checks': checks,
    }
    (OUTPUT / 'verification_KA.json').write_text(json.dumps(report, indent=2)+'\n')
    plot(data, base)
    print(json.dumps(checks, indent=2), flush=True)
    print('Saved K=A maps and three-panel comparison.', flush=True)


if __name__ == '__main__':
    main()
