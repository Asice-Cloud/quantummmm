#!/usr/bin/env python3
"""Plot sqrt(F) with K=2A written as A evaluated at doubled amplitudes.

The factor is absorbed once: A(2E, 2eta, 2tc, 2E0) = 2A(E, eta, tc, E0).
The existing full-state map therefore applies without reintegration. Independent
wavefunction checks use the doubled matrix amplitudes and generator_factor=1.
"""

import hashlib
import json
from pathlib import Path

from analyze_peak_structure import ATOL, MAX_STEP, RTOL, independent_wave
from plot_state_evolution import map_axes, palette
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / 'normalization_comparison'
SOURCE = ROOT / 'state_evolution' / 'endpoint_map.npz'
LEVELS = np.linspace(0., 1., 11)
AMPLITUDES = r"$E_1'=0.02,\quad t_c'=E_0'=0.6\ \mathrm{meV}$"
ETA_RANGE = r"$|t_1|'\in[0.002,0.2]\ \mathrm{meV},\quad |t_1|'/E_1'=|t_1|/E_1$"


def check_equivalence(data, metadata):
    assert metadata['generator'] == 'K=2A'
    assert metadata['E1_map_meV'] == .01
    assert metadata['tc_meV'] == metadata['E0_meV'] == .3
    assert np.allclose(data['f'], np.sqrt(np.maximum(data['F'], 0)), atol=1e-14)
    samples = []
    for row, col in [(0, 20), (50, 60), (100, 200)]:
        tau = float(data['tau'][col])
        ratio = float(10**data['log_ratio'][row])
        wave = independent_wave(tau, .02, .02*ratio, tc=.6, e0=.6,
                                generator_factor=1.)
        sample = {
            'tau': tau,
            'eta_over_E1': ratio,
            'doubled_matrix_E1': .02,
            'doubled_matrix_eta': .02*ratio,
            'F_original_K2A': float(data['F'][row, col]),
            'F_doubled_amplitudes_factor_one': float(wave['F']),
            'B_max_error': float(np.max(abs(wave['b']-data['b'][row, col]))),
            'F_error': float(abs(wave['F']-data['F'][row, col])),
            'sqrt_F_error': float(abs(np.sqrt(wave['F'])-data['f'][row, col])),
        }
        assert sample['B_max_error'] < 2e-8, sample
        assert sample['F_error'] < 2e-8, sample
        samples.append(sample)
    return samples


def draw_model(ax, data, cmap):
    mesh = ax.contourf(data['tau']/100, data['log_ratio'],
                       np.clip(data['f'], 0, 1), levels=LEVELS, cmap=cmap)
    map_axes(ax)
    ax.set_box_aspect(655/778)
    return mesh


def save(fig, stem):
    fig.savefig(OUTPUT / f'{stem}.png', dpi=220)
    fig.savefig(OUTPUT / f'{stem}.pdf')
    plt.close(fig)


def plot(data):
    cmap = palette()
    plt.rcParams.update({
        'font.size': 11, 'axes.titlesize': 12, 'mathtext.fontset': 'stix',
        'axes.linewidth': .8, 'xtick.direction': 'in', 'ytick.direction': 'in',
        'pdf.fonttype': 42,
    })
    fig, ax = plt.subplots(figsize=(7.2, 6.4), layout='constrained')
    mesh = draw_model(ax, data, cmap)
    ax.set_title(AMPLITUDES + '\n' + ETA_RANGE, fontsize=11)
    fig.colorbar(mesh, ax=ax, pad=.025, ticks=np.linspace(0, 1, 6),
                 label=r'$\sqrt{F(6\tau)}$')
    fig.suptitle(r'$\sqrt{F}:\quad K=2A(\mathbf{p})=A(2\mathbf{p})$', fontsize=17)
    fig.supxlabel('Matrix amplitudes absorb the factor 2; the six-stage protocol and duration are unchanged.',
                  fontsize=8)
    save(fig, 'sqrtF_K2A_doubled_matrix_parameters')

    fig, axes = plt.subplots(1, 2, figsize=(12.9, 6.0), layout='constrained')
    paper = np.asarray(Image.open(ROOT / 'paper_alignment' / 'paper_heatmap.png').convert('RGB'))
    axes[0].imshow(paper, origin='upper', extent=(0, 12, -1, 1),
                   aspect='auto', interpolation='nearest')
    map_axes(axes[0])
    axes[0].set_box_aspect(655/778)
    axes[0].set_title('Paper Fig. 1(d)\n' +
                      r'Caption: $E_1=0.01,\quad t_c=E_0=0.3\ \mathrm{meV}$')
    mesh = draw_model(axes[1], data, cmap)
    axes[1].set_title(r'Our $K=2A(\mathbf{p})=A(2\mathbf{p})$' + '\n' + AMPLITUDES)
    axes[1].set_ylabel('')
    fig.colorbar(mesh, ax=axes.tolist(), pad=.018, shrink=.87,
                 ticks=np.linspace(0, 1, 6), label=r'$f=\sqrt{F}$')
    fig.suptitle(r'$f=|\langle\psi_1^+(0)|\psi_1^-(6\tau)\rangle|$' +
                 ': paper and complete-state evolution', fontsize=16)
    fig.supxlabel("Right: primed values are matrix amplitudes after absorbing the factor 2. "
                  'Both panels use the same duration and coupling-ratio axes.', fontsize=9)
    save(fig, 'paper_vs_K2A_sqrtF')


def main():
    OUTPUT.mkdir(exist_ok=True)
    with np.load(SOURCE) as stored:
        data = {key: stored[key] for key in stored.files}
    metadata = json.loads((SOURCE.parent / 'verification.json').read_text())
    samples = check_equivalence(data, metadata)
    report = {
        'identity': 'K=2*A(E,eta,tc,E0)=A(2*E,2*eta,2*tc,2*E0)',
        'absorbed_generator_factor': 1.,
        'matrix_parameters_meV': {'E1': .02, 'eta_range': [.002, .2], 'tc': .6, 'E0': .6},
        'physical_caption_parameters_meV': {'E1': .01, 'eta_range': [.001, .1], 'tc': .3, 'E0': .3},
        'time_scaling': 1.,
        'end_time': '6*tau',
        'tau_range': [float(data['tau'][0]), float(data['tau'][-1])],
        'grid': list(data['f'].shape),
        'source_map': str(SOURCE.relative_to(ROOT.parent)),
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'source_reuse_reason': 'Every point has an identical generator after absorbing the factor 2.',
        'solver': {'method': 'DOP853', 'rtol': RTOL, 'atol': ATOL, 'max_local_step': MAX_STEP},
        'independent_doubled_parameter_wave_checks': samples,
        'paper_image': 'Original PDF crop; not raw fidelity data.',
    }
    (OUTPUT / 'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    plot(data)
    print(json.dumps({'output_directory': str(OUTPUT), 'wave_checks': samples}, indent=2), flush=True)


if __name__ == '__main__':
    main()
