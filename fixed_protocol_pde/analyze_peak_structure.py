#!/usr/bin/env python3
"""Propagate the complete bivector and controlled-duration responses.

The theory retains tau, tc and E0 as design parameters. Numerical slices here
select the Fig. 1 caption's E1=.01, tc=E0=.3, and the existing solver settings.
No slow projection or fitted coefficient enters the complete-state solution.
"""

import json
import os
from pathlib import Path

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/quantumsss-matplotlib')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.signal import find_peaks


ROOT = Path(__file__).resolve().parent
RTOL, ATOL, MAX_STEP = 2e-11, 2e-13, .025
PAIRS = np.triu_indices(5, 1)


def controls(stage, v, tc, e0):
    a = (1-np.cos(np.pi*v))/2
    b = 1-a
    if stage % 3 == 0:
        return tc*a, 0., e0*b, a
    if stage % 3 == 1:
        return tc*b, tc*a, 0., b
    return 0., tc*b, e0*a, 0.


def state_action(z, E, eta, p, r, ed, g):
    """N z, in the same canonical ten-component order as the derivation."""
    d, x1, x2, x3, y1, y2, y3, h1, h2, h3 = np.moveaxis(z, -1, 0)
    q = eta*g
    return np.stack((
        -p*x2-q*y3, r*x2+E*y1-q*h2,
        p*d-r*x1+ed*x3+E*y2-q*h3, -ed*x2+E*y3,
        -E*x1+r*y2+p*h1, -E*x2-r*y1+ed*y3,
        q*d-E*x3-ed*y2-p*h3, -p*y1+ed*h2,
        q*x1-ed*h1+r*h3, q*x2+p*y3-r*h2,
    ), axis=-1)


def readout(b, u=None, h=None):
    d = b[..., 0]
    transverse = b[..., 1:7]
    F = (d*d-d)/2+np.sum(transverse**2, axis=-1)/4
    ell = 1-d*d-np.sum(transverse**2, axis=-1)/2
    result = {'F': F, 'f': np.sqrt(np.maximum(F, 0)), 'd': d, 'ell': ell}
    if u is not None:
        result['F_tau'] = (d-.5)*u[..., 0]+np.sum(transverse*u[..., 1:7], axis=-1)/2
    if h is not None:
        result['F_tautau'] = (u[..., 0]**2+(d-.5)*h[..., 0]
            +np.sum(u[..., 1:7]**2+transverse*h[..., 1:7], axis=-1)/2)
    return result


def propagate(tau, E, eta, tc=.3, e0=.3, response=True, generator_factor=2.):
    """Complete state; optional kappa diagnoses the matrix convention only."""
    tau, E, eta = np.broadcast_arrays(tau, E, eta)
    shape = tau.shape
    tau, E, eta = [x.ravel() for x in (tau, E, eta)]
    blocks = 3 if response else 1
    data = np.zeros((blocks, len(tau), 10))
    data[0, :, 0] = 1
    for stage in range(6):
        def rhs(v, flat):
            z = flat.reshape(data.shape)
            p, r, ed, g = controls(stage, v, tc, e0)
            action = state_action(z, E, eta, p, r, ed, g)
            out = generator_factor*tau[None, :, None]*action
            if response:
                # d_v b_tau = kappa*tau*N b_tau + kappa*N b;
                # d_v b_tautau = kappa*tau*N b_tautau + 2*kappa*N b_tau.
                out[1] += generator_factor*action[0]
                out[2] += 2*generator_factor*action[1]
            return out.ravel()
        sol = solve_ivp(rhs, (0., 1.), data.ravel(), method='DOP853',
                        rtol=RTOL, atol=ATOL, max_step=MAX_STEP, t_eval=[1.])
        if not sol.success:
            raise RuntimeError(sol.message)
        data = sol.y[:, -1].reshape(data.shape)
    fields = [z.reshape(shape+(10,)) for z in data]
    result = readout(*fields)
    result['b'] = fields[0]
    return result


def independent_wave(tau, E, eta, tc=.3, e0=.3, generator_factor=2.):
    """Independent complex wavefunction, with first/second tau derivatives."""
    data = np.zeros((3, 5), dtype=complex)
    initial = np.array([1., -1j, 0., 0., 0.])/np.sqrt(2)
    data[0] = initial
    for stage in range(6):
        def rhs(v, flat):
            z = flat.reshape(3, 5)
            p, r, ed, g = controls(stage, v, tc, e0)
            a = np.stack((E*z[:, 1]+eta*g*z[:, 4], -E*z[:, 0]-p*z[:, 3],
                          r*z[:, 3], p*z[:, 1]-r*z[:, 2]+ed*z[:, 4],
                          -eta*g*z[:, 0]-ed*z[:, 3]), axis=-1)
            out = generator_factor*tau*a
            out[1] += generator_factor*a[0]
            out[2] += 2*generator_factor*a[1]
            return out.ravel()
        sol = solve_ivp(rhs, (0., 1.), data.ravel(), method='DOP853',
                        rtol=RTOL, atol=ATOL, max_step=MAX_STEP, t_eval=[1.])
        if not sol.success:
            raise RuntimeError(sol.message)
        data = sol.y[:, -1].reshape(3, 5)
    amp, amp_tau, amp_tautau = data@initial
    u, v = np.sqrt(2)*data[0].real, -np.sqrt(2)*data[0].imag
    B = np.outer(u, v)-np.outer(v, u)
    return {'b': B[PAIRS], 'F': abs(amp)**2,
            'F_tau': 2*np.real(amp_tau*np.conj(amp)),
            'F_tautau': 2*np.real(amp_tautau*np.conj(amp)+abs(amp_tau)**2)}


def slow_prediction(tau, E, eta):
    """Prior control-defined slow model, retaining all finite slow rotations.

    This comparison assumes tc=E0. It drops coupling to the two fast modes;
    it is not substituted into the complete-state calculation above.
    """
    tau, E, eta = np.broadcast_arrays(tau, E, eta)
    shape = tau.shape
    tau, E, eta = [x.ravel() for x in (tau, E, eta)]
    def gate(v):
        a = (1-np.cos(np.pi*v))/2
        return a, 1-a, np.hypot(a, 1-a)
    I0 = quad(lambda v: .5/gate(v)[2], 0., 1., epsabs=2e-14)[0]
    I1 = quad(lambda v: .5*gate(v)[2], 0., 1., epsabs=2e-14)[0]
    identity = np.broadcast_to(np.eye(3), (len(tau), 3, 3)).copy()
    def rhs(v, flat):
        a, b, rho = gate(v)
        R = flat.reshape(identity.shape)
        out = np.empty_like(R)
        h, k = 2*tau*E*a/rho, 2*tau*eta*b
        out[:, 0] = h[:, None]*R[:, 1]+k[:, None]*R[:, 2]
        out[:, 1] = -h[:, None]*R[:, 0]
        out[:, 2] = -k[:, None]*R[:, 0]
        return out.ravel()
    sol = solve_ivp(rhs, (0., 1.), identity.ravel(), method='DOP853',
                    rtol=RTOL, atol=ATOL, max_step=MAX_STEP)
    if not sol.success:
        raise RuntimeError(sol.message)
    R2 = sol.y[:, -1].reshape(identity.shape)
    def R12(angle):
        result = identity.copy()
        result[:, 0, 0] = result[:, 1, 1] = np.cos(angle)
        result[:, 0, 1], result[:, 1, 0] = np.sin(angle), -np.sin(angle)
        return result
    turn = np.array([[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]])
    U = R12(2*tau*E)@R2@turn@R12(2*tau*(I0*E-I1*eta))
    F = (1-(U@U)[:, 2, 2])**2/4
    return {'F': F.reshape(shape), 'I0': I0, 'I1': I1}


def refine_prominent_peaks(tau_axis, ratios, full):
    """Refine time maxima selected by grid prominence, not all small ripples."""
    candidates = []
    for row, ratio in enumerate(ratios):
        indices, info = find_peaks(full['f'][row], prominence=.15)
        for index, prominence in zip(indices, info['prominences']):
            lo, hi = tau_axis[index-1], tau_axis[index+1]
            if full['F_tau'][row, index-1] > 0 and full['F_tau'][row, index+1] < 0:
                candidates.append((row, ratio, lo, hi, prominence))
    if not candidates:
        return []
    eta = np.array([x[1]*.01 for x in candidates])
    lo = np.array([x[2] for x in candidates])
    hi = np.array([x[3] for x in candidates])
    for _ in range(28):
        mid = (lo+hi)/2
        result = propagate(mid, .01, eta)
        positive = result['F_tau'] > 0
        lo, hi = np.where(positive, mid, lo), np.where(positive, hi, mid)
    tau = (lo+hi)/2
    result = propagate(tau, .01, eta)
    peaks = []
    for i, (row, ratio, _, _, prominence) in enumerate(candidates):
        if result['F_tautau'][i] >= 0:
            continue
        peaks.append({'ratio': float(ratio), 'tau': float(tau[i]),
                      **{key: float(result[key][i]) for key in ('F', 'f', 'd', 'ell', 'F_tau', 'F_tautau')},
                      'grid_prominence': float(prominence)})
    return peaks


def plot_slices(tau_axis, ratios, full, slow_F, peaks):
    fig, axes = plt.subplots(2, 2, figsize=(11., 7.), layout='constrained')
    for row, ax in enumerate(axes.ravel()):
        ax.plot(tau_axis/100, full['f'][row], color='black', label=r'Complete $\mathcal{B}$')
        ax.plot(tau_axis/100, np.sqrt(np.maximum(slow_F[row], 0)), '--', color='#D55E00', label='Slow rotations')
        points = [p for p in peaks if abs(p['ratio']-ratios[row]) < 1e-9]
        ax.scatter([p['tau']/100 for p in points], [p['f'] for p in points],
                   s=15, color='#0072B2', zorder=3, label='Refined prominent peaks')
        ax.set(xlabel=r'$\tau\ (100/\mathrm{meV})$', ylabel=r'$f=\sqrt{F}$',
               title=rf'$|t_1|/E_1={ratios[row]:.5g}$', ylim=(-.03, 1.03), xlim=(0, 12))
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(r'$E_1=0.01$ meV, $t_c=E_0=0.3$ meV; $K=2A$; two cycles')
    fig.savefig(ROOT/'peak_structure.png', dpi=180)
    fig.savefig(ROOT/'peak_structure.pdf')
    plt.close(fig)


def main():
    sample_cases = [(20., .01, .001), (300., .01, .01), (1000., .01, .1)]
    errors = {'b': 0., 'F': 0., 'F_tau': 0., 'F_tautau': 0.}
    for tau, E, eta in sample_cases:
        full = propagate(tau, E, eta)
        wave = independent_wave(tau, E, eta)
        for key in errors:
            errors[key] = max(errors[key], float(np.max(abs(full[key]-wave[key]))))
    if max(errors.values()) > 2e-8:
        raise AssertionError(errors)
    print('Independent wavefunction and duration-response checks passed.', flush=True)

    tau_axis = np.linspace(0., 1200., 241)
    ratios = np.array([.1, 1., 2.344645627, 10.])
    tau, ratio = np.meshgrid(tau_axis, ratios)
    full = propagate(tau, .01, .01*ratio)
    slow = slow_prediction(tau, .01, .01*ratio)
    print('Complete-state and comparison slices propagated.', flush=True)
    peaks = refine_prominent_peaks(tau_axis, ratios, full)
    report = {
        'E1': .01, 'tc': .3, 'E0': .3, 'hbar': 1,
        'generator_convention': 'K=2*A', 'tau_range': [0., 1200.],
        'ratios': ratios.tolist(), 'tau_grid_spacing': 5.,
        'solver': {'method': 'DOP853', 'rtol': RTOL, 'atol': ATOL, 'max_local_step': MAX_STEP},
        'independent_wave_errors': errors,
        'peak_selection': 'Grid prominence in f >= .15, bracketed positive-to-negative F_tau crossing, F_tautau<0; does not enumerate all weak ripples.',
        'prominent_peaks': peaks,
        'slow_comparison': {'I0': slow['I0'], 'I1': slow['I1'],
            'max_F_error_tau_ge_100': float(np.max(abs(full['F']-slow['F'])[tau >= 100.]))},
        'paper_comparison': {
            'axes': 'tau/(100/meV), log10(abs(t1)/E1)', 'color': 'f=sqrt(F)',
            'caption_E1_meV': .01, 'body_E1_meV': .001,
            'scope': 'Model predictions under K=2*A; no quantitative reproduction or parameter fit to the paper image is claimed.'},
    }
    (ROOT/'peak_structure.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez_compressed(ROOT/'peak_structure.npz', tau=tau_axis, ratios=ratios,
                        **full, slow_F=slow['F'])

    plot_slices(tau_axis, ratios, full, slow['F'], peaks)
    print(json.dumps({'independent_wave_errors': errors,
                      'prominent_peak_count': len(peaks),
                      'slow_max_F_error_tau_ge_100': report['slow_comparison']['max_F_error_tau_ge_100']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
