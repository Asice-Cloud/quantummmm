#!/usr/bin/env python3
"""Predict the double-braid stripes from control-defined slow rotations.

No intermediate state, fitted coefficient, or propagator from the full model
enters the prediction. Existing full-model data are loaded only for comparison.
The approximation discards fast/slow mixing; finite E*tau and eta*tau are kept.

Run: .venv/bin/python -B geometric_validation/predict_geometric_stripes.py
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

from verify_bottleneck import PSI_INITIAL, PSI_TARGET, controls, wave_rhs


ROOT = Path(__file__).resolve().parent
TAU, TC = 100., .3
C = 2*TAU
RTOL, ATOL = 2e-11, 2e-13
TURN = np.array([[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]])


def gate(v):
    a = (1-np.cos(np.pi*v))/2
    b = 1-a
    return a, b, np.hypot(a, b)


I0 = quad(lambda v: .5/gate(v)[2], 0, 1, epsabs=2e-14)[0]
I1 = quad(lambda v: .5*gate(v)[2], 0, 1, epsabs=2e-14)[0]


def rotation12(angle):
    angle = np.atleast_1d(angle)
    result = np.broadcast_to(np.eye(3), (len(angle), 3, 3)).copy()
    result[:, 0, 0] = result[:, 1, 1] = np.cos(angle)
    result[:, 0, 1], result[:, 1, 0] = np.sin(angle), -np.sin(angle)
    return result


def slow_generator(v, E, eta):
    a, b, rho = gate(v)
    G = np.zeros(np.broadcast_shapes(np.shape(E), np.shape(eta))+(3, 3))
    G[..., 0, 1], G[..., 1, 0] = C*E*a/rho, -C*E*a/rho
    G[..., 0, 2], G[..., 2, 0] = C*eta*b, -C*eta*b
    return G


def predict(E, eta):
    """Integrate only the second-stage SO(3); other rotations are quadratures."""
    E, eta = np.broadcast_arrays(np.atleast_1d(E), np.atleast_1d(eta))
    E, eta = E.ravel(), eta.ravel()
    shape = (len(E), 3, 3)
    def rhs(v, flat):
        return (slow_generator(v, E, eta)@flat.reshape(shape)).ravel()
    sol = solve_ivp(rhs, (0., 1.), np.broadcast_to(np.eye(3), shape).ravel(),
                    method='DOP853', rtol=RTOL, atol=ATOL, max_step=.04)
    assert sol.success, sol.message
    R2 = sol.y[:, -1].reshape(shape)
    theta1, theta3 = C*(I0*E-I1*eta), C*E
    U = rotation12(theta3)@R2@TURN@rotation12(theta1)
    d = (U@U)[:, 2, 2]
    # Euler chart R2=R12(x) R23(y) R13(z); valid when cos(y)>0.
    x = np.arctan2(R2[:, 0, 1], R2[:, 1, 1])
    y = -np.arcsin(np.clip(R2[:, 2, 1], -1., 1.))
    z = np.arctan2(-R2[:, 2, 0], R2[:, 2, 2])
    phase = theta1+theta3+x-z
    F = (1-d)**2/4
    factorized = np.cos(y)**4*np.cos(phase/2)**4
    return {'F': F, 'd': d, 'U': U, 'R2': R2, 'x': x, 'y': y, 'z': z,
            'alpha': -y, 'phase': phase, 'factorization_error': abs(F-factorized),
            'linear_phase': C*((1+2*I0)*E-(I1+.5)*eta)}


def frame(stage, v):
    """An analytic SO(5) frame determined by the ideal controls, including axes."""
    a, b, rho = gate(v)
    s, t = a/rho, b/rho
    theta_prime = (np.pi/2)*np.sin(np.pi*v)/rho**2
    e1, e2, e3, ea, eb = np.eye(5)
    if stage % 3 == 0:
        n, m = s*e2+t*eb, t*e2-s*eb
        np_, mp = theta_prime*m, -theta_prime*n
        Q = np.column_stack((ea, n, e1, m, e3))
        Qp = np.column_stack((0*ea, np_, 0*e1, mp, 0*e3))
    elif stage % 3 == 1:
        n, m = t*e2-s*e3, s*e2+t*e3
        np_, mp = -theta_prime*m, theta_prime*n
        Q = np.column_stack((ea, n, e1, m, eb))
        Qp = np.column_stack((0*ea, np_, 0*e1, mp, 0*eb))
    else:
        n, m = -t*e3+s*eb, s*e3+t*eb
        np_, mp = theta_prime*m, -theta_prime*n
        Q = np.column_stack((ea, n, e1, e2, m))
        Qp = np.column_stack((0*ea, np_, 0*e1, 0*e2, mp))
    return Q, Qp


def exact_frame_generator(stage, v, E, eta):
    p, r, ed, g = controls(stage, v)
    A = np.zeros((5, 5))
    for i, j, coefficient in ((0, 1, E), (0, 4, eta*g),
                             (1, 3, -p), (2, 3, r), (3, 4, ed)):
        A[i, j], A[j, i] = coefficient, -coefficient
    Q, Qp = frame(stage, v)
    return C*Q.T@A@Q-Q.T@Qp


def checks():
    """Check the frame algebra and angles against independent evolutions."""
    errors = dict(frame_orthogonality=0., frame_derivative=0., frame_determinant=0.,
                  slow_block=0., fast_block=0., frame_wave=0., norm=0., null=0.,
                  euler_rotation=0., euler_phase_integral=0., angle_gradient=0.)
    points = [(0., 0.), (.01, .023), (.02, .023), (.027, .001), (.03, .05)]
    sample_rows = []
    for E, eta in points:
        for stage in range(3):
            for v in np.linspace(.003, .997, 31):
                Q, Qp = frame(stage, v)
                diff = (frame(stage, v+1e-6)[0]-frame(stage, v-1e-6)[0])/2e-6
                a, b, rho = gate(v)
                G = exact_frame_generator(stage, v, E, eta)
                if stage == 0:
                    expected = np.zeros((3, 3))
                    expected[0, 1] = C*(E*b-eta*a*a)/rho
                    expected[1, 0] = -expected[0, 1]
                elif stage == 1:
                    expected = slow_generator(v, E, eta)
                else:
                    expected = np.array([[0., C*E, 0.], [-C*E, 0., 0.], [0., 0., 0.]])
                values = {'frame_orthogonality': np.max(abs(Q.T@Q-np.eye(5))),
                          'frame_determinant': abs(np.linalg.det(Q)-1),
                          'frame_derivative': np.max(abs(diff-Qp)),
                          'slow_block': np.max(abs(G[2:, 2:]-expected)),
                          'fast_block': np.max(abs(G[:2, :2]-C*TC*rho*np.array([[0., 1.], [-1., 0.]])))}
                for key, value in values.items():
                    errors[key] = max(errors[key], float(value))
        case = {'E1': E, 'eta': eta, 'tau': TAU}
        physical = PSI_INITIAL.copy()
        framed = frame(0, 0.)[0].T@physical
        for stage in range(6):
            physical_sol = solve_ivp(lambda v, psi: wave_rhs(stage, v, psi, case),
                (0., 1.), physical, method='DOP853', rtol=8e-13, atol=8e-15)
            frame_sol = solve_ivp(lambda v, psi: exact_frame_generator(stage, v, E, eta)@psi,
                (0., 1.), framed, method='DOP853', rtol=8e-13, atol=8e-15)
            assert physical_sol.success and frame_sol.success
            physical = physical_sol.y[:, -1]
            framed_physical = frame(stage, 1.)[0]@frame_sol.y[:, -1]
            errors['frame_wave'] = max(errors['frame_wave'], float(np.max(abs(physical-framed_physical))))
            framed = frame(stage+1, 0.)[0].T@framed_physical
        errors['norm'] = max(errors['norm'], float(abs(np.vdot(physical, physical)-1)))
        errors['null'] = max(errors['null'], float(abs(np.sum(physical**2))))
        estimate = predict(E, eta)
        # Euler equations also integrate the finite, noncommuting correction.
        def rhs(v, value):
            x, y, z, correction = value[:4]
            a, b, rho = gate(v)
            h, g = C*E*a/rho, C*eta*b
            cx, sx, ty, sy = np.cos(x), np.sin(x), np.tan(y), 1/np.cos(y)
            dx, dy, dz = h-g*cx*ty, g*sx, g*cx*sy
            J = np.array([[g*sx*ty, -g*cx*sy**2, 0.],
                          [g*cx, 0., 0.], [-g*sx*sy, g*cx*sy*ty, 0.]])
            sensitivity = value[4:].reshape(3, 2)
            forcing = np.array([[C*a/rho, -C*b*cx*ty], [0., C*b*sx], [0., C*b*cx*sy]])
            ds = J@sensitivity+forcing
            return np.r_[dx, dy, dz, -g*(cx*(sy+ty)-1), ds.ravel()]
        sol = solve_ivp(rhs, (0., 1.), np.zeros(10), method='DOP853',
                        rtol=2e-12, atol=2e-14, max_step=.025)
        assert sol.success
        x, y, z, correction = sol.y[:4, -1]
        R23 = np.array([[1., 0., 0.], [0., np.cos(y), np.sin(y)], [0., -np.sin(y), np.cos(y)]])
        R13 = np.array([[np.cos(z), 0., np.sin(z)], [0., 1., 0.], [-np.sin(z), 0., np.cos(z)]])
        reconstructed = rotation12(x)[0]@R23@R13
        errors['euler_rotation'] = max(errors['euler_rotation'], float(np.max(abs(reconstructed-estimate['R2'][0]))))
        phase = C*((1+I0)*E-I1*eta)+x-z
        errors['euler_phase_integral'] = max(errors['euler_phase_integral'], float(abs(phase-estimate['linear_phase'][0]-correction)))
        deriv = sol.y[4:, -1].reshape(3, 2)
        phase_grad = C*np.array([1+I0, -I1])+deriv[0]-deriv[2]
        F_grad = (-4*np.cos(y)**3*np.sin(y)*np.cos(phase/2)**4*deriv[1]
                  -2*np.cos(y)**4*np.cos(phase/2)**3*np.sin(phase/2)*phase_grad)
        step = 1e-6
        fd = np.array([(predict(E+step, eta)['F'][0]-predict(E-step, eta)['F'][0])/(2*step),
                       (predict(E, eta+step)['F'][0]-predict(E, eta-step)['F'][0])/(2*step)])
        errors['angle_gradient'] = max(errors['angle_gradient'], float(np.max(abs(F_grad-fd))))
        sample_rows.append({'E': E, 'eta': eta, 'full_F': float(abs(np.vdot(PSI_TARGET, physical))**2),
                            'predicted_F': float(estimate['F'][0]), 'y': float(y),
                            'phase': float(phase), 'phase_correction': float(correction),
                            'gradient_E': float(F_grad[0]), 'gradient_eta': float(F_grad[1])})
    assert max(v for k, v in errors.items() if k != 'angle_gradient') < 2e-8, errors
    assert errors['angle_gradient'] < 2e-5, errors
    return {'errors': errors, 'independent_full_wave_samples': sample_rows}


def refine_peaks(E_axis, eta_axis, F):
    candidates = [(E, eta_axis[row-1], eta_axis[row+1])
                  for col, E in enumerate(E_axis) for row in range(1, len(eta_axis)-1)
                  if F[row, col] >= max(F[row-1, col], F[row+1, col])]
    E = np.array([row[0] for row in candidates])
    lo, hi = np.array([row[1] for row in candidates]), np.array([row[2] for row in candidates])
    ratio = (np.sqrt(5)-1)/2
    for _ in range(30):
        left, right = hi-ratio*(hi-lo), lo+ratio*(hi-lo)
        choose = predict(E, left)['F'] < predict(E, right)['F']
        lo, hi = np.where(choose, left, lo), np.where(choose, hi, right)
    eta = (lo+hi)/2
    values = predict(E, eta)
    return [{'E': float(E[i]), 'eta': float(eta[i]), 'F': float(values['F'][i])}
            for i in range(len(E))]


def refine_phase_lines(E_axis, eta_axis, phase):
    residual = np.sin(phase)
    brackets = [(E, eta_axis[row], eta_axis[row+1])
                for col, E in enumerate(E_axis) for row in range(len(eta_axis)-1)
                if residual[row, col]*residual[row+1, col] < 0]
    E = np.array([row[0] for row in brackets])
    lo, hi = np.array([row[1] for row in brackets]), np.array([row[2] for row in brackets])
    left = np.sin(predict(E, lo)['phase'])
    for _ in range(22):
        mid = (lo+hi)/2
        value = np.sin(predict(E, mid)['phase'])
        replace = left*value > 0
        lo, left = np.where(replace, mid, lo), np.where(replace, value, left)
        hi = np.where(replace, hi, mid)
    eta = (lo+hi)/2
    values = predict(E, eta)
    result = [{'E': float(E[i]), 'eta': float(eta[i]), 'F': float(values['F'][i])}
              for i in range(len(E)) if np.cos(values['phase'][i]) > 0]
    # Include exact roots on the domain edge, notably E=eta=0.
    for row in (0, len(eta_axis)-1):
        for col in np.flatnonzero((abs(residual[row]) < 1e-11) & (np.cos(phase[row]) > 0)):
            result.append({'E': float(E_axis[col]), 'eta': float(eta_axis[row]),
                           'F': float(predict(E_axis[col], eta_axis[row])['F'][0])})
    return result


def compare_branches(predicted, full):
    comparisons = []
    for row in predicted:
        candidates = [point for point in full if abs(point['E']-row['E']) < 1e-10]
        if not candidates:
            continue
        closest = min(candidates, key=lambda point: abs(point['eta']-row['eta']))
        if abs(closest['eta']-row['eta']) > .003:
            continue
        comparisons.append({'E': row['E'], 'predicted_eta': row['eta'],
                            'full_eta': closest['eta'], 'eta_error': row['eta']-closest['eta'],
                            'predicted_F': row['F'],
                            **({'full_F': closest['F'], 'height_error': row['F']-closest['F']}
                               if 'F' in closest else {})})
    return comparisons


def plots(full, model, phases):
    E_axis, eta_axis = full['E_axis'], full['eta_axis']
    extent = [E_axis[0], E_axis[-1], eta_axis[0], eta_axis[-1]]
    fig, axes = plt.subplots(1, 3, figsize=(14., 4.5), layout='constrained')
    for ax, values, title in zip(axes, (full['F'], model['F'], abs(full['F']-model['F'])),
                                ('Full five-mode result', 'Control-defined slow rotations', 'Absolute fidelity error')):
        error = ax is axes[2]
        im = ax.imshow(values, extent=extent, origin='lower', aspect='auto',
                       vmin=0., vmax=.018 if error else 1., cmap='magma' if error else 'viridis',
                       interpolation='nearest')
        ax.set(xlabel=r'$E_1$ (meV)', ylabel=r'$t_1=\eta$ (meV)', title=title)
        fig.colorbar(im, ax=ax, label=r'$|F_{\rm slow}-F|$' if error else r'$F=f^2$')
    for branch in range(-2, 3):
        eta = ((1+2*I0)*E_axis-2*np.pi*branch/C)/(I1+.5)
        mask = (eta >= eta_axis[0]) & (eta <= eta_axis[-1])
        axes[1].plot(E_axis[mask], eta[mask], '--', color='white', lw=.8, alpha=.8)
    axes[1].text(.03, .97, 'Dashed: leading phase alignment', transform=axes[1].transAxes,
                 ha='left', va='top', color='white', fontsize=9)
    fig.savefig(ROOT/'stripe_prediction.png', dpi=210)
    fig.savefig(ROOT/'stripe_prediction.pdf')
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11., 8.), layout='constrained')
    col = np.argmin(abs(E_axis-.01))
    axes[0, 0].plot(eta_axis, full['F'][:, col], color='black', label='Full model')
    axes[0, 0].plot(eta_axis, model['F'][:, col], '--', color='#D55E00', label='Slow prediction')
    axes[0, 0].plot(eta_axis, np.cos(model['linear_phase'][:, col]/2)**4, ':',
                    color='#0072B2', label='Leading phase alone')
    axes[0, 0].set(ylabel=r'$F=f^2$', ylim=(-.02, 1.02), title=r'$E_1=0.01$ meV')
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].plot(eta_axis, full['actual_latitude'][:, col], color='black', label='Full interface latitude')
    axes[0, 1].plot(eta_axis, model['alpha'][:, col], '--', color='#D55E00', label=r'Prediction $\alpha=-y$')
    axes[0, 1].set(ylabel='Latitude (rad)', title='Noncommuting rotations lower the envelope')
    axes[0, 1].legend(fontsize=8)
    exact_phase = np.unwrap(-np.pi*full['azimuth_over_pi'][:, col])
    exact_phase += 2*np.pi*np.round((model['linear_phase'][0, col]-exact_phase[0])/(2*np.pi))
    slow_phase = np.unwrap(model['phase'][:, col])
    slow_phase += 2*np.pi*np.round((model['linear_phase'][0, col]-slow_phase[0])/(2*np.pi))
    axes[1, 0].plot(eta_axis, exact_phase, color='black', label='Full interface phase')
    axes[1, 0].plot(eta_axis, slow_phase, '--', color='#D55E00', label='Slow prediction')
    axes[1, 0].plot(eta_axis, model['linear_phase'][:, col], ':', color='#0072B2', label='Leading quadrature')
    axes[1, 0].set(ylabel=r'$\Phi=-\Delta\varphi$ (rad)', title='Finite phase accumulation')
    axes[1, 0].legend(fontsize=8)
    for E, color in ((.005, '#0072B2'), (.02, '#009E73'), (.025, '#CC79A7')):
        j = np.argmin(abs(E_axis-E))
        axes[1, 1].plot(eta_axis, full['F'][:, j], color=color, label=rf'Full $E_1={E:g}$')
        axes[1, 1].plot(eta_axis, model['F'][:, j], '--', color=color)
    axes[1, 1].set(ylabel=r'$F=f^2$', ylim=(-.02, 1.02), title='Other slices: dashed curves are predictions')
    axes[1, 1].legend(fontsize=8)
    for ax in axes.ravel():
        ax.set(xlabel=r'$t_1=\eta$ (meV)', xlim=(0., .05))
        ax.grid(alpha=.15)
    fig.savefig(ROOT/'stripe_angles.png', dpi=210)
    fig.savefig(ROOT/'stripe_angles.pdf')
    plt.close(fig)


def main():
    validation = checks()
    print('Analytic frame and independent angle checks passed.', flush=True)
    with np.load(ROOT/'parameter_geometry.npz') as saved:
        full = {key: saved[key] for key in saved.files}
    E, eta = np.meshgrid(full['E_axis'], full['eta_axis'])
    shape = E.shape
    flat = predict(E, eta)
    model = {key: value.reshape(shape+value.shape[1:]) for key, value in flat.items()}
    error = abs(model['F']-full['F'])
    phase_error = np.angle(np.exp(1j*(model['phase']+np.pi*full['azimuth_over_pi'])))
    linear_error = np.angle(np.exp(1j*(model['linear_phase']+np.pi*full['azimuth_over_pi'])))
    peaks = refine_peaks(full['E_axis'], full['eta_axis'], model['F'])
    phases = refine_phase_lines(full['E_axis'], full['eta_axis'], model['phase'])
    full_report = json.loads((ROOT/'parameter_geometry.json').read_text())
    peak_matches = compare_branches(peaks, full_report['local_maxima'])
    phase_matches = compare_branches(phases, [p for p in full_report['refined_contours']
                                             if p['quantity'] == 'phase_alignment'])
    selected = []
    for E0 in (.005, .01, .015, .02, .025, .03):
        candidates = [p for p in peak_matches if abs(p['E']-E0) < 1e-10]
        selected.append(max(candidates, key=lambda p: p['predicted_F']))
    maximum = np.unravel_index(np.argmax(error), shape)
    report = {'tau': TAU, 'tc': TC, 'point_count': int(E.size),
              'I0': I0, 'I1': I1,
              'leading_stripe_slope': (1+2*I0)/(I1+.5),
              'leading_stripe_spacing_meV': 2*np.pi/(C*(I1+.5)),
              'F_error_max': float(error.max()), 'F_error_rms': float(np.sqrt(np.mean(error**2))),
              'F_error_quantiles': np.quantile(error, [.5, .9, .99]).tolist(),
              'worst_F_point': {'E': float(E[maximum]), 'eta': float(eta[maximum]),
                               'full_F': float(full['F'][maximum]), 'predicted_F': float(model['F'][maximum])},
              'latitude_error_max_rad': float(np.max(abs(model['alpha']-full['actual_latitude']))),
              'phase_error_max_rad': float(np.max(abs(phase_error))),
              'phase_error_rms_rad': float(np.sqrt(np.mean(phase_error**2))),
              'linear_phase_error_max_rad': float(np.max(abs(linear_error))),
              'linear_phase_error_rms_rad': float(np.sqrt(np.mean(linear_error**2))),
              'minimum_cos_y': float(np.min(np.cos(model['y']))),
              'orthogonality_error_max': float(np.max(abs(flat['U'].transpose(0, 2, 1)@flat['U']-np.eye(3)))),
              'factorization_error_max': float(np.max(flat['factorization_error'])),
              'peaks': peaks, 'phase_alignments': phases,
              'peak_comparisons': peak_matches, 'phase_comparisons': phase_matches,
              'selected_best_peaks': selected, 'validation': validation}
    for key, rows in (('peak', peak_matches), ('phase_alignment', phase_matches)):
        errors = np.array([p['eta_error'] for p in rows])
        report[key+'_matched_count'] = len(rows)
        report[key+'_eta_error_max_meV'] = float(np.max(abs(errors)))
        report[key+'_eta_error_rms_meV'] = float(np.sqrt(np.mean(errors**2)))
    (ROOT/'stripe_prediction.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez_compressed(ROOT/'stripe_prediction.npz', E_axis=full['E_axis'], eta_axis=full['eta_axis'], **model)
    plots(full, model, phases)
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('peaks', 'phase_alignments', 'peak_comparisons', 'phase_comparisons', 'validation')}, indent=2), flush=True)


if __name__ == '__main__':
    main()
