#!/usr/bin/env python3
"""Exact geometry and parameter sensitivities for B(3tau) -> B(5tau).

Propagate the state and its parameter derivatives, transport the terminal
target backward, and resolve fidelity/direction derivatives by time segment.
Independent checks use complex-wavefunction propagation and sensitivities,
including finite input/operator swaps at 3tau. No small-parameter or
adiabatic approximation is made.

Run from the project root:
  .venv/bin/python -B geometric_validation/verify_second_cycle_chain.py
"""

import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/quantumsss-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad_vec, solve_ivp

from verify_bottleneck import (
    B_INITIAL, CASES, INDICES, PSI_INITIAL, PSI_TARGET,
    bivector_matrix, bivector_rhs, controls, wave_rhs,
)


ROOT = Path(__file__).resolve().parent
SAMPLE = np.linspace(0, 1, 201)


def commutator(A, B):
    return A@B-B@A


def parameter_generators(stage, v, case):
    _, _, _, gate = controls(stage, v)
    generators = np.zeros((2, 5, 5))
    generators[0, 0, 1], generators[0, 1, 0] = 1., -1.
    generators[1, 0, 4], generators[1, 4, 0] = gate, -gate
    return 2*case['tau']*generators


def sensitivity_rhs(stage, v, values, case):
    rows = values.reshape(3, 10)
    B = bivector_matrix(rows[0])
    sources = parameter_generators(stage, v, case)
    result = np.array([bivector_rhs(stage, v, row, case) for row in rows])
    for j in range(2):
        result[j+1] += commutator(sources[j], B)[INDICES]
    return result.ravel()


def integrate(rhs, initial, backward=False):
    interval = (1., 0.) if backward else (0., 1.)
    solution = solve_ivp(rhs, interval, initial, method='DOP853',
                         rtol=8e-13, atol=8e-15,
                         max_step=.01, dense_output=True)
    assert solution.success, solution.message
    return solution


def fidelity_to(B, target):
    P, Pt = -B@B, -target@target
    return np.trace(Pt@P-target@B)/4


def fidelity_gradient(B, target):
    Pt = -target@target
    return ((B@Pt+Pt@B+target)/2)[INDICES]


def weight_matching_bound(k, target_k):
    k, target_k = np.clip([k, target_k], 0., 1.)
    return (1+np.sqrt((1-k)*(1-target_k))+np.sqrt(k*target_k))**2/4


def saturating_plane(target, k):
    """A legal plane attaining the fixed-P33 overlap bound."""
    e3 = np.eye(5)[2]
    P = -target@target
    target_k = P[2, 2]
    if target_k > 1e-12:
        u = P@e3/np.linalg.norm(P@e3)
    else:
        column = np.argmax(np.diag(P))
        u = P[:, column]/np.linalg.norm(P[:, column])
    v = -target@u
    v[2] = 0.
    v /= np.linalg.norm(v)
    a = u-u[2]*e3
    if np.linalg.norm(a) < 1e-8:
        a = np.eye(5)[np.argmin(abs(v)+10*e3)]
    a -= (a@v)*v
    a[2] = 0.
    a /= np.linalg.norm(a)
    u_new = np.sqrt(1-k)*a+np.sqrt(k)*e3
    return np.outer(u_new, v)-np.outer(v, u_new)


def three_mode_geometry(B, target):
    pairs = (np.array([0, 0, 1]), np.array([1, 2, 2]))
    b, g = B[pairs], target[pairs]
    beta, gamma = np.linalg.norm(b), np.linalg.norm(g)
    cosine = np.clip(b@g/(beta*gamma), -1., 1.)
    angle = np.arccos(cosine)
    sphere_F = (1+cosine)**2/4
    # The density matrix is 3/2-Lipschitz in the bivector Frobenius norm.
    delta_B = 2*np.sqrt(max(0., 1-beta))
    delta_G = 2*np.sqrt(max(0., 1-gamma))
    error_bound = 1.5*(delta_B+delta_G)
    exact_F = fidelity_to(B, target)
    assert abs(sphere_F-exact_F) <= error_bound+5e-11
    P, Pt = -B@B, -target@target
    return {'B3_three_mode_vector': b.tolist(),
            'required_B3_three_mode_vector': g.tolist(),
            'B3_three_mode_norm': float(beta),
            'required_B3_three_mode_norm': float(gamma),
            'angle_over_pi': float(angle/np.pi),
            'three_mode_F': float(sphere_F),
            'absolute_F_error': float(abs(sphere_F-exact_F)),
            'conservative_absolute_error_bound': float(error_bound),
            'plane_position_loss': float(np.linalg.norm(P-Pt)**2/8),
            'oriented_bivector_loss': float(np.linalg.norm(B-target)**2/8)}


def direct_wave_sensitivity(case):
    initial = np.zeros((3, 5), dtype=complex)
    initial[0] = PSI_INITIAL
    def rhs(stage, v, values):
        rows = values.reshape(3, 5)
        result = np.array([wave_rhs(stage, v, row, case) for row in rows])
        _, _, _, gate = controls(stage, v)
        state = rows[0]
        result[1] += 2*case['tau']*np.array(
            [state[1], -state[0], 0., 0., 0.])
        result[2] += 2*case['tau']*gate*np.array(
            [state[4], 0., 0., 0., -state[0]])
        return result.ravel()
    current = initial.ravel()
    for stage in range(5):
        solution = integrate(lambda v, y: rhs(stage, v, y), current)
        current = solution.y[:, -1]
    rows = current.reshape(3, 5)
    psi, derivatives = rows[0], rows[1:]
    a = np.vdot(PSI_TARGET, psi)
    da = derivatives@PSI_TARGET.conj()
    dF = 2*np.real(a.conjugate()*da)
    dd = 2*np.imag(derivatives[:, 0]*psi[1].conjugate()
                  + psi[0]*derivatives[:, 1].conjugate())
    return psi, dF, dd


def check_case(case, reference):
    initial = np.zeros((3, 10))
    initial[0] = B_INITIAL[INDICES]
    current = initial.ravel()
    forward = []
    for stage in range(5):
        solution = integrate(
            lambda v, y: sensitivity_rhs(stage, v, y, case), current)
        forward.append(solution)
        current = solution.y[:, -1]
    terminal = current.reshape(3, 10)
    B5 = bivector_matrix(terminal[0])
    target5 = -B_INITIAL
    grad5 = fidelity_gradient(B5, target5)
    backward = [None]*5
    current = np.concatenate((target5[INDICES], grad5))
    for stage in range(4, -1, -1):
        def rhs(v, y):
            return np.concatenate([bivector_rhs(stage, v, row, case)
                                   for row in y.reshape(2, 10)])
        solution = integrate(rhs, current, backward=True)
        backward[stage] = solution
        current = solution.y[:, -1]

    errors = {'wave_state': 0., 'backward_overlap_constant': 0.,
              'geometric_adjoint': 0., 'target_plane_constraint': 0.,
              'first_stage_F_bound': 0., 'first_stage_d_bound': 0.}
    endpoint_F = fidelity_to(B5, target5)
    segment_terms = []
    for stage in range(5):
        def integrand(v):
            rows = forward[stage].sol(v).reshape(3, 10)
            B = bivector_matrix(rows[0])
            target, adjoint = backward[stage].sol(v).reshape(2, 10)
            pd = -target
            force = np.array([commutator(A, B)[INDICES]
                              for A in parameter_generators(stage, v, case)])
            return np.concatenate((force@adjoint, force@pd))
        integral, quadrature_error = quad_vec(
            integrand, 0, 1, epsabs=2e-10, epsrel=3e-12, limit=600)
        segment_terms.append({'stage': stage+1,
            'F_E': float(integral[0]), 'F_eta': float(integral[1]),
            'd_E': float(integral[2]), 'd_eta': float(integral[3]),
            'estimated_quadrature_error': float(quadrature_error)})
        for v in (0., .19, .53, .87, 1.):
            b = forward[stage].sol(v).reshape(3, 10)[0]
            B = bivector_matrix(b)
            target, adjoint = backward[stage].sol(v).reshape(2, 10)
            G = bivector_matrix(target)
            errors['backward_overlap_constant'] = max(
                errors['backward_overlap_constant'],
                abs(fidelity_to(B, G)-endpoint_F))
            errors['geometric_adjoint'] = max(
                errors['geometric_adjoint'],
                float(np.max(abs(fidelity_gradient(B, G)-adjoint))))
            errors['target_plane_constraint'] = max(
                errors['target_plane_constraint'], float(np.max(abs(G@G@G+G))))

    boundary_rows = [forward[0].sol(0).reshape(3, 10)]
    boundary_rows += [solution.sol(1).reshape(3, 10) for solution in forward]
    B3, B4 = bivector_matrix(boundary_rows[3][0]), bivector_matrix(boundary_rows[4][0])
    G3 = bivector_matrix(backward[3].sol(0).reshape(2, 10)[0])
    G4 = bivector_matrix(backward[4].sol(0).reshape(2, 10)[0])
    P33 = float((-B3@B3)[2, 2])
    radius = np.sqrt(max(0., 1-P33))
    F_bound = (1+radius)**2/4
    required_P33 = float((-G4@G4)[2, 2])
    future_F_bound = weight_matching_bound(P33, required_P33)
    maximizer = saturating_plane(G4, P33)
    errors['weight_matching_bound_saturation'] = float(abs(
        fidelity_to(maximizer, G4)-future_F_bound))
    errors['weight_matching_bound_violation'] = float(max(
        0., endpoint_F-future_F_bound))
    errors['saturating_plane_P33'] = float(abs((-maximizer@maximizer)[2, 2]-P33))
    geometry = three_mode_geometry(B3, G3)
    errors['geometric_loss_identity'] = abs(
        geometry['plane_position_loss']+geometry['oriented_bivector_loss']
        -(1-endpoint_F))
    B4_path = bivector_matrix(forward[3].sol(SAMPLE)[:10].T)
    errors['first_stage_d_bound'] = float(max(
        0., np.max(abs(B4_path[:, 0, 1]))-radius))
    errors['first_stage_F_bound'] = float(max(
        0., max(fidelity_to(B, target5) for B in B4_path)-F_bound))

    sensitivities = {}
    target3, adjoint3 = backward[3].sol(0).reshape(2, 10)
    for j, parameter in enumerate(('E', 'eta')):
        for label, covector, grad in (
                ('F', adjoint3, grad5), ('d', -target3, np.eye(10)[0])):
            preparation = float(covector@boundary_rows[3][j+1])
            first = segment_terms[3][label+'_'+parameter]
            second = segment_terms[4][label+'_'+parameter]
            total = float(grad@terminal[j+1])
            sensitivities[label+'_'+parameter] = {
                'preparation_0_to_3tau': preparation,
                'first_3_to_4tau': first, 'second_4_to_5tau': second,
                'total': total,
            }
            errors['sensitivity_partition_'+label+'_'+parameter] = abs(
                preparation+first+second-total)
            errors['preparation_integral_'+label+'_'+parameter] = abs(
                sum(t[label+'_'+parameter] for t in segment_terms[:3])-preparation)

    psi, wave_dF, wave_dd = direct_wave_sensitivity(case)
    wave_B5 = 2*np.imag(np.outer(psi, psi.conj()))
    errors['wave_state'] = float(np.max(abs(wave_B5-B5)))
    for j, parameter in enumerate(('E', 'eta')):
        errors['wave_F_derivative_'+parameter] = abs(
            wave_dF[j]-sensitivities['F_'+parameter]['total'])
        errors['wave_d_derivative_'+parameter] = abs(
            wave_dd[j]-sensitivities['d_'+parameter]['total'])
    for i in (3, 4, 5):
        expected_psi = reference[case['name']+'_psi'][i*800]
        expected_B = 2*np.imag(np.outer(expected_psi, expected_psi.conj()))
        errors['reference_B'+str(i)] = float(np.max(abs(
            bivector_matrix(boundary_rows[i][0])-expected_B)))
    for key, value in errors.items():
        threshold = 2e-7 if ('derivative' in key or 'partition' in key
                              or 'integral' in key) else 2e-9
        assert value < threshold, (case['name'], key, value)

    E, eta = case['E1'], case['eta']
    J4 = .3*B4[0, 3]+eta*B4[1, 4]
    curvature = -4*((.3**2+eta**2)*B4[0, 1]-.6*eta*B4[3, 4]
                     +E*(.3*B4[1, 3]-eta*B4[0, 4]))
    report = {**case, 'F5': float(endpoint_F), 'd5': float(B5[0, 1]),
        'P33_at_3tau': P33, 'first_stage_radius': float(radius),
        'first_stage_F_upper_bound': float(F_bound),
        'required_P33_at_4tau': required_P33,
        'future_F_weight_matching_upper_bound': float(future_F_bound),
        'unavoidable_weight_mismatch_loss': float(1-future_F_bound),
        'remaining_alignment_loss': float(future_F_bound-endpoint_F),
        'J_at_4tau': float(J4), 'd_velocity_at_4tau': float(-2*J4),
        'd_curvature_at_4tau': float(curvature),
        'B3': B3.tolist(), 'B4': B4.tolist(),
        'required_B3_for_perfect_F5': G3.tolist(),
        'required_B4_for_perfect_F5': G4.tolist(),
        'three_mode_boundary_geometry': geometry,
        'sensitivities': sensitivities, 'segment_sensitivity_terms': segment_terms,
        'errors': errors}
    plot_time, plot_d, plot_F = [], [], []
    for stage in (3, 4):
        B = bivector_matrix(forward[stage].sol(SAMPLE)[:10].T)
        keep = slice(None) if stage == 3 else slice(1, None)
        plot_time.append(stage+SAMPLE[keep])
        plot_d.append(B[keep, 0, 1])
        plot_F.append(np.array([fidelity_to(b, target5) for b in B[keep]]))
    return report, {'t': np.concatenate(plot_time), 'd': np.concatenate(plot_d),
                    'F': np.concatenate(plot_F)}


def finite_swaps(reports, reference):
    names = ('eta_low', 'high_fidelity', 'E_changed')
    lookup = {r['name']: r for r in reports}
    entries = []
    for incoming in names:
        B3 = np.array(lookup[incoming]['B3'])
        for operator in names:
            target = np.array(lookup[operator]['required_B3_for_perfect_F5'])
            F = float(fidelity_to(B3, target))
            d = float(-np.sum(B3*target)/2)
            case = lookup[operator]
            psi = reference[incoming+'_psi'][2400].copy()
            for stage in (3, 4):
                wave = integrate(lambda v, y: wave_rhs(stage, v, y, case), psi)
                psi = wave.y[:, -1]
            F_wave = float(abs(np.vdot(PSI_TARGET, psi))**2)
            d_wave = float(2*np.imag(psi[0]*psi[1].conjugate()))
            error = max(abs(F-F_wave), abs(d-d_wave))
            assert error < 2e-9, (incoming, operator, error)
            entries.append({'incoming_B3': incoming, 'two_stage_operator': operator,
                            'F5': F, 'd5': d, 'wave_error': error})
    return {'names': names, 'entries': entries}


def plot_results(reports, paths, swaps):
    names = swaps['names']
    labels = {'eta_low': 'E=.01, eta=.010',
              'high_fidelity': 'E=.01, eta=.023',
              'E_changed': 'E=.02, eta=.023'}
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), layout='constrained')
    colors = ('#2271b2', '#d55e00', '#009e73')
    for name, color in zip(names, colors):
        path = paths[name]
        axes[0].plot(path['t'], path['d'], color=color, label=labels[name])
        axes[1].plot(path['t'], path['F'], color=color, label=labels[name])
        case = next(r for r in reports if r['name'] == name)
        axes[1].plot([3, 4], [case['first_stage_F_upper_bound']]*2,
                     color=color, ls=':', lw=1.4)
    for ax, title, ylabel in zip(axes[:2],
            ('(a) Orientation through the second cycle',
             '(b) Fidelity and first-stage upper bounds'), ('d', 'F')):
        ax.axvline(4, color='.5', ls='--', lw=1)
        ax.set(xlabel=r'$t/\tau$', ylabel=ylabel, title=title, xlim=(3, 5))
        ax.grid(alpha=.2)
    axes[0].axhline(0, color='.5', lw=.8)
    axes[0].legend(fontsize=8, loc='lower left')
    axes[1].set_ylim(-.02, 1.02)
    F = np.array([[next(e['F5'] for e in swaps['entries']
                 if e['incoming_B3'] == row and e['two_stage_operator'] == col)
                 for col in names] for row in names])
    im = axes[2].imshow(F, vmin=0, vmax=1, cmap='viridis')
    for i in range(3):
        for j in range(3):
            label = f'{F[i,j]:.2e}' if 0 < F[i,j] < .001 else f'{F[i,j]:.3f}'
            axes[2].text(j, i, label, ha='center', va='center',
                         color='black' if F[i,j] > .55 else 'white')
    axes[2].set_xticks(range(3), ['low eta', 'high F', 'changed E'])
    axes[2].set_yticks(range(3), ['low eta', 'high F', 'changed E'])
    axes[2].set(xlabel='Parameters in stages 4 and 5', ylabel='Incoming state at 3tau',
                title='(c) Exact finite input/operator swaps')
    fig.colorbar(im, ax=axes[2], shrink=.8, label=r'$F(5\tau)$')
    fig.savefig(ROOT/'second_cycle_chain.png', dpi=180)
    fig.savefig(ROOT/'second_cycle_chain.pdf')
    plt.close(fig)


def main():
    reports, paths = [], {}
    with np.load(ROOT/'bottleneck.npz') as reference:
        for case in CASES:
            report, path = check_case(case, reference)
            reports.append(report)
            paths[case['name']] = path
            print(case['name'], 'F5=', report['F5'],
                  'F derivatives=', {k: v for k, v in report['sensitivities'].items()
                                     if k.startswith('F_')}, flush=True)
        swaps = finite_swaps(reports, reference)
    result = {'units': 'hbar=1; energies in meV; derivatives per meV',
              'cases': reports, 'finite_swaps_at_3tau': swaps}
    (ROOT/'second_cycle_chain.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    plot_results(reports, paths, swaps)
    print('Maximum recorded check error:',
          max(max(r['errors'].values()) for r in reports))
    print('Maximum finite-swap wave error:',
          max(e['wave_error'] for e in swaps['entries']))


if __name__ == '__main__':
    main()
