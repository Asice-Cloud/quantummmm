#!/usr/bin/env python3
"""Verify exact segment composition using only the oriented-plane state.

Stage 1 uses two complex channels plus the four mode-3 components; stage 2
uses its 4+6 block; stage 3 uses the exact SO(2) x SO(3) action. All six
segments inherit their complete incoming state. The independent reference
is the complex-wavefunction trajectory already saved in bottleneck.npz.

Run from the project root:
  .venv/bin/python -B geometric_validation/verify_full_protocol.py
"""

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp

from verify_bottleneck import B_INITIAL, CASES, controls
from verify_second_stage import (
    block_rhs,
    fields_from_matrix as second_fields,
    matrix_from_fields as second_matrix,
)


ROOT = Path(__file__).resolve().parent
GRID = np.linspace(0, 1, 801)


def first_fields(B):
    d, hab = B[0, 1], B[3, 4]
    za, zb = B[0, 3] + 1j*B[1, 3], B[0, 4] + 1j*B[1, 4]
    cp, cm = za + 1j*zb, za - 1j*zb
    return np.array([d+hab, cp.real, cp.imag, d-hab, cm.real, cm.imag,
                     B[0, 2], B[1, 2], B[2, 3], B[2, 4]])


def first_matrix(values):
    xp, cp = values[..., 0], values[..., 1] + 1j*values[..., 2]
    xm, cm = values[..., 3], values[..., 4] + 1j*values[..., 5]
    za, zb = (cp+cm)/2, (cp-cm)/(2j)
    B = np.zeros(values.shape[:-1] + (5, 5))
    B[..., 0, 1], B[..., 3, 4] = (xp+xm)/2, (xp-xm)/2
    B[..., 0, 3], B[..., 1, 3] = za.real, za.imag
    B[..., 0, 4], B[..., 1, 4] = zb.real, zb.imag
    B[..., 0, 2], B[..., 1, 2] = values[..., 6], values[..., 7]
    B[..., 2, 3], B[..., 2, 4] = values[..., 8], values[..., 9]
    return B - np.swapaxes(B, -1, -2)


def first_rhs(v, values, case):
    p, _, ed, gate = controls(0, v)
    E, q = case['E1'], case['eta']*gate
    derivative = np.zeros(10)
    for sign, start in ((1, 0), (-1, 3)):
        X, x, y = values[start:start+3]
        V, delta = p-sign*q, E+sign*ed
        derivative[start:start+3] = (-V*x, V*X+delta*y, -delta*x)
    x, u, h, j = values[6:]
    derivative[6:] = (E*u-q*j, -E*x+p*h, ed*j-p*u, -ed*h+q*x)
    return 2*case['tau']*derivative


def integrate(rhs, initial):
    solution = solve_ivp(rhs, (0, 1), initial, method='DOP853',
                         rtol=2e-12, atol=2e-14, t_eval=GRID, max_step=.01)
    assert solution.success, solution.message
    return solution.y.T


def third_rotations(tau):
    def rhs(v, values):
        _, r, ed, _ = controls(2, v)
        D = np.array([[0, r, 0], [-r, 0, ed], [0, -ed, 0]])
        return (2*tau*D@values.reshape(3, 3)).ravel()
    return integrate(rhs, np.eye(3).ravel()).reshape(-1, 3, 3)


def third_path(B, case, rotations):
    d = B[0, 1]
    z = B[0, 2:] + 1j*B[1, 2:]
    zs = (np.exp(-2j*case['E1']*case['tau']*GRID)[:, None]
          * np.einsum('i,tji->tj', z, rotations))
    hs = rotations@B[2:, 2:]@np.swapaxes(rotations, -1, -2)
    result = np.zeros((len(GRID), 5, 5))
    result[:, 0, 1], result[:, 1, 0] = d, -d
    result[:, 0, 2:], result[:, 1, 2:] = zs.real, zs.imag
    result[:, 2:, 0], result[:, 2:, 1] = -zs.real, -zs.imag
    result[:, 2:, 2:] = hs
    return result


def generator(stage, v, case):
    p, r, ed, gate = controls(stage, v)
    E, q = case['E1'], case['eta']*gate
    return 2*np.array([[0, E, 0, 0, q], [-E, 0, 0, -p, 0],
                       [0, 0, 0, r, 0], [0, p, -r, 0, ed],
                       [-q, 0, 0, -ed, 0]])


def generator_derivative(stage, v, case):
    gate_derivative = np.pi*np.sin(np.pi*v)/(2*case['tau'])
    if stage % 3 == 0:
        p, q, r, ed = (.3*gate_derivative, case['eta']*gate_derivative,
                       0., -.3*gate_derivative)
    elif stage % 3 == 1:
        p, q, r, ed = (-.3*gate_derivative, -case['eta']*gate_derivative,
                       .3*gate_derivative, 0.)
    else:
        p, q, r, ed = 0., 0., -.3*gate_derivative, .3*gate_derivative
    return 2*np.array([[0, 0, 0, 0, q], [0, 0, 0, -p, 0],
                       [0, 0, 0, r, 0], [0, p, -r, 0, ed],
                       [-q, 0, 0, -ed, 0]])


def commutator(K, B):
    return K@B - B@K


def fidelity(B):
    return (np.sum(B[..., :2, :]**2, axis=(-1, -2))-2*B[..., 0, 1])/4


def symbolic_checks():
    E, p, q, r, ed = sp.symbols('E p q r ed', real=True)
    d, x, y, z, u, v, w, h, j, k = sp.symbols(
        'd x y z u v w h j k', real=True)
    A = sp.Matrix([[0, E, 0, 0, q], [-E, 0, 0, -p, 0],
                   [0, 0, 0, r, 0], [0, p, -r, 0, ed],
                   [-q, 0, 0, -ed, 0]])
    B = sp.Matrix([[0, d, x, y, z], [-d, 0, u, v, w],
                   [-x, -u, 0, h, j], [-y, -v, -h, 0, k],
                   [-z, -w, -j, -k, 0]])
    DB = A*B - B*A
    W = sp.Matrix([x, u, h, j])
    DW = sp.Matrix([DB[0, 2], DB[1, 2], DB[2, 3], DB[2, 4]])
    M = sp.Matrix([[0, E, 0, -q], [-E, 0, p, 0],
                   [0, -p, 0, ed], [q, 0, -ed, 0]])
    assert sp.expand(DW.subs(r, 0)-M*W) == sp.zeros(4, 1)
    assert M+M.T == sp.zeros(4)
    plus_norm = (d+k)**2 + (y-w)**2 + (v+z)**2
    minus_norm = (d-k)**2 + (y+w)**2 + (v-z)**2
    assert sp.expand(plus_norm-minus_norm-4*(d*k-y*w+z*v)) == 0
    assert sp.expand((plus_norm+minus_norm)/2
                     + (W.T*W)[0] - sum(B[a, b]**2
                     for a in range(5) for b in range(a+1, 5))) == 0
    observables = (d*d-d)/2+(x*x+y*y+z*z+u*u+v*v+w*w)/4
    coords = (d, x, y, z, u, v, w, h, j, k)
    pairs = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2),
             (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
    DF = sum(sp.diff(observables, c)*DB[a, b]
             for c, (a, b) in zip(coords, pairs))
    expected = (p*((1-d)*y+h*u-k*w)
                + q*((1-d)*w-k*y-j*x))/2
    assert sp.expand(DF-expected) == 0
    J = sp.Matrix([[0, 1], [-1, 0]])
    D = sp.Matrix([[0, r, 0], [-r, 0, ed], [0, -ed, 0]])
    third = DB.subs({p: 0, q: 0})
    assert third[0, 1] == 0
    assert third[:2, 2:] == E*J*B[:2, 2:]-B[:2, 2:]*D
    assert third[2:, 2:] == D*B[2:, 2:]-B[2:, 2:]*D
    return {'mode3_subsystem': True, 'channel_radius_identities': True,
            'exact_fidelity_flux': True, 'third_stage_blocks': True}


def check_case(case, reference, rotations):
    current = B_INITIAL.copy()
    segments, first_invariants = [], []
    errors = {}
    for stage in range(6):
        if stage % 3 == 0:
            fields = integrate(lambda v, y: first_rhs(v, y, case),
                               first_fields(current))
            path = first_matrix(fields)
            plus = np.sum(fields[:, :3]**2, axis=1)
            minus = np.sum(fields[:, 3:6]**2, axis=1)
            mode3 = np.sum(fields[:, 6:]**2, axis=1)
            errors[f'stage_{stage+1}_channel_norms'] = float(max(
                np.max(abs(plus-plus[0])), np.max(abs(minus-minus[0])),
                np.max(abs(mode3-mode3[0]))))
            errors[f'stage_{stage+1}_channel_radius_identity'] = float(max(
                np.max(abs(plus-(1-mode3))),
                np.max(abs(minus-(1-mode3)))))
            first_invariants.append({'stage': stage+1,
                'P33': float(mode3[0]),
                'channel_radius': float(np.sqrt(plus[0]))})
        elif stage % 3 == 1:
            X, Y = second_fields(current)
            fields = integrate(lambda v, y: block_rhs(v, y, case),
                               np.concatenate((X, Y)))
            path = second_matrix(fields)
        else:
            path = third_path(current, case, rotations)
            errors[f'stage_{stage+1}_F_constant'] = float(
                np.max(abs(fidelity(path)-fidelity(current))))
        errors[f'stage_{stage+1}_incoming_state'] = float(
            np.max(abs(path[0]-current)))
        segments.append(path)
        current = path[-1]

    full = np.concatenate([segments[0]]+[s[1:] for s in segments[1:]])
    psi = reference[case['name']+'_psi']
    B_wave = 2*np.imag(psi[..., :, None]*psi.conj()[..., None, :])
    errors['full_wave_bivector'] = float(np.max(abs(full-B_wave)))
    errors['full_wave_fidelity'] = float(np.max(
        abs(fidelity(full)-reference[case['name']+'_F'])))
    errors['simple_plane_constraint'] = float(np.max(abs(full@full@full+full)))
    errors['norm_constraint'] = float(np.max(
        abs(np.sum(full**2, axis=(1, 2))-2)))
    for stage in range(5):
        B = segments[stage][-1]
        derivatives = []
        for st, v in ((stage, 1.), (stage+1, 0.)):
            K, dK = generator(st, v, case), generator_derivative(st, v, case)
            first = commutator(K, B)
            second = commutator(dK, B)+commutator(K, first)
            derivatives.append((first, second))
        errors[f'join_{stage+1}_first_derivative'] = float(np.max(
            abs(derivatives[0][0]-derivatives[1][0])))
        errors[f'join_{stage+1}_second_derivative'] = float(np.max(
            abs(derivatives[0][1]-derivatives[1][1])))
    assert max(errors.values()) < 2e-9, errors
    boundaries = []
    for j in range(7):
        B = full[j*800]
        P = -B@B
        boundaries.append({'t_over_tau': j, 'd': float(B[0, 1]),
            'ell': float(1-(P[0, 0]+P[1, 1])/2),
            'F': float(fidelity(B)), 'P33': float(P[2, 2])})
    third_motion = [{'stage': stage+1,
        'B_change_Frobenius': float(np.linalg.norm(
            segments[stage][-1]-segments[stage][0])),
        'F_in': float(fidelity(segments[stage][0])),
        'F_out': float(fidelity(segments[stage][-1]))}
        for stage in (2, 5)]
    return {**case, 'boundaries': boundaries,
            'first_stage_invariants': first_invariants,
            'third_stage_motion': third_motion, 'errors': errors}


def main():
    result = {'symbolic': symbolic_checks()}
    rotations = {tau: third_rotations(tau) for tau in {c['tau'] for c in CASES}}
    with np.load(ROOT/'bottleneck.npz') as reference:
        result['cases'] = [check_case(c, reference, rotations[c['tau']])
                           for c in CASES]
    (ROOT/'full_protocol.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    for case in result['cases']:
        print(case['name'], 'F_end=', case['boundaries'][-1]['F'],
              'max_error=', max(case['errors'].values()))
    print('Symbolic checks:', result['symbolic'])


if __name__ == '__main__':
    main()
