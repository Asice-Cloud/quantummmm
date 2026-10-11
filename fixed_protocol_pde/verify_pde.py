#!/usr/bin/env python3
"""Check the exact fixed-protocol parameter/state PDE against wave derivatives.

This numerical example matches the existing parameter work: tau=100, tc=E0=.3,
the same cosine gates, and local max_step=.025. The derivation retains tau, tc,
and E0 as controlled design symbols; these values only select a check instance.
Run: .venv/bin/python -B fixed_protocol_pde/verify_pde.py
"""

import json
import os
from pathlib import Path

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import numpy as np
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parent
TAU, TC, E0, HBAR = 100., .3, .3, 1.
C = 2*TAU/HBAR
PSI0 = np.array([1., -1j, 0., 0., 0.])/np.sqrt(2)
TARGET = PSI0.conj()


def L(i, j):
    result = np.zeros((5, 5))
    result[i, j], result[j, i] = 1., -1.
    return result


B0 = L(0, 1)
BSTAR, PSTAR = -B0, -B0@B0


def controls(stage, v):
    a = (1-np.cos(np.pi*v))/2
    b = 1-a
    if stage % 3 == 0:
        return TC*a, 0., E0*b, a
    if stage % 3 == 1:
        return TC*b, TC*a, 0., b
    return 0., TC*b, E0*a, 0.


def generator(stage, v, E, eta):
    p, r, ed, gate = controls(stage, v)
    K = C*(E*L(0, 1)+eta*gate*L(0, 4)-p*L(1, 3)+r*L(2, 3)+ed*L(3, 4))
    return K, C*L(0, 1), C*gate*L(0, 4)


def comm(a, b):
    return a@b-b@a


def field_rhs(stage, v, flat, E, eta):
    B, XE, Xeta, SE, Seta, HEE, HEeta, Hetaeta, TEeta, TetaE = flat.reshape(10, 5, 5)
    K, KE, Keta = generator(stage, v, E, eta)
    return np.stack((
        comm(K, B),
        comm(K, XE)+KE,
        comm(K, Xeta)+Keta,
        comm(K, SE)+comm(KE, B),
        comm(K, Seta)+comm(Keta, B),
        comm(K, HEE)+2*comm(KE, SE),
        comm(K, HEeta)+comm(KE, Seta)+comm(Keta, SE),
        comm(K, Hetaeta)+2*comm(Keta, Seta),
        comm(K, TEeta)+comm(KE, Xeta),
        comm(K, TetaE)+comm(Keta, XE),
    )).ravel()


def wave_rhs(stage, v, flat, E, eta):
    """Independent component equations for psi and its first/second derivatives."""
    psi, pe, pt, pee, pet, ptt = flat.reshape(6, 5)
    p, r, ed, gate = controls(stage, v)
    q = eta*gate
    def action(z):
        return C*np.array([E*z[1]+q*z[4], -E*z[0]-p*z[3], r*z[3],
                           p*z[1]-r*z[2]+ed*z[4], -q*z[0]-ed*z[3]])
    def action_e(z):
        return C*np.array([z[1], -z[0], 0., 0., 0.])
    def action_eta(z):
        return C*gate*np.array([z[4], 0., 0., 0., -z[0]])
    return np.stack((action(psi), action(pe)+action_e(psi),
                     action(pt)+action_eta(psi), action(pee)+2*action_e(pe),
                     action(pet)+action_e(pt)+action_eta(pe),
                     action(ptt)+2*action_eta(pt))).ravel()


def outer(a, b):
    return a[..., :, None]*b.conj()[..., None, :]


def inner(a, b):
    return np.sum(a*b, axis=(-2, -1))


def verify(E, eta):
    fields, wave = np.zeros((10, 5, 5)), np.zeros((6, 5), complex)
    fields[0], wave[0] = B0, PSI0
    errors = {key: 0. for key in (
        'B_wave', 'parameter_PDE_absolute', 'parameter_PDE_relative',
        'first_derivatives_wave_relative', 'second_derivatives_wave_relative',
        'connection_curvature_absolute', 'connection_curvature_relative',
        'mixed_parameter_PDE_relative', 'F_wave', 'F_gradient_wave_relative',
        'F_hessian_wave_relative', 'metric_wave_relative',
        'state_constraint', 'wave_normalization', 'metric_PSD_violation', 'interface_K',
    )}
    final = None
    for stage in range(6):
        grid = np.linspace(0., 1., 101)
        fs = solve_ivp(lambda v, y: field_rhs(stage, v, y, E, eta), (0., 1.), fields.ravel(),
                       method='DOP853', rtol=2e-11, atol=2e-13, max_step=.025, t_eval=grid)
        ws = solve_ivp(lambda v, y: wave_rhs(stage, v, y, E, eta), (0., 1.), wave.ravel(),
                       method='DOP853', rtol=8e-13, atol=8e-15, max_step=.025, t_eval=grid)
        assert fs.success and ws.success
        B, XE, Xt, SE, St, HEE, HEt, Htt, TEt, TtE = np.moveaxis(fs.y.T.reshape(-1, 10, 5, 5), 1, 0)
        psi, pe, pt, pee, pet, ptt = np.moveaxis(ws.y.T.reshape(-1, 6, 5), 1, 0)
        rho = outer(psi, psi)
        rE, rt = outer(pe, psi)+outer(psi, pe), outer(pt, psi)+outer(psi, pt)
        rEE = outer(pee, psi)+outer(psi, pee)+2*outer(pe, pe)
        rEt = outer(pet, psi)+outer(psi, pet)+outer(pe, pt)+outer(pt, pe)
        rtt = outer(ptt, psi)+outer(psi, ptt)+2*outer(pt, pt)
        S = np.stack((SE, St), axis=1)
        H = np.stack((HEE, HEt, Htt), axis=1)
        wave_S = 2*np.imag(np.stack((rE, rt), axis=1))
        wave_H = 2*np.imag(np.stack((rEE, rEt, rtt), axis=1))
        connection_S = np.stack((comm(XE, B), comm(Xt, B)), axis=1)
        curvature = TEt-TtE-comm(XE, Xt)
        mixed_from_connection = comm(TEt, B)+comm(Xt, SE)
        P = -B@B
        dP = -(S@B[:, None]+B[:, None]@S)
        ddP = np.stack((-(HEE@B+B@HEE+2*SE@SE),
                        -(HEt@B+B@HEt+SE@St+St@SE),
                        -(Htt@B+B@Htt+2*St@St)), axis=1)
        F = (np.trace(PSTAR@P, axis1=-2, axis2=-1)-np.trace(BSTAR@B, axis1=-2, axis2=-1))/4
        grad = (np.trace(PSTAR@dP, axis1=-2, axis2=-1)-np.trace(BSTAR@S, axis1=-2, axis2=-1))/4
        hess = (np.trace(PSTAR@ddP, axis1=-2, axis2=-1)-np.trace(BSTAR@H, axis1=-2, axis2=-1))/4
        amplitude = psi@TARGET.conj()
        aE, at = pe@TARGET.conj(), pt@TARGET.conj()
        wave_grad = 2*np.real(amplitude.conj()[:, None]*np.column_stack((aE, at)))
        wave_hess = 2*np.real(amplitude.conj()[:, None]*np.column_stack((pee@TARGET.conj(), pet@TARGET.conj(), ptt@TARGET.conj()))
                             +np.column_stack((aE.conj()*aE, aE.conj()*at, at.conj()*at)))
        metric = (inner(dP[:, :, None], dP[:, None, :])+inner(S[:, :, None], S[:, None, :]))/8
        wave_der = np.stack((pe, pt), axis=1)
        gram = np.einsum('nli,nmi->nlm', wave_der.conj(), wave_der)
        overlap = np.einsum('nli,ni->nl', wave_der.conj(), psi)
        wave_metric = np.real(gram-overlap[:, :, None]*overlap.conj()[:, None, :])
        values = {
            'B_wave': np.max(abs(B-2*np.imag(rho))),
            'parameter_PDE_absolute': np.max(abs(S-connection_S)),
            'parameter_PDE_relative': np.max(abs(S-connection_S)/(1+abs(S)+abs(connection_S))),
            'first_derivatives_wave_relative': np.max(abs(S-wave_S)/(1+abs(S)+abs(wave_S))),
            'second_derivatives_wave_relative': np.max(abs(H-wave_H)/(1+abs(H)+abs(wave_H))),
            'connection_curvature_absolute': np.max(abs(curvature)),
            'connection_curvature_relative': np.max(abs(curvature)/(1+abs(TEt)+abs(TtE)+abs(comm(XE, Xt)))),
            'mixed_parameter_PDE_relative': np.max(abs(HEt-mixed_from_connection)/(1+abs(HEt)+abs(mixed_from_connection))),
            'F_wave': np.max(abs(F-abs(amplitude)**2)),
            'F_gradient_wave_relative': np.max(abs(grad-wave_grad)/(1+abs(grad)+abs(wave_grad))),
            'F_hessian_wave_relative': np.max(abs(hess-wave_hess)/(1+abs(hess)+abs(wave_hess))),
            'metric_wave_relative': np.max(abs(metric-wave_metric)/(1+abs(metric)+abs(wave_metric))),
            'state_constraint': np.max(abs(B@B@B+B)),
            'wave_normalization': np.max(abs(np.sum(abs(psi)**2, axis=1)-1)),
            'metric_PSD_violation': max(0., -np.min(np.linalg.eigvalsh(metric))),
        }
        if stage < 5:
            for left, right in zip(generator(stage, 1., E, eta), generator(stage+1, 0., E, eta)):
                values['interface_K'] = max(values.get('interface_K', 0.), float(np.max(abs(left-right))))
        for key, value in values.items():
            errors[key] = max(errors[key], float(value))
        fields, wave = fs.y[:, -1].reshape(10, 5, 5), ws.y[:, -1].reshape(6, 5)
        final = {'F': float(F[-1]), 'F_gradient': grad[-1].tolist(), 'F_hessian_EE_Eeta_etaeta': hess[-1].tolist(),
                 'state_parameter_metric': metric[-1].tolist()}
    normalized = [value for key, value in errors.items() if key.endswith('_relative')]
    assert max(normalized) < 2e-7, errors
    assert max(errors[key] for key in ('B_wave', 'F_wave', 'state_constraint', 'wave_normalization', 'interface_K')) < 2e-9, errors
    assert errors['metric_PSD_violation'] < 1e-8, errors
    return {'E1': E, 'eta': eta, 'errors_over_full_six_segments': errors, 'endpoint': final}


def main():
    rows = [verify(E, eta) for E, eta in ((0., 0.), (.01, .023), (.027, .001))]
    report = {'tau': TAU, 'tc': TC, 'E0': E0, 'hbar': HBAR, 'end_time': 6*TAU,
              'local_max_step': .025, 'method': 'DOP853',
              'note': 'Sample verification of exact PDE identities; not a parameter-space error bound.',
              'samples': rows}
    (ROOT/'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    for row in rows:
        print(json.dumps(row, indent=2), flush=True)


if __name__ == '__main__':
    main()
