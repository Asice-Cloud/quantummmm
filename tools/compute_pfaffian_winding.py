#!/usr/bin/env python3
"""Compute winding (phase winding of q(k)) and sign test q(0)*q(pi) for candidates.

Saves results to results/abs_pfaffian_report.json
"""
import os
import json
import numpy as np

from verify_mzm import map_c_to_params


def q_of_k(t, Delta, mu, k):
    # q(k) = -mu - 2 t cos k + 2 i Delta sin k
    return -mu - 2.0 * t * np.cos(k) + 2j * Delta * np.sin(k)


def winding_number(t, Delta, mu, nk=1001):
    ks = np.linspace(0, 2*np.pi, nk)
    q = q_of_k(t, Delta, mu, ks)
    phases = np.unwrap(np.angle(q))
    dw = phases[-1] - phases[0]
    return float(np.round(dw / (2*np.pi)))


def sign_test(t, Delta, mu):
    q0 = q_of_k(t, Delta, mu, 0.0)
    qpi = q_of_k(t, Delta, mu, np.pi)
    val = np.real(q0) * np.real(qpi)
    return float(np.sign(val))


def load_candidates(path):
    with open(path, 'r') as f:
        return json.load(f)


def main(candidates_path='results/abs_scan_multi_candidates.json', top_k=10):
    os.makedirs('results', exist_ok=True)
    cand = load_candidates(candidates_path)
    if not cand:
        print('No candidates found in', candidates_path)
        return
    cand_sorted = sorted(cand, key=lambda x: x.get('min_abs', 1e9))
    top = cand_sorted[:top_k]
    report = []
    for i, c in enumerate(top, start=1):
        cxy = float(c['cxy'])
        cyx = float(c['cyx'])
        cxz = float(c.get('cxz', 0.0))
        # construct pauli mapping assumed baseline used in scans
        c_xx = 0.05
        c_yy = 0.05
        c_xy = complex(cxy)
        c_yx = complex(cyx)
        c_zz = 0.0
        c_z0 = 0.0
        c_0z = 0.0
        t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
        wnum = winding_number(t, Delta, mu)
        sgn = sign_test(t, Delta, mu)
        def cplx(z):
            return {'re': float(np.real(z)), 'im': float(np.imag(z))}
        report.append({'index': i, 'cxy': cxy, 'cyx': cyx, 'cxz': cxz, 't': cplx(t), 'Delta': cplx(Delta), 'mu': float(mu), 'U': float(U), 'winding': wnum, 'sign_test': sgn})

    out = 'results/abs_pfaffian_report.json'
    with open(out, 'w') as f:
        json.dump(report, f, indent=2)
    print('Pfaffian/winding report saved to', out)


if __name__ == '__main__':
    main()
