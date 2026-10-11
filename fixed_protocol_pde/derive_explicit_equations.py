#!/usr/bin/env python3
"""Derive and check the explicit parameter-response coefficients algebraically.

Run: .venv/bin/python -B fixed_protocol_pde/derive_explicit_equations.py
No time evolution, parameter fit, or auxiliary connection field is needed.
"""

import json
from itertools import combinations
from pathlib import Path

import sympy as sp


def main():
    E, eta, p, r, ed, g = sp.symbols('E eta p r ed g', real=True)
    tau = sp.symbols('tau', positive=True)
    state = sp.Matrix(sp.symbols('d x1 x2 x3 y1 y2 y3 h1 h2 h3', real=True))
    pairs = list(combinations(range(5), 2))
    B = sp.zeros(5)
    for value, (i, j) in zip(state, pairs):
        B[i, j], B[j, i] = value, -value
    A = sp.zeros(5)
    for i, j, value in ((0, 1, E), (0, 4, eta*g), (1, 3, -p),
                         (2, 3, r), (3, 4, ed)):
        A[i, j], A[j, i] = value, -value
    commutator = A*B-B*A
    rhs = sp.Matrix([commutator[i, j] for i, j in pairs])
    expected = sp.Matrix([
        [0, 0, -p, 0, 0, 0, -eta*g, 0, 0, 0],
        [0, 0, r, 0, E, 0, 0, 0, -eta*g, 0],
        [p, -r, 0, ed, 0, E, 0, 0, 0, -eta*g],
        [0, 0, -ed, 0, 0, 0, E, 0, 0, 0],
        [0, -E, 0, 0, 0, r, 0, p, 0, 0],
        [0, 0, -E, 0, -r, 0, ed, 0, 0, 0],
        [eta*g, 0, 0, -E, 0, -ed, 0, 0, 0, -p],
        [0, 0, 0, 0, -p, 0, 0, 0, ed, 0],
        [0, eta*g, 0, 0, 0, 0, 0, -ed, 0, r],
        [0, 0, eta*g, 0, 0, 0, p, 0, -r, 0],
    ])
    N = rhs.jacobian(state)
    assert sp.simplify(N-expected) == sp.zeros(10)
    assert N+N.T == sp.zeros(10)
    d, x1, x2, x3, y1, y2, y3, h1, h2, h3 = state
    E_source = sp.Matrix([0, y1, y2, y3, -x1, -x2, -x3, 0, 0, 0])
    eta_source = g*sp.Matrix([-y3, -h2, -h3, 0, 0, 0, d, 0, x1, x2])
    assert N.diff(E)*state == E_source
    assert N.diff(eta)*state == eta_source
    assert N.diff(E, 2) == N.diff(E, eta) == N.diff(eta, 2) == sp.zeros(10)
    local_generator = 2*tau*N
    assert local_generator.diff(E)*state == 2*tau*E_source
    assert local_generator.diff(eta)*state == 2*tau*eta_source
    fidelity = (d*d-d)/2+(x1*x1+x2*x2+x3*x3+y1*y1+y2*y2+y3*y3)/4
    first_channel = (1-d)*x2+h1*y1-h3*y3
    eta_channel = (1-d)*y3-h3*x2-h2*x1
    fidelity_rate = p*first_channel+eta*g*eta_channel
    # Physical time t has factor 2; local v=s/tau has factor 2*tau.
    assert sp.expand((sp.Matrix([fidelity]).jacobian(state)*(2*N*state))[0]-fidelity_rate) == 0
    assert sp.expand((sp.Matrix([fidelity]).jacobian(state)*(local_generator*state))[0]-tau*fidelity_rate) == 0
    assert fidelity_rate.diff(E) == 0
    assert sp.expand(fidelity_rate.diff(eta)-g*eta_channel) == 0
    report = {
        'tau': str(tau), 'controlled_design_parameters': ['tau', 'tc', 'E0'],
        'local_time': 'v=s/tau', 'end_time': '6*tau',
        'time_generator_factor': str(2*tau), 'physical_time_generator_factor': 2,
        'state_order': [str(s) for s in state],
        'N': [[str(N[i, j]) for j in range(10)] for i in range(10)],
        'E_source_before_factor_2_tau': [str(x) for x in E_source],
        'eta_source_before_factor_2_tau': [str(x) for x in eta_source],
        'F': str(fidelity), 'physical_time_F_rate': str(fidelity_rate),
        'physical_time_F_rate_explicit_eta_source': str(g*eta_channel),
        'checks': {'matrix_from_full_commutator': True, 'skew_symmetric': True,
                   'E_source': True, 'eta_source': True, 'zero_second_parameter_derivatives': True,
                   'symbolic_tau_parameter_sources': True, 'normalized_time_fidelity_rate': True,
                   'fidelity_rate': True, 'explicit_fidelity_E_source_zero': True,
                   'explicit_fidelity_eta_source': True},
        'scope': 'Exact mixed time/parameter response equations; not a closed parameter-only PDE for endpoint fidelity.',
    }
    Path(__file__).with_name('explicit_equations.json').write_text(json.dumps(report, indent=2)+'\n')
    print('Exact symbolic checks passed with symbolic tau for the state equations, parameter sources, and fidelity-rate coupling.')


if __name__ == '__main__':
    main()
