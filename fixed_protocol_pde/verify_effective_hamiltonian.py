#!/usr/bin/env python3
"""Audit the literal effective Hamiltonian through Clifford and Fock matrices.

The spinor and Nambu Hamiltonians are built from the original operator products,
not from the existing five-dimensional K. Exact algebra determines K first;
independent time evolutions then check the bivector and parameter responses.
"""

import json
from pathlib import Path

from analyze_peak_structure import ATOL, MAX_STEP, RTOL, controls, readout, state_action
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parent
KEYS = ('E', 'q', 'p', 'r', 'ed')
ACTIVE = (0, 1, 2, 4, 5)  # Six Majoranas: (1, 2, 3, spectator, a, b).
PAIRS = np.triu_indices(5, 1)


def zero(matrix):
    assert matrix.applyfunc(sp.simplify) == sp.zeros(*matrix.shape), matrix


def comm(a, b):
    return a*b-b*a


def block(a, b, c, d):
    return a.row_join(b).col_join(c.row_join(d))


def exact_construction():
    I = sp.I
    eye, nil = sp.eye(2), sp.zeros(2)
    sx, sy, sz = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.diag(1, -1)
    qi, qj, qk = -I*sx, -I*sy, -I*sz
    zero(qi*qj-qk)
    gamma = [block(nil, eye, eye, nil), block(nil, -qi, qi, nil),
             block(nil, -qj, qj, nil), block(nil, -qk, qk, nil),
             block(eye, nil, nil, -eye)]
    for i in range(5):
        zero(gamma[i]-gamma[i].H)
        for j in range(5):
            zero(gamma[i]*gamma[j]+gamma[j]*gamma[i]-2*int(i == j)*sp.eye(4))
    sigma = {(i, j): comm(gamma[i], gamma[j])/4
             for i in range(5) for j in range(i+1, 5)}
    printed_sigma = {(0, 1): block(qi, nil, nil, -qi)/2,
                     (1, 3): block(qj, nil, nil, qj)/2,
                     (0, 4): block(nil, -eye, eye, nil)/2,
                     (2, 3): block(-qi, nil, nil, -qi)/2,
                     (3, 4): block(nil, qk, qk, nil)/2}
    for indices, matrix in printed_sigma.items():
        zero(sigma[indices]-matrix)

    # Literal products in model.md; physical ordering is intentionally retained.
    H_spin = {'E': I*gamma[0]*gamma[1], 'q': -I*gamma[4]*gamma[0],
              'p': I*gamma[3]*gamma[1], 'r': -I*gamma[3]*gamma[2],
              'ed': I*gamma[3]*gamma[4]}
    canonical = {'E': (0, 1, 1), 'q': (0, 4, 1), 'p': (1, 3, -1),
                 'r': (2, 3, 1), 'ed': (3, 4, 1)}
    A, K_spin = {}, {}
    for key in KEYS:
        i, j, sign = canonical[key]
        A[key] = sp.zeros(5)
        A[key][i, j], A[key][j, i] = sign, -sign
        zero(H_spin[key]-2*I*sign*sigma[i, j])
        K_spin[key] = sp.Matrix(5, 5, lambda row, col:
            sp.trace(gamma[col]*I*comm(H_spin[key], gamma[row]))/4)
        zero(K_spin[key]-2*A[key])

    # Independent physical Fock construction: three ordinary complex fermions.
    lowering = sp.Matrix([[0, 1], [0, 0]])
    fermions = []
    for index in range(3):
        fermions.append(sp.kronecker_product(*[
            sz if k < index else lowering if k == index else eye for k in range(3)]))
    for i in range(3):
        for j in range(3):
            zero(fermions[i]*fermions[j]+fermions[j]*fermions[i])
            zero(fermions[i]*fermions[j].H+fermions[j].H*fermions[i]-int(i == j)*sp.eye(8))
    majorana = []
    for c in fermions:
        majorana.extend((c+c.H, -I*(c-c.H)))
    g1, g2, g3, spectator, ga, gb = majorana
    dot = fermions[2]
    # Paper Eq. (1), phi1=0, phi2=-pi, hence t1=q, t2=p, t3=-i*r.
    H_fock = {'E': I*g1*g2, 'q': (dot.H-dot)*g1,
              'p': I*(dot+dot.H)*g2, 'r': -I*(dot+dot.H)*g3,
              'ed': 2*dot.H*dot}
    H_majorana = {'E': I*g1*g2, 'q': -I*gb*g1, 'p': I*ga*g2,
                  'r': -I*ga*g3, 'ed': I*ga*gb}
    nambu = fermions+[c.H for c in fermions]
    # gamma/sqrt(2) = T * (c1,c2,d,c1^dagger,c2^dagger,d^dagger)^T.
    T = sp.zeros(6)
    for j in range(3):
        T[2*j, j] = T[2*j, 3+j] = 1/sp.sqrt(2)
        T[2*j+1, j], T[2*j+1, 3+j] = -I/sp.sqrt(2), I/sp.sqrt(2)
    zero(T*T.H-sp.eye(6))
    H_bdg = {}
    for key in KEYS:
        zero(H_fock[key]-H_majorana[key]-(sp.eye(8) if key == 'ed' else sp.zeros(8)))
        zero(H_fock[key]-H_fock[key].H)
        H_bdg[key] = sp.Matrix(6, 6, lambda row, col:
            sp.trace(nambu[col].H*comm(nambu[row], H_fock[key]))/4)
        for row in range(6):
            zero(comm(nambu[row], H_fock[key])-sum(
                (H_bdg[key][row, col]*nambu[col] for col in range(6)), sp.zeros(8)))
        zero(H_bdg[key]-H_bdg[key].H)
        K6 = sp.Matrix(6, 6, lambda row, col:
            sp.trace(majorana[col]*I*comm(H_fock[key], majorana[row]))/8)
        zero(T*H_bdg[key]*T.H-I*K6)
        zero(K6.extract(ACTIVE, ACTIVE)-K_spin[key])
        zero(K6[3, :])
        zero(K6[:, 3])

    # Check all ten coefficients of the actual PDE implementation, per basis.
    b_symbols = sp.symbols('d x1 x2 x3 y1 y2 y3 h1 h2 h3', real=True)
    B = sp.zeros(5)
    for value, (i, j) in zip(b_symbols, zip(*PAIRS)):
        B[int(i), int(j)], B[int(j), int(i)] = value, -value
    numerical_basis_errors = {}
    for key in KEYS:
        values = {name: int(name == key) for name in KEYS}
        rhs = comm(K_spin[key], B)
        coefficients = sp.Matrix([rhs[int(i), int(j)] for i, j in zip(*PAIRS)]).jacobian(b_symbols)
        implemented = 2*state_action(np.eye(10), values['E'], values['q'],
                                     values['p'], values['r'], values['ed'], 1.).T
        error = float(np.max(abs(np.array(coefficients, dtype=float)-implemented)))
        assert error == 0.
        numerical_basis_errors[key] = error

    E, q, p, r, ed = sp.symbols('E q p r ed', real=True)
    symbols = dict(zip(KEYS, (E, q, p, r, ed)))
    combined_bdg = sum((symbols[key]*H_bdg[key] for key in KEYS), sp.zeros(6))
    combined_K = sum((symbols[key]*K_spin[key] for key in KEYS), sp.zeros(5))
    report = {
        'checks_exact': {
            'quaternion_and_clifford_relations': True,
            'five_printed_Sigma_matrices': True,
            'literal_H_equals_2i_times_signed_Sigma': True,
            'spinor_commutator_gives_K_equals_2A': True,
            'fermion_anticommutation_relations': True,
            'paper_Eq1_equals_Eq2_plus_Ed_identity': True,
            'Nambu_commutator_closes_exactly': True,
            'BdG_to_Majorana_gives_2iA': True,
            'sixth_Majorana_decouples': True,
        },
        'PDE_basis_coefficient_errors': numerical_basis_errors,
        'K': [[str(x) for x in combined_K.row(i)] for i in range(5)],
        'Nambu_H': [[str(x) for x in combined_bdg.row(i)] for i in range(6)],
        'canonical_H': 'E X12 + q X15 - p X24 + r X34 + ed X45; Xij=i Gamma_i Gamma_j',
        'documentation_corrections': ['Sigma24 coupling label: -p',
                                      'Sigma15 coupling label: +q',
                                      'Sigma34 coupling label: +r'],
        'physical_generator_factor': 2,
    }
    convert = lambda matrices: np.array([np.array(matrices[key], dtype=complex) for key in KEYS])
    return report, np.array(gamma, dtype=complex), convert(H_spin), convert(H_bdg), np.array(T, dtype=complex)


def global_controls(zeta, tc=.3, e0=.3):
    """A separate global clock, without using the PDE's stage-control helper."""
    phase = zeta % 3.
    if phase < 1.:
        opening = .5*(1-np.cos(np.pi*phase))
        return tc*opening, 0., e0*(1-opening), opening
    if phase < 2.:
        opening = .5*(1-np.cos(np.pi*(phase-1.)))
        return tc*(1-opening), tc*opening, 0., 1-opening
    opening = .5*(1-np.cos(np.pi*(phase-2.)))
    return 0., tc*(1-opening), e0*opening, 0.


def reconstruct_spinor(U, gamma):
    operators = U[0].conj().T[None]@gamma@U[0][None]
    rotation = np.einsum('iab,jba->ij', operators, gamma).real/4
    derivatives = []
    for derivative in U[1:]:
        operators_derivative = (derivative.conj().T[None]@gamma@U[0][None]
                                +U[0].conj().T[None]@gamma@derivative[None])
        derivatives.append(np.einsum('iab,jba->ij', operators_derivative, gamma).real/4)
    B0 = np.zeros((5,5))
    B0[0,1], B0[1,0] = 1., -1.
    fields = [rotation@B0@rotation.T]
    fields.extend(derivative@B0@rotation.T+rotation@B0@derivative.T for derivative in derivatives)
    return np.array([matrix[PAIRS] for matrix in fields]), rotation


def parameter_fidelity(fields):
    b = fields[0]
    return (b[0]-.5)*fields[1:,0]+np.sum(b[None,1:7]*fields[1:,1:7], axis=-1)/2


def numerical_check(tau, E, eta, gamma, spin_basis, bdg_basis, T):
    checkpoints = np.linspace(0., 6., 121)
    # Complete bivector with E and eta responses, propagated through six stages.
    field = np.zeros((3,10))
    field[0,0] = 1.
    field_history = []
    for stage in range(6):
        def rhs(v, flat):
            z = flat.reshape(3,10)
            p, r, ed, gate = controls(stage, v, .3, .3)
            result = state_action(z, E, eta, p, r, ed, gate)
            result[1] += state_action(z[0], 1., 0., 0., 0., 0., 0.)
            result[2] += state_action(z[0], 0., 1., 0., 0., 0., gate)
            return (2*tau*result).ravel()
        solution = solve_ivp(rhs, (0.,1.), field.ravel(), method='DOP853',
                             rtol=RTOL, atol=ATOL, max_step=MAX_STEP, t_eval=np.linspace(0.,1.,21))
        assert solution.success, solution.message
        values = solution.y.T.reshape(-1,3,10)
        field_history.extend(values if stage == 0 else values[1:])
        field = values[-1]
    field_history = np.array(field_history)

    # One uninterrupted spinor + Nambu solve: Hamiltonians come from operator algebra.
    initial_spin = np.zeros((3,4,4), complex)
    initial_spin[0] = np.eye(4)
    initial_M = np.array([1., -1j, 0., 0., 0., 0.])/np.sqrt(2)
    initial_wave = np.zeros((3,6), complex)
    initial_wave[0] = T.conj().T@initial_M
    target_Nambu = T.conj().T@initial_M.conj()
    initial = np.concatenate((initial_spin.ravel(), initial_wave.ravel()))
    def rhs_independent(zeta, flat):
        U, wave = flat[:48].reshape(3,4,4), flat[48:].reshape(3,6)
        p, r, ed, gate = global_controls(zeta)
        weights = np.array([E, eta*gate, p, r, ed])
        Hs = np.einsum('k,kij->ij', weights, spin_basis)
        Hb = np.einsum('k,kij->ij', weights, bdg_basis)
        dU = Hs@U
        dw = np.einsum('ij,kj->ki', Hb, wave)
        dU[1] += spin_basis[0]@U[0]
        dU[2] += gate*spin_basis[1]@U[0]
        dw[1] += bdg_basis[0]@wave[0]
        dw[2] += gate*bdg_basis[1]@wave[0]
        return -1j*tau*np.concatenate((dU.ravel(), dw.ravel()))
    independent = solve_ivp(rhs_independent, (0.,6.), initial, method='DOP853',
                            rtol=RTOL, atol=ATOL, max_step=MAX_STEP, t_eval=checkpoints)
    assert independent.success, independent.message
    spin_history = independent.y[:48].T.reshape(-1,3,4,4)
    wave_history = independent.y[48:].T.reshape(-1,3,6)
    errors = {key: 0. for key in (
        'B_vs_spinor', 'B_vs_BdG', 'F_vs_spinor', 'F_vs_BdG',
        'parameter_response_vs_spinor_normalized', 'parameter_response_vs_BdG_normalized',
        'F_parameter_gradient_vs_spinor_normalized', 'F_parameter_gradient_vs_BdG_normalized',
        'spinor_unitarity', 'rotation_orthogonality', 'state_constraint', 'spectator_amplitude',
    )}
    endpoints = []
    def update(key, actual, expected, normalized=False):
        error = abs(actual-expected)
        if normalized:
            error = error/(1+abs(actual)+abs(expected))
        errors[key] = max(errors[key], float(np.max(error)))
    for index, (fields, U, wave) in enumerate(zip(field_history, spin_history, wave_history)):
        spin_fields, R = reconstruct_spinor(U, gamma)
        majorana_wave = wave@T.T
        active_wave = majorana_wave[:, ACTIVE]
        psi = active_wave[0]
        B = 2*np.imag(np.outer(psi, psi.conj()))
        bdg_fields = [B[PAIRS]]
        for response in active_wave[1:]:
            B_response = 2*np.imag(np.outer(response,psi.conj())+np.outer(psi,response.conj()))
            bdg_fields.append(B_response[PAIRS])
        bdg_fields = np.array(bdg_fields)
        update('B_vs_spinor', fields[0], spin_fields[0])
        update('B_vs_BdG', fields[0], bdg_fields[0])
        update('parameter_response_vs_spinor_normalized', fields[1:], spin_fields[1:], True)
        update('parameter_response_vs_BdG_normalized', fields[1:], bdg_fields[1:], True)
        amplitude = wave@target_Nambu.conj()
        F_bdg = abs(amplitude[0])**2
        F = readout(fields[0])['F']
        F_spin = readout(spin_fields[0])['F']
        update('F_vs_spinor', F, F_spin)
        update('F_vs_BdG', F, F_bdg)
        gradient_bdg = 2*np.real(amplitude[0].conj()*amplitude[1:])
        update('F_parameter_gradient_vs_spinor_normalized', parameter_fidelity(fields), parameter_fidelity(spin_fields), True)
        update('F_parameter_gradient_vs_BdG_normalized', parameter_fidelity(fields), gradient_bdg, True)
        update('spinor_unitarity', U[0].conj().T@U[0], np.eye(4))
        update('rotation_orthogonality', R.T@R, np.eye(5))
        update('state_constraint', B@B@B+B, np.zeros((5,5)))
        update('spectator_amplitude', majorana_wave[:,3], 0.)
        if index % 20 == 0:
            endpoints.append({'t_over_tau': float(checkpoints[index]),
                              'F_PDE': float(F), 'F_spinor_adjoint': float(F_spin),
                              'F_BdG_overlap': float(F_bdg)})
    assert max(value for key, value in errors.items() if not key.endswith('_normalized')) < 2e-8, errors
    assert max(value for key, value in errors.items() if key.endswith('_normalized')) < 2e-6, errors
    return {'tau': tau, 'E1': E, 'eta': eta, 'checkpoints': len(checkpoints),
            'errors_over_entire_protocol': errors, 'interfaces_and_endpoints': endpoints}


def main():
    report, gamma, spin, bdg, T = exact_construction()
    print('Exact Hamiltonian, Clifford, Fock/BdG, and PDE coefficient checks passed.', flush=True)
    samples = []
    for tau, E, eta in [(100.,0.,0.), (100.,.01,.001), (300.,.01,.01), (1000.,.01,.1)]:
        row = numerical_check(tau, E, eta, gamma, spin, bdg, T)
        samples.append(row)
        print(json.dumps({key: row[key] for key in ('tau','E1','eta','errors_over_entire_protocol')}), flush=True)
    report.update(
        numerical_samples=samples,
        solver={'method': 'DOP853', 'rtol': RTOL, 'atol': ATOL, 'max_normalized_step': MAX_STEP},
        scope='Exact coefficient identities for the stated effective Hamiltonian; numerical trajectory checks at selected parameters. Does not establish reproduction of the published figure.',
        conclusion='K=2A is fixed by the Hamiltonian and Clifford normalization; the PDE signs agree. Three coupling labels in model.md required correction, but the PDE implementation did not.',
    )
    (ROOT/'effective_hamiltonian_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print('Saved effective_hamiltonian_verification.json', flush=True)


if __name__ == '__main__':
    main()
