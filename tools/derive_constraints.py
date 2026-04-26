#!/usr/bin/env python3
"""Derive symbolic algebraic constraints for XXZ, XYZ (mapping), and Free‑fermion condition.

Produces:
- XXZ: elimination of eta giving polynomial relation between A=c_xx=c_yy and B=c_zz.
- Free‑fermion: attempts to express H_spin (Pauli coefficients) as a quadratic fermion Hamiltonian
  and returns the algebraic consistency conditions (matrix equations) among c_{αβ}.

Run:
  python3 tools/derive_constraints.py
"""
from __future__ import annotations
import sympy as sp
from sympy import I


def xxz_constraint():
    eta = sp.symbols('eta')
    A, B = sp.symbols('A B')
    # A = 1/(2 sin eta), B = cot eta / 2
    poly = sp.simplify(4*A**2 - 4*B**2 - 1)
    return sp.Eq(poly, 0), poly


def free_fermion_constraints():
    # Define Pauli matrices
    sx = sp.Matrix([[0,1],[1,0]])
    sy = sp.Matrix([[0,-I],[I,0]])
    sz = sp.Matrix([[1,0],[0,-1]])
    id2 = sp.eye(2)

    # Pauli labels
    labels = ['0','x','y','z']
    paulis = {'0': id2, 'x': sx, 'y': sy, 'z': sz}

    # symbolic Pauli coefficients c_ab (explicit real and imag parts to avoid re() in Groebner)
    c_re = {}
    c_im = {}
    for a in ['0','x','y','z']:
        for b in ['0','x','y','z']:
            name = f'c_{a}{b}'
            cre, cim = sp.symbols(name + '_re ' + name + '_im', real=True)
            c_re[a+b] = cre
            c_im[a+b] = cim

    # build H_spin using c = cre + I*cim
    H = sp.zeros(4,4)
    for k in c_re:
        cre = c_re[k]
        cim = c_im[k]
        sym = cre + I*cim
        a, b = k[0], k[1]
        H += sym * sp.kronecker_product(paulis[a], paulis[b])

    # Build fermionic operators via JW: c1 = sigma^- ⊗ I, c2 = sigma^z ⊗ sigma^-
    s_minus = (sx - I*sy)/2
    c1 = sp.kronecker_product(s_minus, id2)
    c2 = sp.kronecker_product(sz, s_minus)
    c1d = c1.H
    c2d = c2.H

    # parameters for quadratic fermion Hamiltonian
    mu1, mu2 = sp.symbols('mu1 mu2', real=True)
    t_re, t_im = sp.symbols('t_re t_im', real=True)
    delta_re, delta_im = sp.symbols('delta_re delta_im', real=True)
    const = sp.symbols('const', real=True)

    t = t_re + I*t_im
    delta = delta_re + I*delta_im

    # build quadratic H = sum_i mu_i c_i^† c_i + (t c1^† c2 + h.c.) + 1/2(delta c1^† c2^† + h.c.) + const*I
    Hq = sp.zeros(4,4)
    Hq += mu1 * (c1d*c1)
    Hq += mu2 * (c2d*c2)
    Hq += t * (c1d*c2) + sp.conjugate(t) * (c2d*c1)
    Hq += sp.Rational(1,2) * (delta * (c1d*c2d) + sp.conjugate(delta) * (c2*c1))
    Hq += const * sp.eye(4)

    # Form linear equations by equating matrix elements H - Hq = 0
    M = sp.simplify(H - Hq)
    eqs = []
    for i in range(4):
        for j in range(4):
            eqs.append(sp.Eq(sp.simplify(M[i,j]), 0))

    # Unknowns in Hq
    unknowns = [mu1, mu2, t_re, t_im, delta_re, delta_im, const]

    # Convert to real linear system by separating real and imag parts
    re_eqs = []
    for e in eqs:
        expr = sp.simplify(sp.expand_complex(e.lhs))
        re_eqs.append(sp.Eq(sp.re(expr), 0))
        re_eqs.append(sp.Eq(sp.im(expr), 0))

    # Try linear solve for unknowns; coefficients are linear in unknowns
    sol = sp.solve(re_eqs, unknowns, dict=True)
    if sol:
        return {'sol_exists': True, 'solution': sol}
    else:
        # construct polynomial system and attempt elimination of unknowns via Groebner basis
        polys = [sp.simplify(e.lhs) for e in re_eqs]
        # list of all c real/imag symbols
        c_syms = list(c_re.values()) + list(c_im.values())
        try:
            G = sp.groebner(polys, *unknowns, *c_syms, order='lex')
            eliminated = G.eliminate(unknowns)
            # eliminated is a GroebnerBasis or list of polys in c_syms only
            eliminated_polys = list(eliminated)
            return {'sol_exists': False, 'groebner_eliminated': eliminated_polys}
        except Exception as e:
            import traceback
            print('  Groebner elimination raised exception:', e)
            traceback.print_exc()
            # Try alternative variable orderings for Groebner elimination
            try_orders = [
                tuple(list(unknowns) + list(c_syms)),
                tuple(list(c_syms) + list(unknowns)),
            ]
            eliminated_polys = None
            for vars_order in try_orders:
                try:
                    G = sp.groebner(polys, *vars_order, order='lex')
                    # Some SymPy versions lack GroebnerBasis.eliminate; extract basis
                    basis_polys = list(G)
                    # keep only polys that do not contain any of the unknowns
                    eliminated_polys = []
                    for p in basis_polys:
                        if not any(var in p.free_symbols for var in unknowns):
                            eliminated_polys.append(sp.simplify(p))
                    if eliminated_polys:
                        break
                except Exception as e2:
                    print('   groebner with order failed:', e2)

            if eliminated_polys:
                outp = '\n'.join([str(p) for p in eliminated_polys])
                outpath = 'results/free_fermion_constraints.txt'
                try:
                    with open(outpath, 'w') as f:
                        f.write(outp)
                    print('  Wrote eliminated constraints to', outpath)
                except Exception:
                    print('  Failed to write eliminated constraints to', outpath)
                return {'sol_exists': False, 'groebner_eliminated': eliminated_polys, 'outpath': outpath}
            else:
                A, b = sp.linear_eq_to_matrix([e.lhs for e in re_eqs], unknowns)
                rankA = A.rank()
                Ab = A.row_join(b)
                rankAb = Ab.rank()
                constraints = ['groebner elimination failed: ' + str(e)]
                return {'sol_exists': False, 'rankA': rankA, 'rankAb': rankAb, 'constraints': constraints}


def main():
    print('XXZ constraint (eliminated eta):')
    eq, poly = xxz_constraint()
    print('  Algebraic relation:', sp.pretty(eq))
    print('\nFree-fermion attempt: solving H_spin = quadratic fermion H')
    res = free_fermion_constraints()
    if res.get('sol_exists'):
        print('  Solution exists; quadratic parameters in terms of c_..:')
        print(res['solution'])
    else:
        print('  No linear solution found; rank info:')
        print('   rankA=', res.get('rankA'), ' rankAb=', res.get('rankAb'))
        print('   Further polynomial constraints required (see code).')


if __name__ == '__main__':
    main()
