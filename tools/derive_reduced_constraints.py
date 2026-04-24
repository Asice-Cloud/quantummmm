#!/usr/bin/env python3
"""Derive algebraic constraints for reduced two-site case:
 c_xx = c_yy = A, c_zz = B, other c_ab = 0.

Compute:
 - XXZ elimination relation (should reproduce 4 A^2 - 4 B^2 - 1 = 0)
 - Free-fermion polynomial constraints by eliminating quadratic fermion parameters
"""
from __future__ import annotations
import sympy as sp
from sympy import I


def main():
    A, B = sp.symbols('A B')

    # Pauli matrices
    sx = sp.Matrix([[0,1],[1,0]])
    sy = sp.Matrix([[0,-I],[I,0]])
    sz = sp.Matrix([[1,0],[0,-1]])
    id2 = sp.eye(2)

    # Build H_spin for reduced coefficients
    H = A * sp.kronecker_product(sx, sx) + A * sp.kronecker_product(sy, sy) + B * sp.kronecker_product(sz, sz)

    print('H (reduced) = A xx + A yy + B zz')

    # XXZ relation
    poly_xxz = sp.simplify(4*A**2 - 4*B**2 - 1)
    print('\nXXZ elimination relation:')
    print(' 4*A^2 - 4*B^2 - 1 = 0  =>', sp.factor(poly_xxz))

    # Build JW fermions for two sites
    s_minus = (sx - I*sy)/2
    c1 = sp.kronecker_product(s_minus, id2)
    c2 = sp.kronecker_product(sz, s_minus)
    c1d = c1.H
    c2d = c2.H

    # quadratic fermion unknowns
    mu1, mu2 = sp.symbols('mu1 mu2')
    t_re, t_im = sp.symbols('t_re t_im')
    d_re, d_im = sp.symbols('d_re d_im')
    const = sp.symbols('const')

    t = t_re + I*t_im
    delta = d_re + I*d_im

    Hq = sp.zeros(4,4)
    Hq += mu1 * (c1d*c1)
    Hq += mu2 * (c2d*c2)
    Hq += t * (c1d*c2) + sp.conjugate(t) * (c2d*c1)
    Hq += sp.Rational(1,2) * (delta * (c1d*c2d) + sp.conjugate(delta) * (c2*c1))
    Hq += const * sp.eye(4)

    M = sp.simplify(H - Hq)

    # Build real equations
    eqs = []
    for i in range(4):
        for j in range(4):
            expr = sp.simplify(sp.expand_complex(M[i,j]))
            eqs.append(sp.Eq(sp.re(expr), 0))
            eqs.append(sp.Eq(sp.im(expr), 0))

    unknowns = [mu1, mu2, t_re, t_im, d_re, d_im, const]

    # Convert eqs to polynomial form (lhs), skipping trivial True equations
    polys = []
    for e in eqs:
        # e can be an Equality or a BooleanTrue/BooleanFalse from simplification
        if isinstance(e, sp.Equality):
            lhs = sp.simplify(e.lhs)
            if not lhs.equals(0):
                polys.append(lhs)
        else:
            # if equation simplified to False, system inconsistent -> add 1 as nonzero polynomial
            if e == sp.false:
                polys.append(sp.Integer(1))

    print('\nNumber of nontrivial real equations:', len(polys))

    # Attempt Groebner elimination to remove unknowns, leaving polynomials in A,B
    try:
        G = sp.groebner(polys, *unknowns, A, B, order='lex')
        elim = G.eliminate(unknowns)
        print('\nEliminated polynomials in A,B (free-fermion constraints):')
        if len(elim) == 0:
            print('  (no nontrivial polynomial constraints found after elimination)')
        else:
            for p in elim:
                print('  ', sp.factor(p))
    except Exception as e:
        print('\nGroebner elimination failed:', e)


if __name__ == '__main__':
    main()
