#!/usr/bin/env python3
"""Symbolic check for XYZ/8-vertex style R(u).

Parameterize R(u) near u=0 with entries a(u)=A + a1*u, b(u)=B + b1*u,
c(u)=C + c1*u, d(u)=D + d1*u and compute h = P dR/du / rho, expand on Pauli
and report symbolic Delta.
"""
from __future__ import annotations
import sympy as sp
I = sp.I

# symbols
A,B,C,D = sp.symbols('A B C D')
a1,b1,c1,d1 = sp.symbols('a1 b1 c1 d1')
rho = sp.symbols('rho')

# define a(u) etc linearized: value + derivative*u
a = A + a1*sp.symbols('u')
b = B + b1*sp.symbols('u')
c = C + c1*sp.symbols('u')
d = D + d1*sp.symbols('u')

# R(u) matrix in basis 00,01,10,11 with 8-vertex pattern
# [ a  0  0  d
#   0  b  c  0
#   0  c  b  0
#   d  0  0  a ]
R = sp.Matrix([[a,0,0,d],[0,b,c,0],[0,c,b,0],[d,0,0,a]])

# derivative at u=0: coefficients are a1,b1,c1,d1 (since derivative of value + coef*u)
dRdu = sp.Matrix([[a1,0,0,d1],[0,b1,c1,0],[0,c1,b1,0],[d1,0,0,a1]])

# permutation P
P = sp.zeros(4,4)
for a_idx in range(2):
    for b_idx in range(2):
        in_idx = (a_idx<<1) | b_idx
        out_idx = (b_idx<<1) | a_idx
        P[out_idx, in_idx] = 1

# local h
h = sp.simplify(P * dRdu / rho)

# Pauli matrices
sx = sp.Matrix([[0,1],[1,0]])
sy = sp.Matrix([[0,-I],[I,0]])
sz = sp.Matrix([[1,0],[0,-1]])
id2 = sp.eye(2)

def coeff_pauli(H, A, B):
    return sp.simplify((sp.Trace(sp.kronecker_product(A, B) * H) / 4))

basis = [('I', id2), ('X', sx), ('Y', sy), ('Z', sz)]
coeffs = {}
for la, Aop in basis:
    for lb, Bop in basis:
        key = f'c_{la}{lb}'
        coeffs[key] = coeff_pauli(h, Aop, Bop)

Delta = sp.simplify(coeffs['c_XX'] - coeffs['c_YY'] - I*(coeffs['c_XY'] + coeffs['c_YX']))

print('Symbolic Pauli coefficients (nonzero subset):')
for k in sorted(coeffs.keys()):
    v = sp.simplify(coeffs[k])
    if v!=0:
        print(k, ':', v)

print('\nSymbolic Delta =')
sp.pprint(Delta)

print('\nIs Delta identically zero? ->', sp.simplify(Delta)==0)
