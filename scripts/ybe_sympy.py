#!/usr/bin/env python3
"""Symbolic verification of diagonal exponential R solutions for constant YBE.

Usage:
    python3 scripts/ybe_sympy.py

Requires: sympy
"""
from sympy import symbols, cos, sin, I, simplify, expand, factor, pprint

# Symbols (treat as complex symbols)
bx, by, bz = symbols('bx by bz')

# trig shorthand
cx = cos(bx)
sx = sin(bx)
cy = cos(by)
sy = sin(by)
cz = cos(bz)
sz = sin(bz)

# A coefficients (remove global phase)
A0 = cx*cy*cz + I*sx*sy*sz
Ax = I*sx*cy*cz + cx*sy*sz
Ay = I*cx*sy*cz + sx*cy*sz
Az = I*cx*cy*sz + sx*sy*cz

# display A's
print('\nA0 =')
pprint(simplify(A0))
print('\nAx =')
pprint(simplify(Ax))
print('\nAy =')
pprint(simplify(Ay))
print('\nAz =')
pprint(simplify(Az))

# Define E1,E2,E3 in terms of A's
a, b, c_, d = A0, Ax, Ay, Az
E1 = simplify(a*d*(b - c_))
E2 = simplify(b*c_*(a - d))
E3 = simplify(a*b*c_ - a*b*d - a*c_*d + b*c_*d)

print('\n-- Expanded / Factored forms of constraints --')
print('\nE1 (expanded):')
pprint(expand(E1))
print('\nE1 (factored):')
pprint(factor(E1))

print('\nE2 (expanded):')
pprint(expand(E2))
print('\nE2 (factored):')
pprint(factor(E2))

print('\nE3 (expanded):')
pprint(expand(E3))
print('\nE3 (factored):')
pprint(factor(E3))

# Suggest simple derived conditions
print('\n-- Suggested algebraic branches (from factor patterns) --')
print('1) a=b=c=d (nondegenerate continuous family)')
print('   => corresponds to bx==by==bz (mod pi) in trig params')
print("2) any single-term or I+single-term combos: e.g. A0=0 or Ax=0 etc.")
print('   => these correspond to sin/cos special values (0, +/-1) or their complex analogues')

print('\nDone.')
