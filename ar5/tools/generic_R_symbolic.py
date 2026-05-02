#!/usr/bin/env python3
"""Symbolic check: general two-qubit R(0) and dR/du at u=0.

Constructs symbolic R0 and dR0 matrices (free symbols), forms
h = P * dR0 / rho, expands on Pauli basis and prints the symbolic
expression for Delta = c_xx - c_yy - i(c_xy + c_yx).
"""
from __future__ import annotations
import sympy as sp

I = sp.I

# symbols for R0 entries and dR entries (16 entries each)
r = sp.symbols('r0:16')
# r is a flat tuple length 16; reshape
R_syms = sp.Matrix([[r[i*4 + j] for j in range(4)] for i in range(4)])
d = sp.symbols('d0:16')
dR_syms = sp.Matrix([[d[i*4 + j] for j in range(4)] for i in range(4)])
rho = sp.symbols('rho')

# Pauli matrices
sx = sp.Matrix([[0,1],[1,0]])
sy = sp.Matrix([[0,-I],[I,0]])
sz = sp.Matrix([[1,0],[0,-1]])
id2 = sp.eye(2)

# permutation P on two qubits
P = sp.zeros(4,4)
for a in range(2):
    for b in range(2):
        in_idx = (a<<1) | b
        out_idx = (b<<1) | a
        P[out_idx, in_idx] = 1

# local h
h = sp.simplify(P * dR_syms / rho)

def coeff_pauli(H, A, B):
    return sp.simplify((sp.Trace(sp.kronecker_product(A, B) * H) / 4))

basis = [('I', id2), ('X', sx), ('Y', sy), ('Z', sz)]
coeffs = {}
for la, A in basis:
    for lb, B in basis:
        key = f'c_{la}{lb}'
        coeffs[key] = coeff_pauli(h, A, B)

print('Pauli coefficients (symbolic) keys:')
print(sorted(coeffs.keys()))

cx_xx = coeffs['c_XX']
cx_yy = coeffs['c_YY']
cx_xy = coeffs['c_XY']
cx_yx = coeffs['c_YX']

Delta = sp.simplify(cx_xx - cx_yy - I*(cx_xy + cx_yx))

print('\nSymbolic Delta =')
sp.pprint(Delta)

# Determine whether Delta is identically zero
is_zero = sp.simplify(Delta) == 0
print('\nIs Delta identically zero? ->', bool(is_zero))

# If not identically zero, print a short example: pick random numeric substitution
if not is_zero:
    subs = {rho:1}
    # assign small integer values to d entries to show nonzero
    for idx, sym in enumerate(d):
        subs[sym] = idx+1
    # r entries set to zero
    for sym in r:
        subs[sym] = 0
    print('\nExample numeric Delta with d00..d15 = 1..16, rho=1:')
    print(Delta.subs(subs))
