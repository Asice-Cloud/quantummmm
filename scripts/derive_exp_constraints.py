#!/usr/bin/env python3
"""Derive YBE constraints directly from R = exp(i * sum c_mm sigma^m ⊗ sigma^m).
Prints simplified independent equations in terms of b_x,b_y,b_z (and c00).
"""
from sympy import symbols, Matrix, I, cos, sin, simplify, trigsimp
from sympy import expand, factor
from sympy.matrices import kronecker_product

# symbols
c00, bx, by, bz = symbols('c00 bx by bz')

# Pauli matrices (2x2)
sx = Matrix([[0,1],[1,0]])
sy = Matrix([[0,-I],[I,0]])
sz = Matrix([[1,0],[0,-1]])
I2 = Matrix.eye(2)

# Two-site operators
Sx = kronecker_product(sx, sx)
Sy = kronecker_product(sy, sy)
Sz = kronecker_product(sz, sz)
I4 = Matrix.eye(4)

# local H and R
H = c00*I4 + bx*Sx + by*Sy + bz*Sz
R = (I*H).exp()  # matrix exponential

# Build R12 and R23 as 8x8 matrices (three sites)
I2_3 = Matrix.eye(2)
R12 = kronecker_product(R, I2_3)
R23 = kronecker_product(I2_3, R)

# Compute commutator difference
B = simplify(R12*R23*R12 - R23*R12*R23)

# collect nonzero independent entries
entries = []
for e in B:
    ee = simplify(trigsimp(expand(e)))
    if ee != 0:
        entries.append(ee)

# deduplicate (up to sign or factor) by factoring
uniq = []
for e in entries:
    fe = factor(e)
    if not any(simplify(fe - u) == 0 for u in uniq):
        uniq.append(fe)

print('\nNumber of nonzero constraint expressions:', len(uniq))
for i,u in enumerate(uniq,1):
    print(f'Constraint {i}:')
    print(u)
    print()

# Try to simplify to minimal generating set via gcd-like approach: remove multiples
minimal = []
for u in uniq:
    is_mul = False
    for v in uniq:
        if u == v:
            continue
        # if u is factorable by v (i.e., u = v * q), check if simplify(u/v) is independent
        q = simplify(u / v)
        if q!=0 and q.is_commutative:
            # crude check: if ratio simplifies to expression without bx/by/bz, consider dependent
            if q.free_symbols <= {c00}:
                is_mul = True
                break
    if not is_mul:
        minimal.append(u)

print('Minimal set (heuristic):')
for i,u in enumerate(minimal,1):
    print(f'  {i}:', u)

print('\nDone.')
