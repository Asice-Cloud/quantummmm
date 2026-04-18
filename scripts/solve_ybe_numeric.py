#!/usr/bin/env python3
"""Numerically solve YBE for complex bx,by,bz using sympy.nsolve.
Find solutions to selected independent complex equations (entries of B==0).
"""
from sympy import symbols, Matrix, I, exp, cos, sin
from sympy import kronecker_product, simplify
from sympy import N
from sympy import nsolve
import sympy as sp
import random

# symbols
c00, bx, by, bz = symbols('c00 bx by bz')

# Pauli matrices
sx = Matrix([[0,1],[1,0]])
sy = Matrix([[0,-I],[I,0]])
sz = Matrix([[1,0],[0,-1]])
I2 = Matrix.eye(2)

Sx = kronecker_product(sx, sx)
Sy = kronecker_product(sy, sy)
Sz = kronecker_product(sz, sz)
I4 = Matrix.eye(4)

# local H and R
H = c00*I4 + bx*Sx + by*Sy + bz*Sz
R = (I*H).exp()

# three-site
I2_ = Matrix.eye(2)
R12 = kronecker_product(R, I2_)
R23 = kronecker_product(I2_, R)
B = simplify(R12*R23*R12 - R23*R12*R23)

# choose three complex equations (matrix entries) expected independent
# We'll pick B[0,1], B[0,2], B[1,2]
eqs = [B[0,1], B[0,2], B[1,2]]

# set c00=0 for simplicity (global phase)
subs = {c00: 0}

# prepare numeric functions
eqs_sub = [sp.simplify(e.subs({c00:0})) for e in eqs]

# convert to lambdas? nsolve accepts sympy expressions directly
sols = []
found = []

print('Attempting numeric nsolve for complex bx,by,bz (multiple starts)')
for attempt in range(12):
    # random complex initial guess
    guess = [complex(random.uniform(-1,1), random.uniform(-1,1)) for _ in range(3)]
    try:
        s = nsolve(eqs_sub, (bx, by, bz), guess, tol=1e-14, maxsteps=200)
        s = [complex(N(si)) for si in s]
        # normalize solutions to principal branch (map imaginary small to 0)
        def near_eq(a,b,eps=1e-8):
            return abs(a-b) < eps
        unique = True
        for u in found:
            if all(near_eq(u[i], s[i]) for i in range(3)):
                unique = False
                break
        if unique:
            found.append(s)
            print('Found solution:', s)
    except Exception as e:
        #print('attempt failed', e)
        pass

if not found:
    print('No solutions found with these equation choices and attempts.')
else:
    with open('scripts/ybe_numeric_solutions.txt','w') as f:
        for s in found:
            f.write(','.join(repr(x) for x in s)+"\n")
    print('Wrote scripts/ybe_numeric_solutions.txt with', len(found), 'solutions')

print('Done')
