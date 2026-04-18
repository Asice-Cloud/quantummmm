#!/usr/bin/env python3
"""Enumerate minimal atomic polynomial factor sets (branches) for exponential ansatz.
Reads derive_exp_constraints outputs by re-running the factorization, substitutes
u=exp(I*bx), v=exp(I*by), w=exp(I*bz), clears denominators and uses Groebner
bases to test if a set of factor polynomials implies all constraint polynomials.
Writes minimal branches to scripts/exp_minimal_branches.txt.
"""
import sympy as sp
from sympy import symbols, Matrix, I, exp
from sympy.matrices import kronecker_product
from sympy import simplify, factor_list, expand

# symbols
c00, bx, by, bz = symbols('c00 bx by bz')

# Pauli matrices
sx = Matrix([[0,1],[1,0]])
sy = Matrix([[0,-I],[I,0]])
sz = Matrix([[1,0],[0,-1]])
I2 = Matrix.eye(2)

# Two-site ops
Sx = kronecker_product(sx, sx)
Sy = kronecker_product(sy, sy)
Sz = kronecker_product(sz, sz)
I4 = Matrix.eye(4)

# build R and three-site B
H = c00*I4 + bx*Sx + by*Sy + bz*Sz
R = (I*H).exp()
I2_3 = Matrix.eye(2)
R12 = kronecker_product(R, I2_3)
R23 = kronecker_product(I2_3, R)
B = simplify(R12*R23*R12 - R23*R12*R23)

# collect unique nonzero entries
entries = []
for e in B:
    ee = simplify(sp.trigsimp(sp.expand(e)))
    if ee != 0:
        entries.append(ee)

uniq = []
for e in entries:
    fe = sp.factor(e)
    if not any(sp.simplify(fe - u) == 0 for u in uniq):
        uniq.append(fe)

print('Found', len(uniq), 'unique constraint expressions')

# substitute exp(I*bx)->u etc
u,v,w = sp.symbols('u v w')
# substitute exponentials by algebraic vars; set exp(I*c00)=1 (global phase)
# build a robust substitution map: map exp(k*I*bx) -> u**k for a range of k
subs_map = {}
# map powers of exp(I*c00) to 1 (global phase)
for k in range(-6, 7):
    subs_map[sp.exp(k*I*c00)] = 1
# map powers of exp(I*bx), exp(I*by), exp(I*bz) to u**k, v**k, w**k
for k in range(-6, 7):
    subs_map[sp.exp(k*I*bx)] = u**k
    subs_map[sp.exp(k*I*by)] = v**k
    subs_map[sp.exp(k*I*bz)] = w**k

poly_exprs = []
for e in uniq:
    esub = sp.simplify(e.xreplace(subs_map))
    # clear negative powers by multiplying with u^3 v^3 w^3
    esub = sp.simplify(expand(esub * u**3 * v**3 * w**3))
    # convert to polynomial over algebraic extension
    # ensure no exp(I*c00) remains
    esub = esub.xreplace({sp.exp(I*c00): 1})
    poly = sp.Poly(sp.expand(esub), u, v, w, domain='QQ')
    poly_exprs.append(poly.as_expr())

# Extract irreducible factors from each uniq expression (after substitution and clearing)
factors = []
for expr in poly_exprs:
    fl = factor_list(expr, u, v, w)
    # factor_list returns (coeff, [(factor, exp), ...])
    coeff, facs = fl
    for fac, powr in facs:
        # ensure factor has no exp(I*c00)
        fac = fac.xreplace({sp.exp(I*c00): 1})
        factors.append(sp.Poly(sp.expand(fac), u, v, w, domain='QQ').as_expr())

# unique factors
uniq_factors = []
for f in factors:
    if not any(sp.simplify(f - g) == 0 for g in uniq_factors):
        uniq_factors.append(f)

print('Extracted', len(uniq_factors), 'unique atomic factors')

# enumerate minimal combinations that imply all poly_exprs
from itertools import combinations
valid_combos = []
for r in range(1, min(5, len(uniq_factors))+1):
    for combo in combinations(range(len(uniq_factors)), r):
        pols = [sp.expand(uniq_factors[i]) for i in combo]
        # compute Groebner basis over QQ (rational coefficients)
        G = sp.groebner(pols, u, v, w, order='lex', domain='QQ')
        ok = True
        for expr in poly_exprs:
            rem = G.reduce(sp.Poly(expr, u, v, w, domain='QQ'))[1]
            if sp.simplify(rem.as_expr()) != 0:
                ok = False
                break
        if ok:
            valid_combos.append(tuple(combo))

# minimize combos
minimal = []
for c in valid_combos:
    is_min = True
    for d in valid_combos:
        if set(d) < set(c):
            is_min = False
            break
    if is_min:
        minimal.append(c)

# write human-readable results
with open('scripts/exp_minimal_branches.txt','w') as f:
    f.write('Minimal atomic factor index sets (refers to factors listed below)\n')
    f.write('\nFactors (index: factor_expr):\n')
    for i,fac in enumerate(uniq_factors):
        f.write(f'{i}: {sp.sstr(fac)}\n')
    f.write('\nMinimal combos (indices):\n')
    for c in minimal:
        f.write(str(c) + '\n')

print('Wrote scripts/exp_minimal_branches.txt with', len(minimal), 'minimal combos')
print('Done')
