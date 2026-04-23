#!/usr/bin/env python3
"""Symbolic factorization of constant YBE for R=exp(i H_P)
with H_P = b_x S_x + b_y S_y + b_z S_z (diagonal two-site generators).
Compares factorization in variables u=e^{2 i b_x}, v=e^{2 i b_y}, w=e^{2 i b_z}.
"""
import sympy as sp

# symbols
bx, by, bz = sp.symbols('bx by bz')
I = sp.I

# Pauli
sx = sp.Matrix([[0,1],[1,0]])
sy = sp.Matrix([[0,-I],[I,0]])
sz = sp.Matrix([[1,0],[0,-1]])
I2 = sp.eye(2)

# Two-site S_alpha
Sx = sp.kronecker_product(sx, sx)
Sy = sp.kronecker_product(sy, sy)
Sz = sp.kronecker_product(sz, sz)

# H_P and R (drop global c00)
Hp = bx*Sx + by*Sy + bz*Sz
R = sp.exp(I*Hp)

# Embed on three sites using kron3
def kron3(A,B,C):
    return sp.kronecker_product(A, sp.kronecker_product(B, C))

R12 = kron3(R, I2, I2)
R23 = kron3(I2, R, I2)
# build R13 by acting R on sites 1 and 3: perform permutation with swap of middle space
# construct R13 by explicit action: sum_{a,b,c,d} R_{ab,cd} |a><c|⊗I⊗|b><d|
# but easier: build local operators for each S_alpha on 1 and 3
Sx13 = kron3(sx, I2, sx)
Sy13 = kron3(sy, I2, sy)
Sz13 = kron3(sz, I2, sz)
Hp13 = bx*Sx13 + by*Sy13 + bz*Sz13
R13 = sp.exp(I*Hp13)

# Compute residual
Res = sp.simplify(R12*R23*R12 - R23*R12*R23)

# Collect nonzero matrix entries and factor them in variables u,v,w
u = sp.symbols('u')
v = sp.symbols('v')
w = sp.symbols('w')
# substitution u=exp(2 i bx) etc
subs = {sp.exp(2*I*bx): u, sp.exp(2*I*by): v, sp.exp(2*I*bz): w,
        sp.cos(bx): (u+1/u)/2, sp.cos(by): (v+1/v)/2, sp.cos(bz): (w+1/w)/2,
        sp.sin(bx): (u-1/u)/(2*I), sp.sin(by): (v-1/v)/(2*I), sp.sin(bz): (w-1/w)/(2*I)}

exprs = []
for i in range(8):
    for j in range(8):
        e = sp.simplify(Res[i,j])
        if e!=0:
            ef = sp.simplify(sp.factor(sp.together(e.rewrite(sp.exp)).subs(subs)))
            exprs.append(ef)

# uniquify
exprs_u = []
for e in exprs:
    if not any(sp.simplify(e - ee)==0 for ee in exprs_u):
        exprs_u.append(e)

print('Found', len(exprs_u), 'distinct nonzero factorized expressions.\n')
for k,e in enumerate(exprs_u[:20],1):
    print('Expr',k,':')
    print(sp.factor(e))
    print('---\n')

# Also compute simple factorization for a representative entry (0,0)
rep = sp.simplify(Res[0,0])
print('\nRepresentative (0,0) entry factorized (trig form):')
print(sp.factor(sp.simplify(sp.re(sp.expand_trig(rep)))))

print('\nDone.')
