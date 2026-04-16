#!/usr/bin/env python3
"""Attempt a Groebner-based inclusion test between the YBE residual ideal and candidate factor ideal.

This script:
- builds the YBE residual polynomials in variables A,B,C by substituting exp(iJx)->A etc.
- defines a candidate factor set F_k (A-B, A+B, permutations, and symmetric E-polynomials)
- computes Groebner bases and tests polynomial remainders to check ideal containments.
"""
from sympy import symbols, I, Matrix, simplify, expand, factor, exp
from sympy import kronecker_product as kron
from sympy import groebner, Poly

# build residual expressions with Jx,Jy,Jz, then substitute exp(I*Jx)->A etc.
Jx, Jy, Jz, c00 = symbols('Jx Jy Jz c00')
sigma0 = Matrix([[1,0],[0,1]])
sigma1 = Matrix([[0,1],[1,0]])
sigma2 = Matrix([[0,-I],[I,0]])
sigma3 = Matrix([[1,0],[0,-1]])
K = 0 * kron(sigma0, sigma0) + Jx * kron(sigma1, sigma1) + Jy * kron(sigma2, sigma2) + Jz * kron(sigma3, sigma3)
R = (I * K).exp()
I2 = sigma0
r12 = kron(R, I2)
r23 = kron(I2, R)
P23 = Matrix([[0]*8 for _ in range(8)])
for a in range(2):
    for b in range(2):
        for c in range(2):
            src = a*4 + b*2 + c
            dst = a*4 + c*2 + b
            P23[dst, src] = 1
r13 = P23 * r12 * P23
expr = simplify(r12*r13*r23 - r23*r13*r12)

# variables for algebraic polynomials
A, B, C = symbols('A B C')
subs_map = {exp(I*Jx): A, exp(I*Jy): B, exp(I*Jz): C}

res_polys = []
for i in range(8):
    for j in range(8):
        e = expr[i, j]
        if e != 0:
            # substitute exponentials
            fe = simplify(expand(e.xreplace(subs_map)))
            # clear any rational factors by turning into Poly and extracting numerator
            p = Poly(fe.as_numer_denom()[0], A, B, C)
            # drop zero poly
            if p.total_degree() >= 0 and any(coeff != 0 for coeff in p.coeffs()):
                res_polys.append(p.as_expr())

# deduplicate
res_polys_unique = []
for p in res_polys:
    if not any(simplify(p - q) == 0 for q in res_polys_unique):
        res_polys_unique.append(p)

print('Residual polynomials count:', len(res_polys_unique))

# Candidate factor set (from factorization observations)
F = []
# linear phase factors
F += [A - B, A + B, B - C, B + C, C - A, C + A]
# squared/exponential symmetric polynomials (E_x = A^2 etc)
Ex = A**2; Ey = B**2; Ez = C**2
F += [Ex*Ey + 1,
      Ex*Ey - Ex*Ez + Ey*Ez - 1,
      Ex*Ey + Ex*Ez - Ey*Ez - 1,
      Ex*Ez + 1,
      Ey*Ez + 1]

vars = (A, B, C)

GB_F = groebner(F, *vars, order='lex')
GB_res = groebner(res_polys_unique, *vars, order='lex')

out_lines = []
out_lines.append('Groebner basis of candidate factor ideal (F):\n')
out_lines.append(str(GB_F) + '\n\n')
out_lines.append('Groebner basis of residual ideal (Res):\n')
out_lines.append(str(GB_res) + '\n\n')

# Test 1: each residual polynomial reduces to 0 mod GB_F (i.e., Res ⊆ ideal(F))
res_in_F = True
for p in res_polys_unique:
    r = Poly(p, *vars).rem(*GB_F.polys)
    out_lines.append('Remainder of residual poly modulo GB_F: ' + str(r.as_expr()) + '\n')
    if r.as_expr() != 0:
        res_in_F = False

out_lines.append('\nAll residuals in ideal(F)? ' + str(res_in_F) + '\n\n')

# Test 2: product of factors reduces to 0 mod GB_res (i.e., product ∈ ideal(Res))
prodF = 1
for f in F:
    prodF = simplify(prodF * f)
rprod = Poly(prodF, *vars).rem(*GB_res.polys)
out_lines.append('Remainder of product(F) modulo GB_res: ' + str(rprod.as_expr()) + '\n')
out_lines.append('\nProduct in ideal(Res)? ' + str(rprod.as_expr() == 0) + '\n')

with open('scripts/groebner_proof_out.txt', 'w') as f:
    f.writelines(out_lines)

print('Wrote scripts/groebner_proof_out.txt')
