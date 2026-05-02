#!/usr/bin/env python3
"""Symbolically expand XXZ local h = P dR/du|_0 / sin(eta) on Pauli basis and compute BdG params."""
import sympy as sp
I = sp.I

# symbols
eta = sp.symbols('eta')

# Pauli matrices
sx = sp.Matrix([[0,1],[1,0]])
sy = sp.Matrix([[0,-I],[I,0]])
sz = sp.Matrix([[1,0],[0,-1]])
id2 = sp.eye(2)

def R_xxz_sym(u, eta):
    s = sp.sin
    a = s(u + eta)
    b = s(u)
    c = s(eta)
    R = sp.zeros(4,4)
    R[0,0] = a
    R[1,1] = b
    R[1,2] = c
    R[2,1] = c
    R[2,2] = b
    R[3,3] = a
    return R

u0 = 0
R0 = R_xxz_sym(u0, eta)
# derivative w.r.t u
da = sp.cos(u0 + eta)
db = sp.cos(u0)
dc = 0
dR0 = sp.zeros(4,4)
dR0[0,0] = da
dR0[1,1] = db
dR0[1,2] = dc
dR0[2,1] = dc
dR0[2,2] = db
dR0[3,3] = da

# permutation operator P
P = sp.zeros(4,4)
for a in range(2):
    for b in range(2):
        in_idx = (a<<1) | b
        out_idx = (b<<1) | a
        P[out_idx, in_idx] = 1

rho = sp.sin(eta)
h = P * dR0 / rho

# Pauli basis
basis = [('I', id2), ('X', sx), ('Y', sy), ('Z', sz)]

def coeff_pauli(H, A, B):
    # coeff = Tr[(A⊗B) H] / 4
    return sp.simplify((sp.Trace(sp.kronecker_product(A, B) * H) / 4))

coeffs = {}
for la, A in basis:
    for lb, B in basis:
        key = f'c_{la}{lb}'
        coeffs[key] = sp.simplify(coeff_pauli(h, A, B))

print('Pauli coefficients (symbolic):')
for k in sorted(coeffs.keys()):
    print(k, ':', sp.simplify(coeffs[k]))

# Build BdG mapping (simplified): use mappings from myabs.md
cx_xx = coeffs['c_XX']
cx_yy = coeffs['c_YY']
cx_xy = coeffs['c_XY']
cx_yx = coeffs['c_YX']
cz0 = coeffs['c_ZI']
c0z = coeffs['c_IZ']
c00 = coeffs['c_II']

Delta = sp.simplify(cx_xx - cx_yy - I*(cx_xy + cx_yx))
t_par = sp.simplify(cx_xx + cx_yy + I*(cx_xy - cx_yx))
mu = sp.simplify(-2*(cz0 + c0z))

print('\nDerived BdG params (symbolic):')
print('Delta =', sp.simplify(Delta))
print('t =', sp.simplify(t_par))
print('mu =', sp.simplify(mu))

# Check whether coefficients lie in allowed subspace from Groebner (i.e., which are zero)
zero_list = [k for k,v in coeffs.items() if sp.simplify(v)==0]
nonzero_list = [k for k,v in coeffs.items() if sp.simplify(v)!=0]
print('\nZero coefficients:', zero_list)
print('\nNonzero coefficients:', nonzero_list)

print('\nSimplified Delta (factored):')
print(sp.factor(Delta))
