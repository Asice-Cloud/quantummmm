#!/usr/bin/env python3
"""Solve diagonal CYBE ([H12,H13]+[H12,H23]+[H13,H23]=0) for
H = c_xx σ^x⊗σ^x + c_yy σ^y⊗σ^y + c_zz σ^z⊗σ^z, then build a Kitaev chain
from the quadratic part (set c_zz=0 when needed) and diagonalize BdG.
"""
import sympy as sp
import numpy as np
from numpy.linalg import eigh

# SymPy setup
cxx, cyy, czz = sp.symbols('cxx cyy czz', real=True)
# Pauli matrices
sx = sp.Matrix([[0,1],[1,0]])
sy = sp.Matrix([[0,-sp.I],[sp.I,0]])
sz = sp.Matrix([[1,0],[0,-1]])
I2 = sp.eye(2)

# Two-site operators
XX = sp.kronecker_product(sx, sx)
YY = sp.kronecker_product(sy, sy)
ZZ = sp.kronecker_product(sz, sz)

# Three-site space: build H12, H13, H23 in 2^3=8-dim
def embed12(M):
    return sp.kronecker_product(M, I2)
def embed23(M):
    return sp.kronecker_product(I2, M)
def embed13(M):
    # act on sites 1 and 3: M on 1,3 with identity on 2
    return sp.kronecker_product(M, I2)
# For embed13 we need ordering 1⊗2⊗3; create with explicit kron: (A⊗I) on (1,3) is
# construct using indices: embed13 = sum_{a,b,c,d} M_{ac,bd} |a><c|⊗I⊗|b><d|
# simpler: use sympy's kronecker_product with permutations: (sx on site1)⊗I⊗(sx on site3)

def kron3(A,B,C):
    return sp.kronecker_product(A, sp.kronecker_product(B, C))

H12 = cxx * kron3(sx, sx, I2) + cyy * kron3(sy, sy, I2) + czz * kron3(sz, sz, I2)
H23 = cxx * kron3(I2, sx, sx) + cyy * kron3(I2, sy, sy) + czz * kron3(I2, sz, sz)
H13 = cxx * kron3(sx, I2, sx) + cyy * kron3(sy, I2, sy) + czz * kron3(sz, I2, sz)

# Classical/Yang-Baxter-like constraint (first nontrivial order)
Comm = (H12*H13 - H13*H12) + (H12*H23 - H23*H12) + (H13*H23 - H23*H13)

# Extract independent polynomial equations from Comm == 0
eqs = []
for i in range(8):
    for j in range(8):
        entry = sp.simplify(Comm[i,j])
        if entry != 0:
            eqs.append(sp.Eq(sp.simplify(entry), 0))

# Reduce to a minimal independent set
eqs = sp.poly(eqs[0].lhs, cxx, cyy, czz).coeffs() if eqs else []
# Instead of the above awkward reduction, form the set of unique symbolic expressions
exprs = []
for i in range(8):
    for j in range(8):
        e = sp.simplify(Comm[i,j])
        if e!=0:
            exprs.append(sp.expand(e))
exprs = list({sp.simplify(e):None for e in exprs}.keys())
# Solve exprs == 0
sol = sp.solve(exprs, (cxx, cyy, czz), dict=True)
print('SymPy solutions (exact):')
print(sol)

# If solve fails to find parametric families, compute polynomial GCD relations
# Convert expressions to polynomial ideal and compute Groebner basis to find relations
G = sp.groebner(exprs, cxx, cyy, czz, order='lex')
print('\nGroebner basis:')
print(G)

# Inspect basis for simple algebraic relations
basis = list(G)

# If basis gives relations, show them
if basis:
    print('\nRelations from Groebner basis:')
    for b in basis:
        print(sp.factor(b))

# pick representative numeric solutions from relations
# If basis contains factor cxx-cyy etc, choose sample values
# We'll try to find simple families: e.g., cxx=cyy, czz free, or others.

# For numeric BdG tests we prefer czz = 0 (no interactions). Find simple relations.
# Try to solve for czz in terms of cxx,cyy if possible
czz_expr = None
try:
    sol_czz = sp.solve(sp.Eq(basis[0],0), czz)
    if sol_czz:
        czz_expr = sol_czz[0]
        print('\nSolved czz in terms of cxx,cyy:', czz_expr)
except Exception:
    pass

# Choose representative point(s)
candidates = []
# Prefer sample with czz=0
candidates.append({cxx:1.0, cyy:0.5, czz:0.0})
# Also try symmetric cxx=cyy
candidates.append({cxx:1.0, cyy:1.0, czz:0.0})

# Function to build BdG and compute spectrum for given cxx,cyy,czz (we use czz only to set mu)
def bdg_spectrum(params, L=40):
    cx = float(params[cxx]); cy = float(params[cyy]); cz = float(params[czz])
    t = cx + cy
    Delta = cx - cy
    # chemical potential from czz (using earlier mapping mu_site = 2(c_z0+c_0z) - 4 c_zz,
    # with no z0 terms we take mu = - ( -4 czz )? here use mu = 4*cz (sign conventions vary)
    mu = 0.0 if abs(cz) < 1e-12 else -4.0*cz

    L = int(L)
    # single-particle h and pairing Delta_mat (LxL)
    h = np.zeros((L,L), dtype=float)
    D = np.zeros((L,L), dtype=complex)
    for i in range(L-1):
        h[i,i+1] = -t
        h[i+1,i] = -t
        D[i,i+1] = Delta
        D[i+1,i] = -Delta
    for i in range(L):
        h[i,i] = -mu

    # BdG matrix
    top = np.concatenate([h, D], axis=1)
    bottom = np.concatenate([-D.conjugate(), -h.T], axis=1)
    Hbdg = np.concatenate([top, bottom], axis=0)
    Hbdg = (Hbdg + Hbdg.conj().T)/2.0
    vals, vecs = eigh(Hbdg)
    vals = np.sort(np.real(vals))
    return vals

# Run spectra for candidates and print small eigenvalues
for idx, cand in enumerate(candidates):
    vals = bdg_spectrum(cand, L=60)
    # select smallest positive eigenvalues
    small = np.abs(vals)
    small_sorted = np.sort(small)[:8]
    print(f"\nCandidate {idx+1}: cxx={cand[cxx]}, cyy={cand[cyy]}, czz={cand[czz]}")
    print('smallest |E| values:', np.round(small_sorted, 6))

print('\nDone.')
