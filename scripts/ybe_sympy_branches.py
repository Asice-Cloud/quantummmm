#!/usr/bin/env python3
"""Enumerate algebraic branches of diagonal exponential R that satisfy constant YBE.

This script symbolically constructs A0,Ax,Ay,Az from bx,by,bz (complex allowed),
forms the YBE constraints E1,E2,E3 in abstract variables a,b,c,d and then
enumerates minimal combinations of atomic conditions (zeros or equalities)
that imply the constraints. Outputs branches both in (a,b,c,d) form and as
the corresponding trig conditions in bx,by,bz.

Usage:
    python3 scripts/ybe_sympy_branches.py

Requires: sympy
"""
from sympy import symbols, cos, sin, I, simplify
from itertools import combinations

# trig symbols
bx, by, bz = symbols('bx by bz')

# trig shorthand
cx = cos(bx); sx = sin(bx)
cy = cos(by); sy = sin(by)
cz = cos(bz); sz = sin(bz)

# A coefficients (remove global phase)
A0 = simplify(cx*cy*cz + I*sx*sy*sz)
Ax = simplify(I*sx*cy*cz + cx*sy*sz)
Ay = simplify(I*cx*sy*cz + sx*cy*sz)
Az = simplify(I*cx*cy*sz + sx*sy*cz)

print('\nComputed A coefficients (in terms of bx,by,bz):\n')
print('A0 =', A0)
print('Ax =', Ax)
print('Ay =', Ay)
print('Az =', Az)

# Abstract vars for algebraic classification
a, b, c, d = symbols('a b c d')

# Abstract constraints (these are equivalent to the YBE in the 4-dim basis)
E1 = a*d*(b - c)
E2 = b*c*(a - d)
E3 = a*b*c - a*b*d - a*c*d + b*c*d

atoms = ["a=0", "b=0", "c=0", "d=0", "b=c", "a=d"]

def apply_conditions(expr, conds):
    """Apply a set of atomic conditions to an expression and simplify."""
    subs = {}
    for con in conds:
        if con == 'a=0':
            subs[a] = 0
        elif con == 'b=0':
            subs[b] = 0
        elif con == 'c=0':
            subs[c] = 0
        elif con == 'd=0':
            subs[d] = 0
        elif con == 'b=c':
            subs[c] = b
        elif con == 'a=d':
            subs[d] = a
    return simplify(expr.subs(subs))

valid_branches = []

# enumerate all non-empty combinations of atomic conditions
for r in range(1, len(atoms)+1):
    for combo in combinations(atoms, r):
        # apply to abstract equations
        e1 = apply_conditions(E1, combo)
        e2 = apply_conditions(E2, combo)
        e3 = apply_conditions(E3, combo)
        if simplify(e1) == 0 and simplify(e2) == 0 and simplify(e3) == 0:
            valid_branches.append(tuple(sorted(combo)))
print(f'\nFound {len(valid_branches)} valid branches (including non-minimal supersets).')
# minimize branches (remove supersets)
minimal = []
for bset in valid_branches:
    is_min = True
    for other in valid_branches:
        if set(other) < set(bset):
            is_min = False
            break
    if is_min and bset not in minimal:
        minimal.append(bset)

print('\n\nDiscovered algebraic branches (minimal sets of atomic conditions):\n')
for i, br in enumerate(sorted(minimal), 1):
    print(f'Branch {i}:', br)
    # show what these conditions mean when applied to A-coeffs
    conds = list(br)
    # map abstract a,b,c,d to A0,Ax,Ay,Az
    subs_map = {a: A0, b: Ax, c: Ay, d: Az}
    e1_trig = apply_conditions(E1.subs(subs_map), conds)
    e2_trig = apply_conditions(E2.subs(subs_map), conds)
    e3_trig = apply_conditions(E3.subs(subs_map), conds)
    print('  -> E1 after substitution:', simplify(e1_trig))
    print('  -> E2 after substitution:', simplify(e2_trig))
    print('  -> E3 after substitution:', simplify(e3_trig))
    # also print simplified forms of the imposed conditions in trig terms
    print('  -> imposed (in A-coords):')
    for con in conds:
        if con == 'a=0':
            print('     A0 = 0  =>', simplify(A0))
        elif con == 'b=0':
            print('     Ax = 0  =>', simplify(Ax))
        elif con == 'c=0':
            print('     Ay = 0  =>', simplify(Ay))
        elif con == 'd=0':
            print('     Az = 0  =>', simplify(Az))
        elif con == 'b=c':
            print('     Ax = Ay  =>', simplify(Ax - Ay))
        elif con == 'a=d':
            print('     A0 = Az  =>', simplify(A0 - Az))
    print()

print('\nNote: the branch list includes the nondegenerate family a=b=c=d as the condition {"b=c","a=d","a!=0"},\n      and many degenerate branches where one or more A_* vanish or equal each other.')

print('\nDone.')
