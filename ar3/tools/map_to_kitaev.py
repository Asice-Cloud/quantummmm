#!/usr/bin/env python3
"""Map derived Majorana bilinear coefficients to Kitaev-like parameters.

Outputs symbolic expressions for t_eff, Delta_eff, mu_eff in terms of c_xx,c_yy,c_xy,c_yx (and c00).
"""
from sympy import symbols, I, simplify

# define symbols
c_xx, c_yy, c_xy, c_yx, c00 = symbols('c_xx c_yy c_xy c_yx c00')

# coefficients extracted from derive_mapping.py
hop1 = c_xx + I*c_xy - I*c_yx + c_yy
hop2 = c_xx - I*c_xy + I*c_yx + c_yy
pair2 = c_xx - c_yy + I*(c_xy + c_yx)
pair1 = simplify(pair2.conjugate())

# Interpret t, Delta, mu (conventions: H_kin = -t (c1d c2 + h.c.), H_pair = Delta c1 c2 + h.c.)
# Here hop1 corresponds to coeff of c1d*c2, hop2 to c1*c2d (its hermitian conjugate)
hop_avg = simplify((hop1 + hop2)/2)
t_eff = simplify(-hop_avg)
Delta_eff = simplify(pair2)
# On-site mu not generated in this diagonal Pauli subspace (set to symbolic c00 if provided)
mu_eff = simplify(0 + c00)

out = f"t_eff = {simplify(t_eff)}\nDelta_eff = {simplify(Delta_eff)}\nmu_eff = {simplify(mu_eff)}\n"

print('Kitaev-like mapping (conventions: H_kin = -t(c1d c2 + h.c.), H_pair = Delta c1 c2 + h.c.)')
print(out)

# write to results
import os
os.makedirs('results', exist_ok=True)
with open('results/kitaev_mapping.txt','w') as f:
    f.write('Kitaev-like mapping (conventions: H_kin = -t(c1d c2 + h.c.), H_pair = Delta c1 c2 + h.c.)\n')
    f.write(out)

print('Wrote results/kitaev_mapping.txt')
