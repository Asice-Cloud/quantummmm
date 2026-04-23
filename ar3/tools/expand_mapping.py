#!/usr/bin/env python3
"""Produce an extended mapping from Pauli coefficients c_{mu,nu}
to Kitaev/BdG parameters and other operator classes.

Writes results/extended_mapping.txt with symbolic expressions and brief notes.
"""
from sympy import symbols
import os

os.makedirs('results', exist_ok=True)

# define common c coefficients (listed in R22.md)
names = ['c_xx','c_yy','c_xy','c_yx','c_x0','c_y0','c_z0','c_0x','c_0y','c_0z','c_00','c_zz','c_xz','c_yz','c_zx','c_zy']
cs = symbols(' '.join(names))
vals = dict(zip(names, cs))

# Define mapped quantities (per R22 conventions)
t = vals['c_xx'] + vals['c_yy']
ReDelta = vals['c_xx'] - vals['c_yy']
ImDelta = vals['c_xy'] + vals['c_yx']
mu_site = 2*(vals['c_z0'] + vals['c_0z']) - 4*vals['c_zz']
U = 4*vals['c_zz']
E_bond = vals['c_00'] - vals['c_z0'] - vals['c_0z'] + vals['c_zz']

lines = []
lines.append('Extended mapping (conventions follow R22.md)')
lines.append('')
lines.append('Symbols: ' + ', '.join(names))
lines.append('')
lines.append('Kitaev/BdG parameters:')
lines.append(f'  t (real symmetric hopping) = c_xx + c_yy')
lines.append(f'  Re(Delta) = c_xx - c_yy')
lines.append(f'  Im(Delta) = c_xy + c_yx  (pairing imaginary part)')
lines.append(f'  => Delta (complex) = (c_xx - c_yy) + i*(c_xy + c_yx)')
lines.append(f'  mu_site (linear density term per site) = 2*(c_z0 + c_0z) - 4*c_zz')
lines.append(f'  U (nearest-neighbor density-density) = 4*c_zz')
lines.append(f'  E_bond (per-bond constant) = c_00 - c_z0 - c_0z + c_zz')
lines.append('')
lines.append('Chiral/antisymmetric contributions:')
lines.append('  Imaginary/antisymmetric hopping component ~ c_xy - c_yx (generates i*(c1d c2 - c2d c1) like terms)')
lines.append('  Imaginary/antisymmetric pairing component already included in Im(Delta) via c_xy+c_yx (symmetric imaginary)')
lines.append('')
lines.append('String / nonlocal contributions (require Jordan-Wigner string S_j):')
lines.append('  These arise from single-site x/y terms multiplied by neighbor z: c_x0,c_y0,c_0x,c_0y and c_xz,c_yz,c_zx,c_zy')
lines.append(r'  They map to nonlocal operators S_j (c_j + c_j^\\dagger) (2n_{j+1}-1) etc and are NOT captured by simple local (t,Delta,mu,U) parameters.')
lines.append('')
lines.append('Mapping matrix M (R22 summary) for p=[t,Delta,mu_site,U,E_bond]^T and c vector order [c_xx,c_yy,c_xy,c_yx,c_zz,c_z0,c_0z,c_00]^T:')
lines.append('  M = [[1,1,0,0,0,0,0,0], [1,-1,0,0,0,0,0,0], [0,0,0,0,-4,2,2,0], [0,0,0,0,4,0,0,0], [0,0,0,0,1,-1,-1,1]]')
lines.append('  (see R22.md for derivation and notes)')
lines.append('')
lines.append('Notes:')
lines.append(r' - Sign conventions for t depend on Hamiltonian sign: if H_kin = -t(c1^\\dagger c2 + h.c.) then the scalar t above corresponds to + (c_xx + c_yy) in the expansion; some documents absorb a minus sign into t. Be consistent.')
lines.append(' - Delta here is complex; both Re and Im parts come from c_xx,c_yy and c_xy,c_yx respectively.')
lines.append(r' - mu_site excludes any global I\\otimes I constant c_00 which contributes to bond energy E_bond; include c_00 if you want a chemical potential shift.')
lines.append('')

with open('results/extended_mapping.txt','w') as f:
    f.write('\n'.join(lines))

print('Wrote results/extended_mapping.txt')
