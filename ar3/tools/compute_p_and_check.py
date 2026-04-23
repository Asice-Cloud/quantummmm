#!/usr/bin/env python3
"""Compute p = (t,Delta,U,mu,C_perbond) from c_{mu,nu} in ybe.md and classify ABS vs MZM.

Uses the mapping documented in ybe.md:
  t = c_xx + c_yy + i(c_xy - c_yx)
  Delta = c_xx - c_yy - i(c_xy + c_yx)
  U = 4*c_zz
  mu = 4*c_zz - 2*(c_z0 + c_0z)
  C_perbond = c_zz - (c_z0 + c_0z) + c_00

Classification heuristic (quadratic case): if |Delta|>0 and |mu| < 2|t| => topological (MZM),
else likely non-topological (possible ABS). Outputs numeric values and recommendation.
"""
import re
from pathlib import Path
import numpy as np


def parse_c_from_ybe(path):
    text = Path(path).read_text(encoding='utf-8')
    # look for c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z, c_00
    keys = ['c_xx','c_yy','c_xy','c_yx','c_zz','c_z0','c_0z','c_00']
    vals = {}
    for k in keys:
        # try patterns like c_xx = 0.5 or c_xx: 0.5
        m = re.search(r'%s\s*[:=]\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)' % k, text)
        if m:
            vals[k] = float(m.group(1))
    return vals


def compute_p(c):
    # default 0 if missing
    get = lambda k: c.get(k, 0.0)
    c_xx = get('c_xx')
    c_yy = get('c_yy')
    c_xy = get('c_xy')
    c_yx = get('c_yx')
    c_zz = get('c_zz')
    c_z0 = get('c_z0')
    c_0z = get('c_0z')
    c_00 = get('c_00')

    t = c_xx + c_yy + 1j*(c_xy - c_yx)
    Delta = c_xx - c_yy - 1j*(c_xy + c_yx)
    U = 4.0 * c_zz
    mu = 4.0 * c_zz - 2.0*(c_z0 + c_0z)
    C_perbond = c_zz - (c_z0 + c_0z) + c_00
    return dict(t=t, Delta=Delta, U=U, mu=mu, C_perbond=C_perbond)


def classify(p):
    t = p['t']
    Delta = p['Delta']
    mu = p['mu']
    # simple Kitaev chain criterion (quadratic): |mu| < 2|t| and Delta != 0
    if abs(Delta) > 1e-12 and abs(mu) < 2.0 * abs(t):
        return 'topological (supports MZM under quadratic approx)'
    else:
        return 'non-topological (likely ABS or needs further check)'


def main():
    path = Path(__file__).resolve().parents[1] / 'ybe.md'
    c = parse_c_from_ybe(path)
    if not c:
        print('No coefficients parsed from ybe.md; using defaults (from tools/abs_r_verify).')
        c = {'c_xx':1.0,'c_yy':1.0,'c_xy':0.0,'c_yx':0.0,'c_zz':0.0,'c_z0':0.0,'c_0z':0.0,'c_00':0.0}
    p = compute_p(c)
    cls = classify(p)

    print('Parsed c values:')
    for k in sorted(c.keys()):
        print('  %-6s = %g' % (k, c[k]))
    print('\nComputed p:')
    print('  t      =', p['t'])
    print('  Delta  =', p['Delta'])
    print('  U      =', p['U'])
    print('  mu     =', p['mu'])
    print('  C_perb =', p['C_perbond'])
    print('\nClassification (heuristic):', cls)

    if cls.startswith('non-topological'):
        print('\nRecommendation:')
        print(' - If you expect MZM, reduce/zero terms that break quadratic form (c_zz, c_x0, c_0x, etc).')
        print(' - Or build finite-chain BdG using (t,Delta,mu) and check LDOS or Pfaffian/winding number.')


if __name__ == "__main__":
    main()
