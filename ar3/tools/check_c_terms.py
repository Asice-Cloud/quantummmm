#!/usr/bin/env python3
"""Check ybe.md for coefficients that produce four-fermion or JW-string terms.

Outputs a list of problematic c_{mu,nu} and minimal suggestion (set to 0).
"""
import re
from pathlib import Path

path = Path(__file__).resolve().parents[1] / 'ybe.md'

text = ''
try:
    text = path.read_text(encoding='utf-8')
except Exception as e:
    print('Could not read ybe.md:', e)
    raise SystemExit(1)

# Keys that trigger four-fermion interactions
four_fermion = ['c_zz']
# Keys that trigger JW string / nonlocal terms
string_terms = ['c_x0','c_0x','c_y0','c_0y','c_xz','c_zx','c_yz','c_zy']

pattern = re.compile(r'c_([a-z0-9]+)\s*[:=]\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)')
found = {m.group(1): float(m.group(2)) for m in pattern.finditer(text)}

def report_group(keys):
    res = []
    for k in keys:
        kk = k.split('c_')[-1]
        val = found.get(kk)
        if val is None:
            # try full key
            val = found.get(k)
        res.append((k, val))
    return res

ff = report_group(four_fermion)
st = report_group(string_terms)

print('Parsed coefficients (from ybe.md) for problematic keys:')
any_prob = False
for k,v in ff:
    print(' -', k, ':', v)
    if v is not None and abs(v) > 1e-12:
        any_prob = True

for k,v in st:
    print(' -', k, ':', v)
    if v is not None and abs(v) > 1e-12:
        any_prob = True

if not any_prob:
    print('\nNo significant four-fermion or JW-string coefficients detected. Model is strictly quadratic under JW (within parsing limits).')
else:
    print('\nProblematic nonzero coefficients detected.')
    print('Minimal suggestion: set the following coefficients to 0 (or numerically ≈0):')
    for k,v in ff+st:
        if v is not None and abs(v) > 1e-12:
            print('  -', k, ' (current =', v, ') -> set to 0')
    print('\nAfter making these coefficients zero, re-run the mapping p = M·C and use BdG criteria (t,Δ,μ) to classify ABS vs MZM.')
