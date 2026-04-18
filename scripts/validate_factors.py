#!/usr/bin/env python3
import numpy as np

I = 1j

sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-I],[I,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

def kron(a,b):
    return np.kron(a,b)

Sx = kron(sx,sx)
Sy = kron(sy,sy)
Sz = kron(sz,sz)
I4 = np.eye(4, dtype=complex)

def R_from_params(c00,bx,by,bz):
    # R = e^{i c00} prod_alpha (cos b + i sin b S_alpha)
    Rx = np.cos(bx)*I4 + I*np.sin(bx)*Sx
    Ry = np.cos(by)*I4 + I*np.sin(by)*Sy
    Rz = np.cos(bz)*I4 + I*np.sin(bz)*Sz
    R = np.exp(I*c00) * Rx.dot(Ry).dot(Rz)
    return R

def three_site_B(R):
    I2 = np.eye(2, dtype=complex)
    R12 = np.kron(R, I2)
    R23 = np.kron(I2, R)
    B = R12.dot(R23).dot(R12) - R23.dot(R12).dot(R23)
    return B

def norm(B):
    return np.linalg.norm(B)

cases = []
# Case A: bx = by (sin(bx-b y)=0)
cases.append((0.1, 0.5, 0.5, 0.2))
# Case B: bx = -by (sin(bx+by)=0)
cases.append((0.2, 0.4, -0.4, 0.3))
# Case C: F_-: Z = (XY -1)/(X+Y)
bx = 0.2; by = 0.5
X = np.exp(2j*bx); Y = np.exp(2j*by)
Z = (X*Y - 1)/(X + Y)
bz = -0.5j * np.log(Z)
cases.append((0.0, bx, by, bz))
# Case D: F_+: Z = (1 - XY)/(X+Y)
bx = 0.3; by = 0.6
X = np.exp(2j*bx); Y = np.exp(2j*by)
Z = (1 - X*Y)/(X + Y)
bz = -0.5j * np.log(Z)
cases.append((0.0, bx, by, bz))
# Case E: random non-solution
cases.append((0.0, 0.11, 0.27, 0.93))

print('Running validation samples:')
for i,(c00,bx,by,bz) in enumerate(cases,1):
    R = R_from_params(c00,bx,by,bz)
    B = three_site_B(R)
    n = norm(B)
    print(f'Case {i}: c00={c00}, bx={bx}, by={by}, bz={bz} -> ||B|| = {n:.3e}')

print('\nInterpretation: near-zero norm indicates YBE satisfied (within numeric precision).')

# Quick assertions
print('\nAsserting near-zero for first four (solutions) and non-zero for last:')
for i,(c00,bx,by,bz) in enumerate(cases,1):
    R = R_from_params(c00,bx,by,bz)
    B = three_site_B(R)
    n = norm(B)
    if i <= 4:
        assert n < 1e-8, f"Case {i} expected solution but residual {n}"
    else:
        assert n > 1e-6, f"Case {i} expected non-solution but residual {n}"
print('All checks passed.')
