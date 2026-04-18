#!/usr/bin/env python3
"""Generate representative samples (equality family, F_- and F_+ branches),
compute three-site YBE residual norms and print/save results.
"""
import numpy as np
I = 1j

# Pauli matrices
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-I],[I,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

Sx = np.kron(sx, sx)
Sy = np.kron(sy, sy)
Sz = np.kron(sz, sz)
I4 = np.eye(4, dtype=complex)

def R_from_params(c00,bx,by,bz):
    Rx = np.cos(bx)*I4 + I*np.sin(bx)*Sx
    Ry = np.cos(by)*I4 + I*np.sin(by)*Sy
    Rz = np.cos(bz)*I4 + I*np.sin(bz)*Sz
    return np.exp(I*c00) * Rx.dot(Ry).dot(Rz)

def three_B(R):
    R12 = np.kron(R, I2)
    R23 = np.kron(I2, R)
    return R12.dot(R23).dot(R12) - R23.dot(R12).dot(R23)

def norm(B):
    return np.linalg.norm(B)

out_lines = []

# Sample 1: equality family bx=by=bz
c00 = 0.0
bx = by = bz = 0.5
R = R_from_params(c00,bx,by,bz)
B = three_B(R)
out_lines.append(('equality', c00, bx, by, bz, norm(B)))

# Sample 2: F_- branch. pick bx,by then Z=(XY-1)/(X+Y)
bx = 0.2; by = 0.5
X = np.exp(2j*bx); Y = np.exp(2j*by)
Z = (X*Y - 1)/(X + Y)
# choose principal log; allow complex bz
bz = -0.5j * np.log(Z)
R = R_from_params(0.0, bx, by, bz)
B = three_B(R)
out_lines.append(('F_minus', 0.0, bx, by, bz, norm(B)))

# Sample 3: F_+ branch. pick bx,by then Z=(1-XY)/(X+Y)
bx = 0.3; by = 0.6
X = np.exp(2j*bx); Y = np.exp(2j*by)
Z = (1 - X*Y)/(X + Y)
bz = -0.5j * np.log(Z)
R = R_from_params(0.0, bx, by, bz)
B = three_B(R)
out_lines.append(('F_plus', 0.0, bx, by, bz, norm(B)))

# Sample 4: random non-solution
c00 = 0.0
bx, by, bz = 0.11, 0.27, 0.93
R = R_from_params(c00,bx,by,bz)
B = three_B(R)
out_lines.append(('random', c00, bx, by, bz, norm(B)))

# Print and save to file
print('Validation results (name, c00, bx, by, bz, ||B||):')
for t,c00,bx,by,bz,n in out_lines:
    print(f'{t}: c00={c00}, bx={bx}, by={by}, bz={bz} -> ||B|| = {n:.6e}')

with open('scripts/ybe_exp_validation.txt','w') as f:
    f.write('Validation results (name, c00, bx, by, bz, ||B||):\n')
    for t,c00,bx,by,bz,n in out_lines:
        f.write(f'{t}: c00={c00}, bx={bx}, by={by}, bz={bz} -> ||B|| = {n:.6e}\n')

print('\nWrote scripts/ybe_exp_validation.txt')
