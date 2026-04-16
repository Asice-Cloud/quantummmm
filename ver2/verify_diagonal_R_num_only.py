import numpy as np
from math import pi

sx = np.array([[0,1],[1,0]],dtype=complex)
sy = np.array([[0,-1j],[1j,0]],dtype=complex)
sz = np.array([[1,0],[0,-1]],dtype=complex)
I4 = np.eye(4,dtype=complex)
A=np.kron(sx,sx); B=np.kron(sy,sy); C=np.kron(sz,sz)

def R_from_J(Jx,Jy,Jz, phase=0.0):
    Rx = np.cos(Jx)*I4 + 1j*np.sin(Jx)*A
    Ry = np.cos(Jy)*I4 + 1j*np.sin(Jy)*B
    Rz = np.cos(Jz)*I4 + 1j*np.sin(Jz)*C
    return np.exp(1j*phase) * (Rx @ Ry @ Rz)

# check function
def ybe(R):
    I2 = np.eye(2,dtype=complex)
    S = np.zeros((4,4),dtype=complex)
    for i in range(2):
        for j in range(2):
            S[i*2+j,j*2+i]=1
    R12 = np.kron(R,I2)
    R23 = np.kron(I2,R)
    R13 = np.kron(S,I2) @ R23 @ np.kron(S,I2)
    lhs = R12 @ R13 @ R23
    rhs = R23 @ R13 @ R12
    return np.allclose(lhs,rhs,atol=1e-9)

# isotropic checks
print('Isotropic checks (Jx=Jy=Jz):')
for J in [0.0, pi/4, pi/2, 3*pi/4, pi]:
    print(' J=', J, ' ->', ybe(R_from_J(J,J,J)))

# coarse grid 0..pi/2
grid = np.linspace(0, pi/2, 9)
sols = []
for Jx in grid:
    for Jy in grid:
        for Jz in grid:
            if ybe(R_from_J(Jx,Jy,Jz)):
                sols.append((Jx,Jy,Jz))
print('\nFound', len(sols), 'solutions on coarse grid:')
for s in sols[:50]:
    print(s)
