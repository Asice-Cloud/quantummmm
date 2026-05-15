import numpy as np
from scipy.optimize import linprog
import math, os
I2 = np.array([[1,0],[0,1]], dtype=complex)
SX = np.array([[0,1],[1,0]], dtype=complex)
SY = np.array([[0,-1j],[1j,0]], dtype=complex)
SZ = np.array([[1,0],[0,-1]], dtype=complex)

def kron(a,b): return np.kron(a,b)

def h4(g1,g2,g3,g4,tc=1.0):
    return tc * (-(g1*g2)*kron(SZ,I2) - (g1*g3)*kron(SY,SX) + (g1*g4)*kron(SY,SY) - (g2*g3)*kron(SX,SX) - (g2*g4)*kron(SX,SY) - (g3*g4)*kron(I2,SZ))

def path_segment(k,u):
    if k==1: return u,0.0,1.0-u,1.0
    if k==2: return 1.0-u,u,0.0,1.0
    if k==3: return 0.0,1.0-u,u,1.0
    raise ValueError('bad')

N_per=200
dpoints=[]
for k in [1,2,3]:
    us = np.linspace(0.0,1.0,N_per,endpoint=False)
    for u in us:
        g1,g2,g3,g4 = path_segment(k,u)
        H4 = h4(g1,g2,g3,g4,tc=1.0)
        psi01 = np.array([0,1,0,0], dtype=complex)
        psi10 = np.array([0,0,1,0], dtype=complex)
        Heff = np.array([[psi01.conj().T@H4@psi01, psi01.conj().T@H4@psi10],
                         [psi10.conj().T@H4@psi01, psi10.conj().T@H4@psi10]], dtype=complex)
        d0 = 0.5*np.trace(Heff)
        dx = 0.5*np.trace(Heff @ SX)
        dy = 0.5*np.trace(Heff @ SY)
        dz = 0.5*np.trace(Heff @ SZ)
        dpoints.append([float(dx.real), float(dy.real), float(dz.real)])

dpoints=np.array(dpoints)
min_norm = np.min(np.linalg.norm(dpoints, axis=1))
# LP check
N = dpoints.shape[0]
c = np.zeros(N)
A_eq = np.vstack([dpoints.T, np.ones(N)])
b_eq = np.zeros(4)
b_eq[-1] = 1.0
bounds = [(0.0, None) for _ in range(N)]
res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
print('min|d|=', min_norm)
print('origin_in_convex_hull:', bool(res.success))
