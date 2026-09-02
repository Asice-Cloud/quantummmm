import numpy as np
def q2c(q):  # q=(q0,q1,q2,q3)
    q0,q1,q2,q3=q
    return np.array([[q0+1j*q1, q2+1j*q3],[-q2+1j*q3, q0-1j*q1]],complex)
def Hq(M):  # 2x2 quaternion matrix -> 4x4 complex (block)
    out=np.zeros((4,4),complex)
    for r in range(2):
        for c in range(2):
            out[2*r:2*r+2,2*c:2*c+2]=q2c(M[r][c])
    return out
Z=(0,0,0,0); I=(1,0,0,0); IP=(0,1,0,0); JP=(0,0,1,0); KP=(0,0,0,1)
G1=Hq([[Z,I],[I,Z]])
G2=Hq([[Z,(0,0,-1,0)],[(0,0,1,0),Z]])   # [[0,-j],[j,0]]
G3=Hq([[Z,(0,0,0,-1)],[(0,0,0,1),Z]])   # [[0,-k],[k,0]]
G4=Hq([[Z,(0,0,-1,0)],[(0,0,1,0),Z]]) is None
# careful rebuild
G1=Hq([[Z,I],[I,Z]])
G2=Hq([[Z,(0,0,-1,0)],[(0,0,1,0),Z]])
G3=Hq([[Z,(0,0,0,-1)],[(0,0,0,1),Z]])
G4=Hq([[Z,(0,0,0,-1)],[(0,0,0,1),Z]])  # [[0,-k],[k,0]]
G5=Hq([[I,Z],[Z,(-1,0,0,0)]])          # diag(1,-1)
Gam=[G1,G2,G3,G4,G5]
mx=0
for i in range(5):
    for j in range(5):
        mx=max(mx,np.max(np.abs(Gam[i]@Gam[j]+Gam[j]@Gam[i]-2*(i==j)*np.eye(4))))
print("model Gamma complex Clifford max err:",mx)
# JW
I2=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
g1=np.kron(sx,I2); g2=np.kron(-sy,I2); ga=np.kron(sz,sx); gb=np.kron(sz,-sy); g3=ga@gb@g1@g2
JW=[g1,g2,g3,ga,gb]
# Are the two reps equivalent? Compare X_ij = i g_i g_j spectra & trace of all products
def prof(gs):
    return sorted(round(np.trace(np.linalg.matrix_power(1j*gs[i]@gs[j],2)).real,6) for i in range(5) for j in range(i+1,5))
print("JW  profile:",prof(JW))
print("MOD profile:",prof(Gam))
# try find V via common Cartan: diagonalize X12 & X34 & X45 in each rep
import numpy.linalg as la
def V_from_pair(gs):
    A=1j*gs[0]@gs[1]; B=1j*gs[2]@gs[3]
    # find ordered eigenbases
    va=la.eigh(A)[1]; vb=la.eigh(B)[1]
    # crude: V maps eigenbasis of A in JW -> eigenbasis of A in MOD aligned by commuting B sector? 
    return None
# Test equivalence via invariant: the center-like element gamma1..gamma5 (chirality)
def chi(gs):
    return np.linalg.matrix_power(np.linalg.multi_dot(gs),1)
cJW=chi(JW); cMOD=chi(Gam)
print("JW  chirality trace/offdiag:",np.trace(cJW),np.max(np.abs(cJW-cJW.conj().T)))
print("MOD chirality trace/offdiag:",np.trace(cMOD),np.max(np.abs(cMOD-cMOD.conj().T)))
print("JW  chirality^2:",(cJW@cJW)[0,0])
print("MOD chirality^2:",(cMOD@cMOD)[0,0])
