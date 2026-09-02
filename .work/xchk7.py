import numpy as np
# null-space solve for V: V @ JW_i = MOD_i @ V  for i=0..4  (16 unknowns, 80 eqs)
def q2c(q):
    q0,q1,q2,q3=q
    return np.array([[q0+1j*q1, q2+1j*q3],[-q2+1j*q3, q0-1j*q1]],complex)
def Hq(M):
    out=np.zeros((4,4),complex)
    for r in range(2):
        for c in range(2):
            out[2*r:2*r+2,2*c:2*c+2]=q2c(M[r][c])
    return out
Z=(0,0,0,0); I=(1,0,0,0); IP=(0,1,0,0); JP=(0,0,1,0); KP=(0,0,0,1)
G1=Hq([[Z,I],[I,Z]])
G2=Hq([[Z,(0,0,-1,0)],[(0,0,1,0),Z]])  # [[0,-j],[j,0]]
G3=Hq([[Z,(0,0,0,-1)],[(0,0,0,1),Z]])  # [[0,-k],[k,0]]
G4=Hq([[Z,(0,0,0,-1)],[(0,0,0,1),Z]])  # [[0,-k],[k,0]]  (model Gamma4 = [[0,-k],[k,0]])
G5=Hq([[I,Z],[Z,(-1,0,0,0)]])          # diag(1,-1)
MOD=[G1,G2,G3,G4,G5]
# Verify the complexified Gammas anticommute (sanity that embedding is a homomorphism)
mx=0
for i in range(5):
    for j in range(5):
        mx=max(mx,np.max(np.abs(MOD[i]@MOD[j]+MOD[j]@MOD[i]-2*(i==j)*np.eye(4))))
print("MOD complex Clifford err:",mx)
I2=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
g1=np.kron(sx,I2); g2=np.kron(-sy,I2); ga=np.kron(sz,sx); gb=np.kron(sz,-sy); g3=ga@gb@g1@g2
JW=[g1,g2,g3,ga,gb]
# Build linear system for V (complex, 16 unknowns as vec(V))
Aeq=[]
for i in range(5):
    # V JW_i - MOD_i V = 0  => (JW_i^T kron I - I kron MOD_i) vec(V)=0?  careful: V A - B V =0
    # vec(B V) = (I kron B) vec(V); vec(V A)= (A^T kron I) vec(V). So (A^T kron I - I kron B) vec V =0
    A=JW[i]; B=MOD[i]
    M=np.kron(A.T,np.eye(4)) - np.kron(np.eye(4),B)
    Aeq.append(M)
M=np.vstack(Aeq)  # 80 x 16
# solve nullspace
u,s,vh=np.linalg.svd(M)
print("singular values (should have small ones):",np.round(s,6))
ns=[vh[k] for k in range(len(s)) if s[k]<1e-8]
print("null dim:",len(ns))
if len(ns)>=1:
    V=ns[0].reshape(4,4)
    V=V/np.linalg.norm(V)
    # check unitary
    print("V unitary err:",np.max(np.abs(V@V.conj().T-np.eye(4))))
    # check V JW_i = MOD_i V
    for i in range(5):
        print(" i=%d resid:"%i,np.max(np.abs(V@JW[i]-MOD[i]@V)))
