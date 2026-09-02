import numpy as np
def q2c(q):
    q0,q1,q2,q3=q
    return np.array([[q0+1j*q1, q2+1j*q3],[-q2+1j*q3, q0-1j*q1]],complex)
def Hq(M):
    out=np.zeros((4,4),complex)
    for r in range(2):
        for c in range(2):
            out[2*r:2*r+2,2*c:2*c+2]=q2c(M[r][c])
    return out
Z=(0,0,0,0); I=(1,0,0,0); IM=(0,-1,0,0); JM=(0,0,-1,0); KM=(0,0,0,-1)
G1=Hq([[Z,I],[I,Z]])
G2=Hq([[Z,IM],[(-1)*np.array(IM,dtype=float)[0],Z]]) if False else None
# clean: [[0,-i],[i,0]]
G2=Hq([[Z,(0,-1,0,0)],[(0,1,0,0),Z]])
G3=Hq([[Z,(0,0,-1,0)],[(0,0,1,0),Z]])   # [[0,-j],[j,0]]
G4=Hq([[Z,(0,0,0,-1)],[(0,0,0,1),Z]])   # [[0,-k],[k,0]]
G5=Hq([[I,Z],[Z,(-1,0,0,0)]])
MOD=[G1,G2,G3,G4,G5]
mx=0
for i in range(5):
    for j in range(5):
        mx=max(mx,np.max(np.abs(MOD[i]@MOD[j]+MOD[j]@MOD[i]-2*(i==j)*np.eye(4))))
print("MOD complex Clifford err (want 0):",mx)
I2=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
g1=np.kron(sx,I2); g2=np.kron(-sy,I2); ga=np.kron(sz,sx); gb=np.kron(sz,-sy); g3=ga@gb@g1@g2
JW=[g1,g2,g3,ga,gb]
Aeq=[]
for i in range(5):
    M=np.kron(JW[i].T,np.eye(4))-np.kron(np.eye(4),MOD[i])
    Aeq.append(M)
M=np.vstack(Aeq); u,s,vh=np.linalg.svd(M)
ns=[vh[k] for k in range(len(s)) if s[k]<1e-7]
print("null dim:",len(ns), "(want 1)")
if len(ns)>=1:
    V=ns[0].reshape(4,4); V=V/np.linalg.norm(V)
    print("V unitary err:",np.max(np.abs(V@V.conj().T-np.eye(4))))
    # compare H_EM terms: JW physical gamma mapping gamma1=g1,gamma2=g2,gamma3=g3,gammaA=ga,gammaB=gb
    # model Gamma index: we must figure mapping. Try each permutation? Instead compare the OPERATOR algebra:
    # JW X12=i g1g2 ; if JW~MOD then i JW[i]JW[j] ~ i MOD[i]MOD[j]. Just verify V JW_i V^dag = MOD_i:
    ok=True
    for i in range(5):
        r=np.max(np.abs(V@JW[i]@V.conj().T-MOD[i]))
        if r>1e-6: ok=False
    print("V JW_i Vdag = MOD_i ?",ok)
