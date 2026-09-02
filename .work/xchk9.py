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
Z=(0,0,0,0); I=(1,0,0,0)
G1=Hq([[Z,I],[I,Z]])
G2=Hq([[Z,(0,-1,0,0)],[(0,1,0,0),Z]])
G3=Hq([[Z,(0,0,-1,0)],[(0,0,1,0),Z]])
G4=Hq([[Z,(0,0,0,-1)],[(0,0,0,1),Z]])
G5=Hq([[I,Z],[Z,(-1,0,0,0)]])
MOD=[G1,G2,G3,G4,G5]
I2=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
g1=np.kron(sx,I2); g2=np.kron(-sy,I2); ga=np.kron(sz,sx); gb=np.kron(sz,-sy); g3=ga@gb@g1@g2
JW=[g1,g2,g3,ga,gb]
Aeq=[]
for i in range(5):
    Aeq.append(np.kron(JW[i].T,np.eye(4))-np.kron(np.eye(4),MOD[i]))
M=np.vstack(Aeq); u,s,vh=np.linalg.svd(M)
V=vh[np.argmin(s)].reshape(4,4)
# unitarize via polar
V=V/np.sqrt(V@V.conj().T)
# polar: V = U H ; unitary part U = V (Vdag V)^{-1/2}
def sqrtm_inv(A):
    w,v=np.linalg.eigh(A); return (v*np.sqrt(w))@v.conj().T
U=V@np.linalg.inv(sqrtm_inv(V.conj().T@V))
print("U unitary err:",np.max(np.abs(U@U.conj().T-np.eye(4))))
ok=True
for i in range(5):
    r=np.max(np.abs(U@JW[i]@U.conj().T-MOD[i]))
    if r>1e-8: ok=False; print(" term i=%d mismatch"%i,r)
print("V JW_i Vdag = MOD_i ?",ok)
# Now map physical H_EM operator: physical gamma mapping in JW = (g1,g2,g3,ga,gb) for (gamma1,gamma2,gamma3,gammaA,gammaB)
# H_EM = i Ed gammaA gammaB + i E1 g1g2 + i|t2| gammaA gamma2 - i|t1| gammaB gamma1 - i|t3| gammaA gamma3
# In MOD basis gamma_i->MOD_i ; need mapping which MOD index = gammaA? 
# Determine: compute X_ij=i g_i g_j JW vs MOD and compare term by term to deduce index permutation:
def find(idxj, Xjw):
    # Xjw should equal i MOD[a]MOD[b] for the pair (a,b) that corresponds
    cands=[]
    for a in range(5):
        for b in range(a+1,5):
            Xmod=1j*MOD[a]@MOD[b]
            if np.max(np.abs(Xjw-Xmod))<1e-8: cands.append((a+1,b+1))
            Xmodm=1j*MOD[b]@MOD[a]
            if np.max(np.abs(Xjw-Xmodm))<1e-8: cands.append((b+1,a+1))
    return cands
names=['g1g2','g2g3','g3g4','g4g5']  # index into JW by physical: X12,X23,...
# physical pairs of interest with JW indices 0,1,2,3,4 = gamma1,gamma2,gamma3,gammaA,gammaB
for lbl,i,j in [('X12',0,1),('X23',1,2),('X1A',0,3),('X2A',1,3),('X1B',0,4),('X3A',2,3),('XAB',3,4)]:
    Xjw=1j*JW[i]@JW[j]
    c=find(0,Xjw)
    print(lbl,"-> MOD pair",c)
