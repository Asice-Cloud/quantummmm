import numpy as np
I=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
def kron(a,b): return np.kron(a,b)
# qubit1=code (left), qubit2=ancilla (right)
g1 = kron(sx,I)
g2 = kron(-sy,I)          # so i*g1*g2 = sz(x)I (code occupancy)
ga = kron(sz,sx)
gb = kron(sz,-sy)
g3 = ga@gb@g1@g2          # chirality: anticommutes with all, squares +1 (to verify)
g = [g1,g2,g3,ga,gb]      # index: 0=gamma1,1=gamma2,2=gamma3,3=gammaA,4=gammaB
names=['g1','g2','g3','ga','gb']
# verify Clifford
err=0
for i in range(5):
    for j in range(5):
        M=g[i]@g[j]+g[j]@g[i]-2*(i==j)*np.eye(4)
        err=max(err,np.abs(M).max())
    hdiff=np.abs(g[i]-g[i].conj().T).max(); err=max(err,hdiff)
print("Clifford+hermitian max err:",err)
# i g1 g2 should be 2n1-1 acting on qubit1 only
print("i*g1*g2 diag:",np.round(np.diag(1j*g[0]@g[1]),4))
print("i*ga*gb diag:",np.round(np.diag(1j*g[3]@g[4]),4))
