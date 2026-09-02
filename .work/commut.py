import numpy as np
I=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
def kron(a,b): return np.kron(a,b)
g1=kron(sx,I); g2=kron(-sy,I); ga=kron(sz,sx); gb=kron(sz,-sy); g3=ga@gb@g1@g2
G=[g1,g2,g3,ga,gb]
# term generators in H_EM (model.md sec1.1):
# E1 term: i g1 g2 ; t1 term: -i gb g1 = -i G[4]G[0]; Ed: i ga gb; t2: i ga g2; t3: -i ga g3
X12=1j*G[0]@G[1]; Xb1=-1j*G[4]@G[0]
print("E1 gen X12=i g1 g2\n",np.round(X12,4))
print("t1 gen Xb1=-i gb g1\n",np.round(Xb1,4))
print("Hermitian:",np.allclose(X12,X12.conj().T),np.allclose(Xb1,Xb1.conj().T))
print("tr(X12^2)/4:",np.trace(X12@X12)/4)  # should be -1 (i*... )^2... check
print("anticommute {X12,Xb1}=X12Xb1+Xb1X12:\n",np.round(X12@Xb1+Xb1@X12,4))
C=X12@Xb1-Xb1@X12
print("[X12,Xb1]:\n",np.round(C,4))
# code-subspace restriction: states of gamma1,gamma2 pair = n1 (code), ancilla pair = n2
# X12 diagonalizes in n1: i g1 g2 = kron(sz?,...). compute eigenvalues/vectors
w,v=np.linalg.eigh(X12)
print("eig X12:",np.round(w,4),"-> +1 acts on n1=0 (psi1+ code), -1 on n1=1? map")
# Actually check: i g1 g2 = kron(sz,I) per our JW?  g1=kron(sx,I),g2=kron(-sy,I): g1 g2=kron(sx*(-sy),I)=kron(-i sz,I)... i*g1g2=kron(sz,I)
print("X12 == kron(sz,I):",np.allclose(X12,kron(sz,I)))
print("Xb1 == ?", np.round(Xb1,3))
# t1 term structure: -i gb g1 ; gb=kron(sz,-sy),g1=kron(sx,I) -> gb g1=kron(sz*sx,-sy)=kron(i sy? , ...) sz*sx=i sy? no: sz*sx = [[0,1],[-1,0]] = -i sy. 
# gb g1 = kron(-i sy, -sy)?? sz*sx=-i*sy; times... let's compute: sz@sx=np.array([[0,1],[-1,0]])=-1j*sy
# so gb g1 = kron(-1j sy, -sy); Xb1=-i gb g1 = -i*kron(-i sy,-sy)= kron(-i*-i sy, -sy)? = kron((-1)*sy,-sy)= kron(-sy,-sy)?? compute numerically
# it should couple n1<->ancilla (block off-diag in code): off-block means mixes code n1 flip with ancilla
Bmat=np.kron(sz,sy)  # reference
print("Xb1 resembles kron(sy,sy)?",np.allclose(Xb1,kron(sy,sy)))
print("Xb1 resembles kron(sy,I)?",np.allclose(Xb1,kron(sy,I)))
# Effective code 2x2: X12 restricted to +1 eigenspace (n1=0) = identity; in n1=0 block X12=+1 (diagonal z-type). 
# Show X12 acts as sigma_z (diag +1,-1) on code index n1, ancilla untouched:  i g1g2 = kron(sz,I) confirmed -> in |n1,n2> basis it is sz on code.
# Xb1 = kron(?,?) let's print as 2x2 blocks over code index
print("Xb1 2x2 block over n1:\n",[[np.round(Xb1[0:2,0:2],3),np.round(Xb1[0:2,2:4],3)],[np.round(Xb1[2:4,0:2],3),np.round(Xb1[2:4,2:4],3)]])
