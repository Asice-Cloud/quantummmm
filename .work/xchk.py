import numpy as np
I2=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
g1=np.kron(sx,I2); g2=np.kron(-sy,I2); ga=np.kron(sz,sx); gb=np.kron(sz,-sy); g3=ga@gb@g1@g2
JW=[g1,g2,g3,ga,gb]
# Gamma's as complex, from model.md, but correct Clifford: build via a known faithful rep:
# Use the 4x4 rep: Gamma_i = i-th of [sx,I;...] ... just reuse JW but permuted? The point: are the 
# model Gamma's (as 2x2 quaternion) unitarily equivalent to some JW ordering? Test on a canonical:
# Standard 4x4 rep of Cl(5) with explicit gamma1..5 (complex) e.g. from Fock:
A=np.array([[0,1,0,0],[-1,0,0,0],[0,0,0,-1],[0,0,1,0]],complex)  # i? 
# skip - use JW directly; question is fidelity equivalence between quaternion U and JW, 
# which report already established to 1e-9. Here just do a pragmatic check:
# evolve quaternion U via model's OWN blocks in pure quaternion arithmetic and compare
# eigenvalues of total propagator with JW.

# --- quaternion helpers (Hamilton product), q=(q0,q1,q2,q3) ---
def qmul(a,b):
    a0,a1,a2,a3=a; b0,b1,b2,b3=b
    return (a0*b0-a1*b1-a2*b2-a3*b3,
            a0*b1+a1*b0+a2*b3-a3*b2,
            a0*b2-a1*b3+a2*b0+a3*b1,
            a0*b3+a1*b2-a2*b1+a3*b0)
def qadd(a,b): return tuple(x+y for x,y in zip(a,b))
def qscale(s,a): return tuple(s*x for x in a)
def qconj(a): return (a[0],-a[1],-a[2],-a[3])
Iq=(1,0,0,0); ip=(0,1,0,0); jp=(0,0,1,0); kp=(0,0,0,1)
# 2x2 quaternion mat mult C=A@B (each entry quat)
def mm(A,B):
    return [[qadd(qmul(A[0][0],B[0][0]),qmul(A[0][1],B[1][0])), qadd(qmul(A[0][0],B[0][1]),qmul(A[0][1],B[1][1]))],
            [qadd(qmul(A[1][0],B[0][0]),qmul(A[1][1],B[1][0])), qadd(qmul(A[1][0],B[0][1]),qmul(A[1][1],B[1][1]))]]
def mscale(s,M): return [[qscale(s,x) for x in row] for row in M]
def madd(M,N): return [[qadd(a,b) for a,b in zip(r1,r2)] for r1,r2 in zip(M,N)]
def ident(): return [[Iq,(0,0,0,0)],[(0,0,0,0),Iq]]
def pure(v):  # vector->pure quaternion
    return (0,)+tuple(v)
# model.md blocks (with factor: U_dot=KU, K=sum h_ij Sigma_ij, Sigma entries given in table)
# From Sigma table: Sigma12=diag(ip/2,-ip/2) i.e. A-block entry for E1:
# We need K = h12*Sigma12 + h24*Sigma24 + ...  with h from H_EM coefficients.
# H_EM = E1 X12 + |t2| X24 -|t1| X15 -|t3| X34 + Ed X45 ; X_ij=2 Sigma_ij ; K=sum h_ij Sigma_ij
# => h12=E1, h24=|t2|, h15=-|t1|, h34=-|t3|, h45=Ed
# Sigma12 = diag(ip/2, -ip/2); Sigma24=diag(jp/2, jp/2); Sigma34=diag(-ip/2,-ip/2)? from table:
# Sigma34(-|t3|)= [[-i/2,0],[0,-i/2]]
# Sigma15(-|t1|)= [[0,-1/2],[1/2,0]]  (real offdiag)
# Sigma45(Ed)= [[0,k/2],[k/2,0]]
def K_from_params(E1,t2,t1,t3,Ed):
    # K = E1*S12 + t2*S24 + (-t1)*S15 + (-t3)*S34 + Ed*S45  (t's are |t| magnitudes already signed as in H_EM)
    A11=qadd(qscale(E1/2,ip), qscale(t2/2,jp)); # E1 ip/2 + t2 jp/2 ... but S34 also -ip/2 on diag both
    A11=qadd(qscale(E1/2,ip), qscale(t2/2,jp)); A11=qadd(A11, qscale(-t3/2,ip))
    D22=qadd(qscale(-E1/2,ip), qscale(t2/2,jp)); D22=qadd(D22, qscale(-t3/2,ip))
    B12=( -t1/2,0,0,0)  # Sigma15 offdiag -1/2 real
    B12=qadd(B12, qscale(Ed/2,kp))
    C21=( t1/2,0,0,0)
    C21=qadd(C21, qscale(Ed/2,kp))
    return [[A11,B12],[C21,D22]]
# Note: this should match model.md explicit A,B,C,D (they absorb different factor 2; check by unitarity after integrate)
print("done defs")
