import numpy as np, sys
sys.path.insert(0,'.work')
from xchk2 import qmul,qadd,qscale,mm,mscale,madd,ident,Kblk,quat_to_complex
from dyn import H_EM,G
# constant H: E1=0.01,t2=0.3 (fm=1 mid), t1=0.05, t3=0, Ed=0 ; evolve time T
def jw_const(T,E1,t2,t1,t3,Ed,n=4000):
    H=H_EM(1.0,3.0,E1,t1,t2,Ed,0)  # t inside step; fm(1)=0.5? No - need truly constant: build directly
    # Actually H_EM for phase 0 at s=1.5 of tau=3: fm=0.5(1-cos(pi/2))=0.5, fp=0.5. Not constant. Build H manually:
    I=np.eye(2,dtype=complex)
    # reuse: construct with tc=2*t2 so fm gives t2... too hacky. Build directly:
    I2=np.eye(2,dtype=complex)
    H=1j*(Ed*(G[3]@G[4])+E1*(G[0]@G[1])+abs(t2)*(G[3]@G[1])-abs(t1)*(G[4]@G[0])-abs(t3)*(G[3]@G[2]))
    # exact exp via eigendecomposition (H hermitian)
    w,V=np.linalg.eigh(H)
    U=(V*np.exp(-1j*w*T))@V.conj().T
    return U
def quat_const(T,E1,t2,t1,t3,Ed,factor,n=4000):
    K=Kblk(E1,t2,t1,t3,Ed)
    U=ident(); dt=T/n
    for k in range(n):
        K2=K
        k1=mscale(factor,K2); k1=mm(k1,U)
        k2=mm(mscale(factor,K2), madd(U,mscale(dt/2,k1)))
        k3=mm(mscale(factor,K2), madd(U,mscale(dt/2,k2)))
        k4=mm(mscale(factor,K2), madd(U,mscale(dt,k3)))
        U=madd(U,mscale(dt/6,madd(k1,madd(mscale(2,k2),madd(mscale(2,k3),k4)))))
    return quat_to_complex(U)
T=300.0; E1=0.01; t2=0.3; t1=0.0; t3=0.0; Ed=0.0
for factor in [1,2]:
    Uq=quat_const(T,E1,t2,t1,t3,Ed,factor)
    Uj=jw_const(T,E1,t2,t1,t3,Ed)
    print("factor=%d: max |Uq-Uj| ="%factor, np.max(np.abs(Uq-Uj)))
    print("  eig JW :",np.round(np.sort(np.linalg.eigvals(Uj)),4))
    print("  eig Q  :",np.round(np.sort(np.linalg.eigvals(Uq)),4))
