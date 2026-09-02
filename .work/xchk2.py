import numpy as np
def qmul(a,b):
    a0,a1,a2,a3=a; b0,b1,b2,b3=b
    return (a0*b0-a1*b1-a2*b2-a3*b3,
            a0*b1+a1*b0+a2*b3-a3*b2,
            a0*b2-a1*b3+a2*b0+a3*b1,
            a0*b3+a1*b2-a2*b1+a3*b0)
def qadd(a,b): return tuple(x+y for x,y in zip(a,b))
def qscale(s,a): return tuple(s*x for x in a)
Iq=(1,0,0,0); ip=(0,1,0,0); jp=(0,0,1,0); kp=(0,0,0,1); Z=(0,0,0,0)
def mm(A,B):
    return [[qadd(qmul(A[0][0],B[0][0]),qmul(A[0][1],B[1][0])), qadd(qmul(A[0][0],B[0][1]),qmul(A[0][1],B[1][1]))],
            [qadd(qmul(A[1][0],B[0][0]),qmul(A[1][1],B[1][0])), qadd(qmul(A[1][0],B[0][1]),qmul(A[1][1],B[1][1]))]]
def mscale(s,M): return [[qscale(s,x) for x in row] for row in M]
def madd(M,N): return [[qadd(a,b) for a,b in zip(r1,r2)] for r1,r2 in zip(M,N)]
def ident(): return [[Iq,Z],[Z,Iq]]
def Kblk(E1,t2,t1,t3,Ed):
    t2m=abs(t2); t1m=abs(t1); t3m=abs(t3)
    A=qadd(qscale((E1+t3m)/2,ip), qscale(t2m/2,jp))
    D=qadd(qscale((-E1+t3m)/2,ip), qscale(t2m/2,jp))
    B=qadd(qscale(t1m/2,Iq), qscale(Ed/2,kp))
    C=qadd(qscale(-t1m/2,Iq), qscale(Ed/2,kp))
    return [[A,B],[C,D]]
def integrate_quat(tau,E1,t1c,tc,E0,phase,n=1000):
    def fm(s): return 0.5*(1-np.cos(np.pi*s/tau))
    def fp(s): return 0.5*(1+np.cos(np.pi*s/tau))
    def blocks_at(s,phase):
        if phase==0: t2=tc*fm(s); Ed=E0*fp(s); t3=0.0; t1=t1c*fm(s)
        elif phase==1: t3=tc*fm(s); t2=tc*fp(s); Ed=0.0; t1=t1c*fp(s)
        else: t3=tc*fp(s); Ed=E0*fm(s); t2=0.0; t1=0.0
        return Kblk(E1,t2,t1,t3,Ed)
    U=ident(); dt=tau/n
    for k in range(n):
        s=k*dt
        K=blocks_at(s,phase); K2=blocks_at(s+dt/2,phase); K4=blocks_at(s+dt,phase)
        k1=mm(K,U)
        k2=mm(K2,madd(U,mscale(dt/2,k1)))
        k3=mm(K2,madd(U,mscale(dt/2,k2)))
        k4=mm(K4,madd(U,mscale(dt,k3)))
        U=madd(U,mscale(dt/6, madd(k1,madd(mscale(2,k2),madd(mscale(2,k3),k4)))))
    return U
def quat_to_complex(M):
    out=np.zeros((4,4),complex)
    for r in range(2):
        for c in range(2):
            q=M[r][c]
            out[2*r:2*r+2,2*c:2*c+2]=np.array([[q[0]+1j*q[1],q[2]+1j*q[3]],[-q[2]+1j*q[3],q[0]-1j*q[1]]])
    return out
if __name__=="__main__":
    U=integrate_quat(300,0.0,0.0,0.3,0.3,0)
    Uc=quat_to_complex(U)
    print("step1 quat U unitary err:",np.max(np.abs(Uc@Uc.conj().T-np.eye(4))))
    print("eigs:",np.round(np.linalg.eigvals(Uc),5))
