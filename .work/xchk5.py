import numpy as np, sys
sys.path.insert(0,'.work')
from xchk2 import qmul,qadd,qscale,mm,mscale,madd,ident,Kblk,quat_to_complex
from dyn import step_prop
tc=0.3; E0=0.3
def quat_step(tau,E1,t1c,tc,E0,phase,factor=2.0,n=1500):
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
        k1=mm(mscale(factor,K),U)
        k2=mm(mscale(factor,K2),madd(U,mscale(dt/2,k1)))
        k3=mm(mscale(factor,K2),madd(U,mscale(dt/2,k2)))
        k4=mm(mscale(factor,K4),madd(U,mscale(dt,k3)))
        U=madd(U,mscale(dt/6,madd(k1,madd(mscale(2,k2),madd(mscale(2,k3),k4)))))
    return quat_to_complex(U)
def quat_total(tau,E1,t1c,factor=2.0):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2):
        U=quat_step(tau,E1,t1c,tc,E0,ph,factor)@U
    return U
def jw_total(tau,E1,t1c):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2):
        U=step_prop(tau,E1,t1c,tc,E0,ph,n=3000)@U
    return U
tau=6*100.0
print("protocol double-braid, EIGENVALUE cross-check (U_dot=2KU vs JW):")
for E1,t1c in [(0,0),(0.01,0.001),(0.01,0.3),(0.001,0.01),(0.005,0.05)]:
    a=np.sort(np.linalg.eigvals(jw_total(tau,E1,t1c)))
    b=np.sort(np.linalg.eigvals(quat_total(tau,E1,t1c)))
    # phase-align
    d=min(np.max(np.abs(a-b)), np.max(np.abs(a+b)))
    print(" E1=%.3f t1=%.3f : max eig-diff(module sign)=%.2e"%(E1,t1c,d))
