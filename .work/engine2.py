import numpy as np
from scipy.integrate import solve_ivp
I=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
def kron(a,b): return np.kron(a,b)
g1=kron(sx,I); g2=kron(-sy,I); ga=kron(sz,sx); gb=kron(sz,-sy); g3=ga@gb@g1@g2
G=[g1,g2,g3,ga,gb]
def H_EM(t,tau,E1,t1c,tc,E0,step_phase):
    def fm(x): return 0.5*(1-np.cos(np.pi*x/tau))
    def fp(x): return 0.5*(1+np.cos(np.pi*x/tau))
    if step_phase==0:
        t2=tc*fm(t); Ed=E0*fp(t); t3=0.0; t1=t1c*fm(t)
    elif step_phase==1:
        t3=tc*fm(t); t2=tc*fp(t); Ed=0.0; t1=t1c*fp(t)
    else:
        t3=tc*fp(t); Ed=E0*fm(t); t2=0.0; t1=0.0
    return 1j*( Ed*(G[3]@G[4]) + E1*(G[0]@G[1]) + abs(t2)*(G[3]@G[1]) - abs(t1)*(G[4]@G[0]) - abs(t3)*(G[3]@G[2]) )
def step_prop(tau,E1,t1c,tc,E0,phase,rtol=1e-10,atol=1e-12):
    # integrate dU/ds = -i H U , U(0)=I, s in [0,tau]
    def rhs(s,y):
        U=y.view(complex).reshape(4,4)
        H=H_EM(s,tau,E1,t1c,tc,E0,phase)
        return (-1j*H@U).reshape(-1).view(float)
    y0=np.eye(4,dtype=complex).reshape(-1).view(float)
    sol=solve_ivp(rhs,[0,tau],y0,method='DOP853',rtol=rtol,atol=atol,max_step=tau/40.0)
    U=sol.y[:,-1].view(complex).reshape(4,4)
    return U
def full_prop(tau,E1,t1c,tc,E0,rtol=1e-10,atol=1e-12):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2):
        U=step_prop(tau,E1,t1c,tc,E0,ph,rtol,atol)@U
    return U
def fid(U):
    return 0.5*sum(abs(U[1*2+a,0*2+b])**2 for a in range(2) for b in range(2))
if __name__=='__main__':
    tc=0.3;E0=0.3
    for tau_plot in [3,10]:
        tau=tau_plot*100.0
        print("MZM dbl tau_plot=%g fid=%.6f"%(tau_plot,fid(full_prop(tau,0,0,tc,E0))))
    U1=step_prop(300,0,0,tc,E0,0);U2=step_prop(300,0,0,tc,E0,1);U3=step_prop(300,0,0,tc,E0,2)
    print("MZM single fid=%.6f"%(fid(U3@U2@U1)))
    # compare to old engine at a non-MZM point
    print("E1=0.01 r=1 tau=800 fid=%.6f"%fid(full_prop(800,0.01,0.01,tc,E0)))
    print("E1=0.01 r=0.1 tau=800 fid=%.6f"%fid(full_prop(800,0.01,0.001,tc,E0)))
    print("E1=0.001 r=1 tau=1000 fid=%.6f"%fid(full_prop(1000,0.001,0.001,tc,E0)))
