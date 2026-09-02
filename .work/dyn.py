import numpy as np
I=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
def kron(a,b): return np.kron(a,b)
g1=kron(sx,I); g2=kron(-sy,I); ga=kron(sz,sx); gb=kron(sz,-sy); g3=ga@gb@g1@g2
G=[g1,g2,g3,ga,gb]  # 0,1,2 = gamma1,gamma2,gamma3 ; 3,4 = gammaA,gammaB
def H_EM(t,tau,E1,t1c,tc,E0,step_phase=0):
    # local time within step: s in [0,tau]
    s = t
    def fm(x): return 0.5*(1-np.cos(np.pi*x/tau))
    def fp(x): return 0.5*(1+np.cos(np.pi*x/tau))
    # step_phase: 0=step1,1=step2,2=step3
    if step_phase==0:
        t2=tc*fm(s); Ed=E0*fp(s); t3=0.0; t1=t1c*fm(s)
    elif step_phase==1:
        t3=tc*fm(s); t2=tc*fp(s); Ed=0.0; t1=t1c*fp(s)
    else:
        t3=tc*fp(s); Ed=E0*fm(s); t2=0.0; t1=0.0
    H = 1j*( Ed*(G[3]@G[4]) + E1*(G[0]@G[1]) + abs(t2)*(G[3]@G[1]) - abs(t1)*(G[4]@G[0]) - abs(t3)*(G[3]@G[2]) )
    return H
def step_prop(tau, E1, t1c, tc, E0, phase, n=2000):
    s0=0.0; dt=tau/n
    U=np.eye(4,dtype=complex)
    for k in range(n):
        s=s0+k*dt
        def H(x): return H_EM(x,tau,E1,t1c,tc,E0,phase)
        k1=-1j*H(s)@U
        k2=-1j*H(s+dt/2)@(U+dt/2*k1)
        k3=-1j*H(s+dt/2)@(U+dt/2*k2)
        k4=-1j*H(s+dt)@(U+dt*k3)
        U=U+dt/6*(k1+2*k2+2*k3+k4)
    return U
def full_prop(tau,E1,t1c,tc,E0):
    # double braid = steps 1,2,3, then 1,2,3
    U1=step_prop(tau,E1,t1c,tc,E0,0)
    U2=step_prop(tau,E1,t1c,tc,E0,1)
    U3=step_prop(tau,E1,t1c,tc,E0,2)
    U=np.eye(4,dtype=complex)
    for Uu in (U1,U2,U3,U1,U2,U3): U=Uu@U
    return U
def fidelity(tau,E1,t1c,tc,E0):
    U=full_prop(tau,E1,t1c,tc,E0)
    # code qubit = qubit1 (left factor). n1: i g1 g2 = -1 when occupied (see diag). Use report formula:
    # fid = 1/2 * sum_{a,b} |<1,a|U|0,b>|^2   where code label = qubit1 (rows/cols index = 2*n1+nd? need care)
    # Our row order: index = i1*2 + i2 with i1=qubit1(n1), i2=qubit2(nd). 
    # "code=0" means qubit1 in state |0> -> i1=0 ; code=1 -> i1=1.
    fid=0.0
    for b in range(2):
        col0 = 0*2+b   # |code=0, nd=b>
        for a in range(2):
            row1 = 1*2+a   # <code=1, nd=a|
            fid += abs(U[row1,col0])**2
    return 0.5*fid
# ---- MZM anchor tests ----
tc=0.3; E0=0.3
for tau_plot in [3,10]:
    tau=tau_plot*100.0   # plot->physical: 1 plot = 100/meV
    print("double braid MZM tau_plot=%g fid=% .5f"%(tau_plot, fidelity(tau,0.0,0.0,tc,E0)))
    # single braid: steps 1,2,3 once
    U1=step_prop(tau,0,0,tc,E0,0);U2=step_prop(tau,0,0,tc,E0,1);U3=step_prop(tau,0,0,tc,E0,2)
    U=U3@U2@U1
    fid=0.0
    for b in range(2):
        for a in range(2): fid+=abs(U[1*2+a,0*2+b])**2
    print("single braid MZM tau_plot=%g fid=% .5f"%(tau_plot,0.5*fid))
print("--- extra anchors ---")
tc=0.3; E0=0.3
# no coupling at all: tc=E0=0, E1=t1=0 -> nothing braided -> fid 0
print("no-coupling fid=", fidelity(3*100,0,0,0,0))
# pure E1: tc=E0=t1=0, E1=0.01 -> no braiding, only phase? E1 couples g1g2 (within wire) - acts on code but not between code/ancilla... check
print("pure-E1 fid=", fidelity(3*100,0.01,0,0,0))
# t1 only (E1=0, tc=E0=0): t1 couples gamma_b - gamma1 (ancilla-code) but t3=0, tc=0 means t2,t3 never turn on -> check
print("pure-t1 fid=", fidelity(3*100,0,0.1,0,0))
