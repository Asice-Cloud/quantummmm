import numpy as np, sys
sys.path.insert(0,'.work')
from dyn import step_prop
tc=0.3; E0=0.3
def full_prop(tau,E1,t1c,n=1500):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2):
        U=step_prop(tau,E1,t1c,tc,E0,ph,n=n)@U
    return U
E1=0.01
# oscillation frequency of W vs tau: estimate as function of t1/E1
# W(tau) ~ oscillation; frequency in tau_plot units should scale ~ sqrt(E1^2+t1^2)*tau_phys_scale
taus=np.linspace(0.3,12,80)
for r in [0.1,0.3,1.0,3.0,10.0]:
    t1c=r*E1
    W=np.array([0.5*sum(abs(U[1*2+a,0*2+b])**2 for a in range(2) for b in range(2))
                for U in (full_prop(t*100,E1,t1c) for t in taus)])
    # count zero crossings of (W-mean) as oscillation rate proxy
    w=(W-W.mean())
    zc=np.sum((w[1:]*w[:-1])<0)
    print("r=%.1f t1=%.3f  Wmin=%.4f Wmax=%.4f  zc=%d over dtau range %.1f-%.1f"%(r,t1c,W.min(),W.max(),zc,taus[0],taus[-1]))
