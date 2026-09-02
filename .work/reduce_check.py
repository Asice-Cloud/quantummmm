import numpy as np, sys
from scipy.linalg import expm as _sexpm
sys.path.insert(0,'.work')
from dyn import step_prop, full_prop, fidelity
tc=0.3; E0=0.3

# Reduced 2x2 code-space model per paper: during window where both E1 and t1 active,
# H_reduced ~ a(t) sz + b(t) sy  (E1->sigma_z longitudinal, t1->sigma_y transverse)
# Compare the three-regime shape from full JW engine vs this reduced model.
def reduced_prop(tau,E1,t1c,n=2000):
    sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
    def fm(x): return 0.5*(1-np.cos(np.pi*x/tau))
    def fp(x): return 0.5*(1+np.cos(np.pi*x/tau))
    # step phase 0: E1 on full, t1 ramps fm; step phase 1: t1 on fp... approximate reduced H
    # Use effective E1_eff,t1_eff as piecewise constants acting through the double braid
    U=np.eye(2,dtype=complex)
    # phases: E1 amplitude / t1 amplitude per step (from dyn step table)
    phases=[(1.0,0.5),(1.0,1.0),(0.0,0.0)]*2
    dt=tau/n
    for (ampE,amt1) in phases:
        # full double braid windows only accumulate in windows with E1 (steps 1&2)
        if ampE==0: continue
        for k in range(n):
            Hr = (ampE*E1)*sz + (amt1*t1c)*sy   # noncommuting sz,sy
            Uexp = np.eye(2,dtype=complex)
            # use small-step exp
            Uexp = _sexpm(-1j*Hr*dt)
            U = Uexp@U
    return U

# fidelity of reduced: fidelity of NOT-gate(geometric, geometric part ignored) - compare oscillation freq
E1=0.01
for r in [0.1,1.0,10.0]:
    t1c=r*E1
    print("ratio",r)
    for tau_plot in [1.0,2.0,4.0,6.0,8.0,10.0]:
        Wfull=fidelity(tau_plot*100,E1,t1c,tc,E0)
        Ured=reduced_prop(tau_plot*100,E1,t1c)
        # reduced model (no geometric NOT) -- weight = chance of staying? use sin^2 measure
        # We just report |Ured[0,0]|^2 oscillation to compare with full W envelope
        print("  tau=%.0f Wfull=%.4f  |Ured_00|^2=%.4f"%(tau_plot,Wfull,abs(Ured[0,0])**2))
