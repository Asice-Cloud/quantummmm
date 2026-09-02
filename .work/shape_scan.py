import numpy as np, sys
sys.path.insert(0,'.work')
from dyn import step_prop
tc=0.3; E0=0.3
def full_prop(tau,E1,t1c,n=1200):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2):
        U=step_prop(tau,E1,t1c,tc,E0,ph,n=n)@U
    return U
def fid(U):
    return 0.5*sum(abs(U[1*2+a,0*2+b])**2 for a in range(2) for b in range(2))
taus=np.logspace(-0.3,1.5,40)   # tau_phys 0.5 .. 31
print("tau_phys | MZM | E1=.01 t1=.001 | E1=.01 t1=.01 | E1=.01 t1=.1")
for t in taus:
    row=[fid(full_prop(t,0,0)),fid(full_prop(t,.01,.001)),fid(full_prop(t,.01,.01)),fid(full_prop(t,.01,.1))]
    print("%7.2f %10.4f %16.4f %16.4f %14.4f"%tuple([t]+row))
