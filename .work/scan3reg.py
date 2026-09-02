import numpy as np, sys
sys.path.insert(0,'.work')
from dyn import step_prop
tc=0.3; E0=0.3
def full_prop(tau,E1,t1c,n=3000):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2):
        U=step_prop(tau,E1,t1c,tc,E0,ph,n=n)@U
    return U
def fid(U):
    f=0.0
    for b in range(2):
        for a in range(2): f+=abs(U[1*2+a,0*2+b])**2
    return 0.5*f
E1=0.01  # meV (caption value)
ratios=[0.1,1.0,10.0]
taus=np.linspace(0.2,12,60)
import json
res={}
for r in ratios:
    t1c=r*E1
    Ws=[fid(full_prop(t*100,E1,t1c)) for t in taus]
    res[r]=[round(float(w),6) for w in Ws]
    print("ratio t1/E1=%.1f :"%r, " W[tau=0.5]=%.4f  W[max]=%.4f@tau=%.1f  W[end]=%.4f"%(
        Ws[3], max(Ws), taus[np.argmax(Ws)], Ws[-1]))
# also tau->0 limit check (geometric dominance): very small tau
for r in ratios:
    t1c=r*E1
    W0=fid(full_prop(0.05*100,E1,t1c,n=500))
    print("ratio=%.1f tau_plot=0.05 -> W=%.4f"%(r,W0))
np.save('.work/scan3reg.npy',np.array([taus]+[res[k] for k in ratios]))
print("saved")
