import numpy as np, sys
sys.path.insert(0,'.work')
from xchk2 import integrate_quat, quat_to_complex
from dyn import step_prop
tc=0.3; E0=0.3; tau=6*100.0
def jw_total(tau_phys,E1,t1c):
    U1=step_prop(tau_phys,E1,t1c,tc,E0,0,n=2000)
    U2=step_prop(tau_phys,E1,t1c,tc,E0,1,n=2000)
    U3=step_prop(tau_phys,E1,t1c,tc,E0,2,n=2000)
    U=np.eye(4,dtype=complex)
    for Uu in (U1,U2,U3,U1,U2,U3): U=Uu@U
    return U
def quat_total(tau_phys,E1,t1c):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2):
        U=quat_to_complex(integrate_quat(tau_phys,E1,t1c,tc,E0,ph,n=1500))@U
    return U
for E1,t1c in [(0,0),(0.01,0.3),(0.01,0.001),(0.001,0.01)]:
    a=np.sort(np.linalg.eigvals(jw_total(tau,E1,t1c)))
    b=np.sort(np.linalg.eigvals(quat_total(tau,E1,t1c)))
    # match phases: compare |eig| and relative ang diff
    da=np.sort(np.abs(np.diff(np.angle(a))))
    print("E1=%g t1=%g : JW eigs=%s"% (E1,t1c,np.round(a,4)))
    print("              quat eigs=%s"% (np.round(b,4)))
    print("              |JW|-|quat| max:",np.max(np.abs(np.abs(a)-np.abs(b))))
