import numpy as np, sys
sys.path.insert(0,'.work')
from dyn import step_prop
tc=0.3; E0=0.3
def full_prop(tau,E1,t1c,n=2500):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2):
        U=step_prop(tau,E1,t1c,tc,E0,ph,n=n)@U
    return U
def fid(U):
    return 0.5*sum(abs(U[1*2+a,0*2+b])**2 for a in range(2) for b in range(2))
# adiabatic window: tau_phys from 150 to 1400, check MZM baseline ~1
taus=np.linspace(150,1400,64)
for E1name,E1 in [("E1=0.001meV",0.001),("E1=0.01meV",0.01)]:
    print("=====",E1name)
    for r in [0.1,1.0,10.0]:
        t1c=r*E1
        Ws=np.array([fid(full_prop(t,E1,t1c)) for t in taus])
        # oscillation frequency: FFT of detrended W over tau
        w=Ws-Ws.mean()
        spec=np.abs(np.fft.rfft(w*np.hanning(len(w))))
        freqs=np.fft.rfftfreq(len(w),d=(taus[1]-taus[0]))
        fpk=freqs[np.argmax(spec[1:])+1]
        print("  r=%-4.1f t1=%.4f  W[min]=%.4f W[max]=%.4f W[tau=150]=%.4f W[tau=1400]=%.4f  peak_freq(1/tau)=%.5f"%(r,t1c,Ws.min(),Ws.max(),Ws[0],Ws[-1],fpk))
