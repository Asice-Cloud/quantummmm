import numpy as np, sys, json
sys.path.insert(0,'.work')
from dyn import step_prop, fidelity
tc=0.3; E0=0.3
def full_prop(tau,E1,t1c,n):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2):
        U=step_prop(tau,E1,t1c,tc,E0,ph,n=n)@U
    return U
def fid(U):
    return 0.5*sum(abs(U[1*2+a,0*2+b])**2 for a in range(2) for b in range(2))
# convergence check: MZM double braid at tau=1000, n=600 vs n=2400
Ua=full_prop(1000.0,0,0,2400); Ub=full_prop(1000.0,0,0,600)
print("CONV MZM fid(n2400)=%.6f fid(n600)=%.6f"%(fid(Ua),fid(Ub)))
# adiabatic window
taus=np.linspace(200,1400,49)
N=1200
def curve(E1,t1c):
    return np.array([fid(full_prop(t,E1,t1c,N)) for t in taus])
def describe(W,E1,t1c):
    # zero crossings of (W-0.5) => count oscillations of fidelity through midlevel
    zc=0
    for i in range(1,len(W)):
        if (W[i]-0.5)*(W[i-1]-0.5)<0: zc+=1
    # damping: compare |W-max| after first revival vs max
    return dict(Wmin=float(W.min()),Wmax=float(W.max()),
                W0=float(W[0]),Wend=float(W[-1]),
                halfcross=zc,
                Wsample=[round(float(W[i]),4) for i in (0,8,16,24,32,40,48)])
rows=[]
for E1name,E1 in [("0.01",0.01),("0.001",0.001)]:
    for tag,r in [("pureE1 t1=0",0.0),("r=0.1",0.1),("r=1",1.0),("r=3",3.0),("r=10",10.0)]:
        t1c=r*E1
        W=curve(E1,t1c)
        d=describe(W,E1,t1c)
        d.update(E1=E1name,t1=E1*E1,tag=tag)
        rows.append(d)
        print("E1=%-5s %-10s t1=%.5f Wmin=%.4f Wmax=%.4f W0=%.4f Wend=%.4f halfcross=%2d Wsamp=%s"%(
            E1name,tag,t1c,d['Wmin'],d['Wmax'],d['W0'],d['Wend'],d['halfcross'],d['Wsample']))
json.dump(rows,open('.work/regime_scan.json','w'),indent=1)
print("saved")
