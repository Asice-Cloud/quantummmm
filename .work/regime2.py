import numpy as np, sys, json, time
sys.path.insert(0,'.work')
from engine2 import full_prop, fid
tc=0.3; E0=0.3
def curve(E1,t1c,taus):
    out=[]
    for t in taus:
        out.append(fid(full_prop(t,E1,t1c,tc,E0)))
    return np.array(out)
taus=np.linspace(200,1400,61)
# quick timing on one curve
t0=time.time(); W=curve(0.01,0.01,taus); dt=time.time()-t0
print("timing one curve(61 pts): %.1f s, W[0]=%.4f W[-1]=%.4f"%(dt,W[0],W[-1]))
def stats(W,E1,t1c,tag):
    zc=sum(1 for i in range(1,len(W)) if (W[i]-0.5)*(W[i-1]-0.5)<0)
    # first maximum location in tau (coarse)
    iarg=int(np.argmax(W))
    return dict(tag=tag,E1=E1,t1=t1c,tau_first_max=float(taus[iarg]),
                Wmax=float(W.max()),Wmin=float(W.min()),
                Wend=float(W[-1]),halfcross=zc)
rows=[]
for E1name,E1 in [("0.01",0.01),("0.001",0.001)]:
    for tag,r in [("pure-E1",0.0),("r=0.1",0.1),("r=0.3",0.3),("r=1",1.0),("r=3",3.0),("r=10",10.0)]:
        t1c=r*E1
        W=curve(E1,t1c,taus)
        s=stats(W,E1,t1c,tag)
        rows.append(s)
        # print coarse W trace every ~10 pts
        tr=" ".join("%.2f"%v for v in W[::10])
        print("E1=%-5s %-8s t1=%.5f Wmax=%.4f Wmin=%.4f Wend=%.4f hc=%2d firstmax_tau=%.0f | trace(60pt stride): %s"%(
            E1name,tag,t1c,s['Wmax'],s['Wmin'],s['Wend'],s['halfcross'],s['tau_first_max'],tr))
json.dump(rows,open('.work/regime2.json','w'),indent=1)
np.save('.work/regime2_taus.npy',taus)
print("saved")
