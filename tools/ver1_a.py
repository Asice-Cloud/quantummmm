import sys
sys.path.insert(0,'tools')
import verify_from_R as vr
import numpy as np
from numpy.linalg import norm
import math

def analyze(delta=0.015, N=600):
    us, du, H4s = vr.compute_Hs_from_R(delta, N=N)
    H2_list, ds = vr.build_H2_list_from_H4s(H4s)
    dt = du
    rank = np.linalg.matrix_rank(ds)
    mags = np.linalg.norm(ds, axis=1)
    I0 = float(np.sum(mags) * dt)
    alpha_target = math.pi / I0 if I0>0 else float('nan')
    e_vals, e_vecs, overlaps = vr.eigenbasis_series(H2_list)
    Nn = len(H2_list)
    comm_norms = []
    for j in range(Nn):
        if 0 < j < Nn-1:
            dV = (e_vecs[j+1] - e_vecs[j-1])/(2*dt)
        elif j==0:
            dV = (e_vecs[1] - e_vecs[0])/dt
        else:
            dV = (e_vecs[-1] - e_vecs[-2])/dt
        Vj = e_vecs[j]
        Aj = 1j * (Vj.conj().T @ dV)
        Dj = np.diag(e_vals[j])
        cn = norm(Dj @ Aj - Aj @ Dj, ord='fro')
        an = norm(Aj, ord='fro')
        dn = norm(Dj, ord='fro')
        denom = dn*an if dn*an>0 else 1.0
        comm_norms.append(cn/denom)
    comm_norms = np.array(comm_norms)
    print('d-vector rank =', rank)
    print('I0 =', I0, 'alpha_target =', alpha_target)
    print('comm_norms: min, mean, max =', comm_norms.min(), comm_norms.mean(), comm_norms.max())
    return dict(rank=rank, I0=I0, alpha_target=alpha_target, comm_norms=comm_norms)

res = analyze(0.015, N=600)
print('\nSUMMARY: rank=', res['rank'], 'alpha_target=', res['alpha_target'])
