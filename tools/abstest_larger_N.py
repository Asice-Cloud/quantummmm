python3 - <<'PY'
import os,sys
ROOT = os.path.abspath('.')
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import numpy as np
from numpy.linalg import eigh
# reuse build_midpair_chain from tools.run_midpair_demo
from tools.run_midpair_demo import build_midpair_chain

N = 400
Delta_center = 0.5
center_width = 80
H, t_hop, Delta_bond, mu_site = build_midpair_chain(N, t_val=1.0, Delta_center=Delta_center, center_width=center_width, mu=0.0)
vals, vecs = eigh(H)
# save spectrum plot
os.makedirs('results', exist_ok=True)
try:
    import matplotlib.pyplot as plt
    plt.figure(figsize=(6,4))
    plt.plot(np.sort(vals), '.', markersize=3)
    plt.title(f'Mid-pair spectrum N={N}')
    plt.tight_layout()
    plt.savefig(f'results/xyz_midpair_N{N}_spectrum.png')
    plt.close()
except Exception:
    pass

# lowest modes center/edge fractions
idx = np.argsort(np.abs(vals))
center_idxs = [i for i,v in enumerate(Delta_bond) if abs(v)>1e-12]
if center_idxs:
    start, end = center_idxs[0], center_idxs[-1]
else:
    start, end = N//2 - center_width//2, N//2 + center_width//2
center = list(range(start, end+1))
left_edge = list(range(0, start - 20 if start - 20>0 else 0))
right_edge = list(range(end+1 + 20 if end+1+20 < N else end+1, N))

out_lines = []
for k in range(6):
    i = idx[k]
    vec = vecs[:,i]
    Nloc = vec.shape[0]//2
    u = vec[:Nloc]
    v = vec[Nloc:]
    dens = np.abs(u)**2 + np.abs(v)**2
    total = dens.sum()
    c = dens[center].sum()/total
    le = dens[left_edge].sum()/total if left_edge else 0.0
    re = dens[right_edge].sum()/total if right_edge else 0.0
    out_lines.append((k, float(vals[i]), c, le, re))

# save density and report
vec0 = vecs[:, idx[0]]
Nloc = vec0.shape[0]//2
u = vec0[:Nloc]
v = vec0[Nloc:]
dens0 = np.abs(u)**2 + np.abs(v)**2
np.savetxt(f'results/xyz_midpair_N{N}_density0.txt', dens0)

with open(f'results/xyz_midpair_N{N}_report.txt','w') as f:
    f.write(f'N={N}, center_width={center_width}, Delta_center={Delta_center}\n')
    f.write(f'center indices start={start} end={end}\n')
    f.write('\nLowest energies (abs-sorted first 6):\n')
    for row in out_lines:
        f.write(f'mode{row[0]} energy={row[1]:.6e} center_frac={row[2]:.6e} left_frac={row[3]:.6e} right_frac={row[4]:.6e}\n')

print(f'Saved results/xyz_midpair_N{N}_spectrum.png, results/xyz_midpair_N{N}_density0.txt, results/xyz_midpair_N{N}_report.txt')
PY