#!/usr/bin/env python3
"""Plot Majorana zero-mode counts for grid solutions.

Reads `scripts/ybe_grid_results_coarse.csv`, selects rows with small
`residual_norm`, computes zero-mode counts using
`scripts/majorana_jw_mapping.py`'s `analyze_point`, and saves per-`Jz`
scatter plots and a 3D scatter of (Jx,Jy,Jz) colored by zero-mode count.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import majorana_jw_mapping as mjm
except Exception as e:
    print('Failed to import majorana_jw_mapping:', e)
    raise


CSV = Path('scripts/ybe_grid_results_coarse.csv')
OUT_DIR = Path('scripts')
OUT_DIR.mkdir(exist_ok=True)

tol = 1e-8
points = []
with CSV.open() as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            Jx = float(row['Jx'])
            Jy = float(row['Jy'])
            Jz = float(row['Jz'])
            res = float(row.get('residual_norm', row.get('residual', 1e9)))
        except Exception:
            continue
        if abs(res) <= tol:
            points.append((Jx, Jy, Jz))

if not points:
    print('No points under tolerance found in', CSV)
    raise SystemExit(1)

data = []
for Jx, Jy, Jz in points:
    info = mjm.analyze_point(Jx, Jy, Jz)
    data.append((Jx, Jy, Jz, info['zero_modes']))

data = np.array(data)
Jx = data[:, 0].astype(float)
Jy = data[:, 1].astype(float)
Jz = data[:, 2].astype(float)
Zms = data[:, 3].astype(int)

# 3D scatter
fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(Jx, Jy, Jz, c=Zms, cmap='viridis', s=30)
ax.set_xlabel('Jx')
ax.set_ylabel('Jy')
ax.set_zlabel('Jz')
plt.colorbar(sc, label='zero_modes')
plt.tight_layout()
p3 = OUT_DIR / 'majorana_zero_modes_3d.png'
fig.savefig(p3, dpi=200)
plt.close(fig)

# Per-Jz 2D plots
unique_Jz = np.unique(Jz)
for idx, vz in enumerate(unique_Jz):
    mask = np.isclose(Jz, vz)
    if np.sum(mask) == 0:
        continue
    fig, ax = plt.subplots(figsize=(5, 4))
    sc = ax.scatter(Jx[mask], Jy[mask], c=Zms[mask], cmap='plasma', s=40)
    ax.set_title(f'Jz = {vz:.6f}')
    ax.set_xlabel('Jx')
    ax.set_ylabel('Jy')
    plt.colorbar(sc, label='zero_modes')
    plt.tight_layout()
    outp = OUT_DIR / f'majorana_zero_modes_Jz_{idx}.png'
    fig.savefig(outp, dpi=200)
    plt.close(fig)

print('Saved plots to', OUT_DIR)
