#!/usr/bin/env python3
"""Plot LDOS and Majorana-similarity for validated top candidates.

Reads `results/abs_validate_top.json` produced by `tools/abs_validate_top.py`.
Generates per-candidate LDOS plots for each N and a summary Maj-sim plot.
"""
import json
import os
import numpy as np
import matplotlib.pyplot as plt

RES = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results')
IN_PATH = os.path.join(RES, 'abs_validate_top.json')

with open(IN_PATH, 'r') as f:
    data = json.load(f)

os.makedirs(RES, exist_ok=True)

for i, cand in enumerate(data, start=1):
    per_N = cand.get('per_N', {})
    Ns = sorted(int(k) for k in per_N.keys())

    # Plot LDOS for each N in a single figure (stacked subplots)
    fig, axes = plt.subplots(len(Ns), 1, figsize=(8, 2.5*len(Ns)), sharex=True)
    if len(Ns) == 1:
        axes = [axes]

    for ax, N in zip(axes, Ns):
        entry = per_N[str(N)]
        site_weights = np.array(entry.get('site_weights', []))
        ax.plot(site_weights, '-', lw=1)
        ax.axvline(entry.get('site', 0), color='C1', ls='--', label=f"max site {entry.get('site')}")
        ax.set_ylabel(f'N={N}')
        ax.legend()

    axes[-1].set_xlabel('Site index')
    fig.suptitle(f"Top-{i} LDOS (cxy={cand.get('cxy')}, cyx={cand.get('cyx')}, cxz={cand.get('cxz')})")
    out = os.path.join(RES, f'abs_top{i}_ldos.png')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out)
    plt.close(fig)

    # Plot Maj-sim (1 - maj_norm) vs N
    maj_sims = [per_N[str(N)]['maj_sim'] for N in Ns]
    fig = plt.figure(figsize=(6,3))
    plt.plot(Ns, maj_sims, 'o-')
    plt.xlabel('N')
    plt.ylabel('Majorana similarity (1 - norm)')
    plt.title(f'Top-{i} Majorana similarity')
    out2 = os.path.join(RES, f'abs_top{i}_maj_sim.png')
    fig.tight_layout()
    fig.savefig(out2)
    plt.close(fig)

    print(f'Wrote: {out}, {out2}')

print('All plots written to', RES)
