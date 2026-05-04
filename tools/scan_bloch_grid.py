#!/usr/bin/env python3
"""Scan d0 x T for final Bloch rotation angle and save heatmap.

This script imports run_for_T from tools/bloch_rotation_sim and scans
over a grid of hybridization strengths d0 and step durations T_step.
It saves a NumPy file and a PNG heatmap to results/.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Ensure project root is on sys.path so sibling imports work when running
# this script directly (python tools/scan_bloch_grid.py).
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

try:
    from tools.bloch_rotation_sim import run_for_T
except Exception:
    # fallback: try importlib (should work once parent added to sys.path)
    import importlib
    brs = importlib.import_module('tools.bloch_rotation_sim')
    run_for_T = brs.run_for_T


def compute_angle(b_init, b_final):
    # compute angle between two Bloch vectors in degrees
    ni = np.linalg.norm(b_init)
    nf = np.linalg.norm(b_final)
    if ni == 0 or nf == 0:
        return 0.0
    c = np.dot(b_init, b_final) / (ni * nf)
    c = max(-1.0, min(1.0, c))
    angle = np.arccos(c)
    return np.degrees(angle)


def main():
    os.makedirs('results', exist_ok=True)
    # grid (coarse for speed)
    d0_list = np.linspace(0.0, 0.4, 41)
    T_list = [50, 100, 200, 400, 800]
    nT = len(T_list)
    nd = len(d0_list)
    angles = np.zeros((nd, nT))

    for i, d0 in enumerate(d0_list):
        print(f"Scanning d0={d0:.3f} ({i+1}/{nd})")
        for j, T in enumerate(T_list):
            # use moderate resolution for speed
            t, bloch_traj, psi0, psi_traj = run_for_T(T_step=T, n_per_step=300, alpha=0.8, A=1.0, d_profile=(lambda t, d0=d0: d0))
            b_init = bloch_traj[0]
            b_final = bloch_traj[-1]
            angles[i, j] = compute_angle(b_init, b_final)

    out = 'results/bloch_scan_d0T.npz'
    np.savez(out, d0_list=d0_list, T_list=np.array(T_list), angles=angles)
    print('Saved', out)

    # heatmap
    plt.figure(figsize=(6,4))
    # transpose so y axis is T
    extent = [d0_list[0], d0_list[-1], T_list[0], T_list[-1]]
    plt.imshow(angles.T, origin='lower', aspect='auto', extent=extent, cmap='viridis')
    plt.colorbar(label='final rotation angle (deg)')
    plt.xlabel('d0 (hybridization)')
    plt.ylabel('T_step')
    plt.title('Final Bloch rotation angle vs d0 and T')
    outpng = 'results/bloch_scan_heatmap.png'
    plt.tight_layout()
    plt.savefig(outpng)
    print('Saved', outpng)


if __name__ == '__main__':
    main()
