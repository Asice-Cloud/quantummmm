#!/usr/bin/env python3
"""Higher-resolution scan of d0 x T for final Bloch rotation angle.

Uses `tools/bloch_rotation_sim.run_for_T` with constant d_profile=d0.
Saves `results/bloch_scan_highres.npz` and `results/bloch_scan_highres.png`.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# ensure repo root importable
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools.bloch_rotation_sim import run_for_T


def compute_angle(b_init, b_final):
    ni = np.linalg.norm(b_init)
    nf = np.linalg.norm(b_final)
    if ni == 0 or nf == 0:
        return 0.0
    c = np.dot(b_init, b_final) / (ni * nf)
    c = max(-1.0, min(1.0, c))
    return np.degrees(np.arccos(c))


def main():
    os.makedirs('results', exist_ok=True)
    d0_list = np.linspace(0.0, 0.4, 41)
    T_list = np.linspace(50, 1000, 20)

    angles = np.zeros((len(d0_list), len(T_list)))
    for i, d0 in enumerate(d0_list):
        print(f"d0 {i+1}/{len(d0_list)} = {d0:.3f}")
        for j, T in enumerate(T_list):
            t, bloch_traj, psi0, psi_traj = run_for_T(T_step=float(T), n_per_step=300, alpha=0.8, A=1.0, d_profile=(lambda tt, d0=d0: d0))
            b_init = bloch_traj[0]
            b_final = bloch_traj[-1]
            angles[i, j] = compute_angle(b_init, b_final)

    out = 'results/bloch_scan_highres.npz'
    np.savez(out, d0_list=d0_list, T_list=T_list, angles=angles)
    print('Saved', out)

    plt.figure(figsize=(8,5))
    extent = [d0_list[0], d0_list[-1], T_list[0], T_list[-1]]
    plt.imshow(angles.T, origin='lower', aspect='auto', extent=extent, cmap='plasma')
    plt.colorbar(label='final rotation angle (deg)')
    plt.xlabel('d0')
    plt.ylabel('T')
    plt.title('High-res final Bloch rotation angle vs d0 and T')
    outpng = 'results/bloch_scan_highres.png'
    plt.tight_layout()
    plt.savefig(outpng)
    print('Saved', outpng)


if __name__ == '__main__':
    main()
