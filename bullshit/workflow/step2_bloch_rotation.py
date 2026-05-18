"""
Step 2: Bloch rotation on effective two-level subspace.

For each (u, delta), analyze the effective Hamiltonian H_eff = d_0 I + d·σ
and compute energy level splitting.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from config import RESULTS_STEP2, U_LIST_DEFAULT, DELTA_LIST_DEFAULT, VERBOSE, SAVE_FIGURES, FIG_DPI, FIG_FORMAT
from step1_eight_vertex import eight_vertex_H4, process_u_delta_point

def analyze_bloch_rotation(u, delta):
    """
    Analyze effective two-level subspace for given (u, delta).
    
    Returns effective energy splitting E_± = d_0 ± |d|
    """
    result = process_u_delta_point(u, delta)
    d = result['d']
    d0 = result['d0']
    d_mag = np.linalg.norm(d)
    
    E_plus = d0 + d_mag
    E_minus = d0 - d_mag
    
    return {
        'u': u,
        'delta': delta,
        'd': d,
        'd0': d0,
        '|d|': d_mag,
        'E_+': E_plus,
        'E_-': E_minus,
        'splitting': E_plus - E_minus,
    }

def run(u_list=None, delta_list=None):
    """
    Run Step 2 analysis over all (u, delta).
    """
    if u_list is None:
        u_list = U_LIST_DEFAULT
    if delta_list is None:
        delta_list = DELTA_LIST_DEFAULT
    
    if VERBOSE:
        print(f"[Step 2] Analyzing Bloch rotation for {len(u_list)} x {len(delta_list)} points...")
    
    results = {}
    for u in u_list:
        for delta in delta_list:
            key = (u, delta)
            results[key] = analyze_bloch_rotation(u, delta)
    
    # Save
    output_file = RESULTS_STEP2 / "bloch_rotation_data.npz"
    np.savez(output_file, **{
        'u_list': u_list,
        'delta_list': delta_list,
        'results': results,
    })
    
    if VERBOSE:
        print(f"[Step 2] Results saved to {output_file}")
    
    # Plot: |d| as function of delta for selected u values
    if SAVE_FIGURES and len(u_list) > 0 and len(delta_list) > 1:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        
        # |d| vs delta for fixed u values
        for u in u_list[::max(1, len(u_list)//3)]:  # Sample 3 u values
            d_mags = [results[(u, d)]['|d|'] for d in delta_list]
            axes[0].plot(delta_list, d_mags, 'o-', label=f'u={u:.2f}')
        axes[0].set_xlabel('δ')
        axes[0].set_ylabel('|d|')
        axes[0].set_title('Bloch vector magnitude')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Energy splitting vs delta for fixed u
        for u in u_list[::max(1, len(u_list)//3)]:
            splittings = [results[(u, d)]['splitting'] for d in delta_list]
            axes[1].plot(delta_list, splittings, 's-', label=f'u={u:.2f}')
        axes[1].set_xlabel('δ')
        axes[1].set_ylabel('E₊ - E₋')
        axes[1].set_title('Energy splitting')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        
        fig.tight_layout()
        fig_file = RESULTS_STEP2 / f"bloch_rotation.{FIG_FORMAT}"
        plt.savefig(fig_file, dpi=FIG_DPI)
        if VERBOSE:
            print(f"[Step 2] Figure saved to {fig_file}")
        plt.close()
    
    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Step 2: Bloch rotation analysis')
    parser.add_argument('--u-list', type=str, default=None, help='Comma-separated u values')
    parser.add_argument('--delta-list', type=str, default=None, help='Comma-separated delta values')
    
    args = parser.parse_args()
    
    u_list = U_LIST_DEFAULT if args.u_list is None else [float(x) for x in args.u_list.split(',')]
    delta_list = DELTA_LIST_DEFAULT if args.delta_list is None else [float(x) for x in args.delta_list.split(',')]
    
    results = run(u_list, delta_list)
