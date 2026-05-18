"""
Step 3: Full-chain BdG Hamiltonian and edge localization.

For each (u, delta), map to Kitaev parameters and compute:
- Bulk gap
- Topological criterion |μ| < 2|t|
- Open-chain spectrum E_0(L) for various chain lengths L
- LDOS (local density of states)
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from config import (RESULTS_STEP3, U_LIST_DEFAULT, DELTA_LIST_DEFAULT, 
                    L_LIST_DEFAULT, E_MIN_LDOS, E_MAX_LDOS, NE_LDOS, 
                    NK_BULK_DEFAULT, ETA_LDOS_DEFAULT, VERBOSE, SAVE_FIGURES, FIG_DPI, FIG_FORMAT)
from utils import (build_bdg_chain, diagonalize_bdg, compute_bulk_gap, 
                   topo_criterion, compute_edge_weight, compute_ldos)
from step1_eight_vertex import process_u_delta_point

def analyze_full_chain(u, delta, L_list=None):
    """
    Full-chain BdG analysis for given (u, delta).
    
    Returns:
        dict with bulk_gap, topo_info, E_0(L), edge_weights, etc.
    """
    if L_list is None:
        L_list = L_LIST_DEFAULT
    
    # Get Kitaev parameters from eight-vertex
    pv = process_u_delta_point(u, delta)
    t = pv['t']
    Delta = pv['Delta']
    mu = pv['mu']
    
    # Bulk gap (k-space)
    gap = compute_bulk_gap(t, Delta, mu, nk=NK_BULK_DEFAULT)
    is_topo = topo_criterion(t, Delta, mu, gap_threshold=1e-12)
    
    # Open-chain E_0 and edge weights for various L
    E_0_data = {}
    for L in L_list:
        H_bdg = build_bdg_chain(L, t, Delta, mu, open_boundary=True)
        E, V = diagonalize_bdg(H_bdg)
        E_0 = float(np.min(np.abs(E)))
        
        # Compute edge weight for zero-mode-like state
        if len(V) > 0:
            v = V[:, 0]  # Ground or lowest-energy state
            left, right = compute_edge_weight(v, L, edge_sites=max(1, L//10))
            edge_wt = left + right
        else:
            edge_wt = 0.0
        
        E_0_data[L] = {
            'E_0': E_0,
            'edge_weight': edge_wt,
            'spectrum': E,
        }
    
    # LDOS for longest chain (for visualization)
    L_max = max(L_list)
    H_bdg = build_bdg_chain(L_max, t, Delta, mu, open_boundary=True)
    energies = np.linspace(E_MIN_LDOS, E_MAX_LDOS, NE_LDOS)
    ldos_es, ldos, E_spec = compute_ldos(H_bdg, energies, eta=ETA_LDOS_DEFAULT)
    
    return {
        'u': u,
        'delta': delta,
        't': t,
        'Delta': Delta,
        'mu': mu,
        'bulk_gap': gap,
        'is_topological': is_topo,
        'E_0_data': E_0_data,
        'ldos': ldos,
        'ldos_energies': ldos_es,
        'L_max': L_max,
    }

def run(u_list=None, delta_list=None, L_list=None):
    """
    Run Step 3 for all (u, delta) pairs.
    """
    if u_list is None:
        u_list = U_LIST_DEFAULT
    if delta_list is None:
        delta_list = DELTA_LIST_DEFAULT
    if L_list is None:
        L_list = L_LIST_DEFAULT
    
    if VERBOSE:
        print(f"[Step 3] Full-chain BdG analysis for {len(u_list)} x {len(delta_list)} points, L={L_list}...")
    
    results = {}
    for u in u_list:
        for delta in delta_list:
            key = (u, delta)
            results[key] = analyze_full_chain(u, delta, L_list)
            
            if VERBOSE:
                res = results[key]
                print(f"  u={u:.3f}, δ={delta:.3f}: gap={res['bulk_gap']:.4f}, "
                      f"topo={res['is_topological']}, E_0(L={L_list[-1]})={res['E_0_data'][L_list[-1]]['E_0']:.4e}")
    
    # Save
    output_file = RESULTS_STEP3 / "full_chain_data.npz"
    np.savez(output_file, **{
        'u_list': u_list,
        'delta_list': delta_list,
        'L_list': L_list,
        'results': results,
    })
    
    if VERBOSE:
        print(f"[Step 3] Results saved to {output_file}")
    
    # Plot: E_0(L) for selected (u, delta)
    if SAVE_FIGURES and len(u_list) > 0 and len(delta_list) > 0:
        fig, axes = plt.subplots(len(delta_list), len(u_list), figsize=(12, 8))
        if len(u_list) == 1:
            axes = axes.reshape(-1, 1)
        if len(delta_list) == 1:
            axes = axes.reshape(1, -1)
        
        for i, delta in enumerate(delta_list):
            for j, u in enumerate(u_list):
                ax = axes[i, j]
                res = results[(u, delta)]
                L_vals = sorted(res['E_0_data'].keys())
                E_0_vals = [res['E_0_data'][L]['E_0'] for L in L_vals]
                ax.semilogy(L_vals, E_0_vals, 'o-', markersize=6)
                ax.set_xlabel('L')
                ax.set_ylabel('E₀(L)')
                title = f'u={u:.2f}, δ={delta:.3f}'
                if res['is_topological']:
                    title += ' [Topo]'
                ax.set_title(title)
                ax.grid(True, alpha=0.3, which='both')
        
        fig.tight_layout()
        fig_file = RESULTS_STEP3 / f"E0_vs_L.{FIG_FORMAT}"
        plt.savefig(fig_file, dpi=FIG_DPI)
        if VERBOSE:
            print(f"[Step 3] E₀(L) plot saved to {fig_file}")
        plt.close()
    
    # Plot: LDOS heatmap for selected point
    if SAVE_FIGURES and len(u_list) > 0 and len(delta_list) > 0:
        u_idx, delta_idx = 0, 0  # First point
        u, delta = u_list[u_idx], delta_list[delta_idx]
        res = results[(u, delta)]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        im = ax.imshow(res['ldos'], aspect='auto', origin='lower', 
                       extent=[res['ldos_energies'][0], res['ldos_energies'][-1], 0, res['L_max']],
                       cmap='hot', interpolation='nearest')
        ax.set_xlabel('E')
        ax.set_ylabel('Site index')
        ax.set_title(f'LDOS: u={u:.2f}, δ={delta:.3f}')
        plt.colorbar(im, ax=ax, label='LDOS')
        fig.tight_layout()
        
        fig_file = RESULTS_STEP3 / f"ldos_u{u:.2f}_d{delta:.3f}.{FIG_FORMAT}"
        plt.savefig(fig_file, dpi=FIG_DPI)
        if VERBOSE:
            print(f"[Step 3] LDOS plot saved to {fig_file}")
        plt.close()
    
    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Step 3: Full-chain BdG')
    parser.add_argument('--u-list', type=str, default=None, help='Comma-separated u values')
    parser.add_argument('--delta-list', type=str, default=None, help='Comma-separated delta values')
    parser.add_argument('--L-list', type=str, default=None, help='Comma-separated chain lengths')
    
    args = parser.parse_args()
    
    u_list = U_LIST_DEFAULT if args.u_list is None else [float(x) for x in args.u_list.split(',')]
    delta_list = DELTA_LIST_DEFAULT if args.delta_list is None else [float(x) for x in args.delta_list.split(',')]
    L_list = L_LIST_DEFAULT if args.L_list is None else [int(x) for x in args.L_list.split(',')]
    
    results = run(u_list, delta_list, L_list)
