"""
Step 4: Comprehensive MZM/ABS diagnosis and summary.

Compare subspace properties (d vector, Bloch rotation) with full-chain properties
(bulk gap, edge localization, LDOS) to provide unified MZM/ABS classification.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from config import (RESULTS_STEP4, U_LIST_DEFAULT, DELTA_LIST_DEFAULT, 
                    L_LIST_DEFAULT, VERBOSE, SAVE_FIGURES, FIG_DPI, FIG_FORMAT)
from step1_eight_vertex import process_u_delta_point
from step2_bloch_rotation import analyze_bloch_rotation
from step3_full_chain_bdg import analyze_full_chain

def diagnose_mzm_abs(u, delta, L_list=None):
    """
    Unified diagnosis: combine subspace and full-chain evidence.
    
    Returns:
        dict with subspace_type, global_type, evidence
    """
    if L_list is None:
        L_list = L_LIST_DEFAULT
    
    # Layer 1: Eight-vertex + Pauli
    pv = process_u_delta_point(u, delta)
    
    # Layer 2: Subspace Bloch rotation
    bloch = analyze_bloch_rotation(u, delta)
    d = bloch['d']
    d_z = d[2]
    
    # Layer 3: Full-chain BdG
    full = analyze_full_chain(u, delta, L_list)
    gap = full['bulk_gap']
    is_topo = full['is_topological']
    
    # Subspace classification
    if np.abs(d_z) < 1e-6:
        subspace_type = "MZM-like (δ≈0, planar trajectory)"
    else:
        subspace_type = f"ABS-like (δ≠0, d_z={d_z:.4f})"
    
    # Full-chain classification
    L_max = L_list[-1]
    E_0 = full['E_0_data'][L_max]['E_0']
    
    # MZM signature: exponential decay of E_0(L)
    if len(L_list) > 1:
        E_vals = [full['E_0_data'][L]['E_0'] for L in sorted(L_list)]
        L_vals = np.array(sorted(L_list), dtype=float)
        # Fit log(E_0) ~ -α*L to detect exponential decay
        log_E = np.log(np.array(E_vals) + 1e-10)
        if len(L_vals) >= 2:
            alpha = -np.polyfit(L_vals, log_E, 1)[0]
            has_exponential_decay = alpha > 0.01
        else:
            has_exponential_decay = False
    else:
        alpha = 0.0
        has_exponential_decay = False
    
    if is_topo and has_exponential_decay:
        global_type = "Majorana Zero Mode (MZM)"
        evidence = f"Topological gap={gap:.4f}, exponential decay α={alpha:.4f}, E_0(L={L_max})={E_0:.4e}"
    elif is_topo and E_0 < 0.01:
        global_type = "Topological edge state (likely MZM-like)"
        evidence = f"Topological gap={gap:.4f}, small E_0={E_0:.4e}"
    elif E_0 > 0.01:
        global_type = "Andreev Bound State (ABS)"
        evidence = f"Significant finite E_0(L={L_max})={E_0:.4e}"
    else:
        global_type = "Unclear / gapped"
        evidence = f"Gap={gap:.4f}, E_0={E_0:.4e}"
    
    return {
        'u': u,
        'delta': delta,
        'subspace_type': subspace_type,
        'global_type': global_type,
        'evidence': evidence,
        'd_z': d_z,
        'bulk_gap': gap,
        'E_0': E_0,
        'exponential_decay': has_exponential_decay,
        'alpha': alpha,
    }

def run(u_list=None, delta_list=None, L_list=None):
    """
    Run comprehensive MZM/ABS diagnosis.
    """
    if u_list is None:
        u_list = U_LIST_DEFAULT
    if delta_list is None:
        delta_list = DELTA_LIST_DEFAULT
    if L_list is None:
        L_list = L_LIST_DEFAULT
    
    if VERBOSE:
        print(f"[Step 4] Comprehensive MZM/ABS diagnosis for {len(u_list)} x {len(delta_list)} points...")
    
    results = {}
    summary_lines = [
        "=" * 100,
        "COMPREHENSIVE MZM/ABS DIAGNOSIS",
        "=" * 100,
        "",
    ]
    
    for u in u_list:
        for delta in delta_list:
            key = (u, delta)
            results[key] = diagnose_mzm_abs(u, delta, L_list)
            res = results[key]
            
            summary_lines.append(f"u={u:.3f}, δ={delta:.3f}")
            summary_lines.append(f"  Subspace: {res['subspace_type']}")
            summary_lines.append(f"  Global:   {res['global_type']}")
            summary_lines.append(f"  Evidence: {res['evidence']}")
            summary_lines.append("")
            
            if VERBOSE:
                print(f"  u={u:.3f}, δ={delta:.3f}: {res['global_type']}")
    
    # Save summary
    summary_file = RESULTS_STEP4 / "mzm_abs_summary.txt"
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    if VERBOSE:
        print(f"[Step 4] Summary saved to {summary_file}")
    
    # Plot: Summary table visualization
    if SAVE_FIGURES and len(u_list) > 0 and len(delta_list) > 0:
        fig, ax = plt.subplots(figsize=(12, max(4, 0.6*len(delta_list)*len(u_list))))
        ax.axis('tight')
        ax.axis('off')
        
        # Build table data
        table_data = [['u', 'δ', 'Subspace Type', 'Global Type', 'E₀', 'Gap']]
        colors = []
        header_colors = [['lightgray']*6]
        
        for u in u_list:
            row_colors = []
            for delta in delta_list:
                res = results[(u, delta)]
                subspace_short = "MZM-like" if "MZM-like" in res['subspace_type'] else "ABS-like"
                global_short = res['global_type'].split('(')[0].strip()
                
                table_data.append([
                    f"{u:.2f}",
                    f"{delta:.3f}",
                    subspace_short,
                    global_short,
                    f"{res['E_0']:.2e}",
                    f"{res['bulk_gap']:.4f}",
                ])
                
                # Color code
                if "MZM" in res['global_type']:
                    row_colors.append('lightgreen')
                elif "ABS" in res['global_type']:
                    row_colors.append('lightyellow')
                else:
                    row_colors.append('lightcoral')
            
            if row_colors:
                colors.append(row_colors)
        
        table = ax.table(cellText=table_data[1:], colLabels=table_data[0],
                        cellLoc='center', loc='center',
                        colWidths=[0.08, 0.08, 0.22, 0.25, 0.15, 0.12])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)
        
        # Color cells
        for i in range(len(table_data[0])):
            table[(0, i)].set_facecolor('lightgray')
        for i, row_colors in enumerate(colors):
            for j, color in enumerate(row_colors):
                table[(i+1, j)].set_facecolor(color)
        
        fig.tight_layout()
        fig_file = RESULTS_STEP4 / f"mzm_abs_summary.{FIG_FORMAT}"
        plt.savefig(fig_file, dpi=FIG_DPI, bbox_inches='tight')
        if VERBOSE:
            print(f"[Step 4] Summary plot saved to {fig_file}")
        plt.close()
    
    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Step 4: MZM/ABS comprehensive diagnosis')
    parser.add_argument('--u-list', type=str, default=None, help='Comma-separated u values')
    parser.add_argument('--delta-list', type=str, default=None, help='Comma-separated delta values')
    parser.add_argument('--L-list', type=str, default=None, help='Comma-separated chain lengths')
    
    args = parser.parse_args()
    
    u_list = U_LIST_DEFAULT if args.u_list is None else [float(x) for x in args.u_list.split(',')]
    delta_list = DELTA_LIST_DEFAULT if args.delta_list is None else [float(x) for x in args.delta_list.split(',')]
    L_list = L_LIST_DEFAULT if args.L_list is None else [int(x) for x in args.L_list.split(',')]
    
    results = run(u_list, delta_list, L_list)
