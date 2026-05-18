"""
Step 1: Eight-vertex model Pauli expansion and Kitaev mapping.

For each (u, delta), compute:
- Local Hamiltonian H_4(u, delta)
- Pauli coefficients h_{αβ}
- Effective d vector
- Mapped Kitaev parameters (t, Δ, μ)
"""

import numpy as np
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from config import RESULTS_STEP1, U_LIST_DEFAULT, DELTA_LIST_DEFAULT, VERBOSE
from utils import I2, X, Y, Z, kron, pauli_expand, map_to_kitaev, extract_d_vector

def eight_vertex_H4(u, delta):
    """
    Construct the 4x4 eight-vertex Hamiltonian.
    
    H_4(u, delta) = cos(u) * X⊗X 
                  + (sin u / 2) * (Y⊗X - X⊗Y)
                  + (delta / 2) * (Z⊗I - I⊗Z)
    """
    XX = kron(X, X)
    YX = kron(Y, X)
    XY = kron(X, Y)
    ZI = kron(Z, I2)
    IZ = kron(I2, Z)
    
    H = np.cos(u) * XX + (np.sin(u) / 2.0) * (YX - XY) + (delta / 2.0) * (ZI - IZ)
    return H

def process_u_delta_point(u, delta):
    """
    Process a single (u, delta) point.
    
    Returns:
        dict with 'u', 'delta', 'h', 'd', 'd0', 't', 'Delta', 'mu'
    """
    H4 = eight_vertex_H4(u, delta)
    h = pauli_expand(H4)
    d, d0 = extract_d_vector(h)
    t, Delta, mu = map_to_kitaev(h)
    
    return {
        'u': u,
        'delta': delta,
        'h': h,
        'd': d,
        'd0': d0,
        't': t,
        'Delta': Delta,
        'mu': mu,
        '|d|': np.linalg.norm(d),
    }

def run(u_list=None, delta_list=None):
    """
    Run Step 1 for all (u, delta) pairs.
    """
    if u_list is None:
        u_list = U_LIST_DEFAULT
    if delta_list is None:
        delta_list = DELTA_LIST_DEFAULT
    
    if VERBOSE:
        print(f"[Step 1] Processing {len(u_list)} x {len(delta_list)} = {len(u_list)*len(delta_list)} points...")
    
    # Storage
    results = {}
    
    for u in u_list:
        for delta in delta_list:
            key = (u, delta)
            results[key] = process_u_delta_point(u, delta)
            
            if VERBOSE:
                d = results[key]['d']
                d0 = results[key]['d0']
                t = results[key]['t']
                Delta = results[key]['Delta']
                mu = results[key]['mu']
                print(f"  u={u:.3f}, δ={delta:.3f}: |d|={np.linalg.norm(d):.4f}, d_z={d[2]:.4f}, "
                      f"t={t:.4f}, Δ={Delta:.4f}, μ={mu:.4f}")
    
    # Save
    output_file = RESULTS_STEP1 / "eight_vertex_data.npz"
    np.savez(output_file, **{
        'u_list': u_list,
        'delta_list': delta_list,
        'results': results,
    })
    
    if VERBOSE:
        print(f"[Step 1] Results saved to {output_file}")
    
    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Step 1: Eight-vertex model')
    parser.add_argument('--u-list', type=str, default=None, help='Comma-separated u values')
    parser.add_argument('--delta-list', type=str, default=None, help='Comma-separated delta values')
    
    args = parser.parse_args()
    
    u_list = U_LIST_DEFAULT if args.u_list is None else [float(x) for x in args.u_list.split(',')]
    delta_list = DELTA_LIST_DEFAULT if args.delta_list is None else [float(x) for x in args.delta_list.split(',')]
    
    results = run(u_list, delta_list)
