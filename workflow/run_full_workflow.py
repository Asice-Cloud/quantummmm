#!/usr/bin/env python3
"""
Master workflow script: Orchestrates Steps 1-4 with default parameters.

Run this to execute the complete workflow from eight-vertex model through
subspace Bloch rotation to full-chain BdG and unified MZM/ABS diagnosis.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from config import VERBOSE, OUTPUT_DIR
from step1_eight_vertex import run as run_step1
from step2_bloch_rotation import run as run_step2
from step3_full_chain_bdg import run as run_step3
from step4_mzm_abs_criteria import run as run_step4

def main():
    print("\n" + "="*80)
    print("COMPLETE WORKFLOW: Eight-Vertex → Bloch → BdG → MZM/ABS Diagnosis")
    print("="*80 + "\n")
    
    if VERBOSE:
        print(f"Output directory: {OUTPUT_DIR}\n")
    
    # Step 1: Eight-vertex model
    print("\n[=== STEP 1 ===] Eight-Vertex Model Pauli Expansion")
    print("-" * 80)
    results1 = run_step1()
    
    # Step 2: Bloch rotation
    print("\n[=== STEP 2 ===] Bloch Rotation on Effective Subspace")
    print("-" * 80)
    results2 = run_step2()
    
    # Step 3: Full-chain BdG
    print("\n[=== STEP 3 ===] Full-Chain BdG Hamiltonian and Edge Localization")
    print("-" * 80)
    results3 = run_step3()
    
    # Step 4: MZM/ABS diagnosis
    print("\n[=== STEP 4 ===] Comprehensive MZM/ABS Diagnosis")
    print("-" * 80)
    results4 = run_step4()
    
    print("\n" + "="*80)
    print("WORKFLOW COMPLETE")
    print("="*80)
    print(f"\nAll results saved to: {OUTPUT_DIR}")
    print("\nKey outputs:")
    print("  - step1_eight_vertex/eight_vertex_data.npz")
    print("  - step2_bloch_rotation/bloch_rotation_data.npz")
    print("  - step3_full_chain/full_chain_data.npz & E0_vs_L.png")
    print("  - step4_summary/mzm_abs_summary.txt & mzm_abs_summary.png")
    print("\n")

if __name__ == '__main__':
    main()
