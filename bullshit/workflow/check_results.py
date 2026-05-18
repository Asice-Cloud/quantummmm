#!/usr/bin/env python3
"""
Quick utility to check and summarize workflow results.
"""

from pathlib import Path
import sys

def main():
    results_root = Path(__file__).parent.parent / "results"
    workflow_dir = results_root / "workflow"
    
    print("\n" + "="*80)
    print("WORKFLOW RESULTS SUMMARY")
    print("="*80 + "\n")
    
    if not workflow_dir.exists():
        print("❌ No results found. Run 'python run_full_workflow.py' first.\n")
        return
    
    # Check each step
    steps = [
        ("Step 1", "step1_eight_vertex/eight_vertex_data.npz"),
        ("Step 2", "step2_bloch_rotation/bloch_rotation_data.npz"),
        ("Step 3", "step3_full_chain/full_chain_data.npz"),
        ("Step 4", "step4_summary/mzm_abs_summary.txt"),
    ]
    
    for step_name, rel_path in steps:
        full_path = workflow_dir / rel_path
        if full_path.exists():
            size_mb = full_path.stat().st_size / (1024**2)
            print(f"✓ {step_name}: {full_path.name} ({size_mb:.2f} MB)")
        else:
            print(f"✗ {step_name}: Missing")
    
    # Check figures
    print("\nGenerated Figures:")
    for png_file in sorted(workflow_dir.glob("*/*.png")):
        rel = png_file.relative_to(workflow_dir)
        print(f"  • {rel}")
    
    # Show summary
    summary_file = workflow_dir / "step4_summary/mzm_abs_summary.txt"
    if summary_file.exists():
        print(f"\nDiagnosis Summary (from {summary_file.name}):")
        print("-"*80)
        with open(summary_file) as f:
            lines = f.readlines()
            for line in lines[5:25]:  # Show first 20 diagnostic lines
                print(line.rstrip())
        if len(lines) > 25:
            print(f"... ({len(lines)-25} more lines)")
    
    print("\n" + "="*80)
    print("To view results in detail:")
    print(f"  cat {summary_file}")
    print(f"\nTo re-run with different parameters:")
    print(f"  Edit config.py, then run run_full_workflow.py")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
