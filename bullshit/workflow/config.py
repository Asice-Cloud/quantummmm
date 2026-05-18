"""
Global configuration for the complete workflow.

Define all parameter ranges, paths, and constants here.
"""

import numpy as np
from pathlib import Path

# Repository root
REPO_ROOT = Path(__file__).parent.parent

# Output directory
OUTPUT_DIR = REPO_ROOT / "results" / "workflow"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============ Eight-Vertex Model Parameters ============
# Path segments (u, delta)
U_LIST_DEFAULT = [0.0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]  # 5 points
DELTA_LIST_DEFAULT = [0.0, 0.015, 0.1]  # 3 points

# Eight-vertex specific
T_STEP_DEFAULT = 100.0
N_PER_STEP_DEFAULT = 300

# ============ Full-chain BdG Parameters ============
# Chain lengths for edge_localization scan
L_LIST_DEFAULT = [40, 80, 160, 320]

# Bulk gap sampling
NK_BULK_DEFAULT = 801

# ============ Bloch Rotation Parameters ============
# Broadening for LDOS (Lorentzian eta)
ETA_LDOS_DEFAULT = 1e-3

# Energy range for LDOS visualization
E_MIN_LDOS = -0.05
E_MAX_LDOS = 0.05
NE_LDOS = 501

# ============ Output and Verbosity ============
VERBOSE = True
SAVE_FIGURES = True
FIG_DPI = 150
FIG_FORMAT = "png"

# ============ Directories ============
RESULTS_STEP1 = OUTPUT_DIR / "step1_eight_vertex"
RESULTS_STEP2 = OUTPUT_DIR / "step2_bloch_rotation"
RESULTS_STEP3 = OUTPUT_DIR / "step3_full_chain"
RESULTS_STEP4 = OUTPUT_DIR / "step4_summary"

for d in [RESULTS_STEP1, RESULTS_STEP2, RESULTS_STEP3, RESULTS_STEP4]:
    d.mkdir(parents=True, exist_ok=True)

print(f"[Config] Output directory: {OUTPUT_DIR}")
