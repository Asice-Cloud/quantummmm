"""Paper parameters and common settings for reproducing Chen et al. FIGs.

Units: set energy unit such that Delta (superconducting gap) = 1.0.
Time units in the paper are often quoted as 100/Delta; with Delta=1 that
corresponds to 100.
"""
DELTA = 1.0
T_UNIT = 100.0  # 100 / DELTA

# Model energy parameters (in units of DELTA)
T0 = T_UNIT
t0 = 10.0 * DELTA
UR = 2.0 * DELTA
mu0 = -2.0 * t0
VD = 0.1 * t0

# Chain geometry
L = 100
QD_WIDTH = 10

# Figure specific parameter sets
FIG2_TS = [400.0, 450.0, 500.0]  # T values used in Fig.2 (in units where DELTA=1)
FIG3_T = 100.0
FIG3_N_PER_STEP = 300

# Fig4 / Fig5 magnetic modulation params (in units of DELTA)
FIG4_VX0 = 0.59 * DELTA
FIG4_VX1 = 0.02 * DELTA
FIG5_VX1_OPTIONS = [0.01 * DELTA, 0.03 * DELTA]

OUTDIR = 'results'
