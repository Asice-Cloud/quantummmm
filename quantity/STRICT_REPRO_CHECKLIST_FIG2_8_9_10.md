# Strict Reproduction Checklist (Fig.2, Fig.8, Fig.9, Fig.10)

## 0) Goal and Claim Boundary

This checklist is for **strict paper-level reproduction** of PRB 111, 205411 figures.

Pass condition is not just "looks similar". We require:
- Same physical observables on each panel axis.
- Same control variables (B, tau, Ed, disorder case).
- Quantitative agreement on key features within tolerances.

## 1) Priority Order

1. Fig.2 (uniform finite-length baseline calibration)
2. Fig.8 (main claim: inhomogeneous potential enhances braiding fidelity)
3. Fig.10 (disorder + fidelity + Ed scan robustness)
4. Fig.9 (6-case conductance classification consistency)

## 2) Mandatory Observable Mapping

Do not replace paper observables with surrogate quantities in strict mode.

- Paper braiding panel: use fidelity/overlap vs tau.
- Paper conductance panel: use differential conductance map G(E,B).
- Paper spectrum panel: eigenenergy spectrum vs B or Ed as in panel.
- Effective fit panel: same extracted effective model quantities (E1, t1) shown against B.

Current non-strict surrogate (kernel N) can be computed additionally, but must not replace fidelity in strict panels.

## 3) Figure-by-Figure Checklist

## Fig.2 (baseline)

Panels to reproduce:
- (a) uniform nanowire schematic with Delta(x), V(x)
- (b) spectrum vs B + local estimator eta(B)
- (c) extracted E1(B) and effective t1(B)
- (d)-(f) braiding fidelity vs tau for 3 marked B values
- (g)-(i) dot-nanowire spectrum vs Ed at those 3 B values

Strict checks:
- B window and 3 marker fields match caption values.
- In (d)-(f), effective-model dashed and nanowire solid must overlap trend-wise and numerically.
- In (g)-(i), crossing/anticrossing pattern and splitting scale match.

## Fig.8 (main claim)

Panels to reproduce:
- (a),(d) schematic + spectrum for uniform and inhomogeneous cases
- (b),(e) G(E,B) conductance maps for both cases
- (c),(f) braiding fidelity vs tau at selected B=1.9 T (paper triangle)

Strict checks:
- Both (b) and (e) exhibit ZBP-like structure.
- Fidelity ordering must match claim: inhomogeneous case has stronger fidelity retention than uniform at long tau.

## Fig.9 (six scenarios)

Panels to reproduce:
- 6 conductance maps for the six structure classes in caption.

Strict checks:
- White dashed transition marker placement consistent with each panel.
- Each panel must preserve paper-distinct morphology (not six nearly identical maps).

## Fig.10 (disorder robustness)

Panels to reproduce:
- (a) disorder schematic
- (b) spectrum vs B + eta(B)
- (c) conductance map G(E,B)
- (d)-(f) braiding fidelity vs tau at three B markers
- (g)-(i) spectrum vs Ed at corresponding B markers

Strict checks:
- Disorder case must still produce panel-dependent fidelity behavior.
- (g)-(i) use Ed scan (not delta scan).

## 4) Quantitative Acceptance Criteria

Use these as pass/fail gates per figure.

- Feature position error (peaks/crossings/transition B): <= 5%
- Fidelity min/max/mean error in each panel: <= 10%
- Conductance map structural similarity (SSIM or equivalent): >= 0.85
- Grid refinement stability (time step and B/E grids doubled):
  - qualitative ordering unchanged
  - key metrics drift <= 5%

## 5) Internal Consistency Tests (must pass before paper comparison)

- Unitarity check: ||U^\dagger U - I||_F < 1e-8 per representative run.
- Time-step convergence for fidelity at fixed parameters.
- Sign-convention sensitivity check for Majorana couplings (document chosen convention).
- Initial-state definition fixed and consistent with paper setup.

## 6) Existing Code to Reuse

- `quantity/reproduce_effective_braiding_majorana.py`
  - 6 Majorana effective Hamiltonian time evolution and fidelity core.
- `quantity/scan_braiding_enhanced.py`
  - E1 sweep and fidelity scan utilities.
- `quantity/reproduce_fig8_fig9_fig10_pauli.py`
  - Useful plotting structure only; strict mode must replace surrogate observables.

## 7) Immediate Next Implementation Steps

1. Build `strict_repro_fig2.py`:
   - exact panel mapping for (a)-(i)
   - export panel data tables (CSV/NPZ) for auditing
2. Add automated metric evaluator `strict_repro_metrics.py`:
   - feature extraction + threshold checks
3. Reuse same pipeline for Fig.8 and Fig.10, then Fig.9.

## 8) Reporting Template

For each figure, report:
- Inputs: all physical and numerical parameters.
- Outputs: generated panel file + panel-wise metrics.
- Verdict: PASS / PARTIAL / FAIL with reasons.
- Delta from paper: what mismatches remain and likely cause.
