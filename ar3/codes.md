Project scripts and outputs — summary

This file lists the analysis scripts we created for the diagonal-K YBE study, what each does, and the main outputs produced during the session.

- scripts/convert_factors_trig.py
  - Purpose: Convert representative exponential factor expressions (e^{i·}) into cos/sin form without full expand, dedupe factor list.
  - Outputs: scripts/converted_factors.txt (representative converted factor expressions), scripts/converted_factors_out.txt (auxiliary output).

- scripts/grid_search_ybe.py
  - Purpose: Numeric coarse grid search over (Jx,Jy,Jz) (π/4 steps) for full YBE residual norm.
  - Outputs: scripts/ybe_grid_results_coarse.csv (grid points with residual norms).

- scripts/solve_full_ybe_diag.py
  - Purpose: Build R = exp(i K) for K diagonal (c00,I; Jx XX; Jy YY; Jz ZZ) and factor residual matrix entries symbolically.
  - Outputs: scripts/full_ybe_diag_out.txt (symbolic factorization output and nonzero entry listing).

- scripts/symbolic_solve_theta.py
  - Purpose: Symbolic factorization of representative symmetric exponential polynomial (U) and derive the B-equation branch form.
  - Outputs: scripts/symbolic_solve_theta_out.txt (factorization / B expression).

- scripts/verify_discrete_solutions.py
  - Purpose: Evaluate representative exponential polynomials on the discrete θ-grid {0, π/2, π, 3π/2} (equivalently J ∈ nπ/4) and enumerate valid solutions.
  - Outputs: scripts/discrete_theta_solutions.txt.

- scripts/groebner_prove.py
  - Purpose: Use algebraic variables A=exp(i Jx), B=exp(i Jy), C=exp(i Jz) and Groebner bases to test whether the candidate factor ideal F covers the residual ideal Res (completeness test).
  - Outputs: scripts/groebner_proof_out.txt (Groebner bases and remainder checks). Result: Res ⊆ Ideal(F) and product(F) ∈ Ideal(Res) (see file for details and caveats about phase equivalences).

- scripts/representative_spectra.py
  - Purpose: Compute K and R eigenvalues and residual norms for selected representative points to inspect spectral structure.
  - Outputs: scripts/representative_spectra.csv (saved representative eigen/spectrum info printed during runs).

- scripts/majorana_jw_mapping.py
  - Purpose: Construct Jordan–Wigner fermionic operators for the 2-site model, form Majorana operators, extract the antisymmetric Majorana matrix A and compute its spectrum to count zero modes.
  - Outputs: Printed per-point A matrix, eigenvalues, and zero-mode counts (used as input for plotting scripts).

- scripts/plot_majorana_zero_modes.py
  - Purpose: Read the coarse-grid solutions (residual≈0), compute zero-mode counts via `majorana_jw_mapping.analyze_point`, and produce visualizations.
  - Outputs (images): scripts/majorana_zero_modes_3d.png and scripts/majorana_zero_modes_Jz_*.png (per-Jz plane views).

Other workspace changes
- ybe_re.md — updated to include: factorization summary, Groebner result summary and a Visualization section linking the generated PNGs.

Notes and caveats
- Groebner-based conclusions are algebraic (ideal membership). Converting those algebraic results to phase-space (J modulo 2π) requires careful exclusion of spurious algebraic roots and consideration of |A|=|B|=|C|=1; we validated representative points numerically and with discrete grids to reduce false positives.
- Scripts assume the project's Python venv is active (SymPy, NumPy, Matplotlib available). Run scripts via the project's venv Python binary, for example:

  /home/asice-cloud/projects/pyyy/quantumsss/.venv/bin/python scripts/plot_majorana_zero_modes.py

If you want, I can:
- embed the generated PNGs into `ybe_re.md` more formally (figures + captions),
- generate a CSV report of zero-mode counts for all grid points, or
- convert the plots to interactive HTML (Plotly) for browser exploration.

End of summary.
