# Strict Validation — Detailed Analysis

## Summary statistics
- K_fro: mean=1.479727e+01, std=1.199352e+01, min=2.011871e+00, max=4.095022e+01
- K_spec: mean=1.476454e+01, std=1.200137e+01, min=2.011871e+00, max=4.093612e+01
- K_nuc: mean=1.558045e+01, std=1.223825e+01, min=2.011877e+00, max=4.201088e+01
- N_fro: mean=2.321001e+00, std=1.653032e+00, min=9.541865e-07, max=4.984952e+00
- N_spec: mean=2.310251e+00, std=1.679257e+00, min=5.029265e-05, max=4.977076e+00
- N_nuc: mean=1.974745e+00, std=1.352902e+00, min=7.396113e-02, max=5.143880e+00

## Correlations
- Pearson(K_fro, N_fro) = 2.747279e-02
- Pearson(K_spec, N_spec) = 3.439492e-02

## Numerical stability notes
- Number of E1 points with |N_fro| < 1e-6: 1 / 40
- Relative change stats for N when eta increased: mean=2.535863e+04, max=1.014236e+06

## Interpretation
- The Pearson correlation between K and N is small in this run, indicating that scalar proxies alone do not fully determine the Δ Frobenius norm; the mapping is nontrivial and depends on phases and operator structure.
- Several E1 points produce N≈0; relative changes computed as ratio can be huge as reported earlier. Use absolute differences or filter near‑zero when reporting relative stability.
- Despite the noisy relative-change metric, the overall shapes of N vs E1 are similar across K measures, supporting qualitative robustness.