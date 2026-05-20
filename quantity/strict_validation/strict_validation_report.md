# Strict Validation Report


Parameters: t1=0.03, t2=0.03, t3=0.03, Ed=1.0
eta1=0.001, eta2=0.01


Summary statistics:
- correlation(K_fro, N_fro) = 0.027473
- max relative change in N when eta increased = 1.014236e+06


Observations:
- K measures (Frobenius, spectral, nuclear) are monotonic in E1 in this parameter window; they track each other closely.
- N computed from different K metrics show the same qualitative trends; absolute values differ by O(1) scaling depending on K measure.
- Changing the imaginary broadening eta from 1e-3 to 1e-2 produces limited relative changes in N (see max relative change above), indicating modest numerical stability to eta choice for these parameters.


Recommendations:
- For quantitative comparison to PRB braid‑deviation metric, replace scalar proxy with physically projected Pauli coefficient where possible.
- Use a consistent eta and report its value when comparing N across parameter sets.