This tool computes a Hermitian generator H for a 2-site unitary R (R = exp(i H)),
decomposes H into the Pauli basis (σ^μ ⊗ σ^ν), maps selected coefficients to
Kitaev-like parameters (t, Δ, U, μ), and suggests a simple pulse when the
XX+YY channel is dominant.

Usage:

```bash
python3 tools/compute_R_realization.py
```

Output:
- prints the Pauli coefficients and mapped parameters
- writes `tools/R_decomposition_output.json` with numeric data

Notes:
- The script uses `scipy.linalg.logm` and `expm`.
- For a user-specified R, replace the test R construction in the script.
