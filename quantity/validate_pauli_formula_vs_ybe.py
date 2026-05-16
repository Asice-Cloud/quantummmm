#!/usr/bin/env python3
"""Validate Pauli-coefficient non-Abelian quantization formulas vs exact YBE metric.

Compares, under identical parameter scans:
1) Exact (paper-style numerical) metric:
      N_exact = sqrt(Tr(Delta^dagger Delta))
2) Coefficient-only heuristic formula from markdown:
      N_simple ~ sqrt(\int ds dt sum |h_ab(s)|^2 |h_mn(t)|^2 eps_{b m g}^2)
3) Refined 3rd-order kernel formula from markdown:
      ||Delta^(3)||_F^2 ~ 4*2^3 * sum |\int\int K(s,t) h_ab(s) h_mn(t) ds dt|^2 eps^2

Outputs:
- quantity/validate_formula_scan_lambda.png
- quantity/validate_formula_scan_jyratio.png
- quantity/validate_formula_scaling_smallJ.png
- quantity/validate_formula_metrics.txt
"""

import numpy as np
import matplotlib.pyplot as plt
from importlib import util


# Load existing exact numerical pipeline (already used in repo)
spec = util.spec_from_file_location(
    "pauli_repro", "quantity/pauli_tensor_nonabelian_reproduction_script.py"
)
pauli_repro = util.module_from_spec(spec)
spec.loader.exec_module(pauli_repro)


PAULI_KEYS = ["I", "x", "y", "z"]
PAULI_INDEX = {"I": 0, "x": 1, "y": 2, "z": 3}


def levi_civita(i, j, k):
    """Levi-Civita on {1,2,3}; index 0 (identity) gives 0."""
    if i == 0 or j == 0 or k == 0:
        return 0
    if (i, j, k) in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]:
        return 1
    if (i, j, k) in [(1, 3, 2), (3, 2, 1), (2, 1, 3)]:
        return -1
    return 0


def h_coeffs(u, pars):
    """Return dict h_{ab}(u) with keys in {I,x,y,z}^2 based on current model."""
    Jx = pars["Jx"]
    Jy = pars["Jy"]
    Jz = pars["Jz"]
    lam = pars["lambda_abs"]
    dm = pars["delta_mix"]

    out = {}
    out[("x", "x")] = Jx * np.cos(u)
    out[("y", "y")] = Jy * np.sin(u)
    out[("z", "z")] = Jz
    out[("x", "y")] = lam * np.sin(u)
    out[("y", "x")] = lam * np.cos(u)
    out[("x", "z")] = dm * np.sin(2.0 * u)
    out[("z", "y")] = dm * np.cos(2.0 * u)
    return out


def kernel_K(s, t, u, v):
    """K(s,t) from markdown theta-function expression."""
    theta_st = 1.0 if s > t else 0.0
    theta_ts = 1.0 if t > s else 0.0
    return v * theta_st + min(t, v) - u * theta_ts - max(0.0, u - t)


def N_simple_formula(pars, u=1.0, v=0.7, ns=140, nt=180):
    """Coefficient-only heuristic formula (shape predictor)."""
    s_grid = np.linspace(0.0, u, ns)
    t_grid = np.linspace(0.0, u + v, nt)
    ds = s_grid[1] - s_grid[0]
    dt = t_grid[1] - t_grid[0]

    total = 0.0
    for s in s_grid:
        hs = h_coeffs(s, pars)
        for t in t_grid:
            ht = h_coeffs(t, pars)
            inner = 0.0
            for (a, b), h1 in hs.items():
                bidx = PAULI_INDEX[b]
                for (m, n), h2 in ht.items():
                    midx = PAULI_INDEX[m]
                    eps2 = 0
                    for g in [1, 2, 3]:
                        e = levi_civita(bidx, midx, g)
                        eps2 += e * e
                    if eps2 != 0:
                        inner += (abs(h1) ** 2) * (abs(h2) ** 2) * eps2
            total += inner

    return np.sqrt(max(total * ds * dt, 0.0))


def N_kernel_formula(pars, u=1.0, v=0.7, ns=120, nt=160):
    """Refined third-order kernel formula from markdown."""
    s_grid = np.linspace(0.0, u, ns)
    t_grid = np.linspace(0.0, u + v, nt)
    ds = s_grid[1] - s_grid[0]
    dt = t_grid[1] - t_grid[0]

    sumsq = 0.0
    # Sum over tensor channels and Levi-Civita-compatible middle-site pairs.
    for a in PAULI_KEYS:
        for b in PAULI_KEYS:
            bidx = PAULI_INDEX[b]
            for m in PAULI_KEYS:
                midx = PAULI_INDEX[m]
                eps2 = 0
                for g in [1, 2, 3]:
                    e = levi_civita(bidx, midx, g)
                    eps2 += e * e
                if eps2 == 0:
                    continue
                for n in PAULI_KEYS:
                    integ = 0.0 + 0.0j
                    for s in s_grid:
                        hs = h_coeffs(s, pars).get((a, b), 0.0)
                        if hs == 0.0:
                            continue
                        for t in t_grid:
                            ht = h_coeffs(t, pars).get((m, n), 0.0)
                            if ht == 0.0:
                                continue
                            K = kernel_K(s, t, u, v)
                            integ += K * hs * ht
                    integ *= ds * dt
                    sumsq += eps2 * (abs(integ) ** 2)

    # markdown: ||Delta^(3)||_F^2 ~ 4 * 2^3 * sum ...
    val2 = 4.0 * (2.0 ** 3) * sumsq
    return np.sqrt(max(val2, 0.0))


def fit_scale(x, y):
    """Best linear scale y≈c*x in least squares (through origin)."""
    den = float(np.dot(x, x))
    if den <= 1e-30:
        return 0.0
    return float(np.dot(x, y) / den)


def pearsonr(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    x0 = x - x.mean()
    y0 = y - y.mean()
    den = np.linalg.norm(x0) * np.linalg.norm(y0)
    if den <= 1e-30:
        return 0.0
    return float(np.dot(x0, y0) / den)


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def run_scan_lambda(u=1.0, v=0.7):
    lam_scan = np.linspace(0.0, 1.0, 34)
    N_exact, N_s, N_k = [], [], []
    for lam in lam_scan:
        pars = {
            "Jx": 1.0,
            "Jy": 0.7,
            "Jz": 0.25,
            "lambda_abs": float(lam),
            "delta_mix": 0.10,
        }
        n_ex, _ = pauli_repro.compute_N(pars, u=u, v=v, steps=220)
        N_exact.append(n_ex)
        N_s.append(N_simple_formula(pars, u=u, v=v))
        N_k.append(N_kernel_formula(pars, u=u, v=v))
    return lam_scan, np.array(N_exact), np.array(N_s), np.array(N_k)


def run_scan_jyratio(u=1.0, v=0.7):
    r_scan = np.linspace(0.0, 2.0, 34)
    N_exact, N_s, N_k = [], [], []
    for r in r_scan:
        pars = {
            "Jx": 1.0,
            "Jy": float(r),
            "Jz": 0.25,
            "lambda_abs": 0.20,
            "delta_mix": 0.10,
        }
        n_ex, _ = pauli_repro.compute_N(pars, u=u, v=v, steps=220)
        N_exact.append(n_ex)
        N_s.append(N_simple_formula(pars, u=u, v=v))
        N_k.append(N_kernel_formula(pars, u=u, v=v))
    return r_scan, np.array(N_exact), np.array(N_s), np.array(N_k)


def run_scaling_smallJ(u=1.0, v=0.7):
    g_scan = np.linspace(0.15, 1.0, 16)
    N_exact, N_k = [], []
    for g in g_scan:
        pars = {
            "Jx": 1.0 * g,
            "Jy": 0.7 * g,
            "Jz": 0.25 * g,
            "lambda_abs": 0.20 * g,
            "delta_mix": 0.10 * g,
        }
        n_ex, _ = pauli_repro.compute_N(pars, u=u, v=v, steps=220)
        N_exact.append(n_ex)
        N_k.append(N_kernel_formula(pars, u=u, v=v))
    return g_scan, np.array(N_exact), np.array(N_k)


def exact_delta_matrix(pars, u=1.0, v=0.7, steps=260):
    """Return the exact YBE deviation matrix Delta for a parameter set."""
    Ru = pauli_repro.path_ordered_R(u, pars, steps=steps)
    Rv = pauli_repro.path_ordered_R(v, pars, steps=steps)
    Ruv = pauli_repro.path_ordered_R(u + v, pars, steps=steps)
    return pauli_repro.yb_deviation(Ru, Rv, Ruv)


def extract_third_order_exact(base_pars, g_scan=None, u=1.0, v=0.7, steps=260):
    """Extract Delta^(3) by fitting Delta(g) ≈ g^3 * Delta3 for small g.

    Uses a least-squares fit through the origin on the flattened complex matrix.
    """
    if g_scan is None:
        g_scan = np.array([0.02, 0.035, 0.05])

    deltas = []
    for g in g_scan:
        pars = {
            "Jx": base_pars["Jx"] * g,
            "Jy": base_pars["Jy"] * g,
            "Jz": base_pars["Jz"] * g,
            "lambda_abs": base_pars["lambda_abs"] * g,
            "delta_mix": base_pars["delta_mix"] * g,
        }
        deltas.append(exact_delta_matrix(pars, u=u, v=v, steps=steps))

    g3 = g_scan ** 3
    stacked = np.stack([d.reshape(-1) for d in deltas], axis=0)
    denom = np.sum(g3 ** 2)
    if denom <= 1e-30:
        return np.zeros_like(deltas[0]), g_scan, deltas

    # coefficient for each matrix entry: X = sum g^3 Delta(g) / sum g^6
    coeff = np.tensordot(g3, stacked, axes=(0, 0)) / denom
    delta3 = coeff.reshape(deltas[0].shape)
    return delta3, g_scan, deltas


def compare_kernel_vs_exact_third_order(base_pars, u=1.0, v=0.7, steps=260):
    """Compare the kernel formula against the extracted exact Delta^(3)."""
    delta3_exact, g_scan, deltas = extract_third_order_exact(base_pars, u=u, v=v, steps=steps)
    n3_exact = float(np.sqrt(np.real(np.trace(delta3_exact.conj().T @ delta3_exact))))
    n_kernel = float(N_kernel_formula(base_pars, u=u, v=v))

    # Also compute a best-fit scale between kernel estimate and exact third-order norm.
    scale = 0.0 if abs(n_kernel) < 1e-30 else n3_exact / n_kernel
    residual = abs(n3_exact - scale * n_kernel)

    return {
        "delta3_exact": delta3_exact,
        "n3_exact": n3_exact,
        "n_kernel": n_kernel,
        "scale": scale,
        "residual": residual,
        "g_scan": g_scan,
        "deltas": deltas,
    }


def main():
    os = __import__("os")
    os.makedirs("quantity", exist_ok=True)

    u, v = 1.0, 0.7

    # Scan A: lambda_abs
    lam, ex_l, s_l, k_l = run_scan_lambda(u=u, v=v)
    c_s_l = fit_scale(s_l, ex_l)
    c_k_l = fit_scale(k_l, ex_l)
    s_l_fit = c_s_l * s_l
    k_l_fit = c_k_l * k_l

    # Scan B: Jy/Jx
    rr, ex_r, s_r, k_r = run_scan_jyratio(u=u, v=v)
    c_s_r = fit_scale(s_r, ex_r)
    c_k_r = fit_scale(k_r, ex_r)
    s_r_fit = c_s_r * s_r
    k_r_fit = c_k_r * k_r

    # Scan C: small-coupling scaling
    gg, ex_g, k_g = run_scaling_smallJ(u=u, v=v)

    # log-log slope for scaling ex ~ g^p, kernel ~ g^q
    p_ex = np.polyfit(np.log(gg), np.log(np.maximum(ex_g, 1e-18)), 1)[0]
    p_k = np.polyfit(np.log(gg), np.log(np.maximum(k_g, 1e-18)), 1)[0]

    # Plot lambda scan
    plt.figure(figsize=(7, 5))
    plt.plot(lam, ex_l, "k", lw=2, label="Exact YBE metric")
    plt.plot(lam, s_l_fit, "--", lw=2, label="Pauli coeff formula (scaled)")
    plt.plot(lam, k_l_fit, "-.", lw=2, label="Kernel 3rd-order formula (scaled)")
    plt.xlabel(r"$\lambda_{ABS}$")
    plt.ylabel(r"$\mathcal{N}$")
    plt.title("Formula vs Exact: lambda scan")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("quantity/validate_formula_scan_lambda.png", dpi=200)
    plt.close()

    # Plot ratio scan
    plt.figure(figsize=(7, 5))
    plt.plot(rr, ex_r, "k", lw=2, label="Exact YBE metric")
    plt.plot(rr, s_r_fit, "--", lw=2, label="Pauli coeff formula (scaled)")
    plt.plot(rr, k_r_fit, "-.", lw=2, label="Kernel 3rd-order formula (scaled)")
    plt.xlabel(r"$J_y/J_x$")
    plt.ylabel(r"$\mathcal{N}$")
    plt.title("Formula vs Exact: channel-competition scan")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("quantity/validate_formula_scan_jyratio.png", dpi=200)
    plt.close()

    # Plot small-J scaling
    plt.figure(figsize=(7, 5))
    plt.loglog(gg, ex_g, "ko-", lw=2, label=f"Exact (slope~{p_ex:.2f})")
    plt.loglog(gg, k_g, "o-", lw=2, label=f"Kernel formula (slope~{p_k:.2f})")
    plt.xlabel("global coupling scale g")
    plt.ylabel(r"$\mathcal{N}$")
    plt.title("Small-coupling scaling check")
    plt.grid(alpha=0.3, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig("quantity/validate_formula_scaling_smallJ.png", dpi=200)
    plt.close()

    # Direct third-order comparison: exact Delta^(3) vs kernel formula
    lam3_scan = np.linspace(0.0, 1.0, 9)
    n3_exact_list = []
    n3_kernel_list = []
    n3_scale_list = []
    n3_residual_list = []
    for lam in lam3_scan:
        base = {
            "Jx": 1.0,
            "Jy": 0.7,
            "Jz": 0.25,
            "lambda_abs": float(lam),
            "delta_mix": 0.10,
        }
        cmpres = compare_kernel_vs_exact_third_order(base, u=u, v=v, steps=120)
        n3_exact_list.append(cmpres["n3_exact"])
        n3_kernel_list.append(cmpres["n_kernel"])
        n3_scale_list.append(cmpres["scale"])
        n3_residual_list.append(cmpres["residual"])

    n3_exact_list = np.array(n3_exact_list)
    n3_kernel_list = np.array(n3_kernel_list)
    n3_scale_list = np.array(n3_scale_list)
    n3_residual_list = np.array(n3_residual_list)

    n3_scale_fit = fit_scale(n3_kernel_list, n3_exact_list)
    n3_kernel_fit = n3_scale_fit * n3_kernel_list

    plt.figure(figsize=(7, 5))
    plt.plot(lam3_scan, n3_exact_list, "k", lw=2, marker="o", label=r"Exact $\|\Delta^{(3)}\|_F$")
    plt.plot(lam3_scan, n3_kernel_fit, "--", lw=2, marker="s", label=r"Kernel formula (scaled)")
    plt.xlabel(r"$\lambda_{ABS}$")
    plt.ylabel(r"$\|\Delta^{(3)}\|_F$")
    plt.title("Third-order exact vs kernel formula")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("quantity/validate_third_order_kernel_vs_exact.png", dpi=200)
    plt.close()

    # Metrics
    txt = []
    txt.append("Validation of Pauli-coefficient formulas vs exact YBE metric\n")
    txt.append(f"u={u}, v={v}\n\n")

    txt.append("[lambda scan]\n")
    txt.append(f"pearson(simple_scaled, exact) = {pearsonr(s_l_fit, ex_l):.6f}\n")
    txt.append(f"pearson(kernel_scaled, exact) = {pearsonr(k_l_fit, ex_l):.6f}\n")
    txt.append(f"rmse(simple_scaled, exact) = {rmse(s_l_fit, ex_l):.6f}\n")
    txt.append(f"rmse(kernel_scaled, exact) = {rmse(k_l_fit, ex_l):.6f}\n")
    txt.append(f"scale_simple = {c_s_l:.6f}, scale_kernel = {c_k_l:.6f}\n\n")

    txt.append("[Jy/Jx scan]\n")
    txt.append(f"pearson(simple_scaled, exact) = {pearsonr(s_r_fit, ex_r):.6f}\n")
    txt.append(f"pearson(kernel_scaled, exact) = {pearsonr(k_r_fit, ex_r):.6f}\n")
    txt.append(f"rmse(simple_scaled, exact) = {rmse(s_r_fit, ex_r):.6f}\n")
    txt.append(f"rmse(kernel_scaled, exact) = {rmse(k_r_fit, ex_r):.6f}\n")
    txt.append(f"scale_simple = {c_s_r:.6f}, scale_kernel = {c_k_r:.6f}\n\n")

    txt.append("[small-g scaling]\n")
    txt.append(f"slope_exact ~ {p_ex:.6f}\n")
    txt.append(f"slope_kernel_formula ~ {p_k:.6f}\n")

    txt.append("\n[third-order exact vs kernel scan]\n")
    txt.append(f"pearson(kernel_fit, exact_delta3_norm) = {pearsonr(n3_kernel_fit, n3_exact_list):.6f}\n")
    txt.append(f"rmse(kernel_fit, exact_delta3_norm) = {rmse(n3_kernel_fit, n3_exact_list):.6f}\n")
    txt.append(f"scale_fit = {n3_scale_fit:.6f}\n")
    txt.append(f"mean_residual_after_per_point_scaling = {float(np.mean(n3_residual_list)):.6f}\n")

    with open("quantity/validate_formula_metrics.txt", "w", encoding="utf-8") as f:
        f.writelines(txt)

    print("Saved:")
    print("  quantity/validate_formula_scan_lambda.png")
    print("  quantity/validate_formula_scan_jyratio.png")
    print("  quantity/validate_formula_scaling_smallJ.png")
    print("  quantity/validate_third_order_kernel_vs_exact.png")
    print("  quantity/validate_formula_metrics.txt")


if __name__ == "__main__":
    main()
