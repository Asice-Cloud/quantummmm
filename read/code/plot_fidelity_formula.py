"""Evaluate the explicit scalar formula and compare with PRB 111 Fig. 1(d).

The plotted heatmap is evaluated by nested-integral quadrature (see analytic_continuation.md), not
by solving the original Schrödinger equation. Independent ODE samples only
check the formula implementation. No parameters or time scales are fitted.
"""

import json
import math
import os
from pathlib import Path
import subprocess
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/quantumsss_fidelity_mpl")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from PIL import Image
from scipy.integrate import cumulative_simpson, solve_ivp
from scipy.signal import find_peaks
from scipy.special import ellipe, ellipk

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE.parent / "figures" / "fidelity_formula"
DATA = HERE.parent / "data"
K = ellipk(.5) / np.pi
J = ellipe(.5) / np.pi
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
B = (np.eye(2) + 1j * X) / np.sqrt(2)


def formula(ratio, lam, samples=2049, order=9):
    """Uniform ordered-integral formula, C=I1-I3+... through the chosen odd order.

    I_m(s) = integral_0^s h(u) conjugate(I_{m-1}(u)) du, I_0=1,
    h=kappa exp(2i lambda S). This is nested quadrature of explicit kernels.
    Even coefficients are retained to reconstruct the normalized functions C,S.
    """
    lam = np.atleast_1d(lam).astype(float)
    s = np.linspace(0., 1., samples)
    a = (1 - np.cos(np.pi * s)) / 2
    b = 1 - a
    q = np.hypot(a, b)
    g = a / q
    ap = np.pi * np.sin(np.pi * s) / 2
    gp = ap * b / q**3
    rho = np.hypot(ratio * b, g)
    kappa = ratio * (b * gp + ap * g) / (2 * rho**2)
    action = cumulative_simpson(rho, x=s, initial=0.)
    kernel = kappa[None, :] * np.exp(2j * lam[:, None] * action)
    coefficient = np.ones_like(kernel)
    odd = np.zeros(len(lam), dtype=complex)
    even = np.ones(len(lam), dtype=complex)
    fifth = None
    for m in range(1, order + 1):
        coefficient = cumulative_simpson(
            kernel * coefficient.conj(), x=s, axis=1, initial=0.)
        if m % 2:
            odd += (-1)**((m - 1)//2) * coefficient[:, -1]
        else:
            even += (-1)**(m//2) * coefficient[:, -1]
        if m == 5:
            fifth = odd.copy()
    phase = lam * (action[-1] + 1 - ratio*J + K)
    def readout(c):
        return 4*(1 - abs(c)**2)*(np.exp(-1j*phase)*c).real**2
    f = readout(odd)
    beta = -1j * np.exp(-1j*lam*action[-1]) * odd
    nu = ratio*lam
    zphase = np.exp(-1j*lam*action[-1])
    # Convert the same ordered-integral series back to the normalized C,S convention.
    c_r7 = np.exp(.5j*nu)*(zphase*even + zphase.conj()*odd.conj())/np.sqrt(2)
    s_r7 = 1j*np.exp(.5j*nu)*(zphase*odd - zphase.conj()*even.conj())/np.sqrt(2)
    beta_r7 = -(np.exp(-.5j*nu)*s_r7 + 1j*np.exp(.5j*nu)*c_r7.conj())/np.sqrt(2)
    return {
        "f": f, "f5": readout(fifth), "p": abs(beta)**2,
        "visibility": 4*abs(beta)**2*(1-abs(beta)**2),
        "beta": beta, "C": c_r7, "S": s_r7,
        "r7_beta_error": float(np.max(abs(beta-beta_r7))),
        "unitarity_residual": float(np.max(abs(abs(c_r7)**2+abs(s_r7)**2-1))),
    }


def slow_ode(ratio, lam):
    """Independent check at a few points; never used for the heatmap."""
    def rhs(s, flat):
        a = (1-np.cos(np.pi*s))/2
        b = 1-a
        g = a/np.hypot(a, b)
        return (-1j*lam*(ratio*b*Z+g*Y) @ flat.reshape(2, 2)).ravel()
    sol = solve_ivp(rhs, (0., 1.), np.eye(2, dtype=complex).ravel(),
                    method="DOP853", rtol=2e-12, atol=2e-14)
    if not sol.success:
        raise RuntimeError(sol.message)
    d2 = sol.y[:, -1].reshape(2, 2)
    beta = (B.conj().T @ d2)[0, 1]
    delta = lam*(1+K-ratio*J)
    f = 4*(1-abs(beta)**2)*np.imag(np.exp(-1j*delta)*beta)**2
    return float(f), beta


def paper_image():
    """Render the actual PDF page and crop the panel, without rescaling axes.

    Crop positions are measured on the 200-dpi page; larger rendering changes
    image resolution only. Digitization is auxiliary and not original data.
    """
    target = Path("/tmp/quantumsss_prb111_page2")
    if not target.with_suffix(".png").exists():
        subprocess.run(["pdftoppm", "-f", "2", "-l", "2", "-singlefile",
                        "-r", "300", "-png", str(ROOT/"arts/prb111.pdf"),
                        str(target)], check=True)
    page = Image.open(target.with_suffix(".png")).convert("RGB")
    scale = page.width/1712
    def crop(box):
        return np.array(page.crop(tuple(round(v*scale) for v in box)))
    panel = crop((439, 365, 821, 665))
    Image.fromarray(panel).save(OUT/"paper_fig1d.png")
    heatmap = crop((482, 392, 753, 620))
    colorbar = crop((756, 393, 761, 619))
    palette = np.median(colorbar, axis=1)
    cmap = ListedColormap(palette[::-1]/255.)
    return heatmap, cmap


def main():
    start = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    paper, cmap = paper_image()
    lam = np.linspace(0., 12., 601)
    # Resolve the dense upper fringes without oversampling the broad lower lobes.
    logratios = np.unique(np.round(np.r_[np.linspace(-1., 1., 161),
                                       np.linspace(.5, 1., 161)], 12))
    ratios = 10**logratios
    data_path = DATA/"formula_grid.npz"
    report = {"parameters": {"E1_meV": .01, "tc_meV": .3, "E0_meV": .3,
                             "cycles": 2, "horizontal_axis": "tau/(100/meV) = lambda",
                             "vertical_axis": "log10(|t1|max/E1)"},
              "formula": "Uniform ordered integrals through I9; normalized scalar continuation functions",
              "quadrature_samples": 2049, "fit_parameters": 0,
              "checks": [], "curves": []}
    f = np.empty((len(ratios), len(lam)))
    f5, visibility = np.empty_like(f), np.empty_like(f)
    norm_error = 0.
    cache = {}
    if data_path.exists():
        with np.load(data_path) as saved:
            assert np.array_equal(saved["lam"], lam)
            cache = {round(float(y), 12): (saved["f"][i].copy(), saved["f5"][i].copy(),
                                           saved["visibility"][i].copy())
                     for i, y in enumerate(saved["logratios"])}
            norm_error = float(saved["unitarity_residual"])
        print(f"Reusing {len(cache)} completed formula rows", flush=True)
    for i, (ratio, y) in enumerate(zip(ratios, logratios)):
        if y in cache:
            f[i], f5[i], visibility[i] = cache[y]
            continue
        result = formula(ratio, lam)
        f[i], f5[i], visibility[i] = result["f"], result["f5"], result["visibility"]
        norm_error = max(norm_error, result["unitarity_residual"])
        if i % 20 == 0:
            print(f"Formula grid {i+1}/{len(ratios)}; elapsed {time.monotonic()-start:.1f}s", flush=True)
    np.savez_compressed(data_path, lam=lam, logratios=logratios, f=f, f5=f5,
                        visibility=visibility, unitarity_residual=norm_error)
    report["grid_f_range"] = [float(f.min()), float(f.max())]
    report["grid_max_f9_minus_f5"] = float(np.max(abs(f-f5)))
    report["grid_unitarity_residual"] = norm_error
    a0 = np.pi/4
    report["analytic_f9_truncation_bound"] = (
        8*a0**11/math.factorial(11)/(1-a0*a0/(12*13)))
    report["analytic_f5_truncation_bound"] = (
        8*a0**7/math.factorial(7)/(1-a0*a0/(8*9)))

    for ratio in (.1, 1., (1+2*K)/(.5+J), 10.):
        probes = np.array([0., .15, 1., 3., 7., 12.])
        coarse = formula(ratio, probes, samples=2049)
        fine = formula(ratio, probes, samples=4097, order=13)
        for i, value in enumerate(probes):
            ref, beta_ref = slow_ode(ratio, value)
            report["checks"].append({"ratio": ratio, "lambda": float(value),
                "f9_2049_vs_f13_4097": float(abs(coarse["f"][i]-fine["f"][i])),
                "f9_vs_independent_slow_ODE": float(abs(coarse["f"][i]-ref)),
                "beta9_vs_independent_slow_ODE": float(abs(coarse["beta"][i]-beta_ref)),
                "r7_conversion_error": coarse["r7_beta_error"]})
    report["max_probe_f_error"] = max(row["f9_vs_independent_slow_ODE"] for row in report["checks"])
    report["max_probe_refinement_error"] = max(row["f9_2049_vs_f13_4097"] for row in report["checks"])
    if report["max_probe_f_error"] > 2e-6:
        raise RuntimeError(f"Formula check failed: {report['max_probe_f_error']}")

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.1), constrained_layout=True)
    axes[0].imshow(paper, extent=(0, 12, -1, 1), aspect="auto", interpolation="nearest")
    axes[0].set_title("Paper Fig. 1(d): original PDF colors")
    image = axes[1].pcolormesh(lam, logratios, f, shading="nearest", vmin=0, vmax=1, cmap=cmap)
    axes[1].set_title("Our scalar formula: ordered integrals through $I_9$")
    for ax in axes:
        ax.set(xlabel=r"$\tau\ (100/\mathrm{meV})$", ylabel=r"$\log_{10}(|t_1|_{\max}/E_1)$",
               xticks=np.arange(0, 13, 2), yticks=np.linspace(-1, 1, 5),
               xlim=(0, 12), ylim=(-1, 1))
    fig.colorbar(image, ax=axes[1], label=r"$f=|\langle\psi_1^+(0)|\psi_1^-(6\tau)\rangle|$", shrink=.9)
    fig.suptitle(r"Same labeled axes and parameters: $E_1=0.01$, $t_c=E_0=0.3$ meV; no fitted rescaling", fontsize=13)
    fig.savefig(OUT/"paper_vs_formula.png", dpi=190)
    fig.savefig(OUT/"paper_vs_formula.pdf")
    plt.close(fig)

    # Quantitative checks use full SO(5) arrays with the same conventions.
    from full_model import ensure_full_model_cuts
    cached_path = ensure_full_model_cuts(DATA)
    cached_report = json.loads((DATA/"full_model_cuts.json").read_text())
    assert cached_report["parameters"] == {"E1": .01, "tc": .3, "E0": .3, "K": "2A", "cycles": 2}
    with np.load(cached_path) as cached:
        assert np.allclose(cached["tau"], lam*100)
        for ratio in (.1, 1., 2.3446, 10.):
            result = formula(ratio, lam)
            exact = cached[f"f_ratio_{ratio}"]
            use = lam >= .5
            peaks, _ = find_peaks(result["f"], prominence=.1)
            report["curves"].append({"ratio": ratio,
                "formula_peak_x_prominence_gt_0p1": lam[peaks].tolist(),
                "slow_vs_full_max_abs_x_ge_0p5": float(np.max(abs(result["f"][use]-exact[use]))),
                "slow_vs_full_rmse_x_ge_0p5": float(np.sqrt(np.mean((result["f"][use]-exact[use])**2)))})

    # Check cached full curves with a fresh evaluation of their defining series.
    from full_model import protocol_function, readout
    full_probes = np.array([150., 300., 600.])
    rotation, panels = protocol_function(full_probes, .001, order=30, local_phase=1.5)
    with np.load(cached_path) as saved:
        delta = readout(rotation)[0]-np.interp(full_probes, saved["tau"], saved["f_ratio_0.1"])
    report["full_matrix_cache_fresh_check"] = {"ratio": .1, "tau": full_probes.tolist(),
                                              "max_abs": float(np.max(abs(delta))), "panels": panels}
    assert np.max(abs(delta)) < 2e-8
    report["full_curves_source"] = str(cached_path.relative_to(ROOT))
    report["elapsed_seconds"] = time.monotonic()-start
    report["limitations"] = [
        "Paper panel is a qualitative comparison; no raw paper values are used.",
        "Quantitative reference is the full SO(5) model with the same conventions.",
        "The plotted scalar formula is exact only within the low-energy projection apart from stated truncation/quadrature errors.",
        "Full model curves use K=2A and the original gated t1 protocol.",
        "Caption E1=.01 is used; body text also states .001. Neither is silently adjusted.",
    ]
    (DATA/"verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({key: report[key] for key in ("grid_max_f9_minus_f5", "max_probe_f_error", "max_probe_refinement_error", "full_matrix_cache_fresh_check", "curves", "elapsed_seconds")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
