#!/usr/bin/env python3
"""Full finite-angle SO(5) fidelity and its spectral decomposition.

No expansion in E1, t1 or accumulated angles is used. The time integration
uses fourth-order commutator-free group steps with analytic 5x5 skew-matrix
exponentials. One identical three-step cycle is propagated and squared.
"""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/quantumsss-matplotlib")
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm


def controls(stage, v, tc=.3, e0=.3):
    rise = (1 - np.cos(np.pi * v)) / 2
    fall = 1 - rise
    if stage == 0:
        return np.array([tc * rise, 0., e0 * fall, rise])
    if stage == 1:
        return np.array([tc * fall, tc * rise, 0., fall])
    return np.array([0., tc * fall, e0 * rise, 0.])


def skew(e, eta, t2, t3, ed):
    out = np.zeros((5, 5))
    for i, j, value in [(0, 1, e), (0, 4, eta), (3, 1, t2),
                        (3, 2, -t3), (3, 4, ed)]:
        out[i, j], out[j, i] = value, -value
    return out


def apply_skew(v, e, eta, t2, t3, ed):
    # All coefficients have shape (batch,), v has shape (batch,5,5).
    e, eta, t2, t3, ed = [a[:, None] for a in (e, eta, t2, t3, ed)]
    out = np.empty_like(v)
    out[:, 0] = e * v[:, 1] + eta * v[:, 4]
    out[:, 1] = -e * v[:, 0] - t2 * v[:, 3]
    out[:, 2] = t3 * v[:, 3]
    out[:, 3] = t2 * v[:, 1] - t3 * v[:, 2] + ed * v[:, 4]
    out[:, 4] = -eta * v[:, 0] - ed * v[:, 3]
    return out


def exp_apply(v, e, eta, t2, t3, ed):
    """Apply exp(B) via Cayley-Hamilton; eigenvalues of B are 0,+/-iw1,+/-iw2."""
    s = e*e + eta*eta + t2*t2 + t3*t3 + ed*ed
    p = (e*e + eta*eta)*t3*t3 + (e*ed - eta*t2)**2
    q1 = (s + np.sqrt(np.maximum(s*s - 4*p, 0))) / 2
    q2 = np.divide(p, q1, out=np.zeros_like(p), where=q1 > 0)
    gap = q1 - q2
    w1, w2 = np.sqrt(q1), np.sqrt(q2)
    sinc1, sinc2 = np.sinc(w1 / np.pi), np.sinc(w2 / np.pi)
    cosc1 = .5 * np.sinc(w1 / (2*np.pi))**2
    cosc2 = .5 * np.sinc(w2 / (2*np.pi))**2
    safe = gap > 1e-10 * np.maximum(s, 1e-300)
    a3 = np.divide(sinc2 - sinc1, gap, out=np.zeros_like(gap), where=safe)
    a4 = np.divide(cosc2 - cosc1, gap, out=np.zeros_like(gap), where=safe)
    small = q1 < 1e-3
    qs = q1 + q2
    qs2 = q1*q1 + q1*q2 + q2*q2
    qs3 = q1**3 + q1*q1*q2 + q1*q2*q2 + q2**3
    a3[small] = (1/6 - qs/120 + qs2/5040 - qs3/362880)[small]
    a4[small] = (1/24 - qs/720 + qs2/40320 - qs3/3628800)[small]
    repeated = ~safe & ~small
    if np.any(repeated):
        w = w1[repeated]
        a3[repeated] = (np.sin(w) - w*np.cos(w)) / (2*w**3)
        a4[repeated] = (1 - np.cos(w) - .5*w*np.sin(w)) / w**4
    a1 = sinc1 + q1*a3
    a2 = cosc1 + q1*a4
    out = v.copy()
    power = v
    for coeff in (a1, a2, a3, a4):
        power = apply_skew(power, e, eta, t2, t3, ed)
        out += coeff[:, None, None] * power
    return out


def cycle(tau, e1, eta, steps=400, tc=.3, e0=.3, progress=False):
    tau, e1, eta = np.broadcast_arrays(tau, e1, eta)
    shape = tau.shape
    tau, e1, eta = [a.ravel() for a in (tau, e1, eta)]
    result = np.broadcast_to(np.eye(5), (tau.size, 5, 5)).copy()
    scale = 2*tau/steps
    node1, node2 = .5 - np.sqrt(3)/6, .5 + np.sqrt(3)/6
    early, late = (3 + 2*np.sqrt(3))/12, (3 - 2*np.sqrt(3))/12
    for stage in range(3):
        for step in range(steps):
            c1 = controls(stage, (step + node1)/steps, tc, e0)
            c2 = controls(stage, (step + node2)/steps, tc, e0)
            for combined in (early*c1 + late*c2, late*c1 + early*c2):
                t2, t3, ed, gate = combined
                result = exp_apply(result, .5*scale*e1, scale*eta*gate,
                                   scale*t2, scale*t3, scale*ed)
        if progress:
            print(f"Completed cycle segment {stage+1}/3", flush=True)
    return result.reshape(shape + (5, 5))


def direct_fidelity(c):
    r = c @ c
    return ((r[..., 0, 0] - r[..., 1, 1])**2
            + (r[..., 0, 1] + r[..., 1, 0])**2) / 4


def spectral_fidelity(c):
    d = (c + np.swapaxes(c, -1, -2)) / 2
    eigenvalues, vectors = np.linalg.eigh(d)
    # Each cosine occurs twice; the remaining eigenvalue is 1.
    cosine = np.stack([eigenvalues[..., :2].mean(-1),
                       eigenvalues[..., 2:4].mean(-1)], axis=-1)
    angles = np.arccos(np.clip(cosine, -1, 1))
    weights = []
    for pair in (slice(0, 2), slice(2, 4)):
        v1, v2 = vectors[..., 0, pair], vectors[..., 1, pair]
        weights.append(np.sum(v2*v2 - v1*v1, axis=-1)
                       + 2j*np.sum(v1*v2, axis=-1))
    weights = np.stack(weights, axis=-1)
    amplitude = np.sum(weights * np.sin(angles)**2, axis=-1)
    return np.abs(amplitude)**2, angles, weights


def adaptive_cycle(tau, e1, eta):
    y = np.eye(5).ravel()
    for stage in range(3):
        def rhs(v, flat):
            t2, t3, ed, gate = controls(stage, v)
            return (2*tau*skew(e1, eta*gate, t2, t3, ed)
                    @ flat.reshape(5, 5)).ravel()
        sol = solve_ivp(rhs, (0, 1), y, method="DOP853", rtol=2e-12, atol=2e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
    return y.reshape(5, 5)


def parameter_gradients(tau, e1, eta):
    """Exact tangent equations at the specified finite parameters, not at zero."""
    from verify_so5_fidelity import LE, LT
    current = np.zeros((4, 5, 5))
    current[0] = np.eye(5)
    for stage in range(3):
        def rhs(v, flat):
            mats = flat.reshape(4, 5, 5)
            t2, t3, ed, gate = controls(stage, v)
            physical = skew(e1, eta*gate, t2, t3, ed)
            out = 2*tau*physical @ mats
            out[1] += 2*tau*LE @ mats[0]
            out[2] += 2*tau*gate*LT @ mats[0]
            out[3] += 2*physical @ mats[0]
            return out.ravel()
        sol = solve_ivp(rhs, (0, 1), current.ravel(), method="DOP853",
                        rtol=2e-12, atol=2e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        current = sol.y[:, -1].reshape(4, 5, 5)
    c = current[0]
    r = c @ c
    cm = np.array([1, -1j, 0, 0, 0])/np.sqrt(2)
    amp = cm.T @ r @ cm
    derivatives = []
    for dc in current[1:]:
        dr = dc @ c + c @ dc
        da = cm.T @ dr @ cm
        derivatives.append(float(2*np.real(amp.conjugate()*da)))
    return float(abs(amp)**2), np.array(derivatives), c


def bivector_fidelity(tau, e1, eta):
    """Propagate the oriented plane directly in the ten-dimensional adjoint."""
    indices = np.triu_indices(5, 1)
    initial = np.zeros(10)
    initial[0] = 1

    def matrix(values):
        out = np.zeros((5, 5))
        out[indices] = values
        return out-out.T

    current = initial
    for stage in range(6):
        def rhs(v, values):
            t2, t3, ed, gate = controls(stage % 3, v)
            k = 2*tau*skew(e1, eta*gate, t2, t3, ed)
            b = matrix(values)
            return (k @ b-b @ k)[indices]
        sol = solve_ivp(rhs, (0, 1), current, method="DOP853",
                        rtol=2e-12, atol=2e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        current = sol.y[:, -1]
    b = matrix(current)
    f = (np.sum(b[:2]**2)-2*b[0, 1])/4
    return float(f)


def analyze():
    samples = []
    for tau, eta in [(100., .001), (100., .01), (100., .02344645627),
                     (500., .02344645627), (1200., .1)]:
        f, gradient, c = parameter_gradients(tau, .01, eta)
        fs, angles, weights = spectral_fidelity(c)
        pieces = weights*np.sin(angles)**2
        samples.append({"tau": tau, "E1": .01, "t1": eta, "F": f,
                        "f": np.sqrt(f), "cycle_angles": angles.tolist(),
                        "weights_real": weights.real.tolist(),
                        "weights_imag": weights.imag.tolist(),
                        "amplitude_terms_real": pieces.real.tolist(),
                        "amplitude_terms_imag": pieces.imag.tolist(),
                        "gradient_E1_t1_tau": gradient.tolist(),
                        "spectral_residual": float(abs(fs-f))})
    tau, e1, eta = 100., .01, .01
    _, gradient, _ = parameter_gradients(tau, e1, eta)
    finite = []
    parameters = [tau, e1, eta]
    for index, epsilon in [(1, 1e-7), (2, 1e-7), (0, 1e-4)]:
        plus, minus = parameters.copy(), parameters.copy()
        plus[index] += epsilon
        minus[index] -= epsilon
        fp = direct_fidelity(adaptive_cycle(*plus))
        fm = direct_fidelity(adaptive_cycle(*minus))
        finite.append(float((fp-fm)/(2*epsilon)))
    assert np.max(np.abs(np.array(finite)-gradient)) < 1e-5
    bivector_value = bivector_fidelity(tau, e1, eta)
    reference_value = direct_fidelity(adaptive_cycle(tau, e1, eta))
    assert abs(bivector_value-reference_value) < 1e-10
    report = {"samples": samples,
              "adjoint_bivector_check": {"F_bivector": bivector_value,
                                         "F_reference": float(reference_value)},
              "gradient_check": {
        "analytic_E1_t1_tau": gradient.tolist(), "finite_difference": finite}}
    output = Path(__file__).with_name("so5_full_parameter_analysis.json")
    output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


def validate():
    rng = np.random.default_rng(235)
    errors = []
    for _ in range(5):
        values = rng.normal(size=5)
        b = skew(*values)
        actual = exp_apply(np.eye(5)[None], *[np.array([x]) for x in values])[0]
        errors.append(float(np.max(np.abs(actual - expm(b)))))
    assert max(errors) < 1e-12
    random = rng.normal(size=(5, 5))
    c = expm(random-random.T)
    fs, _, _ = spectral_fidelity(c)
    assert abs(fs-direct_fidelity(c)) < 1e-12
    cases = [(0., .01, .1), (1., .01, .001), (100., .01, .01),
             (500., .01, .02344645627), (1200., .01, .1)]
    comparisons = []
    for tau, e1, eta in cases:
        ref = adaptive_cycle(tau, e1, eta)
        row = {"tau": tau, "E1": e1, "t1": eta,
               "F_adaptive": float(direct_fidelity(ref)), "steps": []}
        for count in (200, 400, 800):
            result = cycle(tau, e1, eta, steps=count)
            fp = float(direct_fidelity(result))
            row["steps"].append({"count": count, "F": fp,
                                 "F_error": abs(fp-row["F_adaptive"]),
                                 "matrix_error": float(np.max(np.abs(result-ref)))})
        comparisons.append(row)
    for row in comparisons:
        assert row["steps"][-1]["F_error"] < 2e-5
    report = {"exponential_errors": errors, "comparisons": comparisons}
    Path(__file__).with_name("so5_full_validation.json").write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


def plot(steps, nx, ny):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    root = Path(__file__).resolve().parent
    x = np.linspace(0, 12, nx)
    y = np.linspace(-1, 1, ny)
    tau, eta = np.broadcast_arrays(100*x[None, :], .01*10**y[:, None])
    c = cycle(tau, .01, eta, steps=steps, progress=True)
    f_squared, angles, weights = spectral_fidelity(c)
    discrepancy = float(np.max(np.abs(f_squared-direct_fidelity(c))))
    orthogonality = float(np.max(np.abs(np.swapaxes(c, -1, -2) @ c - np.eye(5))))
    assert discrepancy < 2e-7
    assert orthogonality < 2e-7
    amplitude = np.sqrt(np.clip(f_squared, 0, 1))
    np.savez_compressed(root/"so5_full_fidelity_fig1d.npz", x=x, y=y,
                        f=amplitude, F=f_squared, cycle_angles=angles,
                        cycle_weights=weights, steps=steps, E1=.01, tc=.3, E0=.3,
                        spectral_identity_error=discrepancy,
                        orthogonality_error=orthogonality)
    cmap = LinearSegmentedColormap.from_list(
        "paper_blue_yellow", ["#24218c", "#285ac6", "#22a6cf", "#2bbf99", "#a8d94e", "#fff05a"])
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 11}):
        fig, ax = plt.subplots(figsize=(8.4, 5.8), layout="constrained")
        im = ax.pcolormesh(x, y, amplitude, shading="auto", cmap=cmap,
                           vmin=0, vmax=1, rasterized=True)
        ax.set(xlabel=r"$\tau\ (100/\mathrm{meV})$", ylabel=r"$\lg(|t_1|/E_1)$",
               xlim=(0, 12), ylim=(-1, 1))
        ax.set_title("Full SO(5) spectral fidelity", fontsize=14, pad=14)
        cb = fig.colorbar(im, ax=ax, pad=.03)
        cb.set_label(r"$f=|b_1\sin^2\theta_1+b_2\sin^2\theta_2|$")
        fig.supxlabel(r"$E_1=0.01$ meV, $t_c=E_0=0.3$ meV; two gated cycles; $\hbar=1$",
                      fontsize=10)
        for extension in ("png", "pdf"):
            fig.savefig(root/f"so5_full_fidelity_fig1d.{extension}", dpi=220)
        plt.close(fig)

        fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
        for ax, ratio in zip(axes.ravel(), (.1, 1., 2.344645627, 10.)):
            cc = cycle(100*x, .01, .01*ratio, steps=steps)
            ff, aa, ww = spectral_fidelity(cc)
            ax.plot(x, np.sqrt(np.clip(ff, 0, 1)), color="#243d8c", lw=1.2)
            ax.set(title=rf"$|t_1|/E_1={ratio:.5g}$", xlabel=r"$\tau\ (100/\mathrm{meV})$",
                   ylabel=r"$f$", xlim=(0, 12), ylim=(0, 1.02))
            ax.grid(alpha=.2)
        fig.suptitle("Full finite-angle parameter dependence", fontsize=14)
        for extension in ("png", "pdf"):
            fig.savefig(root/f"so5_full_fidelity_slices.{extension}", dpi=220)
        plt.close(fig)
    print(json.dumps({"grid": [nx, ny], "steps_per_segment": steps,
                      "spectral_identity_error": discrepancy,
                      "orthogonality_error": orthogonality,
                      "output": str(root/"so5_full_fidelity_fig1d.png")}, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--steps", type=int, default=800)
    parser.add_argument("--nx", type=int, default=201)
    parser.add_argument("--ny", type=int, default=241)
    args = parser.parse_args()
    if args.validate:
        validate()
    if args.analyze:
        analyze()
    if args.plot:
        plot(args.steps, args.nx, args.ny)
