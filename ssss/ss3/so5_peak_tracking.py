r"""Track braiding-fidelity peaks and compare them with relative SO(5) error.

For each tau and t1/E1 ratio this computes the full two-cycle propagator,
the paper matrix element F_pq, and the relative endpoint element

    D(T) = R_0(T)^T R(T),

where R_0 is the same finite-time protocol with E1=t1=0.  The endpoint
group distance is evaluated from the five SO(5) eigenangles, so it remains a
well-defined principal distance even when scipy.linalg.logm changes branch.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.integrate import solve_ivp

from so5_error_decomposition import E1, XMAX, project_so5, omega

N_TAU = int(os.environ.get("SO5_PEAK_NTAU", "81"))
PROMINENCE = float(os.environ.get("SO5_PEAK_PROMINENCE", "0.01"))
RATIOS = tuple(float(x) for x in os.environ.get("SO5_PEAK_RATIOS", "-1,0,1").split(","))
N_RIDGE = int(os.environ.get("SO5_RIDGE_NY", "11"))


def paper_amplitude(r):
    vm = np.array([1, -1j, 0, 0, 0], dtype=complex) / np.sqrt(2.0)
    vp = np.array([1, 1j, 0, 0, 0], dtype=complex) / np.sqrt(2.0)
    return np.vdot(vp, r @ vm)


def principal_so5_distance(d):
    """Geodesic distance from I in SO(5), from the five eigenangles.

    If eigenvalues are exp(+-i theta_j), the Frobenius norm of the principal
    logarithm divided by sqrt(2) is sqrt(sum_j theta_j**2).  Using the full
    eigenvalue list gives the equivalent expression
    sqrt(0.5*sum_k angle(lambda_k)**2), including the possible eigenvalue 1.
    """
    vals = np.linalg.eigvals(project_so5(d))
    angles = np.angle(vals)
    return float(np.sqrt(0.5 * np.sum(angles * angles)))


def quadratic_peak(x0, y0, x1, y1, x2, y2):
    """Three-point parabolic peak interpolation."""
    den = y0 - 2.0 * y1 + y2
    if abs(den) < 1e-14:
        return x1, y1
    delta = 0.5 * (y0 - y2) / den
    dx = x1 - x0
    return x1 + delta * dx, y1 - 0.25 * (y0 - y2) * delta


def peak_positions(xs, fs):
    idx, _ = find_peaks(fs, prominence=PROMINENCE)
    out = []
    for i in idx:
        if i == 0 or i == len(xs) - 1:
            out.append((xs[i], fs[i]))
        else:
            out.append(quadratic_peak(xs[i - 1], fs[i - 1], xs[i], fs[i], xs[i + 1], fs[i + 1]))
    return np.asarray(out, dtype=float).reshape(-1, 2)


def solve_many_endpoint(tau, params):
    """Integrate several propagators together for one tau value."""
    n = len(params)
    if tau == 0:
        return [np.eye(5) for _ in params]

    def rhs(q, y):
        mats = y.reshape(n, 5, 5)
        out = np.empty_like(mats)
        for k, (e1, t1) in enumerate(params):
            out[k] = omega(q, tau, e1, t1) @ mats[k]
        return out.reshape(-1)

    sol = solve_ivp(
        rhs, (0.0, 6.0 * tau), np.tile(np.eye(5), (n, 1, 1)).reshape(-1),
        method="DOP853", rtol=2e-9, atol=2e-11,
    )
    return [project_so5(x) for x in sol.y[:, -1].reshape(n, 5, 5)]


def scan():
    taus = np.linspace(0.0, XMAX, N_TAU)
    f0 = np.zeros_like(taus)
    d0 = np.zeros_like(taus)
    f = np.zeros((len(RATIOS), len(taus)))
    d = np.zeros_like(f)
    params = [(0.0, 0.0)] + [(E1, E1 * 10.0**ratio) for ratio in RATIOS]
    for k, tau in enumerate(taus):
        mats = solve_many_endpoint(tau, params)
        r0 = mats[0]
        f0[k] = abs(paper_amplitude(r0))
        d0[k] = principal_so5_distance(np.eye(5))
        for j, ratio in enumerate(RATIOS):
            r = mats[j + 1]
            f[j, k] = abs(paper_amplitude(r))
            d[j, k] = principal_so5_distance(r0.T @ r)
        if k % max(1, len(taus) // 10) == 0:
            print(f"scan {k + 1}/{len(taus)}: tau={tau:.3f}", flush=True)
    peaks0 = peak_positions(taus, f0)
    peaks = [peak_positions(taus, row) for row in f]
    np.savez(
        "so5_peak_tracking.npz", tau=taus, ratios=np.asarray(RATIOS),
        F_ideal=f0, F_defect=f, d_so5=d,
        peaks_ideal=peaks0,
        peaks_defect=np.array(peaks, dtype=object), E1=E1,
    )
    return taus, f0, f, d, peaks0, peaks


def plot_results(taus, f0, f, d, peaks0, peaks):
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)
    axes[0].plot(taus, f0, color="k", lw=2, label="ideal: $E_1=t_1=0$")
    for j, ratio in enumerate(RATIOS):
        label = rf"$\log_{{10}}(t_1/E_1)={ratio:g}$"
        axes[0].plot(taus, f[j], label=label)
        p = peaks[j]
        if len(p):
            axes[0].scatter(p[:, 0], p[:, 1], s=16)
    axes[0].set_ylabel(r"$F_{\rm pq}=|\langle v_+|R|v_-\rangle|$")
    axes[0].set_ylim(0, 1.05)
    axes[0].grid(alpha=.25)
    axes[0].legend(fontsize=8, ncol=2)
    for j, ratio in enumerate(RATIOS):
        axes[1].plot(taus, d[j], label=rf"$\log_{{10}}(t_1/E_1)={ratio:g}$")
    axes[1].set_xlabel(r"$\tau/(100\,\mathrm{meV}^{-1})$")
    axes[1].set_ylabel(r"$d_{SO(5)}(R_0^TR)$")
    axes[1].grid(alpha=.25)
    axes[1].legend(fontsize=8, ncol=2)
    fig.suptitle(rf"Peak tracking and relative group error, $E_1={E1:g}$ meV")
    fig.savefig("so5_peak_tracking.png", dpi=220)


def scan_ridge():
    """Track the global highest-fidelity peak as a function of t1/E1."""
    taus = np.linspace(0.0, XMAX, N_TAU)
    ratios = np.linspace(-1.0, 1.0, N_RIDGE)
    f = np.zeros((len(ratios), len(taus)))
    d = np.zeros_like(f)
    params = [(0.0, 0.0)] + [(E1, E1 * 10.0**ratio) for ratio in ratios]
    for k, tau in enumerate(taus):
        mats = solve_many_endpoint(tau, params)
        r0 = mats[0]
        for j in range(len(ratios)):
            r = mats[j + 1]
            f[j, k] = abs(paper_amplitude(r))
            d[j, k] = principal_so5_distance(r0.T @ r)
    tau_peak = np.zeros(len(ratios))
    f_peak = np.zeros(len(ratios))
    d_peak = np.zeros(len(ratios))
    for j in range(len(ratios)):
        i = int(np.argmax(f[j]))
        if 0 < i < len(taus) - 1:
            tau_peak[j], f_peak[j] = quadratic_peak(
                taus[i - 1], f[j, i - 1], taus[i], f[j, i], taus[i + 1], f[j, i + 1]
            )
            d_peak[j] = np.interp(tau_peak[j], taus, d[j])
        else:
            tau_peak[j], f_peak[j], d_peak[j] = taus[i], f[j, i], d[j, i]
    np.savez(
        "so5_peak_ridge.npz", ratio=ratios, tau_peak=tau_peak,
        F_peak=f_peak, d_peak=d_peak, E1=E1,
    )
    fig, ax1 = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
    ax1.plot(ratios, tau_peak, "o-", color="tab:blue")
    ax1.set_xlabel(r"$\log_{10}(t_1/E_1)$")
    ax1.set_ylabel(r"global peak $\tau_*$ /(100 meV$^{-1}$)", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax1.grid(alpha=.25)
    ax2 = ax1.twinx()
    ax2.plot(ratios, f_peak, "s--", color="tab:orange", label=r"$F(\tau_*)$")
    ax2.plot(ratios, d_peak, "^:", color="tab:green", label=r"$d_{SO(5)}(\tau_*)$")
    ax2.set_ylabel(r"peak fidelity / relative group distance")
    ax2.legend(loc="center right", fontsize=8)
    ax1.set_title(rf"Global peak ridge, $E_1={E1:g}$ meV")
    fig.savefig("so5_peak_ridge.png", dpi=220)
    print("\nGlobal peak ridge (ratio, tau_peak, F_peak, d_peak):")
    for row in zip(ratios, tau_peak, f_peak, d_peak):
        print("  %.3f  %.6f  %.8f  %.6f" % row)
    print("wrote so5_peak_ridge.png and so5_peak_ridge.npz")


def main():
    taus, f0, f, d, peaks0, peaks = scan()
    plot_results(taus, f0, f, d, peaks0, peaks)
    scan_ridge()
    print("\nIdeal peaks (tau, F):")
    for x, y in peaks0:
        print(f"  {x:.6f}  {y:.8f}")
    for ratio, p in zip(RATIOS, peaks):
        print(f"ratio log10(t1/E1)={ratio:g} peaks (tau, F):")
        for x, y in p:
            print(f"  {x:.6f}  {y:.8f}")
    print("wrote so5_peak_tracking.png and so5_peak_tracking.npz")


if __name__ == "__main__":
    main()
