r"""Decompose E1,t1 errors relative to the ideal SO(5) path.

The ideal reference is the *same finite-time protocol* with E1=t1=0.
Thus D(q)=R0(q)^T R(q) removes the baseline non-adiabatic error and isolates
the defect-induced error.  For D=exp(Xi), Xi is projected onto

    h = so(3) \oplus so(2)       (block-preserving/internal error)
    m = MZM--QD mixing directions (leakage error).

The script produces endpoint maps versus tau and log10(t1/E1), plus time
traces for representative ratios.  Time is measured in q=t/(100 meV^-1),
matching so5_group_fidelity_map.py.
"""
from pathlib import Path
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.linalg import logm

TC = 0.3                         # meV
E0 = 0.3                         # meV
E1 = float(os.environ.get("SO5_E1", "0.01"))
SO5_NORMALIZATION = 2.0
YMIN, YMAX = -1.0, 1.0
XMAX = 12.0                      # tau/(100 meV^-1)
N_X = int(os.environ.get("SO5_ERR_NX", "31"))
N_Y = int(os.environ.get("SO5_ERR_NY", "31"))
N_PATH = int(os.environ.get("SO5_ERR_NPATH", "401"))


def env(u, rising):
    return (1.0 - np.cos(np.pi * u)) / 2 if rising else (1.0 + np.cos(np.pi * u)) / 2


def protocol_values(q, tau):
    """Return (t2,t3,Ed,g1) in meV for dimensionless q and tau."""
    if tau == 0:
        return 0.0, 0.0, 0.0, 0.0
    z = q % (3.0 * tau)
    seg = min(2, int(z / tau))
    u = (z - seg * tau) / tau
    if seg == 0:
        return TC * env(u, True), 0.0, E0 * env(u, False), env(u, True)
    if seg == 1:
        return TC * env(u, False), TC * env(u, True), 0.0, env(u, False)
    return 0.0, TC * env(u, False), E0 * env(u, True), 0.0


def omega(q, tau, e1, t1):
    """Generator dR/dq = O(q) R, with q=t/(100 meV^-1)."""
    t2, t3, ed, g1 = protocol_values(q, tau)
    s = SO5_NORMALIZATION
    o = np.zeros((5, 5), dtype=float)
    o[0, 1], o[1, 0] = s * e1, -s * e1
    o[3, 1], o[1, 3] = s * t2, -s * t2
    o[4, 0], o[0, 4] = -s * t1 * g1, s * t1 * g1
    o[3, 2], o[2, 3] = -s * t3, s * t3
    o[3, 4], o[4, 3] = s * ed, -s * ed
    return 100.0 * o


def pack_pair(r0, r):
    return np.concatenate((r0.reshape(-1), r.reshape(-1)))


def unpack_pair(y):
    return y[:25].reshape(5, 5), y[25:].reshape(5, 5)


def rhs_pair(q, y, tau, e1, t1):
    r0, r = unpack_pair(y)
    o0 = omega(q, tau, 0.0, 0.0)
    oa = omega(q, tau, e1, t1)
    return pack_pair(o0 @ r0, oa @ r)


def project_so5(a):
    """Nearest proper orthogonal matrix, used only to remove roundoff."""
    u, _, vh = np.linalg.svd(a)
    q = u @ vh
    if np.linalg.det(q) < 0:
        u[:, -1] *= -1.0
        q = u @ vh
    return q


def log_so5(d):
    """Real antisymmetric logarithm of a nearby SO(5) matrix."""
    l = logm(d)
    l = np.real_if_close(l, tol=1000).real
    return 0.5 * (l - l.T)


def split_error(xi):
    """Return h=so(3)+so(2), m=mixing, and total Lie-algebra norms."""
    h = np.zeros_like(xi)
    # MZM block so(3): indices 0,1,2; QD block so(2): indices 3,4.
    for i, j in ((0, 1), (0, 2), (1, 2), (3, 4)):
        h[i, j], h[j, i] = xi[i, j], xi[j, i]
    m = xi - h
    norm = lambda z: np.linalg.norm(z, "fro") / np.sqrt(2.0)
    return norm(h), norm(m), norm(xi)


def solve_endpoint(tau, e1, t1):
    if tau == 0:
        return np.eye(5), np.eye(5)
    sol = solve_ivp(
        lambda q, y: rhs_pair(q, y, tau, e1, t1),
        (0.0, 6.0 * tau),
        pack_pair(np.eye(5), np.eye(5)),
        method="DOP853", rtol=2e-10, atol=2e-12,
    )
    r0, r = unpack_pair(sol.y[:, -1])
    return project_so5(r0), project_so5(r)


def endpoint_error(tau, e1, t1):
    r0, r = solve_endpoint(tau, e1, t1)
    d = project_so5(r0.T @ r)
    return split_error(log_so5(d))


def solve_path(tau, e1, t1, qgrid):
    if tau == 0:
        eye = np.repeat(np.eye(5)[None, :, :], len(qgrid), axis=0)
        return eye, eye
    sol = solve_ivp(
        lambda q, y: rhs_pair(q, y, tau, e1, t1),
        (0.0, float(qgrid[-1])),
        pack_pair(np.eye(5), np.eye(5)),
        t_eval=qgrid, method="DOP853", rtol=2e-10, atol=2e-12,
    )
    r0 = sol.y[:25].T.reshape(-1, 5, 5)
    r = sol.y[25:].T.reshape(-1, 5, 5)
    return np.array([project_so5(x) for x in r0]), np.array([project_so5(x) for x in r])


def make_maps():
    xs = np.linspace(0.0, XMAX, N_X)
    ys = np.linspace(YMIN, YMAX, N_Y)
    eh = np.zeros((N_Y, N_X))
    em = np.zeros_like(eh)
    et = np.zeros_like(eh)
    for ix, tau in enumerate(xs):
        for iy, y in enumerate(ys):
            h, m, total = endpoint_error(tau, E1, E1 * 10.0**y)
            eh[iy, ix], em[iy, ix], et[iy, ix] = h, m, total
        if ix % max(1, N_X // 10) == 0:
            print(f"map {ix + 1}/{N_X}: tau={tau:.3f}", flush=True)
    np.savez(
        "so5_error_decomposition.npz", x=xs, y=ys,
        eps_internal=eh, eps_mix=em, eps_total=et,
        E1=E1, tc=TC, E0=E0, normalization=SO5_NORMALIZATION,
    )
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
    for ax, z, title in zip(
        axes, (eh, em, et),
        (r"$\epsilon_h$: $SO(3)\oplus SO(2)$", r"$\epsilon_m$: MZM--QD mixing", r"$\epsilon_{\rm total}$"),
    ):
        im = ax.pcolormesh(xs, ys, z, shading="auto", cmap="magma")
        ax.set_xlabel(r"$\tau/(100\,\mathrm{meV}^{-1})$")
        ax.set_ylabel(r"$\log_{10}(t_1/E_1)$")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=.9)
    fig.savefig("so5_error_decomposition.png", dpi=220)
    print("wrote so5_error_decomposition.png and .npz")


def make_paths():
    tau = float(os.environ.get("SO5_ERR_TAU", "4.0"))
    qgrid = np.linspace(0.0, 6.0 * tau, N_PATH)
    ratios = (-1.0, 0.0, 1.0)
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.0), sharex=True, constrained_layout=True)
    for y in ratios:
        r0, r = solve_path(tau, E1, E1 * 10.0**y, qgrid)
        vals = np.zeros((len(qgrid), 3))
        for k, (a, b) in enumerate(zip(r0, r)):
            vals[k] = split_error(log_so5(project_so5(a.T @ b)))
        label = rf"$\log_{{10}}(t_1/E_1)={y:g}$"
        axes[0].plot(qgrid, vals[:, 0], label=label)
        axes[1].plot(qgrid, vals[:, 1], label=label)
        axes[2].plot(qgrid, vals[:, 2], label=label)
    axes[0].set_ylabel(r"$\epsilon_h(t)$")
    axes[1].set_ylabel(r"$\epsilon_m(t)$")
    axes[2].set_ylabel(r"$\epsilon_{\rm total}(t)$")
    axes[2].set_xlabel(r"$q=t/(100\,\mathrm{meV}^{-1})$")
    for ax in axes:
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    fig.suptitle(rf"Relative error path, $E_1={E1:g}$ meV, $\tau={tau:g}$")
    fig.savefig("so5_error_paths.png", dpi=220)
    print("wrote so5_error_paths.png")


if __name__ == "__main__":
    make_maps()
    make_paths()
