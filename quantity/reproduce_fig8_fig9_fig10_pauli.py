#!/usr/bin/env python3
"""Generate Fig.8/9/10-like panels using our Pauli tensor kernel estimator model.

This script uses our own Pauli tensor formula and the third-order kernel
estimator from the documentation:
    h_{\\alpha\\beta}(u) -> \\mathcal N_{kernel} \\sim \\sqrt{\\sum |\\int\\int K(s,t) h_{ab}(s) h_{mn}(t) ds dt|^2}
The effective Hamiltonian is the eight-vertex local Pauli tensor model
from our documentation, not the paper's 6-Majorana braiding fidelity method.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm


I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def P(a, b):
    basis = {"I": I2, "x": sx, "y": sy, "z": sz}
    return np.kron(basis[a], basis[b])


def effective_hamiltonian_pauli(u, pars):
    """Our Pauli tensor model formula from the documentation."""
    delta = pars.get("delta", 0.0)
    return (
        np.cos(u) * P("x", "x")
        + 0.5 * np.sin(u) * (P("y", "x") - P("x", "y"))
        + 0.5 * delta * (P("z", "I") - P("I", "z"))
    )


def map_B_to_u(B, Bmax=4.5):
    return 2.0 * np.pi * (B / Bmax)


def path_ordered_R(u_final, pars, steps=220):
    us = np.linspace(0.0, u_final, steps)
    du = us[1] - us[0]
    R = np.eye(4, dtype=complex)
    for u in us:
        H = effective_hamiltonian_pauli(u, pars)
        R = expm(-1j * H * du) @ R
    return R


def embed12(R):
    return np.kron(R, I2)


def embed23(R):
    return np.kron(I2, R)


def yb_deviation(Ru, Rv, Ruv):
    R12u = embed12(Ru)
    R23u = embed23(Ru)
    R12v = embed12(Rv)
    R23v = embed23(Rv)
    R12uv = embed12(Ruv)
    R23uv = embed23(Ruv)
    return R12u @ R23uv @ R12v - R23v @ R12uv @ R23u


def nonabelian_measure(Delta):
    return np.sqrt(np.real(np.trace(Delta.conj().T @ Delta)))


def pauli_spectrum_vs_B(pars, B_grid, return_vecs=False):
    vals = []
    vecs = []
    for B in B_grid:
        u = map_B_to_u(B)
        H = effective_hamiltonian_pauli(u, pars)
        ev, evec = np.linalg.eigh(H)
        idx = np.argsort(np.real(ev))
        ev = np.real(ev[idx])
        evec = evec[:, idx]
        vals.append(ev)
        vecs.append(evec)
    vals = np.array(vals)
    if return_vecs:
        return vals, np.array(vecs)
    return vals


def conductance_map_from_spectrum(Evals, Evecs, E_grid, gamma=0.01, lead_indices=(0,)):
    nB = Evals.shape[0]
    nE = len(E_grid)
    G = np.zeros((nE, nB), dtype=float)
    for j in range(nB):
        Ej = Evals[j]
        Vj = Evecs[j]
        weights = np.zeros(Ej.shape[0], dtype=float)
        for idx in lead_indices:
            weights += np.abs(Vj[idx, :]) ** 2
        weights /= np.maximum(np.sum(weights), 1e-12)
        for n, En in enumerate(Ej):
            G[:, j] += weights[n] * gamma * gamma / ((E_grid - En) ** 2 + gamma * gamma)
    G = G / np.max(G)
    return G

PAULI_KEYS = ["I", "x", "y", "z"]
PAULI_INDEX = {"I": 0, "x": 1, "y": 2, "z": 3}


def levi_civita(i, j, k):
    if i == 0 or j == 0 or k == 0:
        return 0
    if (i, j, k) in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]:
        return 1
    if (i, j, k) in [(1, 3, 2), (3, 2, 1), (2, 1, 3)]:
        return -1
    return 0


def h_coeffs(u, pars):
    delta = pars.get("delta", 0.0)
    return {
        ("x", "x"): np.cos(u),
        ("y", "x"): 0.5 * np.sin(u),
        ("x", "y"): -0.5 * np.sin(u),
        ("z", "I"): 0.5 * delta,
        ("I", "z"): -0.5 * delta,
    }


def kernel_K(s, t, u, v):
    theta_st = 1.0 if s > t else 0.0
    theta_ts = 1.0 if t > s else 0.0
    return v * theta_st + min(t, v) - u * theta_ts - max(0.0, u - t)


def N_kernel_formula(pars, u=1.0, v=0.7, ns=120, nt=160):
    s_grid = np.linspace(0.0, u, ns)
    t_grid = np.linspace(0.0, u + v, nt)
    ds = s_grid[1] - s_grid[0]
    dt = t_grid[1] - t_grid[0]

    sumsq = 0.0
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

    val2 = 4.0 * (2.0 ** 3) * sumsq
    return np.sqrt(max(val2, 0.0))


def compute_N(pars, u=1.0, v=0.7, steps=220):
    n_kernel = N_kernel_formula(pars, u=u, v=v, ns=120, nt=160)
    return n_kernel, None


def plot_fig8_like(out_path):
    B = np.linspace(0.0, 4.5, 160)
    E = np.linspace(-0.3, 0.3, 220)

    pars_uniform = {"delta": 0.10}
    pars_inhom = {"delta": 0.04}

    spec_u, vecs_u = pauli_spectrum_vs_B(pars_uniform, B, return_vecs=True)
    G_u = conductance_map_from_spectrum(spec_u, vecs_u, E, gamma=0.012, lead_indices=(0,))
    spec_i, vecs_i = pauli_spectrum_vs_B(pars_inhom, B, return_vecs=True)
    G_i = conductance_map_from_spectrum(spec_i, vecs_i, E, gamma=0.009, lead_indices=(0,))

    u_scan = np.linspace(0.0, 6.0, 80)
    N_uniform = np.array([compute_N(pars_uniform, u=u_val, v=0.7, steps=140)[0] for u_val in u_scan])
    N_inhom = np.array([compute_N(pars_inhom, u=u_val, v=0.7, steps=140)[0] for u_val in u_scan])

    fig, axs = plt.subplots(2, 3, figsize=(12, 7))

    ax = axs[0, 0]
    for n in range(spec_u.shape[1]):
        ax.plot(B, spec_u[:, n], color="blue", lw=1)
    ax.set_xlim(0, 4.5)
    ax.set_ylim(-0.3, 0.3)
    ax.set_title("(a) Uniform-like Spectrum")
    ax.set_xlabel("B (T)")
    ax.set_ylabel("E (meV)")

    ax = axs[0, 1]
    im1 = ax.imshow(
        G_u,
        origin="lower",
        aspect="auto",
        extent=[B[0], B[-1], E[0], E[-1]],
        cmap="turbo",
        vmin=0,
        vmax=1,
    )
    ax.set_title("(b) Uniform-like zero-bias conductance map")
    ax.set_xlabel("B (T)")
    ax.set_ylabel("E (meV)")

    ax = axs[0, 2]
    ax.plot(u_scan, N_uniform, "b", lw=2)
    ax.set_ylim(0, np.max(N_uniform) * 1.1)
    ax.set_xlim(u_scan[0], u_scan[-1])
    ax.set_title(r"(c) Kernel $\mathcal{N}$ (uniform-like)")
    ax.set_xlabel(r"$u$")
    ax.set_ylabel(r"$\mathcal{N}$")

    ax = axs[1, 0]
    for n in range(spec_i.shape[1]):
        ax.plot(B, spec_i[:, n], color="blue", lw=1)
    ax.set_xlim(0, 4.5)
    ax.set_ylim(-0.3, 0.3)
    ax.set_title("(d) Inhomogeneous-like Spectrum")
    ax.set_xlabel("B (T)")
    ax.set_ylabel("E (meV)")

    ax = axs[1, 1]
    im2 = ax.imshow(
        G_i,
        origin="lower",
        aspect="auto",
        extent=[B[0], B[-1], E[0], E[-1]],
        cmap="turbo",
        vmin=0,
        vmax=1,
    )
    ax.set_title("(e) Inhomogeneous-like zero-bias conductance map")
    ax.set_xlabel("B (T)")
    ax.set_ylabel("E (meV)")

    ax = axs[1, 2]
    ax.plot(u_scan, N_inhom, "b", lw=2)
    ax.set_ylim(0, np.max(N_inhom) * 1.1)
    ax.set_xlim(u_scan[0], u_scan[-1])
    ax.set_title(r"(f) Kernel $\mathcal{N}$ (inhomogeneous-like)")
    ax.set_xlabel(r"$u$")
    ax.set_ylabel(r"$\mathcal{N}$")

    cbar = fig.colorbar(im2, ax=axs[:, 1], fraction=0.025, pad=0.02)
    cbar.set_label("G (2e^2/h), normalized")

    fig.suptitle("Fig.8-like Reproduction (Pauli Tensor Kernel N Model)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def plot_fig9_like(out_path):
    B = np.linspace(0.0, 4.5, 140)
    E = np.linspace(-0.3, 0.3, 220)

    scenarios = [
        {"name": "(a) Uniform", "delta": 0.05},
        {"name": "(b) Smooth distortion", "delta": 0.10},
        {"name": "(c) QD-like", "delta": 0.14},
        {"name": "(d) Disorder-like", "delta": 0.18},
        {"name": "(e) Strong distortion", "delta": 0.22},
        {"name": "(f) Bulk disorder", "delta": 0.28},
    ]

    fig, axs = plt.subplots(2, 3, figsize=(12, 7))
    last_im = None

    for k, sc in enumerate(scenarios):
        pars = {"delta": sc["delta"]}
        spec, vecs = pauli_spectrum_vs_B(pars, B, return_vecs=True)
        G = conductance_map_from_spectrum(spec, vecs, E, gamma=0.01, lead_indices=(0,))
        Nval, _ = compute_N(pars, u=1.0, v=0.7, steps=160)

        i, j = divmod(k, 3)
        ax = axs[i, j]
        last_im = ax.imshow(
            G,
            origin="lower",
            aspect="auto",
            extent=[B[0], B[-1], E[0], E[-1]],
            cmap="turbo",
            vmin=0,
            vmax=1,
        )
        ax.set_title(f"{sc['name']} zero-bias conductance map (N={Nval:.3f})")
        ax.set_xlabel("B (T)")
        ax.set_ylabel("E (meV)")

    cbar = fig.colorbar(last_im, ax=axs.ravel().tolist(), fraction=0.02, pad=0.02)
    cbar.set_label("G_{zb} (arb. units), normalized")

    fig.suptitle("Fig.9-like Reproduction (Pauli Tensor Kernel N Model)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def plot_fig10_like(out_path):
    B = np.linspace(0.0, 4.5, 150)
    E = np.linspace(-0.3, 0.3, 220)

    rng = np.random.default_rng(7)
    x = np.linspace(0, 1.0, 120)
    Vx = 0.3 + 1.3 * rng.random(len(x))

    # use B values reported in the paper Fig.10 (approx.) to match panels
    B_markers = [1.84, 1.93, 2.02]
    delta_values = [0.08, 0.14, 0.22]

    pars_d = {"delta": 0.22}
    spec_d, vecs_d = pauli_spectrum_vs_B(pars_d, B, return_vecs=True)
    G_d = conductance_map_from_spectrum(spec_d, vecs_d, E, gamma=0.012, lead_indices=(0,))
    N_d, _ = compute_N(pars_d, u=1.0, v=0.7, steps=160)

    taus = np.linspace(0.2, 6.0, 60)
    curves = []
    for delta_val in delta_values:
        pars_path = {"delta": delta_val}
        curves.append(np.array([compute_N(pars_path, u=u_val, v=0.7, steps=120)[0] for u_val in taus]))

    delta_grid = np.linspace(0.02, 0.30, 120)
    spectrum_vs_delta = []
    for Bval in B_markers:
        u_fixed = map_B_to_u(Bval)
        spec_vs_delta = []
        for delta_val in delta_grid:
            pars_tmp = {"delta": delta_val}
            Htmp = effective_hamiltonian_pauli(u_fixed, pars_tmp)
            spec_vs_delta.append(np.sort(np.real(np.linalg.eigvalsh(Htmp))))
        spectrum_vs_delta.append(np.array(spec_vs_delta))

    fig, axs = plt.subplots(3, 3, figsize=(15, 12))

    ax = axs[0, 0]
    ax.plot(x, Vx, color="red", lw=1.2)
    ax.hlines(0.25, 0, 1.0, color="purple", lw=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.8)
    ax.set_title("(a) Disorder configuration")
    ax.set_xlabel("L (µm)")
    ax.set_ylabel("E (meV)")

    ax = axs[0, 1]
    for n in range(spec_d.shape[1]):
        ax.plot(B, spec_d[:, n], color="blue", lw=1)
    for Bval, color in zip(B_markers, ["orange", "green", "magenta"]):
        ax.axvline(Bval, color=color, ls="--", lw=1)
    ax.set_xlim(0, 4.5)
    ax.set_ylim(-0.3, 0.3)
    ax.set_title("(b) Spectrum")
    ax.set_xlabel("B (T)")
    ax.set_ylabel("E (meV)")
    ax.legend([f"B={Bval:.2f}" for Bval in B_markers], loc="upper right", frameon=False)

    ax = axs[0, 2]
    im = ax.imshow(
        G_d,
        origin="lower",
        aspect="auto",
        extent=[B[0], B[-1], E[0], E[-1]],
        cmap="turbo",
        vmin=0,
        vmax=1,
    )
    ax.set_title(f"(c) Conductance map (N={N_d:.3f})")
    ax.set_xlabel("B (T)")
    ax.set_ylabel("E (meV)")

    for k in range(3):
        ax = axs[1, k]
        ax.plot(taus, curves[k], "b", lw=2, label=rf"$\delta={delta_values[k]:.2f}$")
        ax.set_xlim(taus[0], taus[-1])
        ax.set_ylim(0, np.max(curves[k]) * 1.1)
        ax.set_title(["(d)", "(e)", "(f)"][k] + r" Braiding $\mathcal{N}(\tau)$")
        ax.set_xlabel(r"$\tau$")
        ax.set_ylabel(r"$\mathcal{N}$")
        ax.legend(frameon=False, loc="upper left")

    for k in range(3):
        ax = axs[2, k]
        spec_vs_delta = spectrum_vs_delta[k]
        for n in range(spec_vs_delta.shape[1]):
            ax.plot(delta_grid, spec_vs_delta[:, n], color="blue", lw=0.8)
        ax.set_xlim(delta_grid[0], delta_grid[-1])
        maxabs = float(np.max(np.abs(spec_vs_delta)))
        ylim = max(0.3, maxabs * 1.05)
        ax.set_ylim(-ylim, ylim)
        ax.set_title(f"{['(g)', '(h)', '(i)'][k]} Spectrum vs $\\delta$, B={B_markers[k]:.2f}")
        ax.set_xlabel(r"$\delta$")
        ax.set_ylabel("E (meV)")

    cbar = fig.colorbar(im, ax=axs[0, 2], fraction=0.046, pad=0.04)
    cbar.set_label("G (2e^2/h), normalized")

    fig.suptitle("Fig.10-like Reproduction (Pauli Tensor Kernel N Model)")
    fig.tight_layout()
    fig.subplots_adjust(top=0.93, hspace=0.35, wspace=0.35)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main():
    os.makedirs("quantity", exist_ok=True)
    out8 = "quantity/fig8_like_pauli_repro.png"
    out9 = "quantity/fig9_like_pauli_repro.png"
    out10 = "quantity/fig10_like_pauli_repro.png"

    print("Generating Fig8-like panel...")
    plot_fig8_like(out8)
    print("Saved", out8)

    print("Generating Fig9-like panel...")
    plot_fig9_like(out9)
    print("Saved", out9)

    print("Generating Fig10-like panel...")
    plot_fig10_like(out10)
    print("Saved", out10)


if __name__ == "__main__":
    main()
