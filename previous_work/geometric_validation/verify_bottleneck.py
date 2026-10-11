#!/usr/bin/env python3
"""Verify the exact oriented-plane bottleneck with full five-mode evolution.

Propagate the complex wavefunction and the ten real bivector components
independently. No small-angle expansion or adiabatic projection is used.
Run from the project root: .venv/bin/python -B geometric_validation/verify_bottleneck.py
"""

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/quantumsss-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parent
INDICES = np.triu_indices(5, 1)
PSI_INITIAL = np.array([1, -1j, 0, 0, 0]) / np.sqrt(2)
PSI_TARGET = PSI_INITIAL.conj()
B_INITIAL = np.zeros((5, 5))
B_INITIAL[0, 1], B_INITIAL[1, 0] = 1, -1
CASES = [
    {"name": "reference", "tau": 100., "E1": 0., "eta": 0.},
    {"name": "eta_low", "tau": 100., "E1": .01, "eta": .01},
    {"name": "high_fidelity", "tau": 100., "E1": .01, "eta": .023},
    {"name": "E_changed", "tau": 100., "E1": .02, "eta": .023},
    {"name": "long_protocol", "tau": 500., "E1": .01,
     "eta": .02344645627},
]


def controls(stage, v):
    rise = (1 - np.cos(np.pi * v)) / 2
    fall = 1 - rise
    if stage % 3 == 0:
        return .3 * rise, 0., .3 * fall, rise
    if stage % 3 == 1:
        return .3 * fall, .3 * rise, 0., fall
    return 0., .3 * fall, .3 * rise, 0.


def wave_rhs(stage, v, state, case):
    t2, t3, ed, gate = controls(stage, v)
    e, eta = case["E1"], case["eta"] * gate
    # Direct complex-state component equations, in units hbar=1.
    return 2 * case["tau"] * np.array([
        e * state[1] + eta * state[4],
        -e * state[0] - t2 * state[3],
        t3 * state[3],
        t2 * state[1] - t3 * state[2] + ed * state[4],
        -eta * state[0] - ed * state[3],
    ])


def bivector_matrix(values):
    matrix = np.zeros(values.shape[:-1] + (5, 5))
    matrix[..., INDICES[0], INDICES[1]] = values
    return matrix - np.swapaxes(matrix, -1, -2)


def bivector_rhs(stage, v, values, case):
    t2, t3, ed, gate = controls(stage, v)
    generator = np.zeros((5, 5))
    for i, j, coefficient in [
        (0, 1, case["E1"]), (0, 4, case["eta"] * gate),
        (1, 3, -t2), (2, 3, t3), (3, 4, ed),
    ]:
        generator[i, j] = coefficient
        generator[j, i] = -coefficient
    generator *= 2 * case["tau"]
    b = bivector_matrix(values)
    return (generator @ b - b @ generator)[INDICES]


def orientation_event(v, state):
    return 2 * np.imag(state[0] * state[1].conjugate())


def propagate(case):
    psi, values = PSI_INITIAL.copy(), B_INITIAL[INDICES].copy()
    time_samples, wave_samples, bivector_samples, crossings = [], [], [], []
    for stage in range(6):
        grid = np.linspace(0, 1, 801)
        wave = solve_ivp(
            lambda v, state: wave_rhs(stage, v, state, case), (0, 1), psi,
            method="DOP853", rtol=8e-13, atol=8e-15, t_eval=grid,
            events=orientation_event, max_step=.01,
        )
        bivector = solve_ivp(
            lambda v, state: bivector_rhs(stage, v, state, case),
            (0, 1), values, method="DOP853", rtol=2e-12, atol=2e-14,
            t_eval=grid, max_step=.01,
        )
        if not wave.success or not bivector.success:
            raise RuntimeError((wave.message, bivector.message))
        keep = slice(None) if stage == 0 else slice(1, None)
        time_samples.append(stage + grid[keep])
        wave_samples.append(wave.y.T[keep])
        bivector_samples.append(bivector.y.T[keep])
        for local_time, state in zip(wave.t_events[0], wave.y_events[0]):
            crossings.append({
                "t_over_tau": float(stage + local_time),
                "ell": float(np.sum(abs(state[2:]) ** 2)),
                "F": float(abs(np.vdot(PSI_TARGET, state)) ** 2),
                "d": float(orientation_event(local_time, state)),
            })
        psi, values = wave.y[:, -1], bivector.y[:, -1]

    time = np.concatenate(time_samples)
    states = np.concatenate(wave_samples)
    b = bivector_matrix(np.concatenate(bivector_samples))
    p = -b @ b
    rho_from_b = (p + 1j * b) / 2
    rho_from_wave = states[..., :, None] * states.conj()[..., None, :]
    d = 2 * np.imag(states[:, 0] * states[:, 1].conjugate())
    ell = np.sum(abs(states[:, 2:]) ** 2, axis=1)
    f_wave = abs(states @ PSI_TARGET.conj()) ** 2
    f_b = (np.sum(b[:, :2] ** 2, axis=(1, 2)) - 2 * b[:, 0, 1]) / 4
    norm_error = float(np.max(abs(np.sum(abs(states) ** 2, axis=1) - 1)))
    null_error = float(np.max(abs(np.sum(states * states, axis=1))))
    report = {
        **case,
        "F_end": float(f_wave[-1]), "f_end": float(np.sqrt(f_wave[-1])),
        "d_end": float(d[-1]), "ell_end": float(ell[-1]),
        "max_sampled_ell": float(ell.max()),
        "crossings": crossings,
        "minimum_ell_at_detected_crossings": min(
            (item["ell"] for item in crossings), default=None
        ),
        "maximum_F_at_detected_crossings": max(
            (item["F"] for item in crossings), default=None
        ),
        "errors": {
            "wave_normalization": norm_error,
            "wave_null_constraint": null_error,
            "state_reconstruction": float(np.max(abs(rho_from_b-rho_from_wave))),
            "wave_bivector_fidelity": float(np.max(abs(f_wave-f_b))),
            "F_equals_one_minus_ell_minus_d_over_two": float(
                np.max(abs(f_wave-(1-ell-d)/2))
            ),
            "bivector_constraint": float(np.max(abs(b @ b @ b + b))),
            "lower_geometric_bound_violation": float(
                max(0., np.max((1-d*d)/2-ell))
            ),
            "upper_geometric_bound_violation": float(
                max(0., np.max(ell-(1-abs(d))))
            ),
            "F_le_quarter_for_d_nonnegative_violation": float(
                max(0., np.max(f_wave[d >= 0]-.25))
            ),
        },
    }
    tolerance = 2e-9
    if max(report["errors"].values()) > tolerance:
        raise AssertionError(report)
    for item in crossings:
        if item["ell"] < .5-tolerance or item["F"] > .25+tolerance:
            raise AssertionError(item)
    if report["F_end"] > .25 and not crossings:
        raise AssertionError("A high-fidelity path must encounter d=0.")
    return report, {"t_over_tau": time, "d": d, "ell": ell,
                    "F": f_wave, "B": b, "psi": states}


def plot_paths(reports, trajectories):
    plt.rcParams.update({"font.size": 10, "mathtext.fontset": "dejavusans"})
    fig, (geometry, history) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    d_axis = np.linspace(-1, 1, 801)
    ell_axis = np.linspace(0, 1, 601)
    d_grid, ell_grid = np.meshgrid(d_axis, ell_axis)
    lower = (1-d_grid*d_grid)/2
    upper = 1-abs(d_grid)
    field = np.ma.masked_where(
        (ell_grid < lower) | (ell_grid > upper), (1-ell_grid-d_grid)/2
    )
    image = geometry.pcolormesh(d_grid, ell_grid, field, cmap="Greys_r",
                               vmin=0, vmax=1, shading="auto", alpha=.75,
                               rasterized=True)
    geometry.plot(d_axis, (1-d_axis*d_axis)/2, color=".35", lw=1)
    geometry.plot(d_axis, 1-abs(d_axis), color=".35", lw=1)
    geometry.axvline(0, color="black", ls=":", lw=1)
    geometry.axhline(.5, color="black", ls=":", lw=1)
    geometry.text(.025, .79, r"$d=0$ requires $\ell\geq 1/2$", fontsize=9)
    colors = ["#999999", "#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    for report, trajectory, color in zip(reports, trajectories, colors):
        reference = report["name"] == "reference"
        label = (rf"$E={report['E1']:g},\ \eta={report['eta']:.5g},\ "
                 rf"\tau={report['tau']:g}$")
        geometry.plot(trajectory["d"], trajectory["ell"], color=color,
                      lw=.9 if reference else 1.25, alpha=.85)
        geometry.scatter(trajectory["d"][-1], trajectory["ell"][-1],
                         marker="s", s=35, color=color, edgecolor="white",
                         linewidth=.4, zorder=5)
        history.plot(trajectory["t_over_tau"], trajectory["F"], color=color,
                     lw=1.1 if reference else 1.5, label=label)
    geometry.scatter([1], [0], marker="o", color="black", s=35, clip_on=False)
    geometry.scatter([-1], [0], marker="*", color="black", s=75, clip_on=False)
    geometry.set(xlabel=r"Signed plane overlap $d=\mathcal{B}_{12}$",
                 ylabel=r"Weight outside $(1,2)$: $\ell$",
                 xlim=(-1.08, 1.08), ylim=(-.03, 1.04),
                 title="Exact geometric domain and full trajectories")
    geometry.text(.65, .045, "initial", fontsize=9)
    geometry.text(-1., .045, "target", fontsize=9)
    fig.colorbar(image, ax=geometry, fraction=.046, pad=.04,
                 label=r"Squared fidelity $F=(1-\ell-d)/2$")
    history.axhline(.25, color="black", ls=":", lw=1)
    for boundary in range(1, 6):
        history.axvline(boundary, color=".8", lw=.5)
    history.set(xlabel=r"Time $t/\tau$", ylabel=r"Squared fidelity $F$",
                xlim=(0, 6), ylim=(-.03, 1.04),
                title="Independent full-wavefunction evolution")
    history.legend(loc="upper left", fontsize=8, framealpha=.9)
    history.grid(alpha=.15)
    fig.tight_layout()
    fig.savefig(ROOT / "bottleneck.png", dpi=220)
    fig.savefig(ROOT / "bottleneck.pdf")
    plt.close(fig)


def main():
    reports, trajectories = [], []
    for case in CASES:
        report, trajectory = propagate(case)
        reports.append(report)
        trajectories.append(trajectory)
        print(f"{case['name']}: F={report['F_end']:.9g}, "
              f"d={report['d_end']:.9g}, ell={report['ell_end']:.9g}", flush=True)
    checks = {
        key: max(report["errors"][key] for report in reports)
        for key in reports[0]["errors"]
    }
    result = {
        "protocol": {"cycles": 2, "segments_per_cycle": 3,
                     "tc_meV": .3, "E0_meV": .3, "hbar": 1,
                     "tau_unit": "meV^-1", "energy_unit": "meV"},
        "scope": "Five specified full finite-parameter paths, not a parameter-space scan.",
        "max_errors": checks, "cases": reports,
    }
    (ROOT / "bottleneck.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    data = {f"{report['name']}_{key}": value
            for report, trajectory in zip(reports, trajectories)
            for key, value in trajectory.items()}
    np.savez_compressed(ROOT / "bottleneck.npz", **data)
    plot_paths(reports, trajectories)
    print("Maximum verification errors:", json.dumps(checks))


if __name__ == "__main__":
    main()
