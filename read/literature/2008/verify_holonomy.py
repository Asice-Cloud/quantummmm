"""Independent checks of arts/2008.pdf using 12 Jordan-Wigner Majoranas.

Run from the repository root:
    MPLCONFIGDIR=/tmp/quantumsss-mpl .venv/bin/python read/literature/2008/verify_holonomy.py
Only NumPy and Matplotlib are required. Units: hbar = Theta = 1.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)


def kron_all(parts):
    out = np.array([[1.0]], dtype=complex)
    for part in parts:
        out = np.kron(out, part)
    return out


# Pair order: A12, A34, A05, B12, B34, B05.
# gamma_first=X and gamma_second=-Y imply i gamma_first gamma_second=Z.
gammas = {}
for j, (qubit, pair) in enumerate(
    [(q, p) for q in ("A", "B") for p in ((1, 2), (3, 4), (0, 5))]
):
    for label, pauli in zip(pair, (X, -Y)):
        gammas[qubit, label] = kron_all([Z] * j + [pauli] + [I2] * (5 - j))


def parity(q, j, k):
    return 1j * gammas[q, j] @ gammas[q, k]


def inter(j, k):
    return 1j * gammas["A", j] @ gammas["B", k]


PA = parity("A", 1, 2) @ parity("A", 3, 4) @ parity("A", 0, 5)
PB = parity("B", 1, 2) @ parity("B", 3, 4) @ parity("B", 0, 5)
PX = parity("B", 0, 2)
PZ = parity("B", 0, 5)
I20 = inter(2, 0)
I11 = inter(1, 1)
O4 = parity("A", 1, 2) @ parity("B", 0, 1)


def indices(n, m, fixed_local_parity=False):
    result = []
    for i in range(64):
        bits = [(i >> (5 - j)) & 1 for j in range(6)]
        if bits[1] != n or bits[4] != m or bits[2] != 0:
            continue
        if sum(bits) % 2:
            continue
        if fixed_local_parity and (sum(bits[:3]) % 2 or sum(bits[3:]) % 2):
            continue
        result.append(i)
    return result


def restrict(matrix, ids):
    return matrix[np.ix_(ids, ids)]


def loop(nodes, per_segment=240, power=1.0):
    points = []
    for a, b in zip(nodes[:-1], nodes[1:]):
        for u in np.linspace(0, 1, per_segment, endpoint=False):
            v = ((1 - np.cos(np.pi * u)) / 2) ** power
            points.append(np.array(a) * (1 - v) + np.array(b) * v)
    return points


NODES = [(0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 1, 0),
         (1, 0, 0), (1, 0, 1), (0, 0, 1)]
ASYMMETRIC_NODES = [
    (0, 0, 0, 1), (0, 0, .8, 1), (0, 1.2, .8, 1),
    (0, 1.2, .8, 0), (1, 1.2, .8, 0), (1, 0, .8, 0),
    (1, 0, 0, 0), (1, 0, 0, 1), (0, 0, 0, 1),
]


def berry_phase(blocks, points):
    # The minus sign converts the overlap-product phase into integral i<u|du>.
    states = []
    min_gap = float("inf")
    for point in points:
        h = sum(p * b for p, b in zip((1, *point), blocks))
        energies, vectors = np.linalg.eigh(h)
        min_gap = min(min_gap, float(energies[1] - energies[0]))
        states.append(vectors[:, 0])
    product = 1 + 0j
    for a, b in zip(states, states[1:] + states[:1]):
        overlap = np.vdot(a, b)
        product *= overlap / abs(overlap)
    return -float(np.angle(product)), min_gap


def phases(protocol, energy=1.0, nodes=NODES, power=1.0):
    result, gaps = [], []
    points = loop(nodes, power=power)
    for n, m in ((0, 0), (0, 1), (1, 0), (1, 1)):
        ids = indices(n, m, protocol == "4gamma")
        if protocol == "4gamma":
            terms = [np.zeros((64, 64)), -PX, -O4, -PZ]
        elif protocol == "asymmetric":
            terms = [-energy * (PA + PB), -PX, -I11, -I20, -PZ]
        else:
            terms = [-energy * (PA + PB), -PX, -(I20 + I11), -PZ]
        phi, gap = berry_phase([restrict(t, ids) for t in terms], points)
        result.append(phi)
        gaps.append(gap)
    # Remove a common phase by comparing with the first logical state.
    relative = np.angle(np.exp(1j * (np.array(result) - result[0])))
    return {"phases_rad": result, "relative_phases_rad": relative.tolist(),
            "minimum_gap_within_accessible_blocks": min(gaps)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    np.testing.assert_allclose(I20 @ I11, O4, atol=1e-14)
    np.testing.assert_allclose((I20 + I11) @ (I20 + I11),
                               2 * np.eye(64) + 2 * O4, atol=1e-14)
    # Exact two-auxiliary-spin representation, valid at finite charging energy.
    # sigma_x=P25^B, sigma_y=-P02^B, sigma_z=P05^B;
    # tau_z=P^A and tau_x=I11. s=(-1)^(n+m).
    for n, m in ((0, 0), (0, 1), (1, 0), (1, 1)):
        ids = indices(n, m)
        sx = restrict(parity("B", 2, 5), ids)
        tx = restrict(I11, ids)
        tz = restrict(PA, ids)
        sy = restrict(-PX, ids)
        sz = restrict(PZ, ids)
        s = (-1) ** (n + m)
        np.testing.assert_allclose(restrict(O4, ids), s * sx, atol=1e-14)
        h_exact = restrict(-.83 * (PA + PB) - .37 * PX
                           - .91 * (I20 + I11) - 1.13 * PZ, ids)
        h_spins = -2 * .83 * tz + .37 * sy - 1.13 * sz
        h_spins -= .91 * tx @ (np.eye(4) + s * sx)
        np.testing.assert_allclose(h_exact, h_spins, atol=1e-14)
    target = np.diag([1, -1j, -1j, 1])
    S = np.diag([1, 1j])
    CZ = np.diag([1, 1, 1, -1])
    np.testing.assert_allclose(np.kron(S, S) @ target, CZ)
    report = {"algebra_exact_two_spin_representation_and_CZ_checks": "passed",
              "odd_odd_weight_at_y_axis_for_E_equals_y": float((1 - 1 / np.sqrt(2)) / 2),
              "4gamma": phases("4gamma"),
              "double_tunneling": {},
              "asymmetric_1.2_0.8": phases("asymmetric", nodes=ASYMMETRIC_NODES),
              "different_amplitudes": phases("double", nodes=[
                  (x * .7, y * 1.4, z * 1.1) for x, y, z in NODES], power=1.6)}
    for e in (.2, .5, 1.0, 2.0, 5.0):
        report["double_tunneling"][str(e)] = phases("double", e)
    expected = np.array([0, -np.pi / 2, -np.pi / 2, 0])
    cases = [report["4gamma"], *report["double_tunneling"].values(),
             report["asymmetric_1.2_0.8"], report["different_amplitudes"]]
    for case in cases:
        error = np.angle(np.exp(1j * (np.array(case["relative_phases_rad"]) - expected)))
        assert np.max(np.abs(error)) < 2e-4, error
    (OUT / "verification.json").write_text(json.dumps(report, indent=2) + "\n")

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4), constrained_layout=True)
    e = np.logspace(-2, 2, 500)
    gap = 2 / (np.sqrt(e * e + 1) + np.maximum(e, 1)) * np.minimum(e * e, 1)
    axes[0].loglog(e, gap)
    axes[0].axvline(1, color="gray", ls="--", lw=1)
    axes[0].set(xlabel=r"$E/\Theta_y$", ylabel=r"$\Delta E(\tau_2)/\Theta_y$",
                title="Paper Eq. (27): full-spectrum gap at the anchor")
    axes[0].grid(alpha=.2)
    energies = [float(k) for k in report["double_tunneling"]]
    for k, label in [(1, r"$\phi_{01}-\phi_{00}$"), (3, r"$\phi_{11}-\phi_{00}$")]:
        values = [report["double_tunneling"][str(v)]["relative_phases_rad"][k] / np.pi
                  for v in energies]
        axes[1].semilogx(energies, values, "o-", label=label)
    axes[1].set(xlabel=r"$E/\Theta$", ylabel=r"Relative geometric phase / $\pi$",
                title="Independent Wilson-loop check")
    axes[1].set_ylim(-.6, .12)
    axes[1].legend()
    axes[1].grid(alpha=.2)
    fig.savefig(OUT / "figures" / "gap_and_holonomy.png", dpi=180)
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
