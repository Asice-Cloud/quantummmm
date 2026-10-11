#!/usr/bin/env python3
"""Exact modular check of output-derivative ranks for the cosine protocol.

The local coordinate is s=pi*v, tau=100*pi inverse meV, E1=.01,
t1=.03, tc=E0=.3. K_s(s)=K0+Kcos*cos(s). At s=pi/2 all jets
are rational. Reduction modulo a prime proves full rational rank when
the modular rank is full; a deficient modular rank alone is only a
lower bound, so stage 1 and 3 also require the structural arguments.
"""
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np

PRIME = 10007
PAIRS = list(combinations(range(5), 2))


def generator(i, j):
    result = np.zeros((5, 5), dtype=int)
    result[i, j], result[j, i] = 1, -1
    return result


def hermitian_basis():
    result = []
    for i in range(5):
        mat = np.zeros((5, 5), dtype=int)
        mat[i, i] = 1
        result.append((mat, False))
    for i, j in PAIRS:
        mat = np.zeros((5, 5), dtype=int)
        mat[i, j] = mat[j, i] = 1
        result.append((mat, False))
    for i, j in PAIRS:
        result.append((generator(i, j), True))
    return result


def adjoint_matrix(k):
    result = np.zeros((25, 25), dtype=np.int64)
    for j, (mat, imaginary) in enumerate(hermitian_basis()):
        derivative = k @ mat - mat @ k
        if imaginary:
            result[15:, j] = [derivative[i, l] for i, l in PAIRS]
        else:
            result[:5, j] = np.diag(derivative)
            result[5:15, j] = [derivative[i, l] for i, l in PAIRS]
    return result % PRIME


def modular_rank(mat):
    mat = mat.copy() % PRIME
    row = 0
    for col in range(mat.shape[1]):
        pivots = np.flatnonzero(mat[row:, col])
        if not len(pivots):
            continue
        pivot = row + pivots[0]
        mat[[row, pivot]] = mat[[pivot, row]]
        mat[row] = mat[row] * pow(int(mat[row, col]), -1, PRIME) % PRIME
        for j in range(row + 1, mat.shape[0]):
            mat[j] = (mat[j] - mat[j, col] * mat[row]) % PRIME
        row += 1
        if row == mat.shape[0]:
            break
    return row


def output_rows(k0, kcos, order=25):
    mjet = np.zeros((order + 1, 25, 25), dtype=np.int64)
    mjet[0] = adjoint_matrix(k0)
    for n in range(1, order + 1, 2):
        # Coefficients of cos(pi/2+delta)=-sin(delta).
        coefficient = ((-1) ** ((n + 1) // 2)
                       * pow(math.factorial(n) % PRIME, -1, PRIME)) % PRIME
        mjet[n] = coefficient * adjoint_matrix(kcos) % PRIME
    hjet = np.zeros((order + 1, 25), dtype=np.int64)
    hjet[0, 0] = hjet[0, 1] = pow(2, -1, PRIME)
    hjet[0, 15] = PRIME - 1
    rows = [hjet[0].copy()]
    for n in range(order):
        remaining = order - n
        next_jet = np.zeros((remaining, 25), dtype=np.int64)
        for k in range(remaining):
            value = (k + 1) * hjet[k + 1]
            for j in range(k + 1):
                value = (value + hjet[j] @ mjet[k - j]) % PRIME
            next_jet[k] = value % PRIME
        hjet = next_jet
        rows.append(hjet[0].copy())
    return np.array(rows)


def main():
    le, lt = generator(0, 1), generator(0, 4)
    l2a, l3a, lab = generator(1, 3), generator(2, 3), generator(3, 4)
    stages = [
        (2*le+3*lt-30*l2a+30*lab, -3*lt+30*l2a+30*lab),
        (2*le+3*lt-30*l2a+30*l3a, 3*lt-30*l2a-30*l3a),
        (2*le+30*l3a+30*lab, 30*l3a-30*lab),
    ]
    ranks = [modular_rank(output_rows(*stage)[:25]) for stage in stages]
    assert ranks == [16, 25, 1]
    report = {"tau_inverse_meV": 100*math.pi, "E1_meV": .01,
              "t1_meV": .03, "tc_meV": .3, "E0_meV": .3,
              "local_phase": "pi/2", "prime": PRIME,
              "stage_output_ranks_mod_prime": ranks,
              "interpretation": "Stage 2 full rank certifies a nonzero rational determinant; "
                                "stage 1 and 3 upper bounds follow from invariant blocks."}
    output = Path(__file__).with_name("so5_scalar_observability.json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
