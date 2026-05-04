"""Run a minimal example: H(u)=a(u)XX + b(u)YY + c(u)XY + d(u)YX.

This script:
- computes two-site eigenvalues/eigenvectors along a path in u
- computes discrete Berry connection A(u)
- maps Pauli coefficients to (t, Delta, mu)
- builds BdG chain and inspects low-energy spectrum and localization

Run as: python3 tools/run_minimal_example.py
"""
import numpy as np
from scipy.linalg import eigh
import sys
from math import isclose

import os
sys.path.insert(0, os.path.dirname(__file__))
import minimal_model as mm


def linear_path(u: float):
    # simple interpolating path: turn on a (XX) term from 0->1
    a = u
    b = 0.0
    c = 0.0
    d = 0.0
    return a, b, c, d


def main():
    # discretize u
    N = 101
    us = np.linspace(0.0, 1.0, N)
    Hs = []
    h_list = []
    mapped = []
    for u in us:
        a, b, c, d = linear_path(u)
        H = mm.H_two_site(a, b, c, d)
        Hs.append(H)
        h = mm.expand_H_on_pauli(H)
        h_list.append(h)
        mapped.append(mm.map_h_to_kitaev(h))

    du = us[1] - us[0]
    As = mm.compute_berry_along_path(Hs, du)

    # measure non-Abelian tendency: average magnitude of off-diagonal A entries
    offdiag_norms = [np.linalg.norm(A - np.diag(np.diag(A))) for A in As]
    avg_offdiag = float(np.mean(offdiag_norms))

    print(f"Avg off-diagonal norm of A along path: {avg_offdiag:.3e}")

    # show one sample mapping at u=1
    params_final = mapped[-1]
    print("Final mapped Kitaev params (u=1):")
    print({k: complex(v) for k, v in params_final.items()})

    # build BdG chain for final params and check zero modes
    L = 60
    t = params_final["t"]
    Delta = params_final["Delta"]
    mu = params_final["mu"]
    H_bdg = mm.build_bdg_chain(L, t, Delta, mu, open_boundary=True)
    E, V = mm.bdg_spectrum(H_bdg)

    # sort by absolute energy
    idx = np.argsort(np.abs(E))
    E_sorted = E[idx]
    print("Lowest 10 |E| (BdG) :", [float(np.abs(x)) for x in E_sorted[:10]])

    zeros = mm.find_zero_modes(E, V, L, tol=1e-6)
    if zeros:
        print(f"Detected {len(zeros)} near-zero BdG modes. Example localization fractions:")
        for m in zeros[:3]:
            print(f"  E={m['E']:.3e}, left_frac={m['left_frac']:.3f}, right_frac={m['right_frac']:.3f}")
    else:
        print("No near-zero BdG modes detected (within tol=1e-6).")

    # quick topological check using (mu, t) criterion (approx): |mu| < 2|t|
    # note: this uses simple uniform-chain criterion assuming our mapping and conventions
    tval = np.abs(t)
    muval = np.abs(mu)
    topo = muval < 2 * tval
    print(f"Topological criterion |mu| < 2|t| ? {topo}  (|mu|={muval:.3e}, |t|={tval:.3e})")


if __name__ == "__main__":
    main()
