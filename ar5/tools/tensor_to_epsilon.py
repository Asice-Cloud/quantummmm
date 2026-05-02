#!/usr/bin/env python3
"""Compute per-Pauli contributions to the lowest BdG mode energy for a candidate.

Workflow:
 - rebuild base Pauli mapping from XXZ local h (same convention as scans)
 - load top candidate from results/scan_all_mixtures_candidates.json (first entry)
 - construct full mapping (base + added combo amplitudes)
 - build BdG H_total and diagonalize to get lowest positive mode psi
 - for each Pauli key, build H_k with unit coefficient and compute F_k = psi^† H_k psi
 - report contributions c_k * F_k and check sum ~ E_low

This is a linearized (first-order) attribution of the observed small eigenvalue to
individual Pauli basis terms.
"""
import json
import numpy as np
from xxz_R_and_H import R_xxz, dRdu, expand_on_pauli
from verify_mzm import map_c_to_params, build_kitaev_bdg


def build_base_mapping(eta=0.6, u0=0.0):
    R0 = R_xxz(u0, eta)
    dR0 = dRdu(u0, eta)
    # permutation P
    P = np.zeros((4, 4), dtype=complex)
    for a in range(2):
        for b in range(2):
            in_idx = (a << 1) | b
            out_idx = (b << 1) | a
            P[out_idx, in_idx] = 1.0
    rho = np.sin(eta)
    h_local = P @ dR0 / rho
    return expand_on_pauli(h_local)


def pauli_keys():
    return ['X_X', 'Y_Y', 'Z_Z', 'X_Y', 'Y_X', 'X_Z', 'Y_Z', 'Z_X', 'Z_Y']


def mapping_to_bdgs(mapping, N=120):
    # reduce mapping to the canonical 7 coefficients used by map_c_to_params
    c_xx = mapping.get('X_X', 0.0)
    c_yy = mapping.get('Y_Y', 0.0)
    c_xy = mapping.get('X_Y', 0.0)
    c_yx = mapping.get('Y_X', 0.0)
    c_zz = mapping.get('Z_Z', 0.0)
    c_z0 = mapping.get('Z_I', 0.0)
    c_0z = mapping.get('I_Z', 0.0)
    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    H = build_kitaev_bdg(N, t, Delta, mu)
    return H, (t, Delta, mu, U)


def main():
    # parameters (match scan defaults)
    eta = 0.6
    u0 = 0.0
    N = 120

    base_mapping = build_base_mapping(eta=eta, u0=u0)

    with open('results/scan_all_mixtures_candidates.json', 'r') as f:
        candidates = json.load(f)
    if len(candidates) == 0:
        print('No candidates found in results/scan_all_mixtures_candidates.json')
        return

    cand = candidates[0]
    combo = cand.get('combo', [])
    v = float(cand.get('v', 0.0))

    mapping_mod = dict(base_mapping)
    for k in combo:
        mapping_mod[k] = mapping_mod.get(k, 0.0) + complex(v)

    H_total, params = mapping_to_bdgs(mapping_mod, N=N)

    # diagonalize
    evals, evecs = np.linalg.eigh(H_total)
    # pick smallest absolute eigenvalue index
    idx_sorted = np.argsort(np.abs(evals))
    idx0 = idx_sorted[0]
    E0 = evals[idx0]
    psi0 = evecs[:, idx0]

    # normalize
    psi0 = psi0 / np.linalg.norm(psi0)

    # Build pure-single-term H_k for each Pauli key (do NOT include base mapping)
    keys = pauli_keys() + ['Z_I', 'I_Z']
    H_k_map = {}
    for k in keys:
        mapping_k = {}
        mapping_k[k] = 1.0
        Hk, _ = mapping_to_bdgs(mapping_k, N=N)
        H_k_map[k] = Hk

    # Also build H_base for comparison
    H_base, _ = mapping_to_bdgs(base_mapping, N=N)

    # Compute F_k = psi^† H_k psi and contribution c_k * F_k
    contributions = {}
    for k in keys:
        Fk = np.vdot(psi0, H_k_map[k] @ psi0)
        ck = mapping_mod.get(k, 0.0)
        contributions[k] = {'c': complex(ck), 'F': complex(Fk), 'contrib': complex(ck) * complex(Fk)}

    # Summaries
    total_from_parts = sum([v['contrib'] for v in contributions.values()])

    print('Top-1 candidate combo =', combo, 'v =', v)
    print('Reconstructed params (t,Delta,mu,U):', params)
    print('Lowest eigenvalue E0 (closest to 0):', E0)
    print('\nPer-Pauli contributions (c_k * <psi|H_k|psi>):')
    for k, v in sorted(contributions.items(), key=lambda x: -abs(x[1]['contrib'])):
        print(f"{k}: c={v['c']:.6g}, F={v['F']:.6g}, contrib={v['contrib']:.6g}")

    print('\nSum of per-key contributions =', total_from_parts)
    exp_diff = np.vdot(psi0, (H_total - H_base) @ psi0)
    print('Expectation on (H_total - H_base) =', exp_diff)
    print('Direct expectation psi^† H_total psi =', np.vdot(psi0, H_total @ psi0))


if __name__ == '__main__':
    main()
