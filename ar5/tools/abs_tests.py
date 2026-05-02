#!/usr/bin/env python3
"""Run two quick tests:
A) Inject small c_XY into Pauli mapping and check spectrum for ABS-like near‑zero localized states.
B) Add a single-site mu well (positive and negative) and check for localized subgap states.

Saves results to results/abs_tests.txt and prints summary.
"""
import os
import numpy as np
from numpy.linalg import eigh
from xxz_R_and_H import R_xxz, dRdu, expand_on_pauli
from verify_mzm import map_c_to_params


def build_bdg_from_arrays(N, t_hop, Delta_bond, mu_site):
    A = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)
    for i in range(N):
        A[i, i] = -mu_site[i]
        if i < N - 1:
            t = t_hop[i]
            A[i, i + 1] = -t
            A[i + 1, i] = -t
            Delta = Delta_bond[i]
            B[i, i + 1] = Delta
            B[i + 1, i] = -Delta
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    Hbdg = np.concatenate((top, bottom), axis=0)
    return Hbdg


def pauli_map_from_h(h_mapping):
    # extract required coefficients with defaults
    c_xx = h_mapping.get('X_X', 0.0)
    c_yy = h_mapping.get('Y_Y', 0.0)
    c_xy = h_mapping.get('X_Y', 0.0)
    c_yx = h_mapping.get('Y_X', 0.0)
    c_zz = h_mapping.get('Z_Z', 0.0)
    c_z0 = h_mapping.get('Z_I', 0.0)
    c_0z = h_mapping.get('I_Z', 0.0)
    return map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)


def run_tests(eta=0.6, u0=0.0, N=200, FORCE_CXY=0.05, MU_WELL=3.0):
    os.makedirs('results', exist_ok=True)
    out_file = 'results/abs_tests.txt'
    with open(out_file, 'w') as f:
        f.write('abs_tests results\n')
        f.write('eta=%s, u0=%s\n' % (eta, u0))

        R0 = R_xxz(u0, eta)
        dR0 = dRdu(u0, eta)

        # build P
        P = np.zeros((4, 4), dtype=complex)
        for a in range(2):
            for b in range(2):
                in_idx = (a << 1) | b
                out_idx = (b << 1) | a
                P[out_idx, in_idx] = 1.0

        rho = np.sin(eta)
        h_local = P @ dR0 / rho
        mapping_h = expand_on_pauli(h_local)

        # baseline parameters
        t_val, Delta_val, mu_val, U = pauli_map_from_h(mapping_h)
        f.write('\nBaseline params from h_local:\n')
        f.write(f't={t_val}, Delta={Delta_val}, mu={mu_val}, U={U}\n')

        # A: inject small c_xy
        mapping_mod = dict(mapping_h)
        mapping_mod['X_Y'] = mapping_mod.get('X_Y', 0.0) + FORCE_CXY
        t_mod, Delta_mod, mu_mod, U_mod = pauli_map_from_h(mapping_mod)
        f.write('\nTest A: injected c_XY=%g --> params:\n' % FORCE_CXY)
        f.write(f't_mod={t_mod}, Delta_mod={Delta_mod}, mu_mod={mu_mod}\n')

        # build chain and diagonalize
        t_hop = np.full(N - 1, t_mod, dtype=complex)
        Delta_bond = np.full(N - 1, Delta_mod, dtype=complex)
        mu_site = np.full(N, mu_mod, dtype=float)
        H = build_bdg_from_arrays(N, t_hop, Delta_bond, mu_site)
        eigs, vecs = eigh(H)
        eigs_sorted = np.sort(np.real_if_close(eigs))
        f.write('Lowest 10 eigenvalues (Test A):\n')
        for v in eigs_sorted[:10]:
            f.write(f'{v:.6e}\n')

        # localization of lowest mode
        idx = np.argsort(np.abs(eigs))[:1]
        v0 = vecs[:, idx[0]]
        # compute site weights (u and v parts)
        u = v0[:N]
        v = v0[N:]
        site_weights = np.abs(u)**2 + np.abs(v)**2
        max_idx = int(np.argmax(site_weights))
        f.write(f'Lowest mode max weight at site {max_idx}, weight={site_weights[max_idx]:.4e}\n')

        # B: single-site mu well (both signs)
        t_hop0 = np.full(N - 1, t_val, dtype=complex)
        Delta_bond0 = np.full(N - 1, Delta_val, dtype=complex)
        mu_site0 = np.full(N, mu_val, dtype=float)
        mid = N // 2
        for sign in (+1, -1):
            mu_site = mu_site0.copy()
            mu_site[mid] = mu_val + sign * MU_WELL
            H2 = build_bdg_from_arrays(N, t_hop0, Delta_bond0, mu_site)
            eigs2, vecs2 = eigh(H2)
            eigs2s = np.sort(np.real_if_close(eigs2))
            f.write(f'\nTest B: mu_mid = mu + {sign*MU_WELL} --> Lowest 10 eigenvalues:\n')
            for v in eigs2s[:10]:
                f.write(f'{v:.6e}\n')
            idxs = np.argsort(np.abs(eigs2))[:3]
            for ii in idxs:
                vec = vecs2[:, ii]
                u = vec[:N]; v = vec[N:]
                weights = np.abs(u)**2 + np.abs(v)**2
                f.write(f'mode idx {ii}: max site {int(np.argmax(weights))}, weight at mid {weights[mid]:.4e}\n')

    print('Tests complete. Results in', out_file)


if __name__ == '__main__':
    run_tests()
