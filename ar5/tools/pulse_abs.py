#!/usr/bin/env python3
"""Construct ABS from R(u)-derived local Pauli h and pulse sequences.
- derive hermitian local h = P dR/du|0 / sin(eta)
- expand into Pauli basis to get c_{mu,nu}
- map c to (t,Delta,mu) per bond/site using verify_mzm.map_c_to_params
- build spatial BdG and compute spectrum
- support static SN interface and simple Floquet average of two pulses
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eigh
from tools.xxz_R_and_H import R_xxz, dRdu, expand_on_pauli
from tools.verify_mzm import map_c_to_params, build_kitaev_bdg

# Note: this script assumes tools/ is on PYTHONPATH when run from project root.


def h_from_xxz(eta):
    R0 = R_xxz(0.0, eta)
    dR0 = dRdu(0.0, eta)
    # permutation P
    P = np.zeros((4,4), dtype=complex)
    for a in range(2):
        for b in range(2):
            in_idx = (a<<1) | b
            out_idx = (b<<1) | a
            P[out_idx, in_idx] = 1.0
    rho = np.sin(eta)
    h_local = P @ dR0 / rho
    mapping = expand_on_pauli(h_local)
    return mapping


def pauli_map_to_params(mapping):
    # extract required coefficients with default 0
    c_xx = mapping.get('X_X', 0.0)
    c_yy = mapping.get('Y_Y', 0.0)
    c_xy = mapping.get('X_Y', 0.0)
    c_yx = mapping.get('Y_X', 0.0)
    c_zz = mapping.get('Z_Z', 0.0)
    c_z0 = mapping.get('Z_I', 0.0)
    c_0z = mapping.get('I_Z', 0.0)
    t, Delta, mu, U = map_c_to_params(c_xx, c_yy, c_xy, c_yx, c_zz, c_z0, c_0z)
    # ensure real scalars for mu and U (they should be real by construction)
    mu_r = float(np.real_if_close(mu))
    U_r = float(np.real_if_close(U))
    return complex(t), complex(Delta), mu_r, U_r


def build_chain_from_params(N, t_left, Delta_left, mu_left, t_right, Delta_right, mu_right, interface_width=2):
    # build arrays for t_hop (N-1), Delta_bond (N-1), mu_site (N)
    t_hop = np.zeros(N - 1, dtype=complex)
    Delta_bond = np.zeros(N - 1, dtype=complex)
    mu_site = np.zeros(N, dtype=float)
    mid = N // 2
    left_end = mid - interface_width//2
    right_start = left_end + interface_width
    # left region bonds
    for i in range(0, left_end):
        t_hop[i] = t_left
        Delta_bond[i] = Delta_left
    # interface bonds (linear interpolation)
    for i in range(left_end, right_start):
        alpha = (i - left_end) / max(1, (right_start - left_end - 1))
        t_hop[i] = (1 - alpha) * t_left + alpha * t_right
        Delta_bond[i] = (1 - alpha) * Delta_left + alpha * Delta_right
    # right region
    for i in range(right_start, N - 1):
        t_hop[i] = t_right
        Delta_bond[i] = Delta_right
    # mu sites: choose left/right
    for i in range(0, left_end+1):
        mu_site[i] = mu_left
    for i in range(left_end+1, N):
        mu_site[i] = mu_right
    # build BdG from arrays (site-wise mu_site, bond-wise t_hop and Delta_bond)
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

    H = build_bdg_from_arrays(N, t_hop, Delta_bond, mu_site)
    return H, t_hop, Delta_bond, mu_site


def plot_spectrum(vals, out_path, title='spectrum'):
    plt.figure(figsize=(6,4))
    plt.plot(np.sort(vals), '.', markersize=3)
    plt.title(title)
    plt.xlabel('level')
    plt.ylabel('Energy')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def run_demo():
    os.makedirs('results', exist_ok=True)
    eta = 0.6
    mapping = h_from_xxz(eta)
    t_val, Delta_val, mu_val, U = pauli_map_to_params(mapping)
    print('Derived params from R_xxz h_local:')
    print('t=', t_val, 'Delta=', Delta_val, 'mu=', mu_val, 'U=', U)

    # Optionally force a nonzero Delta on the left to test ABS formation
    FORCE_DELTA_LEFT = 0.3

    # Create a right-side normal region by zeroing pairing
    t_right = t_val
    Delta_right = 0.0
    mu_right = mu_val
    t_left = t_val
    # use forced Delta if provided, otherwise mapping value
    Delta_left = FORCE_DELTA_LEFT if FORCE_DELTA_LEFT is not None else Delta_val
    mu_left = mu_val

    N = 200
    H_sn, t_hop, Delta_bond, mu_site = None, None, None, None
    # build chain
    H, t_hop, Delta_bond, mu_site = build_chain_from_params(N, t_left, Delta_left, mu_left, t_right, Delta_right, mu_right, interface_width=4)
    vals, vecs = eigh(H)
    plot_spectrum(vals, 'results/pulse_sn_spectrum.png', title='Pulse-derived S-N spectrum')
    print('Saved results/pulse_sn_spectrum.png')

    # Floquet average demo: average h_left and h_right (first-order Magnus)
    # construct H_eff by averaging t and Delta
    t_eff = 0.5 * (t_left + t_right)
    Delta_eff = 0.5 * (Delta_left + Delta_right)
    H_eff, _, _, _ = build_chain_from_params(N, t_eff, Delta_eff, mu_left, t_eff, Delta_eff, mu_right, interface_width=4)
    vals_eff, vecs_eff = eigh(H_eff)
    plot_spectrum(vals_eff, 'results/pulse_floquet_spectrum.png', title='Floquet-avg spectrum')
    print('Saved results/pulse_floquet_spectrum.png')

    # Save a small report
    with open('results/pulse_report.txt','w') as f:
        f.write(f't_left={t_left}\nDelta_left={Delta_left}\nmu_left={mu_left}\n')
        f.write(f't_right={t_right}\nDelta_right={Delta_right}\nmu_right={mu_right}\n')
        f.write('\nLowest energies (SN):\n')
        for v in np.sort(np.real_if_close(vals))[:10]:
            f.write(f'{v:.6e}\n')
        f.write('\nLowest energies (Floquet avg):\n')
        for v in np.sort(np.real_if_close(vals_eff))[:10]:
            f.write(f'{v:.6e}\n')
    print('Report written to results/pulse_report.txt')

if __name__ == '__main__':
    run_demo()
