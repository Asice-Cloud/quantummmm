#!/usr/bin/env python3
"""Compare Pauli-projected H_eff with Frobenius-scaled H_eff from Schur Sigma

Produces CSV and plots comparing operator norms, scale factors, and N indicators.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import csv
import importlib.util
import sys


spec = importlib.util.spec_from_file_location('embed_kitaev', os.path.join(os.path.dirname(__file__), '..', 'tools', 'embed_kitaev.py'))
embed_kitaev = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = embed_kitaev
spec.loader.exec_module(embed_kitaev)
build_bdg = embed_kitaev.build_bdg


def indices_for_sites(sites, L):
    idx = []
    for s in sites:
        idx.append(s)
    for s in sites:
        idx.append(s + L)
    return idx


def schur_sigma_bdG(omega, H, P_idx, Q_idx):
    H_PQ = H[np.ix_(P_idx, Q_idx)]
    H_QP = H[np.ix_(Q_idx, P_idx)]
    H_QQ = H[np.ix_(Q_idx, Q_idx)]
    inv = np.linalg.inv(omega * np.eye(len(Q_idx), dtype=complex) - H_QQ)
    Sigma = H_PQ @ inv @ H_QP
    return Sigma


def embed_op_on_three(op, which):
    I2 = np.eye(2, dtype=complex)
    if which == '12':
        top = np.hstack([op, np.zeros((4,2), dtype=complex)])
        bot = np.hstack([np.zeros((2,4), dtype=complex), I2])
        M = np.vstack([top, bot])
    else:
        top = np.hstack([I2, np.zeros((2,4), dtype=complex)])
        bot = np.hstack([np.zeros((4,2), dtype=complex), op])
        M = np.vstack([top, bot])
    return M


def pauli16_basis():
    I = np.eye(2, dtype=complex)
    sx = np.array([[0,1],[1,0]], dtype=complex)
    sy = np.array([[0,-1j],[1j,0]], dtype=complex)
    sz = np.array([[1,0],[0,-1]], dtype=complex)
    mats = [I, sx, sy, sz]
    basis = []
    labels = []
    for a in range(4):
        for b in range(4):
            basis.append(np.kron(mats[a], mats[b]))
            labels.append((a,b))
    return basis, labels


def sigma_to_pauli_coeffs(Sigma_maj):
    basis, labels = pauli16_basis()
    coeffs = []
    dens = []
    for P in basis:
        denom = np.real(np.trace(P.conj().T @ P))
        c = np.trace(P.conj().T @ Sigma_maj) / (denom + 1e-30)
        coeffs.append(c)
        dens.append(denom)
    return np.array(coeffs), basis, labels


def main(outdir='quantity/pauli_schur_compare'):
    os.makedirs(outdir, exist_ok=True)

    L = 3
    t0 = 1.0
    Delta0 = 0.3
    E1_vals = np.linspace(0.001, 0.2, 20)
    u = 0.5
    v = 0.7
    eta = 1e-3

    P_sites = [0,2]
    Q_sites = [1]
    P_idx = indices_for_sites(P_sites, L)
    Q_idx = indices_for_sites(Q_sites, L)

    # BdG->Majorana map
    M = np.array([
        [1,0,1,0],
        [-1j,0,1j,0],
        [0,1,0,1],
        [0,-1j,0,1j]
    ], dtype=complex)

    sy = np.array([[0,-1j],[1j,0]], dtype=complex)
    sx = np.array([[0,1],[1,0]], dtype=complex)
    B_base = np.kron(sy, sx)

    rows = []
    N_pauli_list = []
    N_fro_scaled_list = []
    scale_list = []
    opdiff_list = []

    for E1 in E1_vals:
        mu = E1 * np.ones(L)
        H = build_bdg(mu, t0 * np.ones(L-1), Delta0 * np.ones(L-1))
        omega = 0.0 + 1j * eta
        Sigma = schur_sigma_bdG(omega, H, P_idx, Q_idx)

        # Frobenius scalar
        K_fro = np.linalg.norm(Sigma, ord='fro')

        # Transform to Majorana basis and project to Pauli basis
        Sigma_maj = M @ Sigma @ np.linalg.inv(M)
        coeffs, basis, labels = sigma_to_pauli_coeffs(Sigma_maj)

        # Construct H_eff_pauli (keep complex, then symmetrize to Hermitian)
        H_eff_pauli = np.zeros((4,4), dtype=complex)
        for c, P in zip(coeffs, basis):
            H_eff_pauli += c * P
        # make Hermitian (dynamics requires Hermitian generator)
        H_eff_pauli = 0.5 * (H_eff_pauli + H_eff_pauli.conj().T)

        # H_eff_fro (simple base operator)
        H_eff_fro = K_fro * B_base

        # compute operator norms and scale
        norm_pauli = np.linalg.norm(H_eff_pauli, ord='fro')
        norm_fro = np.linalg.norm(H_eff_fro, ord='fro')
        scale = norm_pauli / (norm_fro + 1e-30)
        H_eff_fro_scaled = H_eff_fro * scale

        opdiff = np.linalg.norm(H_eff_pauli - H_eff_fro_scaled, ord='fro')

        # exponentiate and embed to 6x6 for three-site sequence
        R12_p = scipy.linalg.expm(-1j * H_eff_pauli * u)
        R12_v_p = scipy.linalg.expm(-1j * H_eff_pauli * v)
        R12_uv_p = scipy.linalg.expm(-1j * H_eff_pauli * (u+v))
        R23_p = scipy.linalg.expm(-1j * H_eff_pauli * u)
        R23_v_p = scipy.linalg.expm(-1j * H_eff_pauli * v)
        R23_uv_p = scipy.linalg.expm(-1j * H_eff_pauli * (u+v))

        R12_u_e = embed_op_on_three(R12_p, '12')
        R12_v_e = embed_op_on_three(R12_v_p, '12')
        R12_uv_e = embed_op_on_three(R12_uv_p, '12')
        R23_u_e = embed_op_on_three(R23_p, '23')
        R23_v_e = embed_op_on_three(R23_v_p, '23')
        R23_uv_e = embed_op_on_three(R23_uv_p, '23')

        LHS_p = R12_u_e @ R23_uv_e @ R12_v_e
        RHS_p = R23_v_e @ R12_uv_e @ R23_u_e
        Delta_p = LHS_p - RHS_p
        N_pauli = np.linalg.norm(Delta_p, ord='fro')

        # For Frobenius-scaled H
        R12_f = scipy.linalg.expm(-1j * (H_eff_fro * scale) * u)
        R12_v_f = scipy.linalg.expm(-1j * (H_eff_fro * scale) * v)
        R12_uv_f = scipy.linalg.expm(-1j * (H_eff_fro * scale) * (u+v))
        R23_f = scipy.linalg.expm(-1j * (H_eff_fro * scale) * u)
        R23_v_f = scipy.linalg.expm(-1j * (H_eff_fro * scale) * v)
        R23_uv_f = scipy.linalg.expm(-1j * (H_eff_fro * scale) * (u+v))

        R12_u_fe = embed_op_on_three(R12_f, '12')
        R12_v_fe = embed_op_on_three(R12_v_f, '12')
        R12_uv_fe = embed_op_on_three(R12_uv_f, '12')
        R23_u_fe = embed_op_on_three(R23_f, '23')
        R23_v_fe = embed_op_on_three(R23_v_f, '23')
        R23_uv_fe = embed_op_on_three(R23_uv_f, '23')

        LHS_f = R12_u_fe @ R23_uv_fe @ R12_v_fe
        RHS_f = R23_v_fe @ R12_uv_fe @ R23_u_fe
        Delta_f = LHS_f - RHS_f
        N_fro_scaled = np.linalg.norm(Delta_f, ord='fro')

        rows.append((float(E1), float(K_fro), float(norm_pauli), float(norm_fro), float(scale), float(opdiff), float(N_pauli), float(N_fro_scaled)))
        N_pauli_list.append(N_pauli)
        N_fro_scaled_list.append(N_fro_scaled)
        scale_list.append(scale)
        opdiff_list.append(opdiff)

    csvp = os.path.join(outdir, 'pauli_schur_compare.csv')
    with open(csvp, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['E1','K_fro','norm_pauli','norm_fro','scale','opdiff','N_pauli','N_fro_scaled'])
        for r in rows:
            w.writerow(r)

    E = np.array([r[0] for r in rows])
    Np = np.array(N_pauli_list)
    Nf = np.array(N_fro_scaled_list)
    scales = np.array(scale_list)

    plt.figure(figsize=(8,4))
    plt.plot(E, Np, 'b-o', label='N_pauli')
    plt.plot(E, Nf, 'r--s', label='N_fro_scaled')
    plt.xlabel('E1')
    plt.ylabel('N')
    plt.legend()
    plt.title('N: Pauli-projected vs Frobenius-scaled')
    plt.tight_layout()
    p1 = os.path.join(outdir, 'N_pauli_vs_fro_scaled.png')
    plt.savefig(p1, dpi=200)
    plt.close()

    plt.figure(figsize=(6,4))
    plt.plot(E, scales, '-o')
    plt.xlabel('E1')
    plt.ylabel('scale (||H_pauli||_F / ||H_fro||_F)')
    plt.title('Operator scale factor')
    plt.tight_layout()
    p2 = os.path.join(outdir, 'pauli_fro_scale.png')
    plt.savefig(p2, dpi=200)
    plt.close()

    print('Saved CSV:', csvp)
    print('Saved plots:', p1, p2)
    return csvp, p1, p2


if __name__ == '__main__':
    main()
