#!/usr/bin/env python3
"""Plot Pauli-projected coefficients in the complex plane at key E1 values

Saves PNGs to quantity/pauli_phase_plots/ and a CSV of coeffs for inspected E1s.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import importlib.util
import sys

ROOT = os.path.dirname(__file__)

spec = importlib.util.spec_from_file_location('embed_kitaev', os.path.join(ROOT, '..', 'tools', 'embed_kitaev.py'))
embed_kitaev = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = embed_kitaev
spec.loader.exec_module(embed_kitaev)
build_bdg = embed_kitaev.build_bdg


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
            labels.append(f'{a}{b}')
    return basis, labels


def sigma_to_pauli_coeffs(Sigma_maj):
    basis, labels = pauli16_basis()
    coeffs = []
    for P in basis:
        denom = np.real(np.trace(P.conj().T @ P))
        c = np.trace(P.conj().T @ Sigma_maj) / (denom + 1e-30)
        coeffs.append(c)
    return np.array(coeffs), labels


def choose_key_e1s(csvpath, k=3):
    data = np.genfromtxt(csvpath, delimiter=',', names=True)
    E = data['E1']
    if 'N_pauli' in data.dtype.names:
        metric = data['N_pauli']
    else:
        metric = data['K_fro']
    # pick max, min, median-index
    idx_max = int(np.nanargmax(metric))
    idx_min = int(np.nanargmin(metric))
    idx_med = int(len(E)//2)
    chosen = sorted(set([idx_max, idx_min, idx_med]))
    return E[chosen]


def main():
    outdir = os.path.join(ROOT, 'pauli_phase_plots')
    os.makedirs(outdir, exist_ok=True)

    csvp = os.path.join(ROOT, 'pauli_schur_compare', 'pauli_schur_compare.csv')
    if not os.path.exists(csvp):
        raise SystemExit('Expected CSV not found: ' + csvp)

    E_keys = choose_key_e1s(csvp)

    # parameters matching compare_pauli_schur
    L = 3
    t0 = 1.0
    Delta0 = 0.3
    u = 0.5
    v = 0.7
    eta = 1e-3

    P_sites = [0,2]
    Q_sites = [1]
    P_idx = []
    for s in P_sites:
        P_idx.append(s)
    for s in P_sites:
        P_idx.append(s + L)
    Q_idx = []
    for s in Q_sites:
        Q_idx.append(s)
    for s in Q_sites:
        Q_idx.append(s + L)

    # BdG->Majorana map same as compare script
    M = np.array([
        [1,0,1,0],
        [-1j,0,1j,0],
        [0,1,0,1],
        [0,-1j,0,1j]
    ], dtype=complex)

    coeffs_out = []
    labels = None

    for E1 in E_keys:
        mu = E1 * np.ones(L)
        H = build_bdg(mu, t0 * np.ones(L-1), Delta0 * np.ones(L-1))
        omega = 0.0 + 1j * eta

        # compute Schur Sigma (reuse local small inv)
        H_PQ = H[np.ix_(P_idx, Q_idx)]
        H_QP = H[np.ix_(Q_idx, P_idx)]
        H_QQ = H[np.ix_(Q_idx, Q_idx)]
        inv = np.linalg.inv(omega * np.eye(len(Q_idx), dtype=complex) - H_QQ)
        Sigma = H_PQ @ inv @ H_QP

        Sigma_maj = M @ Sigma @ np.linalg.inv(M)
        coeffs, labels = sigma_to_pauli_coeffs(Sigma_maj)
        coeffs_out.append((float(E1), coeffs))

    # save coeffs CSV
    coeff_csv = os.path.join(outdir, 'pauli_coeffs_key_E1s.csv')
    with open(coeff_csv, 'w') as fh:
        fh.write('E1,' + ','.join(labels) + '\n')
        for E1, c in coeffs_out:
            row = [str(E1)] + [('%g+%gj' % (np.real(x), np.imag(x))) for x in c]
            fh.write(','.join(row) + '\n')

    # plotting
    for E1, c in coeffs_out:
        mags = np.abs(c)
        order = np.argsort(-mags)
        fig, ax = plt.subplots(figsize=(6,6))
        ax.axhline(0, color='gray', lw=0.5)
        ax.axvline(0, color='gray', lw=0.5)
        sc = ax.scatter(c.real, c.imag, c=mags, cmap='viridis', s=50)
        # annotate top 4 coefficients
        for idx in order[:4]:
            ax.annotate(labels[idx], (c.real[idx], c.imag[idx]), textcoords='offset points', xytext=(5,5))
        ax.set_xlabel('Real')
        ax.set_ylabel('Imag')
        ax.set_title('Pauli coeffs (complex) at E1=%.6f' % E1)
        plt.colorbar(sc, label='|coeff|')
        plt.tight_layout()
        png = os.path.join(outdir, 'pauli_coeffs_E1_%g.png' % E1)
        plt.savefig(png, dpi=200)
        plt.close()

    print('Saved coeff CSV:', coeff_csv)
    print('Saved plots in:', outdir)


if __name__ == '__main__':
    main()
