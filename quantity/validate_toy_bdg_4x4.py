import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import importlib.util
import sys


spec = importlib.util.spec_from_file_location(
    'embed_kitaev', os.path.join(os.path.dirname(__file__), '..', 'tools', 'embed_kitaev.py')
)
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
    H_PP = H[np.ix_(P_idx, P_idx)]
    H_PQ = H[np.ix_(P_idx, Q_idx)]
    H_QP = H[np.ix_(Q_idx, P_idx)]
    H_QQ = H[np.ix_(Q_idx, Q_idx)]
    inv = np.linalg.inv(omega * np.eye(len(Q_idx), dtype=complex) - H_QQ)
    Sigma = H_PQ @ inv @ H_QP
    return H_PP, Sigma


def sigma_majorana_projection(Sigma):
    M = np.array([
        [1, 0, 1, 0],
        [-1j, 0, 1j, 0],
        [0, 1, 0, 1],
        [0, -1j, 0, 1j],
    ], dtype=complex)
    Sigma_maj = M @ Sigma @ np.linalg.inv(M)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    B = np.kron(sy, sx)
    denom = np.real(np.trace(B.conj().T @ B))
    coeff = np.trace(B.conj().T @ Sigma_maj)
    return float(np.abs(coeff) / (denom + 1e-30))


def run_toy_validation(outdir='quantity/toy_bdg_validation'):
    os.makedirs(outdir, exist_ok=True)

    L = 3
    t0 = 1.0
    Delta0 = 0.3
    E1_vals = np.linspace(0.001, 0.2, 20)
    u = 0.5
    v = 0.7
    eta = 1e-3

    P_sites = [0, 2]
    Q_sites = [1]
    P_idx = indices_for_sites(P_sites, L)
    Q_idx = indices_for_sites(Q_sites, L)

    res = []
    for E1 in E1_vals:
        mu = E1 * np.ones(L)
        H = build_bdg(mu, t0 * np.ones(L - 1), Delta0 * np.ones(L - 1))
        omega = 0.0 + 1j * eta
        _, Sigma = schur_sigma_bdG(omega, H, P_idx, Q_idx)

        K_fro = np.linalg.norm(Sigma, ord='fro')
        K_pauli_proj = sigma_majorana_projection(Sigma)

        # same effective braiding indicator construction as the Kitaev validator
        sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
        sx = np.array([[0, 1], [1, 0]], dtype=complex)
        h12_base = np.kron(sy, sx)
        H12 = K_pauli_proj * np.kron(h12_base, np.eye(2, dtype=complex))
        H23 = K_pauli_proj * np.kron(np.eye(2, dtype=complex), h12_base)

        R12_u = scipy.linalg.expm(-1j * H12 * u)
        R23_uv = scipy.linalg.expm(-1j * H23 * (u + v))
        R12_v = scipy.linalg.expm(-1j * H12 * v)
        R23_v = scipy.linalg.expm(-1j * H23 * v)
        R12_uv = scipy.linalg.expm(-1j * H12 * (u + v))
        R23_u = scipy.linalg.expm(-1j * H23 * u)

        Delta = R12_u @ R23_uv @ R12_v - R23_v @ R12_uv @ R23_u
        N_actual = np.linalg.norm(Delta, ord='fro')

        A_lhs = np.linalg.norm(R12_u @ R23_uv @ R12_v, ord='fro')
        A_rhs = np.linalg.norm(R23_v @ R12_uv @ R23_u, ord='fro')
        A_mean = 0.5 * (A_lhs + A_rhs)
        inner = np.real(np.trace((R12_u @ R23_uv @ R12_v).conj().T @ (R23_v @ R12_uv @ R23_u)))
        cos_theta = inner / (A_lhs * A_rhs + 1e-20)
        theta = np.arccos(np.clip(cos_theta, -1, 1))
        N_theory = 2 * A_mean * np.sin(theta / 2)
        rel_err = np.abs(N_theory - N_actual) / (N_actual + 1e-12)

        res.append((E1, K_fro, K_pauli_proj, theta, N_theory, N_actual, rel_err))

    csvp = os.path.join(outdir, 'toy_bdg_4x4_validation.csv')
    with open(csvp, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['E1', 'K_fro', 'K_pauli_proj', 'theta_rad', 'N_theory', 'N_actual', 'rel_err'])
        for row in res:
            w.writerow(row)

    E = np.array([r[0] for r in res])
    Kf = np.array([r[1] for r in res])
    Kp = np.array([r[2] for r in res])
    Nth = np.array([r[4] for r in res])
    Nac = np.array([r[5] for r in res])
    theta_deg = np.degrees(np.array([r[3] for r in res]))

    plt.figure(figsize=(8, 4))
    plt.plot(E, Kf, '-o', label='K_fro')
    plt.plot(E, Kp, '-^', label='K_pauli_proj')
    plt.xlabel('E1')
    plt.ylabel('K')
    plt.title('Toy BdG 4x4: K measures')
    plt.legend()
    plt.tight_layout()
    p1 = os.path.join(outdir, 'toy_K_compare.png')
    plt.savefig(p1, dpi=200)
    plt.close()

    plt.figure(figsize=(8, 4))
    plt.plot(E, Nac, 'r-o', label='N_actual')
    plt.plot(E, Nth, 'b--', label='N_theory')
    plt.xlabel('E1')
    plt.ylabel('N')
    plt.title('Toy BdG 4x4: phase-driven identity')
    plt.legend()
    plt.tight_layout()
    p2 = os.path.join(outdir, 'toy_N_compare.png')
    plt.savefig(p2, dpi=200)
    plt.close()

    plt.figure(figsize=(8, 3))
    plt.plot(E, theta_deg, 'k-o')
    plt.xlabel('E1')
    plt.ylabel('theta (deg)')
    plt.title('Toy BdG 4x4: matrix angle')
    plt.tight_layout()
    p3 = os.path.join(outdir, 'toy_theta.png')
    plt.savefig(p3, dpi=200)
    plt.close()

    kitaev_csv = os.path.join('quantity', 'kitaev_validation', 'kitaev_propagation_validation.csv')
    summary = {}
    if os.path.exists(kitaev_csv):
        kitaev = np.genfromtxt(kitaev_csv, delimiter=',', names=True)
        summary['corr_Kpauli'] = float(np.corrcoef(Kp, kitaev['K_pauli_proj'][:len(Kp)])[0, 1])
        summary['max_abs_diff_N'] = float(np.max(np.abs(Nac - kitaev['N_actual'][:len(Nac)])))

    report = os.path.join(outdir, 'toy_bdg_4x4_report.md')
    with open(report, 'w') as f:
        f.write('# Toy BdG 4x4 Validation\n\n')
        f.write(f'- CSV: {csvp}\n')
        f.write(f'- K comparison: {p1}\n')
        f.write(f'- N comparison: {p2}\n')
        f.write(f'- theta plot: {p3}\n')
        for k, v in summary.items():
            f.write(f'- {k}: {v:.6e}\n')

    print('Saved results to', outdir)
    return csvp, p1, p2, p3, report


def main():
    csvp, p1, p2, p3, report = run_toy_validation()
    print('Saved toy BdG validation:', csvp)
    print('Saved plots:', p1, p2, p3)
    print('Report:', report)


if __name__ == '__main__':
    main()