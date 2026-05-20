import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

import importlib.util
spec = importlib.util.spec_from_file_location('embed_kitaev', os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'tools', 'embed_kitaev.py')))
embed_kitaev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(embed_kitaev)
build_bdg = embed_kitaev.build_bdg

def indices_for_sites(sites, L):
    # return BdG indices (c part then c† part) for given site list
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

def embed_op_on_three(op, which):
    # op is 4x4 acting on sites [0,1] or [1,2]; which='12' or '23'
    # embed into 6x6 ordering [site0,site1,site2] with Nambu ordering (c,c† concatenated)
    I2 = np.eye(2, dtype=complex)
    if which == '12':
        # op ⊕ I2
        top = np.hstack([op, np.zeros((4,2), dtype=complex)])
        bot = np.hstack([np.zeros((2,4), dtype=complex), I2])
        M = np.vstack([top, bot])
    else:
        # I2 ⊕ op
        top = np.hstack([I2, np.zeros((2,4), dtype=complex)])
        bot = np.hstack([np.zeros((4,2), dtype=complex), op])
        M = np.vstack([top, bot])
    return M

def run_validation(outdir='quantity/kitaev_validation', L=40, E1_vals=None):
    os.makedirs(outdir, exist_ok=True)
    if E1_vals is None:
        E1_vals = np.linspace(0.001, 0.2, 20)

    t0 = 1.0
    Delta0 = 0.3
    VD = 2.0

    res = []
    for E1 in E1_vals:
        # build chain with global mu = E1, plus QD at left implemented in build_bdg if needed
        mu0 = E1
        H = build_bdg(mu0 * np.ones(L), t0 * np.ones(L - 1), Delta0 * np.ones(L - 1))

        # pick center three sites
        c = L // 2
        sites = [c - 1, c, c + 1]
        idx3 = indices_for_sites(sites, L)
        H3 = H[np.ix_(idx3, idx3)]  # 6x6 BdG for three sites

        # partition P = sites 0 and 2, Q = site1
        P_sites = [0, 2]
        Q_sites = [1]
        P_idx = indices_for_sites(P_sites, 3)  # indices in 3-site block
        Q_idx = indices_for_sites(Q_sites, 3)

        omega = 0.0 + 1j * 1e-3
        H_PP, Sigma = schur_sigma_bdG(omega, H3, P_idx, Q_idx)
        K = np.linalg.norm(Sigma, ord='fro')

        # Pauli projection: transform Sigma from BdG (c,c†) basis to Majorana basis
        # ordering of Sigma (P-space) is [c_site0, c_site1, c_site0†, c_site1†]
        # build linear map M such that gamma = M @ psi_old, where
        # gamma = [γ0_1, γ0_2, γ1_1, γ1_2] and psi_old = [c0, c1, c0†, c1†]
        M = np.array([
            [1, 0, 1, 0],
            [-1j, 0, 1j, 0],
            [0, 1, 0, 1],
            [0, -1j, 0, 1j]
        ], dtype=complex)
        Minv = np.linalg.inv(M)
        Sigma_maj = M @ Sigma @ Minv

        # project onto Pauli tensor operator in Majorana basis
        sx = np.array([[0,1],[1,0]], dtype=complex)
        sy = np.array([[0,-1j],[1j,0]], dtype=complex)
        B = np.kron(sy, sx)
        denom = np.real(np.trace(B.conj().T @ B))
        # project coefficient (use full complex trace and take magnitude)
        coeff = np.trace(B.conj().T @ Sigma_maj)
        K_pauli_proj = float(np.abs(coeff) / (denom + 1e-30))

        # construct H12 and H23 from H3
        # H3 ordering: [c0,c1,c2, c0†,c1†,c2†]
        idx12 = indices_for_sites([0,1], 3)
        idx23 = indices_for_sites([1,2], 3)
        H12 = H3[np.ix_(idx12, idx12)]
        H23 = H3[np.ix_(idx23, idx23)]

        # exponentiate to get local evolution
        u = 0.5
        v = 0.7
        R12_u = expm(-1j * H12 * u)
        R12_v = expm(-1j * H12 * v)
        R23_u = expm(-1j * H23 * u)
        R23_v = expm(-1j * H23 * v)
        R23_uv = expm(-1j * H23 * (u + v))
        R12_uv = expm(-1j * H12 * (u + v))

        # embed into 6x6
        R12_u_e = embed_op_on_three(R12_u, '12')
        R12_v_e = embed_op_on_three(R12_v, '12')
        R12_uv_e = embed_op_on_three(R12_uv, '12')
        R23_u_e = embed_op_on_three(R23_u, '23')
        R23_v_e = embed_op_on_three(R23_v, '23')
        R23_uv_e = embed_op_on_three(R23_uv, '23')

        LHS = R12_u_e @ R23_uv_e @ R12_v_e
        RHS = R23_v_e @ R12_uv_e @ R23_u_e
        Delta = LHS - RHS
        N_actual = np.linalg.norm(Delta, ord='fro')

        # angle and theory
        A_lhs = np.linalg.norm(LHS, ord='fro')
        A_rhs = np.linalg.norm(RHS, ord='fro')
        A_mean = 0.5 * (A_lhs + A_rhs)
        inner = np.real(np.trace(LHS.conj().T @ RHS))
        cos_theta = inner / (A_lhs * A_rhs + 1e-20)
        theta = np.arccos(np.clip(cos_theta, -1, 1))
        N_theory = 2 * A_mean * np.sin(theta / 2)
        rel_err = np.abs(N_theory - N_actual) / (N_actual + 1e-12)

        res.append((E1, K, K_pauli_proj, A_lhs, A_rhs, theta, N_theory, N_actual, rel_err))

    # save CSV
    import csv
    csvp = os.path.join(outdir, 'kitaev_propagation_validation.csv')
    with open(csvp, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['E1','K_fro','K_pauli_proj','A_lhs','A_rhs','theta_rad','N_theory','N_actual','rel_err'])
        for row in res:
            w.writerow(row)

    # quick plot N vs E1 and theta
    E = np.array([r[0] for r in res])
    # unpack with new K_pauli_proj field
    Nth = np.array([r[6] for r in res])
    Nac = np.array([r[7] for r in res])
    theta_deg = np.degrees(np.array([r[4] for r in res]))
    Kpa_proj = np.array([r[2] for r in res])
    Kfro = np.array([r[1] for r in res])

    plt.figure(figsize=(8,4))
    plt.plot(E, Nac, 'r-o', label='N_actual')
    plt.plot(E, Nth, 'b--', label='N_theory')
    plt.xlabel('E1')
    plt.ylabel('N')
    plt.legend()
    plt.title('Kitaev local propagation: N theory vs actual')
    plt.tight_layout()
    p1 = os.path.join(outdir, 'kitaev_N_compare.png')
    plt.savefig(p1, dpi=200)
    plt.close()

    plt.figure(figsize=(8,3))
    plt.plot(E, theta_deg, 'k-o')
    plt.xlabel('E1')
    plt.ylabel('theta (deg)')
    plt.title('Matrix angle between LHS and RHS')
    plt.tight_layout()
    p2 = os.path.join(outdir, 'kitaev_theta.png')
    plt.savefig(p2, dpi=200)
    plt.close()

    # plot K_fro vs K_pauli_proj
    plt.figure(figsize=(6,4))
    plt.scatter(Kfro, Kpa_proj, s=12)
    mn = min(Kfro.min(), Kpa_proj.min())
    mx = max(Kfro.max(), Kpa_proj.max())
    plt.plot([mn,mx],[mn,mx],'k--',linewidth=0.7)
    plt.xlabel('K_fro')
    plt.ylabel('K_pauli_proj')
    plt.title('K_fro vs K_pauli_proj (Kitaev)')
    plt.tight_layout()
    p3 = os.path.join(outdir, 'kitaev_Kproj_compare.png')
    plt.savefig(p3, dpi=200)
    plt.close()

    # plot K vs E1 (both Frobenius and Pauli-projected)
    plt.figure(figsize=(8,4))
    plt.plot(E, Kfro, '-o', label='K_fro')
    plt.plot(E, Kpa_proj, '-s', label='K_pauli_proj')
    plt.xlabel('E1')
    plt.ylabel('K')
    plt.title('K vs E1 (Kitaev local propagation)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    p4 = os.path.join(outdir, 'kitaev_K_vs_E1.png')
    plt.savefig(p4, dpi=200)
    plt.close()

    print('Saved CSV:', csvp)
    print('Saved plots:', p1, p2, p3)
    return csvp, p1, p2, p3


if __name__ == '__main__':
    run_validation()
