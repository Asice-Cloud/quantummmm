import numpy as np
import matplotlib.pyplot as plt
import os
import importlib.util
import sys

# dynamic import of schur_sigma_majorana
spec = importlib.util.spec_from_file_location('schur_sigma_majorana', os.path.join(os.path.dirname(__file__), 'schur_sigma_majorana.py'))
sch = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = sch
spec.loader.exec_module(sch)
build_A = sch.build_A
H_bdG_from_A = sch.H_bdG_from_A
schur_sigma = sch.schur_sigma

# Build effective K from Sigma: use Re Tr Sigma at omega ~ 0

def compute_K_from_Sigma(Sigma):
    # derive an effective scalar coupling from Sigma: use Frobenius norm
    val = np.linalg.norm(Sigma, ord='fro')
    return float(val)

# Pauli matrices
I = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)

# Two-site Pauli operator: σ^y ⊗ σ^x
h12_base = np.kron(sy, sx)

# Three-site operators
def op12(K):
    return K * np.kron(h12_base, I)

def op23(K):
    return K * np.kron(I, h12_base)


def R_of_H(H, param):
    # exponentiate -i H param
    return scipy.linalg.expm(-1j * H * param)


if __name__ == '__main__':
    import scipy.linalg
    outdir = 'quantity/delta_N_results'
    os.makedirs(outdir, exist_ok=True)

    # fixed hopping params
    t1 = 0.03
    t2 = 0.03
    t3 = 0.03
    Ed = 1.0

    taus = np.linspace(0.1, 2.0, 10)
    E1_vals = np.linspace(0.001, 0.2, 40)  # sweep E1 as mu proxy
    N_vals = []

    # choose u and v for braiding-like sequence
    u = 0.5
    v = 0.7

    for E1 in E1_vals:
        params = {'E1': float(E1), 'E2': 0.0}
        # build full H
        A = build_A(params, t1, t2, t3, Ed)
        H = H_bdG_from_A(A)
        # partition P indices and Q indices same as in schur_sigma_majorana
        P_idx = [1, 3, 5]
        Q_idx = [0, 2, 4]
        omega = 0.0 + 1j * 1e-3
        Sigma = schur_sigma(omega, H, P_idx, Q_idx)
        K = compute_K_from_Sigma(Sigma)
        # Build h_eff operators (3-site Hilbert space dim=8)
        H12 = op12(K)
        H23 = op23(K)
        # exponentiate
        R12_u = scipy.linalg.expm(-1j * H12 * u)
        R23_uv = scipy.linalg.expm(-1j * H23 * (u + v))
        R12_v = scipy.linalg.expm(-1j * H12 * v)
        R23_v = scipy.linalg.expm(-1j * H23 * v)
        R12_uv = scipy.linalg.expm(-1j * H12 * (u + v))
        R23_u = scipy.linalg.expm(-1j * H23 * u)

        # compute Delta
        Delta = R12_u @ R23_uv @ R12_v - R23_v @ R12_uv @ R23_u
        N = np.sqrt(np.real(np.trace(Delta.conj().T @ Delta)))
        N_vals.append(float(N))

    N_vals = np.array(N_vals)

    plt.figure()
    plt.plot(E1_vals, N_vals, '-o')
    plt.xlabel('E1 (mu proxy)')
    plt.ylabel('N = ||Delta||_F')
    plt.title('Non-Abelian indicator N vs E1 (from Schur Sigma)')
    plt.grid(True)
    fimg = os.path.join(outdir, 'N_vs_E1.png')
    plt.tight_layout()
    plt.savefig(fimg, dpi=200)

    with open(os.path.join(outdir, 'N_vs_E1.txt'), 'w') as f:
        for e, n in zip(E1_vals, N_vals):
            f.write(f"{e:.6e} {n:.6e}\n")

    print('Saved:', fimg, os.path.join(outdir, 'N_vs_E1.txt'))
