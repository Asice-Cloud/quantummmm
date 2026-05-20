import numpy as np
import matplotlib.pyplot as plt
import os

# Build Majorana single-particle coupling matrix A (6x6) from reproduce_effective_braiding_majorana conventions
# Indices: 0..5 correspond to [gamma1, gamma2, gamma3, gamma4, gamma_a, gamma_b]


def build_A(params, t1, t2, t3, Ed):
    A = np.zeros((6, 6), dtype=complex)
    # set antisymmetric entries according to build_H terms (coeff before multiplying by i)
    # A[i,j] = -A[j,i]
    E1 = params.get('E1', 0.0)
    E2 = params.get('E2', 0.0)

    def set_pair(i, j, val):
        A[i, j] = val
        A[j, i] = -val

    set_pair(0, 1, E1)   # γ0-γ1 (E1)
    set_pair(2, 3, E2)   # γ2-γ3 (E2)
    set_pair(4, 5, Ed)   # γ4-γ5 (Ed)
    set_pair(4, 1, t2)   # γ4-γ1 (t2)
    set_pair(5, 0, -t1)  # γ5-γ0 (-t1)
    set_pair(4, 2, -t3)  # γ4-γ2 (-t3)
    return A


def H_bdG_from_A(A):
    # choose H_bdG = 1j * A (generator for Majorana quadratic form)
    return 1j * A


def schur_sigma(omega, H, P_idx, Q_idx):
    # H is full single-particle matrix (n x n)
    # partition and compute Sigma_PP = H_PQ (omega - H_QQ)^{-1} H_QP
    H_PP = H[np.ix_(P_idx, P_idx)]
    H_PQ = H[np.ix_(P_idx, Q_idx)]
    H_QP = H[np.ix_(Q_idx, P_idx)]
    H_QQ = H[np.ix_(Q_idx, Q_idx)]
    inv = np.linalg.inv(omega * np.eye(len(Q_idx), dtype=complex) - H_QQ)
    Sigma = H_PQ @ inv @ H_QP
    return Sigma


def run_scan(outdir):
    os.makedirs(outdir, exist_ok=True)
    # choose parameters similar to toy: t1,t2,t3, Ed, params E1,E2
    params = {'E1': 0.001, 'E2': 0.0}
    t1 = 0.03
    t2 = 0.03
    t3 = 0.03
    Ed = 1.0

    A = build_A(params, t1, t2, t3, Ed)
    H = H_bdG_from_A(A)

    # partition: choose P = odd indices (1,3,5) as gamma_B; Q = even indices (0,2,4) as gamma_A
    P_idx = [1, 3, 5]
    Q_idx = [0, 2, 4]

    E_grid = np.linspace(-2.0, 2.0, 1001)
    eta = 1e-3
    sig_vals = []
    for E in E_grid:
        omega = E + 1j * eta
        Sigma = schur_sigma(omega, H, P_idx, Q_idx)
        # for simplicity, record trace and selected diagonal elements
        sig_vals.append(Sigma)
    sig_vals = np.array(sig_vals)  # shape (len(E), 3, 3)

    # pick element (0,0) and trace
    sigma00 = sig_vals[:, 0, 0]
    sigma_trace = np.trace(sig_vals, axis1=1, axis2=2)

    # plot Im Sigma_trace and Re
    plt.figure(figsize=(6, 4))
    plt.plot(E_grid, np.imag(sigma_trace), label='Im Tr Sigma')
    plt.xlabel('Energy E')
    plt.ylabel('Im Tr Sigma')
    plt.title('Im Tr Sigma(omega) -- Majorana 6-mode')
    plt.legend()
    p1 = os.path.join(outdir, 'maj_sigma_im_trace.png')
    plt.tight_layout()
    plt.savefig(p1, dpi=200)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.plot(E_grid, np.real(sigma_trace), label='Re Tr Sigma')
    plt.xlabel('Energy E')
    plt.ylabel('Re Tr Sigma')
    plt.title('Re Tr Sigma(omega) -- Majorana 6-mode')
    plt.legend()
    p2 = os.path.join(outdir, 'maj_sigma_re_trace.png')
    plt.tight_layout()
    plt.savefig(p2, dpi=200)
    plt.close()

    # element (0,0)
    plt.figure(figsize=(6, 4))
    plt.plot(E_grid, np.imag(sigma00), label='Im Sigma_00')
    plt.xlabel('E')
    plt.ylabel('Im Sigma_00')
    plt.title('Im Sigma_00 (P[0],P[0])')
    plt.legend()
    p3 = os.path.join(outdir, 'maj_sigma_im_00.png')
    plt.tight_layout()
    plt.savefig(p3, dpi=200)
    plt.close()

    # summary
    with open(os.path.join(outdir, 'maj_sigma_summary.txt'), 'w') as f:
        f.write('Parameters:\n')
        f.write(f"t1={t1}, t2={t2}, t3={t3}, Ed={Ed}, E1={params['E1']}\n")
        f.write('\nTrace Sigma extremes:\n')
        f.write(f"max Im Tr Sigma = {np.max(np.imag(sigma_trace)):.6e}\n")
        f.write(f"min Im Tr Sigma = {np.min(np.imag(sigma_trace)):.6e}\n")
    return p1, p2, p3, os.path.join(outdir, 'maj_sigma_summary.txt')

if __name__ == '__main__':
    out = 'quantity/sigma_scan_majorana'
    files = run_scan(out)
    print('Saved:', files)
