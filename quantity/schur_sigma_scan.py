import numpy as np
import matplotlib.pyplot as plt
import os

# Simple Schur-complement Sigma(omega) scanner for toy P-Q-P block
# P and Q are 1x1 (scalars) or small dims; here use scalars for clarity

def sigma_scalar(omega, H_PQ, H_QQ, H_QP):
    # omega: complex scalar
    return H_PQ * (1.0 / (omega - H_QQ)) * H_QP


def run_scan(params):
    E_grid = np.linspace(-1.0, 1.0, 1001)
    eta = params.get('eta', 1e-3)
    H_PQ = params['H_PQ']
    H_QQ = params['H_QQ']
    H_QP = params['H_QP']

    Sig = np.array([sigma_scalar(E + 1j * eta, H_PQ, H_QQ, H_QP) for E in E_grid])
    return E_grid, Sig


def make_plots(outdir):
    os.makedirs(outdir, exist_ok=True)
    # parameters: use the toy values from doc: H_PQ = i t , H_QP = -i mu/2
    t = 0.5
    mu = 0.2
    H_PQ = 1j * t
    H_QP = -1j * mu / 2.0

    # Kitaev-point simplification H_QQ = 0
    params_kitaev = {'H_PQ': H_PQ, 'H_QQ': 0.0, 'H_QP': H_QP, 'eta': 1e-3}
    E, Sig_k = run_scan(params_kitaev)

    # Off-point: H_QQ = 0.1 (small finite Q-Q block)
    params_off = {'H_PQ': H_PQ, 'H_QQ': 0.15 + 0.0j, 'H_QP': H_QP, 'eta': 1e-3}
    E, Sig_o = run_scan(params_off)

    # Plot Im Sigma
    plt.figure(figsize=(6, 4))
    plt.plot(E, np.imag(Sig_k), label='Im Sigma (Kitaev H_QQ=0)')
    plt.plot(E, np.imag(Sig_o), label='Im Sigma (H_QQ=0.15)')
    plt.xlabel('Energy E')
    plt.ylabel('Im Sigma(E + i eta)')
    plt.title('Imaginary part of Schur complement (broadening)')
    plt.legend()
    plt.tight_layout()
    p1 = os.path.join(outdir, 'sigma_im.png')
    plt.savefig(p1, dpi=200)
    plt.close()

    # Plot Re Sigma
    plt.figure(figsize=(6, 4))
    plt.plot(E, np.real(Sig_k), label='Re Sigma (Kitaev H_QQ=0)')
    plt.plot(E, np.real(Sig_o), label='Re Sigma (H_QQ=0.15)')
    plt.xlabel('Energy E')
    plt.ylabel('Re Sigma(E + i eta)')
    plt.title('Real part of Schur complement (Lamb shift)')
    plt.legend()
    plt.tight_layout()
    p2 = os.path.join(outdir, 'sigma_re.png')
    plt.savefig(p2, dpi=200)
    plt.close()

    # Save simple numeric summary
    with open(os.path.join(outdir, 'sigma_summary.txt'), 'w') as f:
        f.write('Parameters:\n')
        f.write(f't = {t}, mu = {mu}\n')
        f.write('H_PQ = i t, H_QP = -i mu/2\n')
        f.write('\nKitaev H_QQ = 0 results:\n')
        f.write(f'Max Im Sigma = {np.max(np.imag(Sig_k)):.6e}\n')
        f.write(f'Min Im Sigma = {np.min(np.imag(Sig_k)):.6e}\n')
        f.write('\nOff-point H_QQ = 0.15 results:\n')
        f.write(f'Max Im Sigma = {np.max(np.imag(Sig_o)):.6e}\n')
        f.write(f'Min Im Sigma = {np.min(np.imag(Sig_o)):.6e}\n')

    return p1, p2, os.path.join(outdir, 'sigma_summary.txt')

if __name__ == '__main__':
    out = 'quantity/sigma_scan'
    p1, p2, txt = make_plots(out)
    print('Saved:', p1, p2, txt)
