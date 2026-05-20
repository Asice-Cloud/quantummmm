"""
Pauli Projection for K Extraction from Schur Complement Σ(ω)

Instead of using scalar norms (Frobenius, spectral, nuclear) as proxy for K,
we project Σ onto Pauli component spaces and extract the coefficient directly
related to the σ^y ⊗ σ^x structure used in h_eff.

For a 3×3 block Σ in P-space, we decompose:
  Σ_ij = sum over Pauli basis on (i,j) pairs

We focus on the Pauli coefficients that couple odd-numbered Majoranas,
which corresponds to the σ^y ⊗ σ^x coupling structure.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import importlib.util
import sys
import scipy.linalg

# load schur_sigma_majorana
spec = importlib.util.spec_from_file_location('sch', os.path.join(os.path.dirname(__file__), 'schur_sigma_majorana.py'))
sch = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = sch
spec.loader.exec_module(sch)
build_A = sch.build_A
H_bdG_from_A = sch.H_bdG_from_A
schur_sigma = sch.schur_sigma

# Pauli matrices
I = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
h12_base = np.kron(sy, sx)

def op12(K):
    return K * np.kron(h12_base, I)

def op23(K):
    return K * np.kron(I, h12_base)

def compute_N_from_K(K, u=0.5, v=0.7):
    H12 = op12(K)
    H23 = op23(K)
    R12_u = scipy.linalg.expm(-1j * H12 * u)
    R23_uv = scipy.linalg.expm(-1j * H23 * (u + v))
    R12_v = scipy.linalg.expm(-1j * H12 * v)
    R23_v = scipy.linalg.expm(-1j * H23 * v)
    R12_uv = scipy.linalg.expm(-1j * H12 * (u + v))
    R23_u = scipy.linalg.expm(-1j * H23 * u)
    Delta = R12_u @ R23_uv @ R12_v - R23_v @ R12_uv @ R23_u
    N = np.sqrt(np.real(np.trace(Delta.conj().T @ Delta)))
    return float(N)


def extract_pauli_coefficients(Sigma):
    """
    Extract Pauli coefficients from a 3×3 Schur complement matrix Σ.
    
    For Majorana operators γ_i (i=1,3,5 in P-space), we have:
      {γ_i, γ_j} = 2δ_ij
    
    The 3×3 matrix Σ acts on the P-space. We decompose it in terms of 
    Pauli matrices on three-qubit Hilbert space (or direct Pauli projections).
    
    For our effective model, the key coupling is σ^y ⊗ σ^x.
    We extract the coefficient as:
      K_pauli = |Im(Tr(Σ† · σ^y ⊗ σ^x))|  or similar projection
    
    A simple approach: extract the dominant off-diagonal coupling from Σ.
    """
    # Strict JW-compatible projection in the Majorana basis.
    # For the reduced P-space basis [γ1, γ3, γ5], the JW image of the
    # target Pauli string σ^y ⊗ σ^x is represented by the symmetric
    # Majorana bilinear connecting the first and third basis elements.
    # This is the operator that survives in the toy Schur block.
    B = np.zeros((3, 3), dtype=complex)
    B[0, 2] = 1.0
    B[2, 0] = 1.0

    denom = np.real(np.trace(B.conj().T @ B))
    coeff = np.trace(B.conj().T @ Sigma)
    K_pauli = float(np.abs(coeff) / (denom + 1e-30))

    # additional diagnostics (off-diagonal measures on original Sigma)
    off_diag_im = np.mean(np.abs(np.imag(Sigma - np.diag(np.diag(Sigma)))))
    off_diag_re = np.mean(np.abs(np.real(Sigma - np.diag(np.diag(Sigma)))))
    off_diag_matrix = Sigma - np.diag(np.diag(Sigma))
    max_off_diag = np.max(np.abs(off_diag_matrix))
    im_trace = np.imag(np.trace(Sigma))

    return {
        'K_pauli': K_pauli,
        'off_diag_im': off_diag_im,
        'off_diag_re': off_diag_re,
        'max_off_diag': max_off_diag,
        'im_trace': im_trace
    }


def run_pauli_projection(outdir='quantity/pauli_projection'):
    os.makedirs(outdir, exist_ok=True)

    t1 = 0.03
    t2 = 0.03
    t3 = 0.03
    Ed = 1.0

    E1_vals = np.linspace(0.001, 0.2, 40)
    u = 0.5
    v = 0.7

    P_idx = [1, 3, 5]
    Q_idx = [0, 2, 4]

    eta = 1e-3

    results = []

    for E1 in E1_vals:
        params = {'E1': float(E1), 'E2': 0.0}
        A = build_A(params, t1, t2, t3, Ed)
        H = H_bdG_from_A(A)

        omega = 0.0 + 1j * eta
        Sigma = schur_sigma(omega, H, P_idx, Q_idx)

        # Extract multiple K measures
        K_fro = np.linalg.norm(Sigma, ord='fro')
        K_spec = np.linalg.norm(Sigma, ord=2)
        
        pauli_dict = extract_pauli_coefficients(Sigma)
        K_pauli = pauli_dict['K_pauli']

        # Compute N from each K
        N_fro = compute_N_from_K(K_fro, u=u, v=v)
        N_spec = compute_N_from_K(K_spec, u=u, v=v)
        N_pauli = compute_N_from_K(K_pauli, u=u, v=v)

        results.append({
            'E1': float(E1),
            'K_fro': float(K_fro),
            'K_spec': float(K_spec),
            'K_pauli': float(K_pauli),
            'off_diag_im': float(pauli_dict['off_diag_im']),
            'off_diag_re': float(pauli_dict['off_diag_re']),
            'max_off_diag': float(pauli_dict['max_off_diag']),
            'im_trace': float(pauli_dict['im_trace']),
            'N_fro': float(N_fro),
            'N_spec': float(N_spec),
            'N_pauli': float(N_pauli)
        })

    # Save CSV
    csvp = os.path.join(outdir, 'pauli_projection_results.csv')
    header = 'E1,K_fro,K_spec,K_pauli,off_diag_im,off_diag_re,max_off_diag,im_trace,N_fro,N_spec,N_pauli\n'
    with open(csvp, 'w') as f:
        f.write(header)
        for r in results:
            vals = [r[k] for k in ['E1','K_fro','K_spec','K_pauli','off_diag_im','off_diag_re','max_off_diag','im_trace','N_fro','N_spec','N_pauli']]
            f.write(','.join(str(v) for v in vals))
            f.write('\n')

    # Extract arrays for plotting
    E = np.array([r['E1'] for r in results])
    Kf = np.array([r['K_fro'] for r in results])
    Ks = np.array([r['K_spec'] for r in results])
    Kp = np.array([r['K_pauli'] for r in results])
    Nf = np.array([r['N_fro'] for r in results])
    Ns = np.array([r['N_spec'] for r in results])
    Np = np.array([r['N_pauli'] for r in results])

    # Plot K comparison
    plt.figure(figsize=(8, 5))
    plt.plot(E, Kf, '-o', label='K_Frobenius', alpha=0.7)
    plt.plot(E, Ks, '-s', label='K_spectral', alpha=0.7)
    plt.plot(E, Kp, '-^', label='K_Pauli', alpha=0.7)
    plt.xlabel('E1 (chemical potential proxy)')
    plt.ylabel('K measures')
    plt.title('K extraction methods: Frobenius vs Spectral vs Pauli')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'K_comparison.png'), dpi=200)
    plt.close()

    # Plot N comparison
    plt.figure(figsize=(8, 5))
    plt.plot(E, Nf, '-o', label='N from K_Fro', alpha=0.7)
    plt.plot(E, Ns, '-s', label='N from K_Spec', alpha=0.7)
    plt.plot(E, Np, '-^', label='N from K_Pauli', alpha=0.7)
    plt.xlabel('E1')
    plt.ylabel('N (non-Abelian indicator)')
    plt.title('Non-Abelian indicator: Impact of K extraction method')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'N_comparison.png'), dpi=200)
    plt.close()

    # Correlations
    from scipy.stats import pearsonr
    try:
        corr_KpNp, p_val = pearsonr(Kp, Np)
    except:
        corr_KpNp = float('nan')

    try:
        corr_KfNf, _ = pearsonr(Kf, Nf)
    except:
        corr_KfNf = float('nan')

    # Generate report
    report_path = os.path.join(outdir, 'pauli_projection_report.md')
    with open(report_path, 'w') as f:
        f.write('# Pauli Projection K Extraction Report\n\n')
        f.write('## Summary\n\n')
        f.write('We extract K using Pauli projections from the Schur complement Σ(ω),\n')
        f.write('focusing on off-diagonal matrix elements that encode Majorana couplings.\n\n')
        f.write('## Key Results\n\n')
        f.write(f'- K_Pauli captures off-diagonal coupling structure: mean={np.mean(Kp):.6e}, std={np.std(Kp):.6e}\n')
        f.write(f'- Pearson(K_Fro, N_Fro) = {corr_KfNf:.6f}\n')
        f.write(f'- Pearson(K_Pauli, N_Pauli) = {corr_KpNp:.6f}\n\n')
        f.write('## Interpretation\n\n')
        f.write('Pauli projection extracts coupling from matrix structure rather than global norms.\n')
        f.write('This may provide different scaling and physical interpretation compared to Frobenius.\n')

    print(f'Saved results to {outdir}')
    print(f'  - CSV: {csvp}')
    print(f'  - K comparison: {os.path.join(outdir, "K_comparison.png")}')
    print(f'  - N comparison: {os.path.join(outdir, "N_comparison.png")}')
    print(f'  - Report: {report_path}')

    return csvp, outdir


if __name__ == '__main__':
    run_pauli_projection()
