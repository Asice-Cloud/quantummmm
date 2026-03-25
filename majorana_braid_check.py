import numpy as np
from scipy.linalg import eigh

def build_bdg_two_site(a=0.0,b=0.0,c=0.0,d=0.0,mu0=0.0):
    """Construct 4x4 BdG single-particle Hamiltonian for a 2-site Kitaev chain.
    Parameters follow user's R->Kitaev mapping conventions.
    Returns H_bdg (4x4 complex).
    """
    t = b + c
    Delta = b - c
    mu = 4.0 * d + mu0
    # h matrix (2x2)
    h = np.array([[-mu/2.0, -t],[-t, -mu/2.0]], dtype=complex)
    # pairing matrix (2x2)
    Delta_mat = np.array([[0.0, Delta],[-Delta, 0.0]], dtype=complex)
    # BdG form: [[h, Delta],[ -Delta.conj(), -h.T]]
    top = np.hstack([h, Delta_mat])
    bottom = np.hstack([-Delta_mat.conj(), -h.T])
    H = np.vstack([top, bottom])
    return H, {'t':t,'Delta':Delta,'mu':mu}

def single_particle_spectrum(H):
    """Return eigenvalues (sorted) of Hermitian BdG H."""
    w, v = eigh(H)
    return w, v

def extract_gap_and_zero_modes(H):
    w, v = single_particle_spectrum(H)
    # sort by absolute value
    idx = np.argsort(np.abs(w))
    w_sorted = w[idx]
    return w_sorted

# simple fitting helper
def linear_fit(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    A = np.vstack([x, np.ones_like(x)]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    # compute R^2
    y_pred = m * x + c
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - ss_res/ss_tot if ss_tot != 0 else 1.0
    return m, c, r2

if __name__ == '__main__':
    # quick self-test
    H, params = build_bdg_two_site(b=0.3,c=0.1,d=0.02,mu0=0.0)
    w = extract_gap_and_zero_modes(H)
    print('params:', params)
    print('eigs:', np.round(w,6))
