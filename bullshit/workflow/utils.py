"""
Utility functions: Pauli mapping, BdG construction, diagnostics.
"""

import numpy as np
from scipy.linalg import eigh

# Pauli matrices
I2 = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

def kron(a, b):
    """Kronecker product."""
    return np.kron(a, b)

def pauli_expand(H4):
    """
    Expand 4x4 Hamiltonian in Pauli tensor basis.
    Returns dict h[(a,b)] with a,b in {'I','X','Y','Z'}.
    """
    syms = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
    h = {}
    for a, A in syms.items():
        for b, B in syms.items():
            P = kron(A, B)
            h[(a, b)] = 0.25 * np.trace(P @ H4).real
    return h

def map_to_kitaev(h):
    """
    Map Pauli coefficients h to Kitaev chain parameters (t, Delta, mu).
    
    Returns:
        t, Delta, mu: complex
    """
    c_xx = h.get(('X', 'X'), 0.0)
    c_yy = h.get(('Y', 'Y'), 0.0)
    c_xy = h.get(('X', 'Y'), 0.0)
    c_yx = h.get(('Y', 'X'), 0.0)
    c_zz = h.get(('Z', 'Z'), 0.0)
    c_z0 = h.get(('Z', 'I'), 0.0)
    c_0z = h.get(('I', 'Z'), 0.0)
    
    t = c_xx + c_yy + 1j * (c_xy - c_yx)
    Delta = c_xx - c_yy - 1j * (c_xy + c_yx)
    mu = 4.0 * c_zz - 2.0 * (c_z0 + c_0z)
    
    return complex(t), complex(Delta), complex(mu)

def extract_d_vector(h):
    """
    Extract effective d vector from Pauli coefficients.
    d = (d_x, d_y, d_z) + d_0
    """
    d_x = h.get(('X', 'X'), 0.0) + h.get(('Y', 'Y'), 0.0)
    d_y = -h.get(('X', 'Y'), 0.0) + h.get(('Y', 'X'), 0.0)
    d_z = h.get(('Z', 'I'), 0.0) - h.get(('I', 'Z'), 0.0)
    d_0 = h.get(('I', 'I'), 0.0) - h.get(('Z', 'Z'), 0.0)
    
    return np.array([d_x, d_y, d_z]), d_0

def build_bdg_chain(L, t, Delta, mu, open_boundary=True):
    """
    Build single-particle BdG Hamiltonian for a Kitaev chain.
    
    Returns:
        H_bdg: (2L, 2L) hermitian matrix
    """
    h = np.zeros((L, L), dtype=complex)
    Delta_mat = np.zeros((L, L), dtype=complex)
    
    # Onsite
    for i in range(L):
        h[i, i] = -mu / 2.0
    
    # Hopping and pairing on links
    for i in range(L - 1):
        h[i, i + 1] = -t / 2.0
        h[i + 1, i] = -t.conjugate() / 2.0
        Delta_mat[i, i + 1] = Delta / 2.0
        Delta_mat[i + 1, i] = -Delta / 2.0
    
    # Periodic boundary (optional)
    if not open_boundary:
        h[0, L - 1] = -t / 2.0
        h[L - 1, 0] = -t.conjugate() / 2.0
        Delta_mat[0, L - 1] = Delta / 2.0
        Delta_mat[L - 1, 0] = -Delta / 2.0
    
    # Construct full BdG
    top = np.hstack([h, Delta_mat])
    bottom = np.hstack([-Delta_mat.conjugate(), -h.T])
    H_bdg = np.vstack([top, bottom])
    H_bdg = 0.5 * (H_bdg + H_bdg.conj().T)  # Ensure hermiticity
    
    return H_bdg

def diagonalize_bdg(H_bdg):
    """
    Diagonalize BdG and return eigenvalues, eigenvectors.
    """
    E, V = eigh(H_bdg)
    idx = np.argsort(np.abs(E))  # Sort by absolute value
    return E[idx], V[:, idx]

def compute_edge_weight(v, L, edge_sites=1):
    """
    Compute edge weight (left + right) from a BdG eigenvector.
    
    v: length 2L, first L are particle, next L are hole
    """
    probs = np.abs(v[:L])**2 + np.abs(v[L:])**2
    left = float(np.sum(probs[:edge_sites]))
    right = float(np.sum(probs[-edge_sites:]))
    return left, right

def compute_bulk_gap(t, Delta, mu, nk=801):
    """
    Compute bulk gap from Kitaev parameters (k-space approach).
    """
    ks = np.linspace(0.0, np.pi, nk)
    min_gap = 1e9
    for k in ks:
        eik = np.exp(1j * k)
        h_k = -mu / 2.0 - 0.5 * (t * eik + t.conjugate() * np.conjugate(eik))
        Delta_k = 0.5 * (Delta * eik - Delta * np.conjugate(eik))
        # 2x2 BdG eigenvalues
        H = np.array([[h_k, Delta_k], [Delta_k.conjugate(), -h_k.conjugate()]], dtype=complex)
        vals = np.linalg.eigvalsh(H)
        min_e = np.min(np.abs(vals))
        if min_e < min_gap:
            min_gap = min_e
    return float(min_gap)

def topo_criterion(t, Delta, mu, gap_threshold=1e-12):
    """
    Simple topological criterion: |mu| < 2|t| and gap > threshold.
    """
    is_topo = (np.abs(mu) < 2.0 * np.abs(t)) and (gap_threshold > 0)
    return bool(is_topo)

def compute_ldos(H_bdg, energies, eta=1e-3):
    """
    Compute local density of states.
    
    Returns:
        energies: energy axis
        ldos: (L, len(energies)) matrix
        E: full spectrum
    """
    E, V = eigh(H_bdg)
    idx = np.argsort(E)
    E = E[idx]
    V = V[:, idx]
    
    L = H_bdg.shape[0] // 2
    nE = len(energies)
    ldos = np.zeros((L, nE), dtype=float)
    
    lorentz_norm = eta / np.pi
    for n in range(len(E)):
        En = float(E[n])
        probs = np.abs(V[:L, n])**2 + np.abs(V[L:, n])**2
        lor = lorentz_norm / ((energies - En)**2 + eta**2)
        ldos += np.outer(probs, lor)
    
    return energies, ldos, E
