r"""Minimal two-site Pauli model and chain BdG utilities.

Implements H(u) = a(u) XX + b(u) YY + c(u) XY + d(u) YX on two spins,
expands on Pauli basis, maps Pauli coefficients to Kitaev params (t, Delta, mu),
and builds a 1D BdG chain for numeric tests.

This module is intentionally small and has no external dependencies
other than numpy/scipy which are common in the workspace.
"""
from typing import Dict, Tuple, List
import numpy as np
from scipy.linalg import eigh


def paulis() -> List[np.ndarray]:
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return [I, X, Y, Z]


def H_two_site(a: complex, b: complex, c: complex, d: complex) -> np.ndarray:
    """Return 4x4 two-site Hamiltonian H = a XX + b YY + c XY + d YX.

    Coefficients should be real for H to be Hermitian (but complex allowed).
    """
    I, X, Y, Z = paulis()
    H = a * np.kron(X, X) + b * np.kron(Y, Y) + c * np.kron(X, Y) + d * np.kron(Y, X)
    # enforce hermiticity (numerical safety)
    H = 0.5 * (H + H.conj().T)
    return H


def expand_H_on_pauli(H: np.ndarray) -> Dict[Tuple[int, int], complex]:
    r"""Compute h_{alpha,beta} = 1/4 Tr[(sigma^a \otimes sigma^b) H].

    Returns a dict keyed by (a,b) with a,b in {0,1,2,3} meaning I,X,Y,Z.
    """
    P = paulis()
    h = {}
    for a in range(4):
        for b in range(4):
            M = np.kron(P[a], P[b])
            h[(a, b)] = 0.25 * np.trace(M @ H)
    return h


def map_h_to_kitaev(h: Dict[Tuple[int, int], complex]) -> Dict[str, complex]:
    """Map Pauli coefficients h_{ab} to Kitaev chain parameters (t, Delta, mu).

    Uses the convention discussed in the docs:
      t = h_xx + h_yy + i(h_xy - h_yx)
      Delta = h_xx - h_yy - i(h_xy + h_yx)
      mu = 4 h_zz - 2(h_z0 + h_0z)

    Indices: 0=I,1=X,2=Y,3=Z
    """
    h_xx = h.get((1, 1), 0.0)
    h_yy = h.get((2, 2), 0.0)
    h_xy = h.get((1, 2), 0.0)
    h_yx = h.get((2, 1), 0.0)
    h_zz = h.get((3, 3), 0.0)
    h_z0 = h.get((3, 0), 0.0)
    h_0z = h.get((0, 3), 0.0)

    t = h_xx + h_yy + 1j * (h_xy - h_yx)
    Delta = h_xx - h_yy - 1j * (h_xy + h_yx)
    mu = 4.0 * h_zz - 2.0 * (h_z0 + h_0z)
    return {"t": t, "Delta": Delta, "mu": mu}


def build_bdg_chain(L: int, t: complex, Delta: complex, mu: complex, open_boundary: bool = True) -> np.ndarray:
    r"""Construct the BdG single-particle Hamiltonian (2L x 2L) for a spinless chain.

    Conventions (nearest neighbor):
      H = -mu sum_j c_j^\dagger c_j - t sum_j (c_j^\dagger c_{j+1} + h.c.)
          + Delta sum_j (c_j c_{j+1} + h.c.)

    We place matrix elements with factors 1/2 so that the standard BdG block
    H_BdG = [[h, Delta_mat], [-Delta_mat.conj(), -h.T]] is Hermitian.
    """
    h = np.zeros((L, L), dtype=complex)
    Delta_mat = np.zeros((L, L), dtype=complex)
    # on-site
    for i in range(L):
        h[i, i] = -mu / 2.0
    # nearest-neighbor
    for i in range(L - 1):
        h[i, i + 1] = -t / 2.0
        h[i + 1, i] = -t.conjugate() / 2.0
        # pairing: antisymmetric on indices
        Delta_mat[i, i + 1] = Delta / 2.0
        Delta_mat[i + 1, i] = -Delta / 2.0

    if not open_boundary:
        # periodic boundary
        h[0, L - 1] = -t / 2.0
        h[L - 1, 0] = -t.conjugate() / 2.0
        Delta_mat[0, L - 1] = Delta / 2.0
        Delta_mat[L - 1, 0] = -Delta / 2.0

    top = np.hstack([h, Delta_mat])
    bottom = np.hstack([-Delta_mat.conjugate(), -h.T])
    H_bdg = np.vstack([top, bottom])
    # hermitize numerically
    H_bdg = 0.5 * (H_bdg + H_bdg.conj().T)
    return H_bdg


def bdg_spectrum(H_bdg: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return eigenvalues (sorted) and eigenvectors of BdG matrix."""
    # Use eigh for Hermitian BdG
    E, V = eigh(H_bdg)
    return E, V


def find_zero_modes(E: np.ndarray, V: np.ndarray, L: int, tol: float = 1e-6) -> List[Dict]:
    """Locate near-zero BdG eigenvalues and return localization info.

    For each eigenvector v (length 2L) corresponding to |E|<tol, compute spatial
    density per site and fraction of weight near edges.
    """
    idx = np.where(np.abs(E) < tol)[0]
    modes = []
    for k in idx:
        v = V[:, k]
        u = v[:L]
        vpart = v[L:]
        density = np.abs(u) ** 2 + np.abs(vpart) ** 2
        left_frac = np.sum(density[: max(1, L // 10)])
        right_frac = np.sum(density[-max(1, L // 10) :])
        modes.append({"eig_index": int(k), "E": float(E[k]), "left_frac": float(left_frac), "right_frac": float(right_frac), "density": density})
    return modes


def compute_berry_along_path(Hs: List[np.ndarray], du: float) -> List[np.ndarray]:
    r"""Compute Berry connection matrices A(u) along a discrete path of Hamiltonians.

    Steps:
      - diagonalize H at each point (eigenvectors columns ordered by eigenvalue)
      - fix phases by maximizing overlap with previous step's eigenvectors (diagonal phase gauge)
      - approximate derivative and form A = i V^\dagger dV/du

    Returns list of A matrices (one per interval, length len(Hs)-1).
    """
    Vs = []
    Es = []
    for H in Hs:
        E, V = eigh(H)
        Es.append(E)
        Vs.append(V)

    As = []
    V_prev = Vs[0]
    for i in range(1, len(Vs)):
        V_curr = Vs[i]
        # overlap
        S = V_prev.conjugate().T @ V_curr
        # phase alignment (diagonal gauge) -- works for non-degenerate levels
        phases = np.angle(np.diag(S))
        V_curr = V_curr * np.exp(-1j * phases)[np.newaxis, :]
        dV = (V_curr - V_prev) / du
        A = 1j * V_prev.conjugate().T @ dV
        As.append(A)
        V_prev = V_curr
    return As


if __name__ == "__main__":
    # small smoke test
    H = H_two_site(1.0, 0.0, 0.0, 0.0)
    h = expand_H_on_pauli(H)
    params = map_h_to_kitaev(h)
    print("h pauli components (sample):", {k: float(v) for k, v in h.items() if k in [(1, 1), (2, 2), (3, 3)]})
    print("mapped params:", params)
