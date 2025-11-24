"""
fibonacci_braid.py
A small module to approximate 2x2 target unitaries with Fibonacci anyon braids.

Provides:
 - build_generators(): returns dict of generators {"1","2","1-","2-"} -> 2x2 complex matrices
 - braid_word_to_matrix(word, gens)
 - phase_invariant_distance(U,V)
 - find_braid_approximations(target, max_len=8, keep_top=10)

Usage example:
    from fibonacci_braid import find_braid_approximations
    best = find_braid_approximations(target=H, max_len=10, keep_top=5)
"""
import numpy as np
from collections import deque

# --- constants and F, R definitions ---
phi = (1 + 5**0.5) / 2
phi_inv = 1 / phi

# F-matrix for tau,tau,tau -> tau block (2x2)
F = np.array([[phi_inv, phi**(-0.5)],
              [phi**(-0.5), -phi_inv]], dtype=complex)

# R diag: R^{tau tau}_1 and R^{tau tau}_tau
R1 = np.exp(-4j * np.pi / 5)    # fusion channel = 1
Rtau = np.exp(3j * np.pi / 5)   # fusion channel = tau
Rdiag = np.diag([R1, Rtau])

# Braid generators in the 3-anyon fusion basis
sigma1 = Rdiag.copy()                  # exchange of first two anyons
sigma2 = np.linalg.inv(F) @ Rdiag @ F # sigma2 = F^{-1} R F

GENS = {
    "1": sigma1,
    "2": sigma2,
    "1-": np.linalg.inv(sigma1),
    "2-": np.linalg.inv(sigma2),
}

# --- utilities --- 
def braid_word_to_matrix(word, gens=GENS):
    """
    word: list of tokens, e.g. ['1','2','1-']
    gens: generator dictionary
    returns the 2x2 matrix corresponding to that braid (left-multiplicative order)
    """
    U = np.eye(2, dtype=complex)
    for g in word:
        U = gens[g] @ U
    return U

def phase_invariant_distance(U, V):
    """
    Frobenius norm minimised over a global phase:
    min_phi || U - e^{i phi} V ||_F
    optimal phi=arg(trace(V^\daggerU))
    returns the Frobenius norm (a nonnegative float)
    """
    tr = np.trace(V.conj().T @ U)
    phi_opt = np.angle(tr)
    diff = U - np.exp(1j * phi_opt) * V
    return np.linalg.norm(diff, ord='fro')

def find_braid_approximations(target, max_len=8, keep_top=10):
    """
    BFS exhaustive search up to length max_len over alphabet ["1","2","1-","2-"].
    Prunes by coarse-graining matrices (rounding) to avoid repeat visits.
    Returns a list of dicts: {distance, length, word (string), matrix (nested list)}
    NOTE: complexity ~ O(4^L) (but pruning/coarse-hash reduces duplicates).
    """
    alphabet = ["1","2","1-","2-"]
    inverse = {"1":"1-","1-":"1","2":"2-","2-":"2"}
    seen = {}   # coarse hash -> (word, matrix)
    best = []
    q = deque([[]])   # start with empty word (identity)
    while q:
        w = q.popleft()
        U = braid_word_to_matrix(w)
        # pruning key: round to reduce floating noise (tunable)
        key = tuple(np.round(np.concatenate([U.real.flatten(), U.imag.flatten()]), 12))
        if key in seen:
            # skip if we've seen this coarse matrix before
            pass
        else:
            seen[key] = (w, U)
            d = phase_invariant_distance(U, target)
            best.append((float(d), w, U))
            if len(w) < max_len:
                for g in alphabet:
                    # avoid immediate inverse cancellation like "... 1 1-" 
                    if len(w) > 0 and inverse[g] == w[-1]:
                        continue
                    q.append(w + [g])
    # sort by distance and return top keep_top
    best_sorted = sorted(best, key=lambda x: x[0])
    results = []
    for d, w, U in best_sorted[:keep_top]:
        results.append({
            "distance": float(d),
            "length": len(w),
            "word": " ".join(w) if w else "(identity)",
            "matrix": U.tolist()
        })
    return results

# expose some symbols for convenience
__all__ = ["F","Rdiag","sigma1","sigma2","GENS","braid_word_to_matrix","phase_invariant_distance","find_braid_approximations"]
