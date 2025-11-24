# demo_run.py
import numpy as np
from fibonacci_braid import find_braid_approximations

# target gates (2x2)
X = np.array([[0,1],[1,0]], dtype=complex)
H = (1/np.sqrt(2)) * np.array([[1,1],[1,-1]], dtype=complex)
T = np.array([[1,0],[0,np.exp(1j*np.pi/4)]], dtype=complex)  # T gate

# run searches (example: max length 8)
print("Searching for approximations up to braid length 8...")
best_H = find_braid_approximations(H, max_len=8, keep_top=6)
best_X = find_braid_approximations(X, max_len=8, keep_top=6)

def pretty_print(name, lst):
    print(f"\n=== Best approximations for {name} ===")
    for it in lst:
        print(f"distance={it['distance']:.6e}  length={it['length']:2d}  word={it['word']}")

pretty_print("Hadamard (H)", best_H)
pretty_print("Pauli-X (X)", best_X)
