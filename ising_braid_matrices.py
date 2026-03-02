#!/usr/bin/env python3
import numpy as np

# Ising model R and F data
R_1 = np.exp(-1j * np.pi / 8)
R_psi = np.exp(3j * np.pi / 8)
D = np.diag([R_1, R_psi])
F = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)

# Basis order (a,b,c) with total T=1:
# 0: (1,1,1)  -> |00>
# 1: (psi,1,psi) -> |10>
# 2: (1,psi,psi) -> |01>
# 3: (psi,psi,1) -> |11>

basis = [(1,1,1), ('psi',1,'psi'), (1,'psi','psi'), ('psi','psi',1)]

def diag_R_over_label(label_list):
    # label_list elements are 1 or 'psi'
    mapping = {1: R_1, 'psi': R_psi}
    return np.diag([mapping[l] for l in label_list])

# Build B1, B3, B5 as diagonal in this basis
alist = [b[0] for b in basis]
blist = [b[1] for b in basis]
clist = [b[2] for b in basis]

B1 = diag_R_over_label(alist)
B3 = diag_R_over_label(blist)
B5 = diag_R_over_label(clist)

# Build B2 and B4 by embedding the 2x2 block M = F^{-1} D F
M = np.linalg.inv(F) @ D @ F

# Empirical embedding: mix indices (0,2) and (1,3)
def embed_M_pairs(pairs):
    mat = np.zeros((4,4), dtype=complex)
    for (i,j) in pairs:
        mat[i,i] = M[0,0]
        mat[i,j] = M[0,1]
        mat[j,i] = M[1,0]
        mat[j,j] = M[1,1]
    return mat

# pairs for B2: (0,2) and (1,3)
B2 = embed_M_pairs([(0,2),(1,3)])
# pairs for B4: act similarly on the right side (mixing pairs by symmetry)
B4 = embed_M_pairs([(0,1),(2,3)])  # alternative placement

def print_mat(name, M):
    np.set_printoptions(precision=4, suppress=True)
    print(f"{name} =\n", M)
    print()

def verify_braid(Bi, Bj):
    left = Bi @ Bj @ Bi
    right = Bj @ Bi @ Bj
    return np.allclose(left, right)

if __name__ == '__main__':
    print_mat('B1 (σ1)', B1)
    print_mat('B2 (σ2)', B2)
    print_mat('B3 (σ3)', B3)
    print_mat('B4 (σ4)', B4)
    print_mat('B5 (σ5)', B5)

    print('Verify braid relations:')
    print('σ1σ2σ1 = σ2σ1σ2 ? ->', verify_braid(B1, B2))
    print('σ2σ3σ2 = σ3σ2σ3 ? ->', verify_braid(B2, B3))
    print('σ3σ4σ3 = σ4σ3σ4 ? ->', verify_braid(B3, B4))
    print('σ4σ5σ4 = σ5σ4σ5 ? ->', verify_braid(B4, B5))
