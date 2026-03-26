import numpy as np

def build_micro_BdG(N, a=0.0, b=0.0, c=0.0, d=0.0, mu0=0.0):
    """Construct a simple microscopic BdG for an N-site chain.
    Microscopic choices (default example):
      - normal nearest-neighbor hopping = b
      - p-wave pairing nearest-neighbor = c
      - onsite chemical potential = mu0 + 4*d + a
    Returns 2N x 2N BdG Hamiltonian (particle/hole blocks).
    """
    t_hop = b
    Delta = c
    mu = mu0 + 4.0*d + a
    H = np.zeros((2*N, 2*N), dtype=complex)
    # particle block h
    for i in range(N):
        H[i,i] = -mu/2.0
        H[i+N, i+N] = mu/2.0
    for i in range(N-1):
        j = i+1
        H[i,j] += -t_hop
        H[j,i] += -t_hop
        H[i+N,j+N] += t_hop
        H[j+N,i+N] += t_hop
        # pairing
        H[i, j+N] += Delta
        H[j, i+N] += -Delta
        H[i+N, j] += -Delta
        H[j+N, i] += Delta
    return H

# extract effective parameters by averaging the corresponding matrix elements
def extract_effective_params_from_micro(H):
    # H is 2N x 2N
    N2 = H.shape[0]
    N = N2//2
    # particle block h = H[:N,:N]
    h = H[:N,:N]
    # pairing block Delta = H[:N, N:]
    D = H[:N, N:]
    # effective t: negative of average nearest-neighbor h[i,i+1]
    t_list = []
    delta_list = []
    for i in range(N-1):
        t_list.append(-np.real(h[i, i+1]))
        delta_list.append(np.real(D[i, i+1]))
    t_eff = float(np.mean(t_list)) if len(t_list)>0 else 0.0
    delta_eff = float(np.mean(delta_list)) if len(delta_list)>0 else 0.0
    # effective mu: -2 * average diagonal of h
    mu_eff = float(-2.0 * np.mean(np.real(np.diag(h))))
    return {'t_eff':t_eff, 'Delta_eff':delta_eff, 'mu_eff':mu_eff}

if __name__ == '__main__':
    H = build_micro_BdG(4, a=0.1, b=0.3, c=0.2, d=0.01, mu0=0.05)
    params = extract_effective_params_from_micro(H)
    print('effective params:', params)
