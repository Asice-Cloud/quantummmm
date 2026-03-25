#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.linalg import eigh

def build_kitaev_chain_BdG(N, t, Delta, mu, periodic=False):
    # N sites, Majorana ordering not needed; return 2N x 2N BdG Hamiltonian
    # basis (c_1,...,c_N, c_1^
    H = np.zeros((2*N,2*N), dtype=complex)
    # h_ij terms and pairing Δ_ij
    for i in range(N):
        # onsite
        H[i,i] = -mu/2.0
        H[i+N,i+N] = mu/2.0
    # nearest neighbor
    for i in range(N-1 + (1 if periodic else 0)):
        j = (i+1) % N
        # hopping -t: contributes to particle block h_{i,j} = -t
        H[i,j] += -t
        H[j,i] += -t
        H[i+N,j+N] += t
        H[j+N,i+N] += t
        # pairing Δ: particle->hole block
        H[i,j+N] += Delta
        H[j,i+N] += -Delta
        H[i+N,j] += -Delta
        H[j+N,i] += Delta
    return H

def spectrum_open_chain(N, t, Delta, mu):
    H = build_kitaev_chain_BdG(N,t,Delta,mu,periodic=False)
    w, v = eigh(H)
    return np.sort(np.real_if_close(w))

def gap_kspace(t, Delta, mu):
    # sample k grid and compute min positive eigenvalue of H(k)
    ks = np.linspace(0, 2*np.pi, 801)
    minE = 1e9
    for k in ks:
        h_k = np.array([[-mu - 2*t*np.cos(k), -2j*Delta*np.sin(k)], [2j*Delta*np.sin(k), mu + 2*t*np.cos(k)]], dtype=complex)
        # eigenvalues
        w = np.linalg.eigvals(h_k)
        w = np.real_if_close(w)
        # take absolute values
        emin = np.min(np.abs(w))
        if emin < minE:
            minE = emin
    return float(minE)

if __name__ == '__main__':
    # parameter grid
    t = 1.0
    Delta = 0.5
    Ns = [20, 40]
    mus = np.linspace(-4.0, 4.0, 161)

    results = {'params':{'t':t,'Delta':Delta,'Ns':Ns}, 'data': []}

    gap_list = []
    topo_indicator = []
    for mu in mus:
        gapk = gap_kspace(t, Delta, mu)
        # Kitaev chain topological criterion (for real Delta): |mu| < 2|t|
        topo = 1 if abs(mu) < 2*abs(t) else 0
        gap_list.append(gapk)
        topo_indicator.append(topo)
        # open chain gap and zero-mode measure for representative N
        w20 = spectrum_open_chain(Ns[0], t, Delta, mu)
        # count near-zero modes magnitude < 1e-2
        nz = np.sum(np.abs(w20) < 1e-2)
        results['data'].append({'mu':float(mu),'gap_kspace':float(gapk),'topo':int(topo),'open_zero_count_N20':int(nz)})

    # save JSON
    with open('spectrum_topology_results.json','w') as f:
        json.dump(results, f, indent=2)

    # plot gap vs mu with topo region shading
    plt.figure(figsize=(6,4))
    plt.plot(mus, gap_list, label='bulk gap (k-space)')
    plt.fill_between(mus, 0, max(gap_list), where=np.array(topo_indicator)==1, color='orange', alpha=0.2, label='topological region |mu|<2|t|')
    plt.xlabel('mu')
    plt.ylabel('min gap')
    plt.title('Bulk gap vs mu (t=1, Delta=0.5)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('gap_vs_mu.png', dpi=200)
    plt.close()

    # plot zero-mode count for open chain N=20
    zero_counts = [d['open_zero_count_N20'] for d in results['data']]
    plt.figure(figsize=(6,4))
    plt.plot(mus, zero_counts)
    plt.xlabel('mu')
    plt.ylabel('zero-mode count (|E|<1e-2)')
    plt.title('Open chain (N=20) zero-mode count vs mu')
    plt.tight_layout()
    plt.savefig('zero_modes_N20_vs_mu.png', dpi=200)
    plt.close()

    print('Saved spectrum_topology_results.json and PNGs')
