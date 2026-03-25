#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.linalg import eigh

# reuse builder from majorana_braid_check but here build general chain

def build_kitaev_chain_BdG(N, t, Delta, mu, periodic=False):
    H = np.zeros((2*N,2*N), dtype=complex)
    for i in range(N):
        H[i,i] = -mu/2.0
        H[i+N,i+N] = mu/2.0
    for i in range(N-1 + (1 if periodic else 0)):
        j = (i+1) % N
        H[i,j] += -t
        H[j,i] += -t
        H[i+N,j+N] += t
        H[j+N,i+N] += t
        H[i,j+N] += Delta
        H[j,i+N] += -Delta
        H[i+N,j] += -Delta
        H[j+N,i] += Delta
    return H

# winding number via phase winding of q(k) = -mu - 2t cos k + 2i Delta sin k
def compute_winding(t, Delta, mu, nk=2001):
    ks = np.linspace(0, 2*np.pi, nk)
    q = -mu - 2.0*t*np.cos(ks) + 2.0j*Delta*np.sin(ks)
    phases = np.angle(q)
    # unwrap and compute net change
    phases_un = np.unwrap(phases)
    winding = int(round((phases_un[-1] - phases_un[0]) / (2*np.pi)))
    return winding, ks, q, phases

# Pfaffian sign change test (k=0, pi) using determinant of q(k) (2x2 off-diagonal block)
# For Kitaev chain q(k) is scalar here, so sign of Re[q(0)]*Re[q(pi)] can be used

def pfaffian_sign_test(t, Delta, mu):
    q0 = -mu - 2.0*t*np.cos(0) + 2.0j*Delta*np.sin(0)
    qpi = -mu - 2.0*t*np.cos(np.pi) + 2.0j*Delta*np.sin(np.pi)
    # use sign of real part
    s = np.sign(np.real(q0) * np.real(qpi))
    return int(-1 if s<0 else 1)

# zero-mode scaling
def zero_mode_scaling(t, Delta, mu, Ns=[40,60,80,120,200]):
    Emins = []
    for N in Ns:
        H = build_kitaev_chain_BdG(N,t,Delta,mu,periodic=False)
        w, v = eigh(H)
        w = np.sort(np.real_if_close(w))
        # take smallest positive eigenvalue
        pos = w[w>1e-12]
        if len(pos)==0:
            Emin = 0.0
        else:
            Emin = pos[0]
        Emins.append(float(Emin))
    # fit ln(Emin) = A - N/xi  (exclude zeros)
    Ns_arr = np.array(Ns)
    E_arr = np.array(Emins)
    mask = E_arr>0
    if np.sum(mask) < 2:
        xi = None
        fit = None
    else:
        y = np.log(E_arr[mask])
        x = Ns_arr[mask]
        A, B = np.polyfit(x, y, 1)
        # slope = B? Actually polyfit returns y = m x + c, where m = slope
        m = A
        c = B
        # here y = m x + c = ln E = m N + c => m = -1/xi
        xi = -1.0/m if m!=0 else None
        fit = {'slope':float(m),'intercept':float(c),'xi':float(xi)}
    return {'Ns':Ns, 'Emins':Emins, 'fit':fit}

if __name__ == '__main__':
    t = 1.0
    Delta = 0.5
    mus = np.linspace(-3.0,3.0,61)

    results = {'winding':[], 'pfaffian':[], 'scaling':{}}
    for mu in mus:
        wnum, ks, q, phases = compute_winding(t, Delta, mu, nk=2001)
        ptest = pfaffian_sign_test(t, Delta, mu)
        results['winding'].append({'mu':float(mu),'winding':int(wnum)})
        results['pfaffian'].append({'mu':float(mu),'pf_sign':int(ptest)})
    # save winding map
    with open('winding_pfaffian_results.json','w') as f:
        json.dump(results, f, indent=2)

    # pick representative mu in topological and trivial for scaling
    mu_topo = 0.0
    mu_triv = 3.0
    scaling_topo = zero_mode_scaling(t,Delta,mu_topo)
    scaling_triv = zero_mode_scaling(t,Delta,mu_triv)
    with open('zero_mode_scaling.json','w') as f:
        json.dump({'topo_mu':mu_topo,'triv_mu':mu_triv,'topo':scaling_topo,'triv':scaling_triv}, f, indent=2)

    # plot winding number vs mu
    mus_plot = [d['mu'] for d in results['winding']]
    wvals = [d['winding'] for d in results['winding']]
    plt.figure(figsize=(6,3))
    plt.plot(mus_plot, wvals, '-o', ms=3)
    plt.xlabel('mu')
    plt.ylabel('winding number')
    plt.title('Winding number vs mu')
    plt.tight_layout()
    plt.savefig('winding_vs_mu.png', dpi=200)
    plt.close()

    # plot scaling fits
    import math
    def plot_scaling(scaling, name):
        Ns = scaling['Ns']
        E = scaling['Emins']
        plt.figure(figsize=(6,4))
        plt.semilogy(Ns, E, 'o-')
        if scaling['fit'] is not None:
            xs = np.linspace(min(Ns), max(Ns), 200)
            lnE = scaling['fit']['slope']*xs + scaling['fit']['intercept']
            plt.semilogy(xs, np.exp(lnE), 'r--', label=f"fit xi={scaling['fit']['xi']:.2f}")
            plt.legend()
        plt.xlabel('N')
        plt.ylabel('E_min (open chain)')
        plt.title(name)
        plt.tight_layout()
        plt.savefig(f"zero_mode_scaling_{name.replace(' ','_')}.png", dpi=200)
        plt.close()

    plot_scaling(scaling_topo, 'topological_mu_0.0')
    plot_scaling(scaling_triv, 'trivial_mu_3.0')

    print('Saved winding_pfaffian_results.json, zero_mode_scaling.json and PNGs')
