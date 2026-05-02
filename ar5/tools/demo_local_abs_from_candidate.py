#!/usr/bin/env python3
"""Demo: use an R(u)-mapped candidate and add a local mu well to create a single localized ABS.

Saves results/demo_local_abs_* images and prints energies, maj_sim and peak info.
"""
import os
import json
import numpy as np
import importlib.util

here = os.path.dirname(__file__)
# load build_bdg_from_params from scan_all_mixtures if available
spec = importlib.util.spec_from_file_location('scan_all_mixtures', os.path.join(here, 'scan_all_mixtures.py'))
scan_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scan_mod)
build_bdg_from_params = scan_mod.build_bdg_from_params

def maj_sim_from_vec(vec):
    half = vec.shape[0]//2
    u = vec[:half]
    v = vec[half:]
    num = np.linalg.norm(u - np.conj(v))
    den = np.linalg.norm(u) + np.linalg.norm(v)
    return 1.0 - (num / den) if den>0 else 0.0


def ldos_from_vec(vec):
    half = vec.shape[0]//2
    u = vec[:half]
    v = vec[half:]
    return np.abs(u)**2 + np.abs(v)**2


def run_demo(candidate_idx=0, N=120, delta_mu=-2.5, j0=None):
    os.makedirs('results', exist_ok=True)
    data = json.load(open('results/scan_all_mixtures_validated.json'))
    cand = data[candidate_idx]
    # extract params
    t = complex(cand['t'])
    Delta = complex(cand['Delta'])
    # `cand['mu']` may be a complex serialized as a string; take real part
    try:
        mu0 = float(cand['mu'])
    except Exception:
        mu0 = float(np.real(complex(cand['mu'])))

    if j0 is None:
        j0 = N//2

    # build site-dependent mu
    mu_array = np.ones(N, dtype=float) * mu0
    mu_array[j0] += delta_mu

    # build BdG with site-dependent mu by assembling block matrices
    A = np.zeros((N,N), dtype=complex)
    for i in range(N):
        A[i,i] = -mu_array[i]
    for i in range(N-1):
        A[i,i+1] = -t
        A[i+1,i] = -np.conj(t)
    B = np.zeros((N,N), dtype=complex)
    for i in range(N-1):
        B[i,i+1] = Delta
        B[i+1,i] = -Delta
    H = np.block([[A, B], [-B.conj(), -A.conj()]])

    vals, vecs = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))
    modes = []
    for k in range(2):
        v = vecs[:, idx[k]]
        E = float(vals[idx[k]])
        ldos = ldos_from_vec(v)
        ms = maj_sim_from_vec(v)
        peak_idx = int(np.argmax(ldos))
        modes.append({'E': E, 'maj_sim': float(ms), 'peak_idx': peak_idx, 'max_ldos': float(np.max(ldos))})

    # save combined LDOS plot
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    v1 = vecs[:, idx[0]]
    v2 = vecs[:, idx[1]]
    l1 = ldos_from_vec(v1)
    l2 = ldos_from_vec(v2)
    s = l1 + l2

    x = np.arange(N)
    plt.figure(figsize=(8,3)); plt.plot(x, l1, label='mode1'); plt.plot(x, l2, label='mode2'); plt.plot(x, s, label='sum', linewidth=2); plt.axvline(j0, color='k', linestyle='--'); plt.legend(); plt.title(f"demo_local_abs candidate#{candidate_idx} j0={j0} dmu={delta_mu}"); plt.xlabel('site'); plt.ylabel('LDOS'); plt.tight_layout(); plt.savefig('results/demo_local_abs_sum.png'); plt.close()

    # save majorana components ρA, ρB for mode1
    half = H.shape[0]//2
    u = v1[:half]; vcomp = v1[half:]
    rhoA = np.abs(u + np.conj(vcomp))**2
    rhoB = np.abs(u - np.conj(vcomp))**2
    plt.figure(figsize=(8,3)); plt.plot(x, rhoA, label='rhoA'); plt.plot(x, rhoB, label='rhoB'); plt.axvline(j0, color='k', linestyle='--'); plt.legend(); plt.title('Majorana components (mode1)'); plt.tight_layout(); plt.savefig('results/demo_local_abs_majorana_mode1.png'); plt.close()

    # print summary
    print('Candidate idx:', candidate_idx)
    print('params: t,Delta,mu0 =', t, Delta, mu0)
    for i,m in enumerate(modes, start=1):
        print(f"mode{i}: E={m['E']:.3e}, maj_sim={m['maj_sim']:.3f}, peak_idx={m['peak_idx']}, max_ldos={m['max_ldos']:.4f}")
    print('Saved plots: results/demo_local_abs_sum.png, results/demo_local_abs_majorana_mode1.png')


if __name__ == '__main__':
    run_demo()
