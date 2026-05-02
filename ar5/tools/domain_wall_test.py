#!/usr/bin/env python3
"""Domain-wall test: use top-2 candidates to build left/right parameter regions and analyze low-energy modes.

Saves JSON report and LDOS/Majorana plots to results/.
"""
import os
import json
import numpy as np
import importlib.util

here = os.path.dirname(__file__)

def parse_complex(x):
    try:
        return complex(x)
    except Exception:
        return complex(x.replace(' ', ''))


def build_site_dependent_bdG(N, t_left, Delta_left, mu_left, t_right, Delta_right, mu_right):
    # site-dependent mu per site, bond-dependent t and Delta per bond
    mu = np.zeros(N, dtype=float)
    for i in range(N):
        mu[i] = float(np.real(mu_left) if i < N//2 else np.real(mu_right))
    # bond arrays length N-1
    t_bond = np.zeros(N-1, dtype=complex)
    Delta_bond = np.zeros(N-1, dtype=complex)
    for i in range(N-1):
        if i < N//2 - 1:
            t_bond[i] = t_left
            Delta_bond[i] = Delta_left
        elif i > N//2:
            t_bond[i] = t_right
            Delta_bond[i] = Delta_right
        else:
            # interface bond: average
            t_bond[i] = 0.5*(t_left + t_right)
            Delta_bond[i] = 0.5*(Delta_left + Delta_right)

    A = np.zeros((N,N), dtype=complex)
    for i in range(N):
        A[i,i] = -mu[i]
    for i in range(N-1):
        A[i,i+1] = -t_bond[i]
        A[i+1,i] = -np.conj(t_bond[i])
    B = np.zeros((N,N), dtype=complex)
    for i in range(N-1):
        B[i,i+1] = Delta_bond[i]
        B[i+1,i] = -Delta_bond[i]
    H = np.block([[A, B], [-B.conj(), -A.conj()]])
    return H


def maj_sim(vec):
    half = vec.shape[0]//2
    u = vec[:half]
    v = vec[half:]
    num = np.linalg.norm(u - np.conj(v))
    den = np.linalg.norm(u) + np.linalg.norm(v)
    return float(1.0 - (num / den)) if den>0 else 0.0


def psiA_B(vec):
    half = vec.shape[0]//2
    u = vec[:half]
    v = vec[half:]
    A = u + np.conj(v)
    B = u - np.conj(v)
    return A, B


def overlap_norm(a, b):
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na==0 or nb==0:
        return 0.0
    return float(np.sum(np.abs(a)*np.abs(b)) / (na*nb))


def run(N=240):
    os.makedirs('results', exist_ok=True)
    data = json.load(open('results/scan_all_mixtures_validated.json'))
    data_sorted = sorted(data, key=lambda d: d.get('min_abs', 1e9))
    left = data_sorted[0]
    right = data_sorted[1]
    tL = parse_complex(left['t'])
    dL = parse_complex(left['Delta'])
    muL = parse_complex(left['mu'])
    tR = parse_complex(right['t'])
    dR = parse_complex(right['Delta'])
    muR = parse_complex(right['mu'])

    H = build_site_dependent_bdG(N, tL, dL, muL, tR, dR, muR)
    vals, vecs = np.linalg.eigh(H)
    idx = np.argsort(np.abs(vals))
    modes = []
    for k in range(4):
        i = idx[k]
        v = vecs[:, i]
        E = float(vals[i])
        ms = maj_sim(v)
        Acomp, Bcomp = psiA_B(v)
        rho = np.abs(v[:N])**2 + np.abs(v[N:])**2
        peak = int(np.argmax(rho))
        modes.append({'k': k, 'E': E, 'maj_sim': ms, 'peak': peak})

    # compute pair overlaps for lowest two modes
    v1 = vecs[:, idx[0]]
    v2 = vecs[:, idx[1]]
    A1, B1 = psiA_B(v1)
    A2, B2 = psiA_B(v2)
    S_A1A2 = overlap_norm(A1, A2)
    S_B1B2 = overlap_norm(B1, B2)
    S_A1B1 = overlap_norm(A1, B1)

    report = {
        'left_combo': left.get('combo'), 'right_combo': right.get('combo'),
        'left_params': {'t': str(tL), 'Delta': str(dL), 'mu': str(muL)},
        'right_params': {'t': str(tR), 'Delta': str(dR), 'mu': str(muR)},
        'modes': modes,
        'overlaps': {'S_A1A2': S_A1A2, 'S_B1B2': S_B1B2, 'S_A1B1': S_A1B1}
    }
    json.dump(report, open('results/domain_wall_report.json','w'), indent=2)

    # save LDOS and majorana plots
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # combined LDOS for two lowest modes
    l1 = np.abs(v1[:N])**2 + np.abs(v1[N:])**2
    l2 = np.abs(v2[:N])**2 + np.abs(v2[N:])**2
    s = l1 + l2
    x = np.arange(N)
    plt.figure(figsize=(8,3)); plt.plot(x, l1, label='mode1'); plt.plot(x, l2, label='mode2'); plt.plot(x, s, label='sum', linewidth=2); plt.axvline(N//2, color='k', linestyle='--'); plt.legend(); plt.title('domain wall combined LDOS'); plt.tight_layout(); plt.savefig('results/domain_wall_ldos.png'); plt.close()

    # majorana components for mode1
    rhoA1 = np.abs(A1)**2; rhoB1 = np.abs(B1)**2
    plt.figure(figsize=(8,3)); plt.plot(x, rhoA1, label='rhoA1'); plt.plot(x, rhoB1, label='rhoB1'); plt.axvline(N//2, color='k', linestyle='--'); plt.legend(); plt.title('domain wall Majorana components (mode1)'); plt.tight_layout(); plt.savefig('results/domain_wall_majorana_mode1.png'); plt.close()

    print('Wrote results/domain_wall_report.json, domain_wall_ldos.png, domain_wall_majorana_mode1.png')


if __name__ == '__main__':
    run()
