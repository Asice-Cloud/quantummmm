#!/usr/bin/env python3
"""Detailed ABS diagnostics for selected local potentials h.

For each h: diagonalize H, find lowest two states, plot |u_j|^2 and |v_j|^2,
select Majorana pair by support, compute U_eff fidelity vs ideal B, and save results.
"""
import numpy as np
from scipy.linalg import expm, eigh, norm
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def kron(*mats):
    res = np.array([[1.0]])
    for m in mats:
        res = np.kron(res, m)
    return res


X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def pauli_op(pauli, site, L):
    ops = []
    for i in range(L):
        if i == site:
            ops.append({'x':X,'y':Y,'z':Z,'0':I2}[pauli])
        else:
            ops.append(I2)
    return kron(*ops)


def build_spin_ops(L):
    sx = [pauli_op('x', i, L) for i in range(L)]
    sy = [pauli_op('y', i, L) for i in range(L)]
    sz = [pauli_op('z', i, L) for i in range(L)]
    return sx, sy, sz


def jw_c_op(j, L):
    ops = []
    for k in range(L):
        if k < j:
            ops.append(Z)
        elif k == j:
            ops.append((X - 1j * Y) / 2)
        else:
            ops.append(I2)
    return kron(*ops)


def make_majoranas(L):
    c = [jw_c_op(j, L) for j in range(L)]
    gam = []
    for j in range(L):
        gam.append(c[j] + c[j].conj().T)
        gam.append(-1j*(c[j] - c[j].conj().T))
    return gam


def build_H_total(L, bond_coeffs, onsite_h):
    sx, sy, sz = build_spin_ops(L)
    H = np.zeros((2**L,2**L), dtype=complex)
    for j in range(L-1):
        c = bond_coeffs.get(f'{j}_{j+1}', {})
        if 'c_xx' in c:
            H += c['c_xx'] * (sx[j] @ sx[j+1])
        if 'c_yy' in c:
            H += c['c_yy'] * (sy[j] @ sy[j+1])
        if 'c_xy' in c:
            H += c['c_xy'] * (sx[j] @ sy[j+1])
        if 'c_yx' in c:
            H += c['c_yx'] * (sy[j] @ sx[j+1])
        if 'c_zz' in c:
            H += c['c_zz'] * (sz[j] @ sz[j+1])
    for j,h in onsite_h.items():
        H += h * sz[j]
    return H.real


def select_majorana_pair_by_support(gam, low_evecs):
    weights = []
    for j,g in enumerate(gam):
        Gp = low_evecs.conj().T @ (g @ low_evecs)
        w = norm(Gp)
        weights.append((w,j))
    weights.sort(reverse=True)
    return weights[0][1], weights[1][1], weights


def plot_ldos_uv(u, v, L, out_png, title=''):
    sites = np.arange(1, L+1)
    plt.figure(figsize=(6,4))
    plt.plot(sites, np.abs(u)**2, '-o', label='|u|^2')
    plt.plot(sites, np.abs(v)**2, '-s', label='|v|^2')
    plt.xlabel('Site')
    plt.ylabel('Probability')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def run_for_h(h, L=40):
    bond_coeffs = {}
    for j in range(L-1):
        bond_coeffs[f'{j}_{j+1}'] = {'c_xx':1.1, 'c_yy':1.0, 'c_xy':0.0, 'c_yx':0.0, 'c_zz':0.0}
    onsite_h = {L//2: h}
    H = build_H_total(L, bond_coeffs, onsite_h)
    evals, evecs = eigh(H)
    idx = np.argsort(np.abs(evals))
    low_idx = idx[:2]
    low_evecs = evecs[:, low_idx]
    P = low_evecs @ low_evecs.conj().T

    R = expm(1j * H)
    U_eff = P @ R @ P

    # majoranas built on L sites (same JW ordering)
    gam = make_majoranas(L)
    p, q, weights = select_majorana_pair_by_support(gam, low_evecs)
    B = expm((np.pi/4.0) * (gam[p] @ gam[q]))
    PB = P @ B @ P
    fid = norm(PB - U_eff) / (norm(PB) + 1e-12)

    # LDOS for lowest two states
    Lsites = L
    c_ops = [jw_c_op(j, L) for j in range(L)]

    diagnostics = {'h': h, 'fid': float(fid), 'energies': list(evals[idx[:4]]) , 'pair': (int(p), int(q))}

    outdir = Path('results')
    outdir.mkdir(exist_ok=True)

    for k,i in enumerate(low_idx):
        psi = evecs[:, i]
        u = psi[:L]
        v = psi[L:]
        outpng = outdir / f'ldos_h{h:.3f}_state{k}.png'
        plot_ldos_uv(u, v, L, outpng, title=f'h={h:.3f} state{k} E={evals[i]:.6f}')
        diagnostics[f'state{k}_u_sum'] = float(np.sum(np.abs(u)**2))
        diagnostics[f'state{k}_v_sum'] = float(np.sum(np.abs(v)**2))

    import json
    with open(outdir / f'abs_diag_h{h:.3f}.json','w') as f:
        json.dump(diagnostics, f, indent=2)

    return diagnostics


def main():
    hs = [0.1, 2.1]
    all_diag = {}
    for h in hs:
        print('Running diagnostics for h=', h)
        diag = run_for_h(h, L=40)
        print(' h=', h, ' fid=', diag['fid'], ' energies(lowest)=', diag['energies'][:2])
        all_diag[h] = diag
    print('\nSaved diagnostics and LDOS plots to results/')


if __name__ == '__main__':
    main()
