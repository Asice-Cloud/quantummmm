#!/usr/bin/env python3
"""Compute lowest positive BdG energy vs chain length L and plot decay.

Saves results/ldep.png and prints the numeric table.
"""
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path


def parse_c_from_ybe(path):
    text = Path(path).read_text(encoding='utf-8')
    keys = ['c_xx','c_yy','c_xy','c_yx','c_zz','c_z0','c_0z','c_00']
    vals = {}
    import re
    for k in keys:
        m = re.search(r'%s\s*[:=]\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)' % k, text)
        if m:
            vals[k] = float(m.group(1))
    return vals


def compute_p(c):
    get = lambda k: c.get(k, 0.0)
    c_xx = get('c_xx')
    c_yy = get('c_yy')
    c_xy = get('c_xy')
    c_yx = get('c_yx')
    c_zz = get('c_zz')
    c_z0 = get('c_z0')
    c_0z = get('c_0z')
    c_00 = get('c_00')
    t = c_xx + c_yy + 1j*(c_xy - c_yx)
    Delta = c_xx - c_yy - 1j*(c_xy + c_yx)
    U = 4.0 * c_zz
    mu = 4.0 * c_zz - 2.0*(c_z0 + c_0z)
    C_perbond = c_zz - (c_z0 + c_0z) + c_00
    return dict(t=t, Delta=Delta, U=U, mu=mu, C_perbond=C_perbond)


def build_bdg(L, t, Delta, mu):
    H0 = np.zeros((L,L), dtype=complex)
    for i in range(L-1):
        H0[i, i+1] = -t
        H0[i+1, i] = -t
    H0 += -mu * np.eye(L)
    D = np.zeros((L,L), dtype=complex)
    for i in range(L-1):
        D[i, i+1] = Delta
        D[i+1, i] = -Delta
    top = np.hstack([H0, D])
    bottom = np.hstack([-D.conj(), -H0.T])
    Hbdg = np.vstack([top, bottom])
    return Hbdg


def lowest_energy_for_L(L, t, Delta, mu):
    H = build_bdg(L, t, Delta, mu)
    evals, _ = eigh(H)
    energies = np.sort(np.abs(evals))
    # smallest non-zero (positive) energy
    return energies[0]


def main():
    c = parse_c_from_ybe(Path('ybe.md'))
    if not c:
        c = {'c_xx':1.1,'c_yy':1.0,'c_xy':0.0,'c_yx':0.0,'c_zz':0.0,'c_z0':0.0,'c_0z':0.0,'c_00':0.0}
    p = compute_p(c)
    t = p['t']; Delta = p['Delta']; mu = p['mu']

    Ls = list(range(20, 201, 20))
    Emins = []
    for L in Ls:
        E = lowest_energy_for_L(L, t, Delta, mu)
        Emins.append(E)
        print(f'L={L:3d}  Emin={E:.6e}')

    outdir = Path('results')
    outdir.mkdir(exist_ok=True)

    plt.figure(figsize=(6,4))
    plt.plot(Ls, Emins, 'o-')
    plt.yscale('log')
    plt.xlabel('Chain length L')
    plt.ylabel('Lowest |E| (log scale)')
    plt.title('Finite-size scaling of lowest BdG energy')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outdir / 'ldep.png')
    print('\nSaved plot to', outdir / 'ldep.png')


if __name__ == '__main__':
    main()
