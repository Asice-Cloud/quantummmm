#!/usr/bin/env python3
"""Compute lowest BdG energy vs chain length L for specified h values and plot decay.
Saves results/ldep_h0.1.png, results/ldep_h2.1.png and combined results/ldep_h_combo.png
"""
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_c_from_ybe(path):
    import re
    text = Path(path).read_text(encoding='utf-8')
    keys = ['c_xx','c_yy','c_xy','c_yx','c_zz','c_z0','c_0z','c_00']
    vals = {}
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
    mu = 4.0 * c_zz - 2.0*(c_z0 + c_0z)
    return t, Delta, mu


def build_bdg(L, t, Delta, mu_vec):
    H0 = np.zeros((L,L), dtype=complex)
    for i in range(L-1):
        H0[i, i+1] = -t
        H0[i+1, i] = -t
    for i in range(L):
        H0[i,i] += -mu_vec[i]
    D = np.zeros((L,L), dtype=complex)
    for i in range(L-1):
        D[i, i+1] = Delta
        D[i+1, i] = -Delta
    top = np.hstack([H0, D])
    bottom = np.hstack([-D.conj(), -H0.T])
    return np.vstack([top, bottom])


def run_for_h(h, Ls, t, Delta, mu0, out_prefix):
    Emins = []
    for L in Ls:
        mu_vec = np.zeros(L) + mu0
        mu_vec[L//2] += 2.0 * h
        H = build_bdg(L, t, Delta, mu_vec)
        evals, _ = eigh(H)
        emin = np.min(np.abs(evals))
        Emins.append(emin)
        print(f'h={h:.3f} L={L:3d} Emin={emin:.6e}')

    outdir = Path('results')
    outdir.mkdir(exist_ok=True)
    plt.figure(figsize=(6,4))
    plt.plot(Ls, Emins, 'o-')
    plt.yscale('log')
    plt.xlabel('L')
    plt.ylabel('Lowest |E| (log)')
    plt.title(f'Lowest energy vs L (h={h})')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outdir / f'{out_prefix}.png')
    return np.array(Emins)


def main():
    vals = parse_c_from_ybe(Path('ybe.md'))
    if not vals:
        vals = {'c_xx':1.1,'c_yy':1.0,'c_xy':0.0,'c_yx':0.0,'c_zz':0.0,'c_z0':0.0,'c_0z':0.0,'c_00':0.0}
    t, Delta, mu0 = compute_p(vals)
    Ls = list(range(20, 201, 20))
    hs = [0.1, 2.1]
    results = {}
    for h in hs:
        Emins = run_for_h(h, Ls, t, Delta, mu0, f'ldep_h{h}')
        results[h] = (Ls, Emins)

    # combined plot
    outdir = Path('results')
    plt.figure(figsize=(6,4))
    for h in hs:
        Ls, Emins = results[h]
        plt.plot(Ls, Emins, 'o-', label=f'h={h}')
    plt.yscale('log')
    plt.xlabel('L')
    plt.ylabel('Lowest |E| (log)')
    plt.title('Lowest energy vs L')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outdir / 'ldep_h_combo.png')
    print('Saved combined plot to results/ldep_h_combo.png')


if __name__ == '__main__':
    main()
