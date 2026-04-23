#!/usr/bin/env python3
"""Refined BdG scan around a chosen h center; save CSV and LDOS for best h.
"""
import argparse
import numpy as np
from scipy.linalg import eigh
from pathlib import Path
import matplotlib.pyplot as plt


def build_bdg(L, t, Delta, mu_vec):
    H0 = np.zeros((L, L), dtype=complex)
    for i in range(L-1):
        H0[i, i+1] = -t
        H0[i+1, i] = -t
    for i in range(L):
        H0[i, i] += -mu_vec[i]
    D = np.zeros((L, L), dtype=complex)
    for i in range(L-1):
        D[i, i+1] = Delta
        D[i+1, i] = -Delta
    top = np.hstack([H0, D])
    bottom = np.hstack([-D.conj(), -H0.T])
    return np.vstack([top, bottom])


def parse_ybe_for_params():
    import re
    text = Path('ybe.md').read_text(encoding='utf-8')
    def getc(key):
        m = re.search(r'%s\s*[:=]\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)' % key, text)
        return float(m.group(1)) if m else 0.0
    c = {k: getc(k) for k in ['c_xx','c_yy','c_xy','c_yx','c_zz','c_z0','c_0z','c_00']}
    t = c['c_xx'] + c['c_yy'] + 1j*(c['c_xy'] - c['c_yx'])
    Delta = c['c_xx'] - c['c_yy'] - 1j*(c['c_xy'] + c['c_yx'])
    mu0 = 4.0 * c['c_zz'] - 2.0*(c['c_z0'] + c['c_0z'])
    return t, Delta, mu0


def compute_ldos(evec, L):
    # evec is 2L vector [u; v]
    u = evec[:L]
    v = evec[L:]
    ldos = np.abs(u)**2 + np.abs(v)**2
    return ldos


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--hcenter', type=float, required=True)
    p.add_argument('--delta', type=float, default=0.1)
    p.add_argument('--n', type=int, default=101)
    p.add_argument('--L', type=int, default=120)
    args = p.parse_args()

    t, Delta, mu0 = parse_ybe_for_params()
    L = args.L
    mid = L//2
    hs = np.linspace(args.hcenter - args.delta, args.hcenter + args.delta, args.n)
    out = Path('results')
    out.mkdir(exist_ok=True)
    rows = []
    for h in hs:
        mu_vec = np.zeros(L, dtype=float) + mu0
        mu_vec[mid] += 2.0 * h
        Hbdg = build_bdg(L, t, Delta, mu_vec)
        evals, evecs = eigh(Hbdg)
        idx = np.argsort(np.abs(evals))
        emin = np.min(np.abs(evals))
        # pick eigenvector with smallest |E|
        ev = evecs[:, idx[0]]
        ldos = compute_ldos(ev, L)
        rows.append((h, emin, ldos))
        print(f'h={h:.6f} Emin={emin:.6e}')

    # find best
    emins = np.array([r[1] for r in rows])
    imin = np.argmin(emins)
    hbest, ebest, ldos_best = rows[imin]
    print('\nBest h=', hbest, ' Emin=', ebest)

    # save CSV (without dumping full LDOS per row), and save best LDOS as png
    csvp = out / f'refine_scan_h_{args.hcenter:.3f}_delta_{args.delta:.3f}_L{L}.csv'
    with open(csvp, 'w') as f:
        f.write('h,Emin\n')
        for h, e, _ in rows:
            f.write(f'{h},{e}\n')
    print('Saved summary CSV to', csvp)

    # save best LDOS plot
    plt.figure(figsize=(8,3))
    plt.plot(np.arange(L), ldos_best, '-o', markersize=2)
    plt.xlabel('site')
    plt.ylabel('LDOS (best state)')
    plt.title(f'LDOS best h={hbest:.6f} Emin={ebest:.3e}')
    pngp = out / f'ldos_best_h_{hbest:.6f}_L{L}.png'
    plt.tight_layout()
    plt.savefig(pngp, dpi=150)
    print('Saved LDOS plot to', pngp)


if __name__ == '__main__':
    main()
