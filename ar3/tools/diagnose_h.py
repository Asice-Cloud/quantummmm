#!/usr/bin/env python3
"""Detailed BdG diagnostics for specified local potential values h.
Generates spectrum and LDOS (particle+hole) for lowest states and saves images.
"""
import numpy as np
from pathlib import Path
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
    U = 4.0 * c_zz
    mu = 4.0 * c_zz - 2.0*(c_z0 + c_0z)
    C_perbond = c_zz - (c_z0 + c_0z) + c_00
    return dict(t=t, Delta=Delta, U=U, mu=mu, C_perbond=C_perbond)


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


def plot_ldos_uv(u, v, outpath, title=''):
    L = len(u)
    x = np.arange(1, L+1)
    plt.figure(figsize=(6,3))
    plt.plot(x, np.abs(u)**2, label='|u|^2 (particle)')
    plt.plot(x, np.abs(v)**2, label='|v|^2 (hole)')
    plt.xlabel('Site')
    plt.ylabel('Probability')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def main():
    vals = parse_c_from_ybe(Path('ybe.md'))
    if not vals:
        vals = {'c_xx':1.1,'c_yy':1.0,'c_xy':0.0,'c_yx':0.0,'c_zz':0.0,'c_z0':0.0,'c_0z':0.0,'c_00':0.0}
    p = compute_p(vals)
    t = p['t']; Delta = p['Delta']; mu0 = p['mu']

    L = 40
    hs = [0.1, 2.1]
    outdir = Path('results')
    outdir.mkdir(exist_ok=True)

    for h in hs:
        mu_vec = np.zeros(L, dtype=float) + mu0
        mu_vec[L//2] += 2.0 * h
        Hbdg = build_bdg(L, t, Delta, mu_vec)
        evals, evecs = eigh(Hbdg)
        idx = np.argsort(np.abs(evals))
        print(f'h={h}: lowest energies (abs):', np.round(np.abs(evals[idx[:4]]),6))

        for k in range(2):
            i = idx[k]
            vec = evecs[:, i]
            u = vec[:L]; v = vec[L:]
            plot_ldos_uv(u, v, outdir / f'ldos_h{h}_state{k}.png', title=f'h={h} state{k} E={evals[i]:.6f}')
            print(f'  state{k} energy = {evals[i]:.6f}, sum|u|^2={np.sum(np.abs(u)**2):.6f}, sum|v|^2={np.sum(np.abs(v)**2):.6f}')

    print('Saved LDOS plots to results/')


if __name__ == '__main__':
    main()
