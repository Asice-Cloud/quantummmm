#!/usr/bin/env python3
"""Build finite Kitaev BdG from (t,Delta,mu) parsed from ybe.md and plot spectrum + LDOS.

Saves `results/bdg_spectrum.png` and `results/ldos_lowest.png`.
"""
import re
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_c_from_ybe(path):
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


def build_bdg(L, t, Delta, mu):
    # H0: hopping and chemical potential (spinless)
    H0 = np.zeros((L,L), dtype=complex)
    for i in range(L-1):
        H0[i, i+1] = -t
        H0[i+1, i] = -t
    H0 += -mu * np.eye(L)
    # pairing for nearest-neighbor p-wave: Delta c_j c_{j+1} + h.c.
    D = np.zeros((L,L), dtype=complex)
    for i in range(L-1):
        D[i, i+1] = Delta
        D[i+1, i] = -Delta
    # BdG matrix
    top = np.hstack([H0, D])
    bottom = np.hstack([-D.conj(), -H0.T])
    Hbdg = np.vstack([top, bottom])
    return Hbdg


def plot_spectrum(energies, out):
    plt.figure(figsize=(6,4))
    plt.plot(np.arange(len(energies)), energies, 'o-')
    plt.axhline(0, color='k', lw=0.6)
    plt.xlabel('Index')
    plt.ylabel('Energy')
    plt.title('BdG spectrum (sorted)')
    plt.tight_layout()
    plt.savefig(out)
    plt.close()


def plot_ldos(evecs, energies, L, out, nstates=2):
    # take lowest-energy (absolute) nstates positive energies
    idx = np.argsort(np.abs(energies))[:nstates]
    probs = np.zeros((nstates, L))
    for k,i in enumerate(idx):
        vec = evecs[:, i]
        # particle components are first L entries
        u = vec[:L]
        probs[k,:] = np.abs(u)**2
    plt.figure(figsize=(6,4))
    for k in range(nstates):
        plt.plot(np.arange(1, L+1), probs[k,:], label=f'state {k}')
    plt.xlabel('Site')
    plt.ylabel('Particle probability')
    plt.title('LDOS (particle component) of lowest-energy states')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out)
    plt.close()


def main():
    path = Path(__file__).resolve().parents[1] / 'ybe.md'
    c = parse_c_from_ybe(path)
    if not c:
        c = {'c_xx':1.0,'c_yy':1.0,'c_xy':0.0,'c_yx':0.0,'c_zz':0.0,'c_z0':0.0,'c_0z':0.0,'c_00':0.0}
    p = compute_p(c)
    t = p['t']
    Delta = p['Delta']
    mu = p['mu']

    L = 40
    Hbdg = build_bdg(L, t, Delta, mu)
    evals, evecs = eigh(Hbdg)
    energies = np.sort(evals)

    outdir = Path('results')
    outdir.mkdir(exist_ok=True)
    plot_spectrum(energies, outdir / 'bdg_spectrum.png')
    plot_ldos(evecs, evals, L, outdir / 'ldos_lowest.png', nstates=2)

    print('Saved plots to results/bdg_spectrum.png and results/ldos_lowest.png')


if __name__ == '__main__':
    main()
