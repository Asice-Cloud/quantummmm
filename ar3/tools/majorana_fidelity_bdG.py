#!/usr/bin/env python3
"""Compute Majorana-mode transformation under single-particle R=exp(i Hbdg)
and compare to ideal braid S (gamma1->gamma2, gamma2->-gamma1).

S_ideal = [[0,1],[-1,0]]

Usage: python tools/majorana_fidelity_bdG.py --h 0.1 --L 120
"""
import argparse
from pathlib import Path
import re
import numpy as np
from scipy.linalg import eigh, expm


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


def psi_to_majorana_vec(psi):
    # psi length 2L: [u; v]
    L2 = psi.shape[0]
    L = L2 // 2
    u = psi[:L]
    v = psi[L:]
    m = np.zeros(2*L, dtype=complex)
    for j in range(L):
        m[2*j]   = 0.5*(u[j] + v[j])
        m[2*j+1] = 0.5j*(u[j] - v[j])
    return m


def majorana_vec_to_psi(m):
    L2 = m.shape[0]
    L = L2 // 2
    u = np.zeros(L, dtype=complex)
    v = np.zeros(L, dtype=complex)
    for j in range(L):
        u[j] = m[2*j] - 1j*m[2*j+1]
        v[j] = m[2*j] + 1j*m[2*j+1]
    return np.concatenate([u, v])


def orthonormalize(vecs):
    # Gram-Schmidt with complex inner product, return real-normalized basis
    B = []
    for v in vecs:
        w = v.copy()
        for b in B:
            proj = np.vdot(b, w) * b
            w = w - proj
        norm = np.sqrt(np.real(np.vdot(w, w)))
        if norm < 1e-12:
            continue
        w = w / norm
        B.append(w)
    return B


def compute_fidelity(H, L, h, center=True):
    # mu background from ybe.md
    vals = parse_c_from_ybe(Path('ybe.md'))
    if not vals:
        vals = {'c_xx':1.1,'c_yy':1.0,'c_xy':0.0,'c_yx':0.0,'c_zz':0.0,'c_z0':0.0,'c_0z':0.0,'c_00':0.0}
    t, Delta, mu0 = compute_p(vals)
    mu_vec = np.zeros(L) + mu0
    if center:
        mu_vec[L//2] += 2.0*h
    else:
        mu_vec[0] += 2.0*h
    Hbdg = build_bdg(L, t, Delta, mu_vec)
    evals, evecs = eigh(Hbdg)
    # take two lowest absolute-energy eigenvectors
    idx2 = np.argsort(np.abs(evals))[:2]
    psi1 = evecs[:, idx2[0]]
    psi2 = evecs[:, idx2[1]]

    m1 = psi_to_majorana_vec(psi1)
    m2 = psi_to_majorana_vec(psi2)
    # form real-like orthonormal basis
    basis = orthonormalize([m1, m2])
    if len(basis) < 2:
        # try real/imag parts
        basis = orthonormalize([np.real(m1), np.imag(m1), np.real(m2), np.imag(m2)])
    if len(basis) < 2:
        raise RuntimeError('Could not construct two orthonormal Majorana vectors')
    alpha = basis[0]
    beta = basis[1]

    # single-particle R
    R = expm(1j * Hbdg)

    # transform basis vectors under R
    def transform(mvec):
        psi = majorana_vec_to_psi(mvec)
        psi2 = R.dot(psi)
        return psi_to_majorana_vec(psi2)

    a_t = transform(alpha)
    b_t = transform(beta)

    # real inner products
    M = np.zeros((2,2), dtype=float)
    M[0,0] = float(np.real(np.vdot(alpha, a_t)))
    M[0,1] = float(np.real(np.vdot(alpha, b_t)))
    M[1,0] = float(np.real(np.vdot(beta, a_t)))
    M[1,1] = float(np.real(np.vdot(beta, b_t)))

    S_ideal = np.array([[0.0, 1.0], [-1.0, 0.0]])
    diff = M - S_ideal
    frob = np.linalg.norm(diff)
    norm_ideal = np.linalg.norm(S_ideal)
    fidelity_like = 1.0 - frob / (norm_ideal + 1e-12)
    return dict(M=M, S_ideal=S_ideal, frob=frob, fidelity_like=fidelity_like, evals=evals[:8])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--h', type=float, default=0.1)
    parser.add_argument('--L', type=int, default=120)
    parser.add_argument('--center', action='store_true', help='apply h at chain center (default)')
    args = parser.parse_args()

    res = compute_fidelity(None, args.L, args.h, center=True)
    print(f'h={args.h:.3f} L={args.L} fidelity_like={res["fidelity_like"]:.6f} frob_diff={res["frob"]:.6e}')
    print('M =')
    print(res['M'])
    print('first few BdG eigenvalues:')
    print(res['evals'])


if __name__ == '__main__':
    main()
