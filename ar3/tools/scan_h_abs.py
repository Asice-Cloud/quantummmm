#!/usr/bin/env python3
"""Scan local potential h at the middle site to find parameters producing smallest |E| and report fidelity.
"""
import numpy as np
from scipy.linalg import expm, eigh, norm
from pathlib import Path


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
    return weights[0][1], weights[1][1]


def main():
    # try to detect if model is quadratic and switch to BdG scanning to save memory
    # parse c's from ybe.md
    import re
    from pathlib import Path
    text = Path('ybe.md').read_text(encoding='utf-8')
    def getc(key):
        m = re.search(r'%s\s*[:=]\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)' % key, text)
        return float(m.group(1)) if m else 0.0

    c_zz = getc('c_zz')
    # string-triggering keys
    string_keys = ['c_x0','c_0x','c_y0','c_0y','c_xz','c_zx','c_yz','c_zy']
    has_string = any(abs(getc(k))>1e-12 for k in string_keys)

    if abs(c_zz) < 1e-12 and not has_string:
        # use BdG single-particle scan
        print('Detected quadratic model: using BdG single-particle scan (low memory).')
        # compute p
        c = {k: getc(k) for k in ['c_xx','c_yy','c_xy','c_yx','c_zz','c_z0','c_0z','c_00']}
        t = c['c_xx'] + c['c_yy'] + 1j*(c['c_xy'] - c['c_yx'])
        Delta = c['c_xx'] - c['c_yy'] - 1j*(c['c_xy'] + c['c_yx'])
        mu0 = 4.0 * c['c_zz'] - 2.0*(c['c_z0'] + c['c_0z'])

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

        L = 40
        hs = np.linspace(0.0, 3.0, 31)
        results = []
        for h in hs:
            # map local sigma^z h to local chemical potential change: delta_mu = 2*h
            mu_vec = np.zeros(L, dtype=float) + mu0
            mid = L//2
            mu_vec[mid] += 2.0 * h
            Hbdg = build_bdg(L, t, Delta, mu_vec)
            evals, _ = eigh(Hbdg)
            emin = np.min(np.abs(evals))
            fid = float('nan')
            p,q = -1,-1
            results.append((h, emin, fid, p, q))
            print(f'h={h:.3f}  Emin={emin:.6e}  fid=NaN  (BdG single-particle scan)')
    else:
        # fall back to many-body spin scan (as before)
        L = 5
        bond_coeffs = {}
        for j in range(L-1):
            bond_coeffs[f'{j}_{j+1}'] = {'c_xx':1.1, 'c_yy':1.0, 'c_xy':0.0, 'c_yx':0.0, 'c_zz':0.0}

        hs = np.linspace(0.0, 3.0, 31)
        results = []
        for h in hs:
            onsite_h = {L//2: h}
            H = build_H_total(L, bond_coeffs, onsite_h)
            evals, evecs = eigh(H)
            emin = np.min(np.abs(evals))

            R = expm(1j * H)
            idx = np.argsort(np.abs(evals))
            low_idx = idx[:2]
            low_evecs = evecs[:, low_idx]
            P = low_evecs @ low_evecs.conj().T

            gam = make_majoranas(L)
            p,q = select_majorana_pair_by_support(gam, low_evecs)
            B = expm((np.pi/4.0) * (gam[p] @ gam[q]))
            U_eff = P @ R @ P
            PB = P @ B @ P
            fid = norm(PB - U_eff) / (norm(PB) + 1e-12)

            results.append((h, emin, fid, p, q))
            print(f'h={h:.3f}  Emin={emin:.6e}  fid={fid:.6f}  pair=({p},{q})')

    results = np.array(results, dtype=object)
    imin = np.argmin(results[:,1].astype(float))
    hbest, ebest, fbest, pbest, qbest = results[imin]
    print('\nBest: h=', hbest, ' Emin=', ebest, ' fid=', fbest, ' pair=', (pbest, qbest))

    out = Path('results')
    out.mkdir(exist_ok=True)
    with open(out / 'scan_h_results.csv','w') as f:
        f.write('h,Emin,fid,p,q\n')
        for r in results:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]}\n")
    print('Saved scan results to results/scan_h_results.csv')


if __name__ == '__main__':
    main()
