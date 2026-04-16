import numpy as np

def paulis():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    si = np.eye(2, dtype=complex)
    return si, sx, sy, sz

def swap4():
    P = np.zeros((4,4), dtype=complex)
    for a in range(2):
        for b in range(2):
            i = 2*a + b
            j = 2*b + a
            P[i, j] = 1
    return P

def R_from_J(Jx, Jy, Jz, c00=0.0):
    si, sx, sy, sz = paulis()
    H = Jx * np.kron(sx, sx) + Jy * np.kron(sy, sy) + Jz * np.kron(sz, sz) + c00 * np.eye(4, dtype=complex)
    vals, vecs = np.linalg.eigh(H)
    R = vecs @ np.diag(np.exp(1j * vals)) @ vecs.conj().T
    return R

def hecke_residual(R, q):
    M = (R - q*np.eye(4)) @ (R + (1.0/q)*np.eye(4))
    return np.linalg.norm(M)

def find_scale_for_hecke(R, tol=1e-8):
    # attempt to find s (complex) and q such that (sR - qI)(sR + q^{-1}I) ~ 0
    lam, vecs = np.linalg.eig(R)
    best = {'res': np.inf}
    n = len(lam)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            lam1 = lam[i]
            lam2 = lam[j]
            denom = lam1 * lam2
            if abs(denom) < 1e-15:
                continue
            s2 = -1.0 / denom
            for sign in [1, -1]:
                s = sign * np.sqrt(s2)
                q = s * lam1
                res = hecke_residual(s*R, q)
                # check how many eigenvalues map close to q or -1/q
                mapped = s * lam
                counts = sum([1 for mv in mapped if abs(mv - q) < 1e-6 or abs(mv + 1.0/q) < 1e-6])
                if res < best['res']:
                    best = {'res': res, 's': s, 'q': q, 'counts': counts, 'mapped': mapped}
    return best

def check_PR_variants(R):
    P = swap4()
    results = {}
    for name, M in [('R', R), ('P·R', P @ R), ('R·P', R @ P)]:
        lam, _ = np.linalg.eig(M)
        # simple Hecke-like test via eigenvalue pattern
        uniq = np.unique(np.round(lam, 8))
        results[name] = {'eigenvalues': lam, 'unique': uniq}
    return results

def fit_quadratic(R):
    I = np.eye(4, dtype=complex)
    R2 = R @ R
    # solve vec(R2) = alpha * vec(R) + beta * vec(I)
    A = np.vstack([R.reshape(-1), I.reshape(-1)]).T
    b = R2.reshape(-1)
    # least squares for complex by stacking real and imag
    A2 = np.vstack([np.hstack([A.real, -A.imag]), np.hstack([A.imag, A.real])])
    b2 = np.concatenate([b.real, b.imag])
    x, *_ = np.linalg.lstsq(A2, b2, rcond=None)
    # recover complex alpha,beta from stacked real solution
    alpha = x[0] + 1j * x[len(x)//2] if len(x) > 1 else x[0]
    beta = x[1] + 1j * x[len(x)//2+1] if len(x) > 2 else x[1]
    resid = np.linalg.norm(R2 - (alpha * R + beta * I))
    return {'alpha': alpha, 'beta': beta, 'resid': resid}

def run_diagnostics(tests):
    out = []
    for Jx, Jy, Jz in tests:
        R = R_from_J(Jx, Jy, Jz)
        lam = np.linalg.eigvals(R)
        hecke_direct = False
        best_q = None
        for qcand in lam:
            if abs(hecke_residual(R, qcand)) < 1e-8:
                hecke_direct = True
                best_q = qcand
                break
        scale_best = find_scale_for_hecke(R)
        pr = check_PR_variants(R)
        quad = fit_quadratic(R)
        out.append({'J': (Jx, Jy, Jz), 'lam': lam, 'hecke_direct': hecke_direct, 'best_q': best_q, 'scale_best': scale_best, 'PR': pr, 'quad': quad})
    return out

def main():
    tests = [
        (np.pi/4, np.pi/4, np.pi/4),
        (0.0, 0.0, 0.0),
        (np.pi/2, 0.0, 0.0),
        (0.0, np.pi/2, 0.0),
        (0.0, 0.0, np.pi/2),
        (np.pi/4, 0.0, 0.0),
    ]
    results = run_diagnostics(tests)
    for r in results:
        J = r['J']
        print(f"J = ({J[0]:.6f},{J[1]:.6f},{J[2]:.6f})")
        print(f"  eigenvalues: {[complex(x) for x in r['lam']]}")
        print(f"  hecke_direct: {r['hecke_direct']}, best_q: {r['best_q']}")
        sb = r['scale_best']
        print(f"  scale search best residual: {sb['res']:.3e}, s≈{sb.get('s')}, q≈{sb.get('q')}, matched eigencount={sb.get('counts')}")
        print(f"  quad fit: alpha≈{r['quad']['alpha']}, beta≈{r['quad']['beta']}, resid={r['quad']['resid']:.3e}")
        for name, info in r['PR'].items():
            print(f"  {name} unique eigenvals (rounded): {np.unique(np.round(info['eigenvalues'],6))}")
        print()

if __name__ == '__main__':
    main()
