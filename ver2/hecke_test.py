import numpy as np

def paulis():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    si = np.eye(2, dtype=complex)
    return si, sx, sy, sz

def R_from_J(Jx, Jy, Jz, c00=0.0):
    si, sx, sy, sz = paulis()
    H = Jx * np.kron(sx, sx) + Jy * np.kron(sy, sy) + Jz * np.kron(sz, sz) + c00 * np.eye(4, dtype=complex)
    # Use eigendecomposition for Hermitian H to compute exp(i H)
    vals, vecs = np.linalg.eigh(H)
    R = vecs @ np.diag(np.exp(1j * vals)) @ vecs.conj().T
    return R

def test_hecke(R, tol=1e-8):
    eigs, _ = np.linalg.eig(R)
    # normalize tiny numerical phases
    eigs = np.array([e if abs(abs(e) - 1) < 1e-6 else e for e in eigs])
    for q in eigs:
        if abs(q) < 1e-12:
            continue
        ok = True
        for lam in eigs:
            if not (abs(lam - q) < 1e-6 or abs(lam + 1.0 / q) < 1e-6):
                ok = False
                break
        if ok:
            M = (R - q * np.eye(R.shape[0])) @ (R + (1.0 / q) * np.eye(R.shape[0]))
            resid = np.linalg.norm(M)
            return True, q, resid, eigs
    return False, None, None, eigs

def main():
    tests = [
        (np.pi/4, np.pi/4, np.pi/4),
        (0.0, 0.0, 0.0),
        (np.pi/2, 0.0, 0.0),
        (0.0, np.pi/2, 0.0),
        (0.0, 0.0, np.pi/2),
        (np.pi/4, 0.0, 0.0),
    ]
    for Jx, Jy, Jz in tests:
        R = R_from_J(Jx, Jy, Jz)
        ok, q, resid, eigs = test_hecke(R)
        print(f"J = ({Jx:.6f},{Jy:.6f},{Jz:.6f}) -> Hecke-like: {ok}")
        if ok:
            print(f"  q ≈ {q:.6g}, residual norm = {resid:.3e}")
        print(f"  eigenvalues: {[complex(e) for e in eigs]}")
        print()

if __name__ == '__main__':
    main()
