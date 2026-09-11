"""Direct four-point test of E1,t1 separability.

Tests both the paper matrix-element fidelity and the SO(5) group fidelity:
  Delta_add = F(E,t)-F(E,0)-F(0,t)+F(0,0)
  Delta_prod = F(E,t)-F(E,0)*F(0,t)/F(0,0)
"""
import sys
import numpy as np
from scipy.linalg import polar

sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
from verify_n4_expansion import solve_exact, amplitude

B23 = np.array([[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]])
B = B23 @ B23


def f_pq(R):
    return abs(amplitude(R))


def f_group(R):
    M = R[:3, :3]
    N = R[3:, :3]
    fsurv = np.sum(M*M) / 3.0
    U, _ = polar(M)
    if np.linalg.det(U) < 0:
        u, s, vh = np.linalg.svd(M)
        u[:, -1] *= -1.0
        U = u @ vh
    frot = (1.0 + np.trace(B.T @ U)) / 4.0
    return fsurv * frot


def four_point(E, t, fn):
    vals = [fn(solve_exact(e, q)) for e, q in
            [(E, t), (E, 0.0), (0.0, t), (0.0, 0.0)]]
    f_et, f_e0, f_0t, f_00 = vals
    add = f_et - f_e0 - f_0t + f_00
    prod = f_et - f_e0 * f_0t / f_00
    return vals, add, prod


def main():
    print("Direct four-point non-separability test")
    print("protocol: tau=1 meV^-1, T=6 tau, Omega_ij=2 h_ij")
    for fn, name in [(f_pq, "F_pq"), (f_group, "F_group")]:
        print(f"\n{name}")
        base_E, base_t = 2e-4, 3e-4
        for scale in [1.0, 0.5, 0.25]:
            vals, add, prod = four_point(scale*base_E, scale*base_t, fn)
            print(f"s={scale:4.2f} values={vals}")
            print(f"  Delta_add  = {add:+.12e}, /s^2={add/scale**2:+.12e}")
            print(f"  Delta_prod = {prod:+.12e}, /s^2={prod/scale**2:+.12e}")


if __name__ == '__main__':
    main()
