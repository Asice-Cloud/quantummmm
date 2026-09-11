"""Direct SO(5) four-point test of E1,t1 non-separability.

For each (E,t), compute the full two-cycle propagator and compare
F(E,t) with both the additive and normalized multiplicative ansatzes.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import svd

TC = 0.3
E0 = 0.3
TAU = 1.0
SO5_NORMALIZATION = 2.0


def env(u, rising):
    return (1.0 - np.cos(np.pi * u)) / 2.0 if rising else (1.0 + np.cos(np.pi * u)) / 2.0


def protocol(q):
    z = q % (3.0 * TAU)
    seg = min(2, int(z / TAU))
    u = (z - seg * TAU) / TAU
    if seg == 0:
        return TC * env(u, True), 0.0, E0 * env(u, False), env(u, True)
    if seg == 1:
        return TC * env(u, False), TC * env(u, True), 0.0, env(u, False)
    return 0.0, TC * env(u, False), E0 * env(u, True), 0.0


def omega(q, E, t1):
    t2, t3, ed, g1 = protocol(q)
    s = SO5_NORMALIZATION
    o = np.zeros((5, 5), dtype=float)
    o[0, 1], o[1, 0] = s * E, -s * E
    o[3, 1], o[1, 3] = s * t2, -s * t2
    o[4, 0], o[0, 4] = -s * t1 * g1, s * t1 * g1
    o[3, 2], o[2, 3] = -s * t3, s * t3
    o[3, 4], o[4, 3] = s * ed, -s * ed
    return 100.0 * o


def project_so5(r):
    u, _, vh = svd(r)
    q = u @ vh
    if np.linalg.det(q) < 0:
        u[:, -1] *= -1.0
        q = u @ vh
    return q


def propagator(E, t1):
    sol = solve_ivp(
        lambda q, y: (omega(q, E, t1) @ y.reshape(5, 5)).reshape(-1),
        (0.0, 6.0 * TAU), np.eye(5).reshape(-1), method="DOP853",
        rtol=2e-12, atol=2e-14,
    )
    return project_so5(sol.y[:, -1].reshape(5, 5))


def paper_fidelity(r):
    vm = np.array([1.0, -1j, 0, 0, 0]) / np.sqrt(2.0)
    vp = np.array([1.0, 1j, 0, 0, 0]) / np.sqrt(2.0)
    return abs(np.vdot(vp, r @ vm))


def group_fidelity(r):
    m = r[:3, :3]
    n = r[3:, :3]
    fsurv = np.sum(m * m) / 3.0
    u, _, vh = svd(m)
    q = u @ vh
    if np.linalg.det(q) < 0:
        u[:, -1] *= -1.0
        q = u @ vh
    b23 = np.array([[1.0, 0, 0], [0, 0, -1.0], [0, 1.0, 0]])
    b = b23 @ b23
    frot = (1.0 + np.trace(b.T @ q)) / 4.0
    return fsurv * frot, fsurv, frot, np.linalg.norm(n, "fro")


def report(E, t1):
    points = {(0.0, 0.0): propagator(0.0, 0.0), (E, 0.0): propagator(E, 0.0),
              (0.0, t1): propagator(0.0, t1), (E, t1): propagator(E, t1)}
    print(f"E={E:.6g}, t1={t1:.6g} meV")
    for name, fun in (("F_pq", paper_fidelity), ("F_group", lambda r: group_fidelity(r)[0])):
        vals = [fun(points[key]) for key in [(E, t1), (E, 0.0), (0.0, t1), (0.0, 0.0)]]
        add = vals[0] - vals[1] - vals[2] + vals[3]
        prod = vals[0] - vals[1] * vals[2] / vals[3]
        print(f"  {name}: F(E,t), F(E,0), F(0,t), F(0,0) = {vals}")
        print(f"    additive residual  Delta_oplus = {add:+.12e}")
        print(f"    product residual   Delta_otimes = {prod:+.12e}")
        print(f"    normalized residuals / |E*t| = {add/abs(E*t1):+.6e}, {prod/abs(E*t1):+.6e}")
    fg, fs, fr, leak = group_fidelity(points[(E, t1)])
    print(f"  mixed-point group components: F_group={fg:.12e}, F_surv={fs:.12e}, F_rot_cond={fr:.12e}, ||N||F={leak:.12e}")


def main():
    for scale in (1.0, 0.5, 0.25):
        report(scale * 2e-4, scale * 3e-4)
    r = propagator(2e-4, 3e-4)
    print("SO(5) orthogonality residual =", np.linalg.norm(r.T @ r - np.eye(5)))


if __name__ == "__main__":
    main()
