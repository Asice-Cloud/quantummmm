"""Verify the E1,t1 Dyson/Taylor expansion in n4.md.

The coefficient matrices R_mn are defined by
R(E,t)=sum_{m,n} E**m*t**n*R_mn + O(3).
They are obtained from a single augmented SO(5) variational ODE.
"""
import numpy as np
from scipy.integrate import solve_ivp

TC, E0 = 0.3, 0.3           # meV
TAU = 1.0                   # meV^{-1}; total protocol time is 6*tau
SO5_NORMALIZATION = 2.0

def env(u, rising):
    return (1.0 - np.cos(np.pi*u))/2 if rising else (1.0 + np.cos(np.pi*u))/2

def protocol_values(t):
    # two successive three-step cycles
    q = t % (3.0*TAU)
    seg = min(2, int(q/TAU))
    u = (q-seg*TAU)/TAU
    if seg == 0:
        return TC*env(u, True), 0.0, E0*env(u, False), env(u, True)
    if seg == 1:
        return TC*env(u, False), TC*env(u, True), 0.0, env(u, False)
    return 0.0, TC*env(u, False), E0*env(u, True), 0.0

def generators(t):
    t2, t3, ed, g1 = protocol_values(t)
    s = SO5_NORMALIZATION
    O0 = np.zeros((5,5)); OE = np.zeros((5,5)); Ot = np.zeros((5,5))
    O0[3,1], O0[1,3] = s*t2, -s*t2
    O0[3,2], O0[2,3] = -s*t3, s*t3
    O0[3,4], O0[4,3] = s*ed, -s*ed
    OE[0,1], OE[1,0] = s, -s
    Ot[4,0], Ot[0,4] = -s*g1, s*g1
    return O0, OE, Ot

def unpack(y, count):
    return y.reshape(count,5,5)

def pack(a):
    return a.reshape(-1)

def rhs_coeff(t, y):
    # order: R00,R10,R01,R20,R11,R02
    r = unpack(y, 6)
    O0, OE, Ot = generators(t)
    dr = np.zeros_like(r)
    dr[0] = O0 @ r[0]
    dr[1] = O0 @ r[1] + OE @ r[0]
    dr[2] = O0 @ r[2] + Ot @ r[0]
    dr[3] = O0 @ r[3] + OE @ r[1]
    dr[4] = O0 @ r[4] + OE @ r[2] + Ot @ r[1]
    dr[5] = O0 @ r[5] + Ot @ r[2]
    return pack(dr)

def rhs_exact(t, y, E, tt):
    R = y.reshape(5,5)
    O0, OE, Ot = generators(t)
    return pack((O0 + E*OE + tt*Ot) @ R)

def solve_coeff():
    y0 = np.zeros((6,5,5)); y0[0] = np.eye(5)
    sol = solve_ivp(rhs_coeff, (0,6*TAU), pack(y0), method='DOP853',
                    rtol=2e-12, atol=2e-14)
    return unpack(sol.y[:,-1], 6)

def solve_exact(E, tt):
    sol = solve_ivp(lambda x,y: rhs_exact(x,y,E,tt), (0,6*TAU),
                    pack(np.eye(5)), method='DOP853',
                    rtol=2e-12, atol=2e-14)
    return sol.y[:,-1].reshape(5,5)

def amplitude(R):
    vm = np.array([1,-1j,0,0,0],complex)/np.sqrt(2)
    vp = np.array([1,1j,0,0,0],complex)/np.sqrt(2)
    return np.vdot(vp, R @ vm)

def main():
    c = solve_coeff()
    a = np.array([amplitude(x) for x in c])
    print('Taylor amplitude coefficients a_mn [00,10,01,20,11,02]:')
    for key, val in zip(['00','10','01','20','11','02'], a):
        print(f'a{key} = {val.real:+.12e} {val.imag:+.12e}i')
    print('\\nScaling test: exact amplitude minus quadratic Taylor polynomial')
    for scale in [1.0, 0.5, 0.25, 0.125]:
        E, tt = scale*2e-4, scale*3e-4
        exact = amplitude(solve_exact(E,tt))
        approx = a[0] + E*a[1] + tt*a[2] + E**2*a[3] + E*tt*a[4] + tt**2*a[5]
        err = abs(exact-approx)
        print(f's={scale:5.3f}  |A_exact-A_2|={err:.6e}  err/s^3={err/scale**3:.6e}')
    E, tt = 2e-4, 3e-4
    R = solve_exact(E,tt)
    print('\\nSO(5) check ||R^T R-I||_F =', np.linalg.norm(R.T@R-np.eye(5)))
    print('F_pq exact      =', abs(amplitude(R)))
    A2 = a[0] + E*a[1] + tt*a[2] + E**2*a[3] + E*tt*a[4] + tt**2*a[5]
    print('F_pq quadratic  =', abs(A2))
    print('mixed coefficient a11 =', a[4])

if __name__ == '__main__':
    main()
