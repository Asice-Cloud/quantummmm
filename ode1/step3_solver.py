import numpy as np
from scipy.integrate import solve_ivp
import csv

# Parameters (choose typical paper values)
C_t3 = 2.0
C_Ed = 2.0
C_E1 = 2.0
gamma1 = 1.0
gamma3 = 1.0

# Physical parameters (meV, time units arbitrary)
t_c = 0.3
E0 = 0.3
E1 = 0.01

# protocol time for a single step
tau = 1.0

# gating profiles
pi = np.pi

def f_plus(t):
    x = pi * t / tau
    return 0.5 * (1 + np.cos(x))

def f_minus(t):
    x = pi * t / tau
    return 0.5 * (1 - np.cos(x))

# angular velocity components
def omega1(t):
    return gamma1 * C_t3 * t_c * f_plus(t)

def omega3(t):
    return gamma3 * C_Ed * E0 * f_minus(t)

# ODE for a,b,c
def rhs(t, y):
    a, b, c = y
    w1 = omega1(t)
    w3 = omega3(t)
    # guard cos(b)
    cb = np.cos(b)
    if np.abs(cb) < 1e-8:
        # avoid division by zero; return large derivative to signal stiffness
        dc = np.sign(cb) * 1e6
    else:
        dc = w3 * np.cos(a) / cb
    da = w1 - w3 * np.cos(a) * np.tan(b)
    db = w3 * np.sin(a)
    return [da, db, dc]

# quaternion helpers
def quat_mul(q1, q2):
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    a = a1*a2 - b1*b2 - c1*c2 - d1*d2
    b = a1*b2 + b1*a2 + c1*d2 - d1*c2
    c = a1*c2 - b1*d2 + c1*a2 + d1*b2
    d = a1*d2 + b1*c2 - c1*b2 + d1*a2
    return np.array([a,b,c,d])

def quat_from_axis_angle(axis, theta):
    # axis is 'i','j','k' mapped to unit imaginary
    ca = np.cos(theta/2)
    sa = np.sin(theta/2)
    if axis == 'i':
        return np.array([ca, sa, 0.0, 0.0])
    if axis == 'j':
        return np.array([ca, 0.0, sa, 0.0])
    if axis == 'k':
        return np.array([ca, 0.0, 0.0, sa])
    raise ValueError('axis')

# convert quaternion to 2x2 complex matrix (a + b i + c j + d k)
def quat_to_mat(q):
    a,b,c,d = q
    # mapping: [[a+ib, c+id], [-c+id, a-ib]]
    m00 = a + 1j*b
    m01 = c + 1j*d
    m10 = -c + 1j*d
    m11 = a - 1j*b
    return np.array([[m00, m01],[m10, m11]], dtype=complex)

# integrate over one step (0..tau)
t_span = (0.0, tau)
t_eval = np.linspace(0, tau, 201)
y0 = [0.0, 0.0, 0.0]
sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval, rtol=1e-8, atol=1e-10)

# compute q(t), V(t) matrices, phase
qs = []
Vs = []
phis = []
for ti, ai, bi, ci in zip(sol.t, sol.y[0], sol.y[1], sol.y[2]):
    # quaternion q = exp(a i/2) * exp(b j/2) * exp(c k/2)
    qA = quat_from_axis_angle('i', ai)
    qB = quat_from_axis_angle('j', bi)
    qC = quat_from_axis_angle('k', ci)
    qAB = quat_mul(qA, qB)
    q = quat_mul(qAB, qC)
    qs.append(q)
    Vmat = quat_to_mat(q)
    Vs.append(Vmat)
    phi = C_E1 * E1 * ti
    phis.append(phi)

# save results
out_csv = 'step3_results.csv'
with open(out_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['t','a','b','c','q0','q1','q2','q3','phi'])
    for i, tval in enumerate(sol.t):
        a,b,c = sol.y[0,i], sol.y[1,i], sol.y[2,i]
        q0,q1,q2,q3 = qs[i]
        phi = phis[i]
        writer.writerow([tval, a, b, c, q0, q1, q2, q3, phi])

print('Wrote', out_csv)
print('Final quaternion q(tau)=', qs[-1])
print('Final phase Phi(tau)=', phis[-1])
print('Final V matrix=\n', Vs[-1])

# Optionally plot if matplotlib is available
try:
    import matplotlib.pyplot as plt
    plt.plot(sol.t, [q[0] for q in qs], label='q0')
    plt.plot(sol.t, [q[1] for q in qs], label='q1')
    plt.plot(sol.t, [q[2] for q in qs], label='q2')
    plt.plot(sol.t, [q[3] for q in qs], label='q3')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('quaternion components')
    plt.savefig('step3_q_components.png')
    print('Saved step3_q_components.png')
except Exception:
    pass
