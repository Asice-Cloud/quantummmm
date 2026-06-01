import numpy as np
from scipy.linalg import expm
from scipy.integrate import cumulative_trapezoid
import csv

# load step3_results.csv
import os
csvfile = 'step3_results.csv'
if not os.path.exists(csvfile):
    raise SystemExit('step3_results.csv not found; run step3_solver.py first')

data = np.loadtxt(csvfile, delimiter=',', skiprows=1)
# columns: t,a,b,c,q0,q1,q2,q3,phi
t = data[:,0]
q0 = data[:,4]
q1 = data[:,5]
q2 = data[:,6]
q3 = data[:,7]
phi = data[:,8]

# parameters consistent with step3_solver.py
C_t3 = 2.0
C_Ed = 2.0
gamma1 = 1.0
gamma3 = 1.0
t_c = 0.3
E0 = 0.3
tau = 1.0
pi = np.pi

def f_plus(ti):
    x = pi * ti / tau
    return 0.5*(1+np.cos(x))

def f_minus(ti):
    x = pi * ti / tau
    return 0.5*(1-np.cos(x))

k1 = C_t3 * t_c * np.array([f_plus(tt) for tt in t])  # corresponds to sigma_x coeff
k3 = C_Ed * E0 * np.array([f_minus(tt) for tt in t])  # corresponds to sigma_z coeff

# Build K(t) matrices: K = k1*σx + k3*σz
sigma_x = np.array([[0,1],[1,0]], dtype=complex)
sigma_y = np.array([[0,-1j],[1j,0]], dtype=complex)
sigma_z = np.array([[1,0],[0,-1]], dtype=complex)

# Ω1(t) = ∫_0^t K(t1) dt1 (matrix valued)
# compute integrals for coefficients
int_k1 = cumulative_trapezoid(k1, t, initial=0)
int_k3 = cumulative_trapezoid(k3, t, initial=0)

# Ω2 coefficient (scalar) times σ_y: i * double integral (k3(t1)k1(t2) - k1(t1)k3(t2))
# We'll compute for each final time T: I(T) = ∫_0^T dt1 ∫_0^{t1} dt2 [k3(t1)k1(t2) - k1(t1)k3(t2)]
# Numerically: for each index i (t[i]=T), compute sum over j<=i of k3[i]*sum{k1[0:j]*dt} - k1[i]*sum{k3[0:j]*dt}
dt = t[1]-t[0]
Nk = len(t)
I = np.zeros(Nk)
# precompute cumulative integrals
cum_k1 = cumulative_trapezoid(k1, t, initial=0)
cum_k3 = cumulative_trapezoid(k3, t, initial=0)
for i in range(Nk):
    T = t[i]
    # approximate integral using discrete sums
    s = 0.0
    for m in range(i+1):
        t1 = t[m]
        # inner integral ∫_0^{t1} k1(t2) dt2 = cum_k1[m]
        inner = k3[m]*cum_k1[m] - k1[m]*cum_k3[m]
        s += inner * (t[1]-t[0])
    I[i] = s

# Ω2 = i * I * σ_y

# Now compute Magnus approx V_mag = exp(Ω1 + Ω2)
from numpy.linalg import norm
V_num = []
V_mag = []
errors = []
for i in range(Nk):
    O1 = int_k1[i]*sigma_x + int_k3[i]*sigma_z
    O2 = 1j * I[i] * sigma_y
    Ome = O1 + O2
    Vapprox = expm(Ome)
    V_mag.append(Vapprox)
    # numerical V from quaternion in step3_solver: convert q -> matrix
    q = np.array([q0[i], q1[i], q2[i], q3[i]])
    # quaternion to matrix mapping used earlier
    Vnum = np.array([[q[0]+1j*q[1], q[2]+1j*q[3]],[-q[2]+1j*q[3], q[0]-1j*q[1]]], dtype=complex)
    V_num.append(Vnum)
    err = norm(Vapprox - Vnum)
    errors.append(err)

# save comparison
out = 'magnus_compare.csv'
with open(out, 'w') as f:
    f.write('t,err\n')
    for ti, e in zip(t, errors):
        f.write(f'{ti},{e}\n')

print('Wrote', out)
print('Max error', max(errors))
print('Error at final time', errors[-1])

try:
    import matplotlib.pyplot as plt
    plt.plot(t, errors)
    plt.yscale('log')
    plt.xlabel('t')
    plt.ylabel('||V_mag - V_num||')
    plt.title('Magnus (1+2) vs numerical')
    plt.savefig('magnus_error.png')
    print('Saved magnus_error.png')
except Exception:
    pass
