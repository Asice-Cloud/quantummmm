#!/usr/bin/env python3
"""
Step 2 interaction picture: peel off H0 = E1 * X1, 
propagate the remainder in the rotating frame, compare with full T exp.
"""
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp, cumulative_trapezoid
import csv
import os

# ── 1. Parameters (paper-typical) ──────────────────────────────────
t_c   = 0.3          # coupling peak
E0    = 0.3          # dot level (≈0 in Step 2, but we include it)
E1    = 0.3          # constant left overlap (varied below)
t1_c  = 0.01         # small t1 amplitude for ABS
tau   = 1.0          # one step duration
pi    = np.pi

# ── 2. Gate profiles ───────────────────────────────────────────────
def f_plus(t):
    """rising: (1+cos)/2"""
    x = pi * t / tau
    return 0.5 * (1.0 + np.cos(x))

def f_minus(t):
    """falling: (1-cos)/2"""
    x = pi * t / tau
    return 0.5 * (1.0 - np.cos(x))

def t2_s2(t):
    """Step 2: t2 descending from t_c to 0"""
    return t_c * f_plus(t)

def t3_s2(t):
    """Step 2: t3 ascending from 0 to t_c"""
    return t_c * f_minus(t)

def t1_s2(t):
    """Step 2: small t1 (ABS)"""
    return t1_c * np.ones_like(t)

def Ed_s2(t):
    """Step 2: Ed ≈ 0"""
    return np.zeros_like(t)

# ── 3. Build 5×5 A(t) matrix ───────────────────────────────────────
#  Basis order: (γ1, γ2, γ3, γa, γb)
def build_A(t):
    """so(5) antisymmetric matrix A(t) = -A(t)^T for Step 2."""
    t2 = t2_s2(t)
    t3 = t3_s2(t)
    t1 = t1_s2(t)
    Ed = Ed_s2(t)
    A = np.zeros((5,5), dtype=float)
    # E1 term: i E1 γ1 γ2
    A[0,1] =  2 * E1
    A[1,0] = -2 * E1
    # t2 term: i|t2| γa γ2  (γa is index 3, γ2 is index 1)
    A[1,3] = -2 * t2      # γ2 row, γa col
    A[3,1] =  2 * t2
    # t1 term: -i|t1| γb γ1 = i|t1| γ1 γb  (γ1 idx 0, γb idx 4)
    A[0,4] =  2 * t1
    A[4,0] = -2 * t1
    # t3 term: -i|t3| γa γ3 = i|t3| γ3 γa  (γ3 idx 2, γa idx 3)
    A[2,3] =  2 * t3
    A[3,2] = -2 * t3
    # Ed term: i Ed γa γb  (usually 0 in Step 2)
    A[3,4] =  2 * Ed
    A[4,3] = -2 * Ed
    return A

def build_A0():
    """H0 = constant E1 part only."""
    A0 = np.zeros((5,5), dtype=float)
    A0[0,1] =  2 * E1
    A0[1,0] = -2 * E1
    return A0

def build_A_rem(t):
    """Remainder A_rem = A - A0."""
    A = build_A(t)
    A0 = build_A0()
    return A - A0

# ── 4. Full numerical T exp (direct propagation) ───────────────────
def ode_rhs_full(t, y_flat):
    """y = vec(R) where R is 5×5, dR/dt = A(t) R."""
    A = build_A(t)
    R = y_flat.reshape(5,5)
    dR = A @ R
    return dR.reshape(-1)

def direct_propagator(t_span, t_eval):
    y0 = np.eye(5).reshape(-1)
    sol = solve_ivp(ode_rhs_full, t_span, y0, t_eval=t_eval,
                    rtol=1e-9, atol=1e-12, method='RK45')
    R_list = [sol.y[:,i].reshape(5,5) for i in range(len(sol.t))]
    return sol.t, R_list

# ── 5. Interaction picture propagation ─────────────────────────────
def ode_rhs_IP(t, y_flat):
    """dR_I/dt = A_I(t) R_I where A_I = exp(-A0 t) A_rem(t) exp(A0 t)."""
    A0 = build_A0()
    A_rem = build_A_rem(t)
    # exp(-A0 t) * A_rem * exp(A0 t)
    Rt = expm(-A0 * t)
    Rt_inv = expm(A0 * t)
    A_I = Rt @ A_rem @ Rt_inv
    R_I = y_flat.reshape(5,5)
    dR_I = A_I @ R_I
    return dR_I.reshape(-1)

def IP_propagator(t_span, t_eval):
    y0 = np.eye(5).reshape(-1)
    sol = solve_ivp(ode_rhs_IP, t_span, y0, t_eval=t_eval,
                    rtol=1e-9, atol=1e-12, method='RK45')
    R_I_list = [sol.y[:,i].reshape(5,5) for i in range(len(sol.t))]
    return sol.t, R_I_list

# ── 6. Magnus (order 1+2) on A_I ───────────────────────────────────
def magnus_I(t_eval):
    """Magnus Ω1+Ω2 on the IP Hamiltonian A_I(t)."""
    A0 = build_A0()
    dt = t_eval[1] - t_eval[0]
    N = len(t_eval)

    # precompute A_I(t) for all t
    A_I_all = np.zeros((N, 5, 5))
    for i, ti in enumerate(t_eval):
        Rt = expm(-A0 * ti)
        Rt_inv = expm(A0 * ti)
        A_I_all[i] = Rt @ build_A_rem(ti) @ Rt_inv

    # Ω1(t) = ∫_0^t A_I(τ) dτ
    int_A_I = cumulative_trapezoid(A_I_all, t_eval, axis=0, initial=0)

    # Ω2 commutator integral
    # Ω2(t) = 1/2 ∫_0^t dτ1 ∫_0^τ1 dτ2 [A_I(τ1), A_I(τ2)]
    omega2 = np.zeros((N, 5, 5))
    for i in range(N):
        comm_sum = np.zeros((5,5))
        for m in range(i+1):
            A1 = A_I_all[m]
            # inner integral: cumulative_A_I up to m
            inner = A1 @ int_A_I[m] - int_A_I[m] @ A1
            comm_sum += inner * dt
        omega2[i] = 0.5 * comm_sum

    R_mag_I = []
    for i in range(N):
        Omega = int_A_I[i] + omega2[i]
        R_mag_I.append(expm(Omega))
    return R_mag_I

# ── 7. Run and compare ─────────────────────────────────────────────
def main():
    t_span = (0.0, tau)
    t_eval = np.linspace(0, tau, 201)

    print("=== Step 2 Interaction Picture Analysis ===")
    print(f"Parameters: E1={E1}, t_c={t_c}, tau={tau}")

    # Direct
    t_full, R_full = direct_propagator(t_span, t_eval)

    # IP
    t_ip, R_I = IP_propagator(t_span, t_eval)
    A0 = build_A0()

    # Reconstruct full from IP: R_full_IP = exp(A0 t) * R_I
    R_from_ip = []
    for i, ti in enumerate(t_ip):
        R0 = expm(A0 * ti)
        R_from_ip.append(R0 @ R_I[i])

    # Magnus on IP
    R_mag_I = magnus_I(t_eval)
    R_from_mag = []
    for i, ti in enumerate(t_eval):
        R0 = expm(A0 * ti)
        R_from_mag.append(R0 @ R_mag_I[i])

    # Compare
    err_ip_vs_direct = []
    err_mag_vs_direct = []
    for i in range(len(t_full)):
        err_ip_vs_direct.append(np.linalg.norm(R_full[i] - R_from_ip[i]))
        err_mag_vs_direct.append(np.linalg.norm(R_full[i] - R_from_mag[i]))

    print(f"IP recon max error vs direct:  {max(err_ip_vs_direct):.2e}")
    print(f"IP recon err at final time:    {err_ip_vs_direct[-1]:.2e}")
    print(f"Magnus(1+2) max err vs direct: {max(err_mag_vs_direct):.2e}")
    print(f"Magnus(1+2) err at final time: {err_mag_vs_direct[-1]:.2e}")

    # Check unitarity: R_full should be in SO(5), i.e. det ≈ 1, R R^T = I
    Rf = R_full[-1]
    print(f"det(R_full final) = {np.linalg.det(Rf):.6f}")
    print(f"||R_full R_full^T - I|| = {np.linalg.norm(Rf @ Rf.T - np.eye(5)):.2e}")
    print(f"R_full final =\n{Rf}")

    # Save
    out = 'step2_ip_results.csv'
    with open(out, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['t', 'err_IP', 'err_Magnus'])
        for i, ti in enumerate(t_full):
            writer.writerow([ti, err_ip_vs_direct[i], err_mag_vs_direct[i]])
    print(f'Wrote {out}')

    # ── Plot ───────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        ax = axes[0]
        ax.plot(t_full, [np.linalg.norm(r) for r in err_ip_vs_direct],
                'b-', label='IP recon error')
        ax.plot(t_full, [np.linalg.norm(r) if np.isscalar(r) else r for r in err_mag_vs_direct],
                'r--', label='Magnus(1+2) on IP error')
        ax.set_yscale('log')
        ax.set_xlabel('t')
        ax.set_ylabel('||error||')
        ax.set_title(f'Step 2 IP (E1={E1})')
        ax.legend()
        ax.grid(True, alpha=0.3)

        ax = axes[1]
        # Show A_I(t) norm evolution
        A_I_norms = []
        for ti in t_eval:
            A0_ = build_A0()
            Rt = expm(-A0_ * ti)
            Rt_inv = expm(A0_ * ti)
            A_I = Rt @ build_A_rem(ti) @ Rt_inv
            A_I_norms.append(np.linalg.norm(A_I, 'fro'))
        ax.plot(t_eval, A_I_norms, 'g-', label='||A_I(t)||_F')
        ax.set_xlabel('t')
        ax.set_ylabel('Frobenius norm')
        ax.set_title('Interaction-picture generator norm')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('step2_ip_analysis.png', dpi=120)
        print('Saved step2_ip_analysis.png')
    except Exception as e:
        print(f'Plot skipped: {e}')

    return t_full, R_full, R_from_ip, R_from_mag

if __name__ == '__main__':
    _, _, _, _ = main()
