#!/usr/bin/env python3
"""验证 Nitsch τ₁₂ Riccati ODE 分析：扫描 λ，找使 |q| 最小的值"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

eta = 0.3  # MBS overlap (medium value for clear effect)

# ==== 四元数运算 ====
def qm(p, q):
    pw, px, py, pz = p; qw, qx, qy, qz = q
    return np.array([pw*qw-px*qx-py*qy-pz*qz, pw*qx+px*qw+py*qz-pz*qy,
                     pw*qy-px*qz+py*qw+pz*qx, pw*qz+px*qy-py*qx+pz*qw])

# ==== Riccati ODE for τ₁₂ ====
def riccati_ode(theta, y, lam):
    q = y  # 4-component quaternion [w, x, y_im, z]
    cth = np.cos(theta); sth = np.sin(theta)
    
    # A = i/2*(λ + cosθ), D = i/2*(λ - cosθ)
    # B = i/2*(1-η)*sinθ * j, C = -B†
    # In quaternion: i→(0,1,0,0), j→(0,0,1,0)
    
    # q' = C + Dq - qA - qBq
    # where A = i*a₀, D = i*d₀, B = i*bⱼ*j, C = -i*bⱼ*j
    a0 = 0.5*(lam + cth)
    d0 = 0.5*(lam - cth)
    bj = 0.5*(1 - eta)*sth
    
    # Build quaternions
    A_q = np.array([0, a0, 0, 0])    # i*a₀ → (0, a₀, 0, 0)
    D_q = np.array([0, d0, 0, 0])    # i*d₀
    B_q = np.array([0, 0, bj, 0])    # i*bⱼ*j (i*j = k? no, just separate axes)
    C_q = np.array([0, 0, -bj, 0])   # -i*bⱼ*j
    
    # Compute dq/dθ = C + Dq - qA - qBq
    # Need time parameterization: d/dθ = (dt/dθ) * d/dt
    # Actually, we want dq/dθ directly, and the Hamiltonian is 
    # already in terms of θ. The physical time evolution uses 
    # adiabatic parameterization. We just integrate in θ.
    
    dq = C_q + qm(D_q, q) - qm(q, A_q) - qm(qm(q, B_q), q)
    return dq

# ==== Scan λ values ====
lam_vals = np.linspace(0.0, 0.5, 50)
max_q = []

for lam in lam_vals:
    sol = solve_ivp(
        lambda t, y: riccati_ode(t, y, lam),
        [0, np.pi/2],           # θ from 0 to π/2
        [0.0, 0.0, 0.0, 0.0],  # q(0) = 0
        method='RK45',
        rtol=1e-8, atol=1e-10,
        max_step=0.05
    )
    q_norms = np.sqrt(np.sum(sol.y**2, axis=0))
    max_q.append(np.max(q_norms))

# ==== Optimal λ ====
max_q = np.array(max_q)
opt_idx = np.argmin(max_q)
lam_opt = lam_vals[opt_idx]
lam_theory = 2*eta/np.pi

print(f"η = {eta}")
print(f"Optimal λ (numerical)  = {lam_opt:.4f}")
print(f"λ_theory = 2η/π        = {lam_theory:.4f}")
print(f"Ratio λ_opt/η           = {lam_opt/eta:.4f}")
print(f"Ratio λ_theory/η        = {lam_theory/eta:.4f}")
print(f"Max |q| at λ_opt        = {max_q[opt_idx]:.4f}")

# ==== Plot: |q|max vs λ ====
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

ax1.plot(lam_vals, max_q, 'b-', lw=2)
ax1.axvline(lam_opt, color='r', ls='--', label=f'λ_opt={lam_opt:.3f}')
ax1.axvline(lam_theory, color='g', ls=':', label=f'λ_theory={lam_theory:.3f}')
ax1.set_xlabel('λ', fontsize=12); ax1.set_ylabel('max |q(θ)|', fontsize=12)
ax1.set_title(f'Riccati: η={eta}, optimal λ via min|q|')
ax1.legend(); ax1.grid(True, alpha=0.3)

# ==== Plot: trajectories at optimal vs suboptimal λ ====
for lam, label, ls in [(lam_opt, f'λ={lam_opt:.3f} (opt)', 'b-'),
                         (0.0, 'λ=0', 'r--'),
                         (0.5, 'λ=0.5', 'g--')]:
    sol = solve_ivp(lambda t,y: riccati_ode(t,y,lam), [0,np.pi/2],
                    [0,0,0,0], method='RK45', rtol=1e-8, atol=1e-10, max_step=0.05)
    q_n = np.sqrt(np.sum(sol.y**2, axis=0))
    ax2.plot(sol.t, q_n, ls, lw=2, label=label)
ax2.set_xlabel('θ', fontsize=12); ax2.set_ylabel('|q(θ)|', fontsize=12)
ax2.set_title(f'|q(θ)| trajectories (η={eta})')
ax2.legend(); ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('riccati_verify.png', dpi=150)
print("\nSaved: riccati_verify.png")
