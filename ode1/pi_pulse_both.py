#!/usr/bin/env python3
"""
验证：E1 ≠ 0, t1 ≠ 0 时，π 脉冲同时消除两者
===============================================
理论：{σ_x, σ_y} = {σ_x, σ_z} = 0，[σ_x, σ_x] = 0
→ U_π = -iσ_x 翻转 σ_y 和 σ_z 动力学项，保留 σ_x 编织项

对称协议 + 中点 π 脉冲（在完整 SO(5) 中实现）：
  U_total = U(后) · U_π · U(前)
  = U_π† · U(前)† · U_π · U_π · U(前)
  = U_π† = iσ_x

数值验证：完整 SO(5) RK4 + SO(5) π 脉冲矩阵
"""
import numpy as np; pi = np.pi
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

tc=0.3; E0=0.3

def fp(t,tau): return 0.5*(1+np.cos(pi*t/tau))
def fm(t,tau): return 0.5*(1-np.cos(pi*t/tau))

def A1(t,tau,e,t1c):
    A=np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3]=-2*tc*fm(t,tau); A[3,1]=2*tc*fm(t,tau)
    A[0,4]=2*t1c*fm(t,tau); A[4,0]=-2*t1c*fm(t,tau)
    A[3,4]=2*E0*fp(t,tau); A[4,3]=-2*E0*fp(t,tau); return A
def A2(t,tau,e,t1c):
    A=np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[1,3]=-2*tc*fp(t,tau); A[3,1]=2*tc*fp(t,tau)
    A[2,3]=2*tc*fm(t,tau); A[3,2]=-2*tc*fm(t,tau)
    A[0,4]=2*t1c*fp(t,tau); A[4,0]=-2*t1c*fp(t,tau); return A
def A3(t,tau,e,t1c):
    A=np.zeros((5,5)); A[0,1]=2*e; A[1,0]=-2*e
    A[2,3]=2*tc*fp(t,tau); A[3,2]=-2*tc*fp(t,tau)
    A[3,4]=2*E0*fm(t,tau); A[4,3]=-2*E0*fm(t,tau); return A

# SO(5) π_x 脉冲矩阵
# σ_x = iγ₂γ₃，绕此轴转 π → SU(2): U_π = -iσ_x
# 在 Majorana 算符上：σ_y = iγ₃γ₁ → -iγ₃γ₁，σ_z = iγ₁γ₂ → -iγ₁γ₂
# 这要求 γ₁ → -γ₁（γ₂,γ₃ 不变）
# 验证：H_EM 中 E1∝iγ₁γ₂→翻号，t1∝iγ₁γ_b→翻号，编织项∝iγ_aγ₂,iγ_aγ₃,iγ_aγ_b→不变
R_pi = np.eye(5)
R_pi[0,0] = -1  # γ₁ sign flip

def so5_protocol(tau, e, t1c, n=500, pi_pulse=False):
    """SO(5) 传播，可选中点 π 脉冲"""
    dt = tau / n
    R = np.eye(5)
    A_fns = [A1, A2, A3]

    for step_idx, A_fn in enumerate(A_fns):
        for s in range(n):
            t = s * dt
            A0 = A_fn(t, tau, e, t1c)
            k1 = A0 @ R
            A_h = A_fn(t+0.5*dt, tau, e, t1c)
            k2 = A_h @ (R+0.5*dt*k1)
            k3 = A_h @ (R+0.5*dt*k2)
            A_f = A_fn(t+dt, tau, e, t1c)
            k4 = A_f @ (R+dt*k3)
            R += dt/6*(k1+2*k2+2*k3+k4)

            # π 脉冲：t = 1.5*tau (step 1, s = n/2)
            if pi_pulse and step_idx == 1 and s == n//2:
                R = R_pi @ R

    return R

def fid_so5(R):
    Rd = R @ R  # double braid
    ov = 0.5*(Rd[0,0]+1j*Rd[1,0]+1j*Rd[0,1]-Rd[1,1])
    return np.abs(ov)**2

# ═══════════════════════════════════════════════════════════════
print("="*65)
print("π-Pulse in FULL SO(5): Simultaneous E1 + t1 Elimination")
print("="*65)

test_points = [
    (0.0, 0.0,   "pure MZM"),
    (0.0, 0.01,  "E1=0, t1=0.01"),
    (0.01, 0.0,  "E1=0.01, t1=0"),
    (0.01, 0.01, "E1=t1=0.01"),
    (0.02, 0.02, "E1=t1=0.02"),
    (0.05, 0.05, "E1=t1=0.05"),
    (0.01, 0.05, "E1<<t1"),
    (0.05, 0.01, "E1>>t1"),
]

tau = 50.0
print(f"\n{'Case':20s} {'No π':>10s} {'With π':>10s} {'Theory':>10s}")
print("-"*55)

for e, t1c, label in test_points:
    R_no = so5_protocol(tau, e, t1c, pi_pulse=False)
    R_pi = so5_protocol(tau, e, t1c, pi_pulse=True)
    print(f"{label:20s} {fid_so5(R_no):10.6f} {fid_so5(R_pi):10.6f} {'1.000000':>10s}")

# ═══ 扫 (E1, t1) 热力图 ═══
print("\n\nFull SO(5) (E1, t1) sweep...")
n_pts = 25
E1_r = np.logspace(-3, -1, n_pts)
t1_r = np.logspace(-3, -1, n_pts)

F_no = np.zeros((n_pts, n_pts))
F_pi = np.zeros((n_pts, n_pts))

for i, e in enumerate(E1_r):
    for j, t1c in enumerate(t1_r):
        F_no[i,j] = fid_so5(so5_protocol(tau, e, t1c, pi_pulse=False))
        F_pi[i,j] = fid_so5(so5_protocol(tau, e, t1c, pi_pulse=True))
    if (i+1) % 5 == 0: print(f"  {i+1}/{n_pts}")

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
Ev, Tv = np.meshgrid(E1_r, t1_r)

for ci, (ax, D, title) in enumerate([
    (axes[0], F_no, 'No π pulse (original)'),
    (axes[1], F_pi, 'With π pulse at midpoint'),
    (axes[2], F_pi - 1.0, 'Deviation from 1.0'),
]):
    if ci < 2:
        im = ax.pcolormesh(Ev, Tv, D, shading='auto', cmap='viridis', vmin=0, vmax=1)
        ax.contour(Ev, Tv, D, levels=[0.3,0.5,0.7,0.9], colors='white', linewidths=0.5)
    else:
        vm = max(abs(D.min()), abs(D.max()), 1e-6)
        im = ax.pcolormesh(Ev, Tv, D, shading='auto', cmap='RdBu_r', vmin=-vm, vmax=vm)
    plt.colorbar(im, ax=ax)
    ax.set_xlabel('E₁ (meV)'); ax.set_ylabel('t₁ (meV)')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_title(title, fontsize=12)

fig.suptitle(f'π-Pulse: Simultaneous E₁ + t₁ Elimination (Full SO(5), τ={tau})',
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig('pi_pulse_both_so5.png', dpi=200)
plt.close(fig)
print("  Saved: pi_pulse_both_so5.png")

print(f"\n  Max |fid_pi - 1| = {abs(F_pi-1.0).max():.4f}")
print("\nDone.")

