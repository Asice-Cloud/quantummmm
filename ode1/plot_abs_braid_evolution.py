#!/usr/bin/env python3
"""
ABS braiding 全步骤演化图 (类似 PRB105 / PRB111 Fig 1b)。
t₁=0, E₁≠0, 大 τ 绝热, 6 步双次 swap, 追踪 4 个波函数重叠。
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi
t_c = 0.3
E0  = 0.3

def fp(t, tau): return 0.5*(1.0+np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1.0-np.cos(pi*t/tau))

def make_rhs(tau, E1v, step_fn):
    """返回 ODE 右端函数 dy/dt = A(t) y"""
    def rhs(t, y):
        A = np.zeros((5,5))
        if step_fn == 1:
            A[0,1]=2*E1v; A[1,0]=-2*E1v
            t2=t_c*fm(t,tau); A[1,3]=-2*t2; A[3,1]=2*t2
            Ed=E0*fp(t,tau); A[3,4]=2*Ed; A[4,3]=-2*Ed
        elif step_fn == 2:
            A[0,1]=2*E1v; A[1,0]=-2*E1v
            t2=t_c*fp(t,tau); A[1,3]=-2*t2; A[3,1]=2*t2
            t3=t_c*fm(t,tau); A[2,3]=2*t3; A[3,2]=-2*t3
        elif step_fn == 3:
            A[0,1]=2*E1v; A[1,0]=-2*E1v
            t3=t_c*fp(t,tau); A[2,3]=2*t3; A[3,2]=-2*t3
            Ed=E0*fm(t,tau); A[3,4]=2*Ed; A[4,3]=-2*Ed
        R = y.reshape(5,5)
        return (A @ R).reshape(-1)
    return rhs

def propagate_with_history(tau, E1v, step_fn, n_pts_per_step=80):
    """传播一步并返回所有时间点的 R(t)"""
    rhs = make_rhs(tau, E1v, step_fn)
    y0 = np.eye(5).reshape(-1)
    t_eval = np.linspace(0, tau, n_pts_per_step)
    sol = solve_ivp(rhs, (0, tau), y0, t_eval=t_eval,
                    rtol=1e-10, atol=1e-14, method='RK45')
    return sol.t, sol.y.reshape(5,5,-1), sol.y[:,-1].reshape(5,5)

# ── 波函数重叠计算 ─────────────────────────────────────────────────
# MZM 基态和激发态（在 iγ₁γ₂ 本征基下）:
# |ψ₁⟩ ∝ (γ₁ − iγ₂)/√2  → 系数 c₁=1, c₂=−i, c₃=0
# |ψ₂⟩ ∝ (γ₁ + iγ₂)/√2  → 系数 d₁=1, d₂=+i, d₃=0
# 演化后 γ_i' = Σ_j R_{ij} γ_j
# 演化态 |ψ₁(t)⟩ ∝ γ₁' − iγ₂' = Σ_j (R[0,j] − iR[1,j]) γ_j

def compute_overlaps(R):
    """计算四个重叠 |⟨ψ_a(t)|ψ_b(0)⟩|², a,b∈{1,2}，包含 ancilla 分量归一化"""
    c_init_full = np.array([1.0, -1j, 0.0, 0.0, 0.0])  # |ψ₁(0)⟩
    d_init_full = np.array([1.0,  1j, 0.0, 0.0, 0.0])  # |ψ₂(0)⟩
    
    c_evolved_full = R[0,:] - 1j*R[1,:]      # |ψ₁(t)⟩ (full 5-mode)
    d_evolved_full = R[0,:] + 1j*R[1,:]      # |ψ₂(t)⟩ (full 5-mode)
    
    norm_c = np.sqrt(2.0)  # |c_init|² = 2
    norm_d = np.sqrt(2.0)
    norm_ce = np.sqrt(np.sum(np.abs(c_evolved_full)**2))
    norm_de = np.sqrt(np.sum(np.abs(d_evolved_full)**2))
    
    o11 = np.abs(np.dot(c_init_full.conj(), c_evolved_full) / (norm_c * norm_ce))**2
    o12 = np.abs(np.dot(c_init_full.conj(), d_evolved_full) / (norm_c * norm_de))**2
    o21 = np.abs(np.dot(d_init_full.conj(), c_evolved_full) / (norm_d * norm_ce))**2
    o22 = np.abs(np.dot(d_init_full.conj(), d_evolved_full) / (norm_d * norm_de))**2
    
    return o11, o12, o21, o22

# ── 主计算 ─────────────────────────────────────────────────────────
tau = 50
E1_abs = 0.01   # ABS 情况
E1_mzm = 0.0    # MZM 对照

print(f"=== ABS Braiding Evolution (tau={tau}, E1={E1_abs}, t1=0) ===")

# 累积时间和重叠
all_t_mzm = []
all_o_mzm = []
all_t_abs = []
all_o_abs = []

R_cum_mzm = np.eye(5)
R_cum_abs = np.eye(5)

step_names = ['Step1', 'Step2', 'Step3', 'Step1', 'Step2', 'Step3']
step_fns   = [1, 2, 3, 1, 2, 3]

t_offset_mzm = 0.0
t_offset_abs = 0.0

for s_idx, (sname, sfn) in enumerate(zip(step_names, step_fns)):
    # MZM
    t_m, R_hist_m, R_final_m = propagate_with_history(tau, E1_mzm, sfn, n_pts_per_step=60)
    for i in range(len(t_m)):
        R_now = R_hist_m[:,:,i] @ R_cum_mzm
        o11, o12, o21, o22 = compute_overlaps(R_now)
        all_t_mzm.append(t_offset_mzm + t_m[i])
        all_o_mzm.append([o11, o12, o21, o22])
    R_cum_mzm = R_final_m @ R_cum_mzm
    t_offset_mzm += tau
    
    # ABS
    t_a, R_hist_a, R_final_a = propagate_with_history(tau, E1_abs, sfn, n_pts_per_step=60)
    for i in range(len(t_a)):
        R_now = R_hist_a[:,:,i] @ R_cum_abs
        o11, o12, o21, o22 = compute_overlaps(R_now)
        all_t_abs.append(t_offset_abs + t_a[i])
        all_o_abs.append([o11, o12, o21, o22])
    R_cum_abs = R_final_a @ R_cum_abs
    t_offset_abs += tau
    
    print(f"  {sname}: MZM o11={compute_overlaps(R_cum_mzm)[0]:.4f}, "
          f"ABS o11={compute_overlaps(R_cum_abs)[0]:.4f}")

all_t_mzm = np.array(all_t_mzm)
all_o_mzm = np.array(all_o_mzm)
all_t_abs = np.array(all_t_abs)
all_o_abs = np.array(all_o_abs)

print(f"\nFinal MZM overlaps: o11={all_o_mzm[-1,0]:.6f}, o12={all_o_mzm[-1,1]:.6f}, "
      f"o21={all_o_mzm[-1,2]:.6f}, o22={all_o_mzm[-1,3]:.6f}")
print(f"Final ABS overlaps: o11={all_o_abs[-1,0]:.6f}, o12={all_o_abs[-1,1]:.6f}, "
      f"o21={all_o_abs[-1,2]:.6f}, o22={all_o_abs[-1,3]:.6f}")

# ── 绘图 ───────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

colors = ['#2c7bb6', '#fdae61', '#d7191c', '#1a9641']
labels_mzm = [
    r'$|\langle\psi_1(t)|\psi_1(0)\rangle|$',
    r'$|\langle\psi_1(t)|\psi_2(0)\rangle|$',
    r'$|\langle\psi_2(t)|\psi_1(0)\rangle|$',
    r'$|\langle\psi_2(t)|\psi_2(0)\rangle|$',
]
labels_abs = [
    r'$|\langle\psi_1(t)|\psi_1(0)\rangle|$',
    r'$|\langle\psi_1(t)|\psi_2(0)\rangle|$',
    r'$|\langle\psi_2(t)|\psi_1(0)\rangle|$',
    r'$|\langle\psi_2(t)|\psi_2(0)\rangle|$',
]

# ── 左图: MZM 对照 ──────────────────────────────────────────────
ax_mzm = axes[0]
for j in range(4):
    ls = '--' if j == 3 else '-'
    ax_mzm.plot(all_t_mzm/tau, all_o_mzm[:,j], color=colors[j], ls=ls, lw=1.8,
                label=labels_mzm[j])

# 标注步边界
for k in range(1, 6):
    ax_mzm.axvline(x=k, color='gray', ls=':', lw=0.8, alpha=0.5)
ax_mzm.text(1.5, 0.02, 'Step1', ha='center', fontsize=8, color='gray')
ax_mzm.text(2.5, 0.02, 'Step2', ha='center', fontsize=8, color='gray')
ax_mzm.text(3.5, 0.02, 'Step3', ha='center', fontsize=8, color='gray')
ax_mzm.text(4.5, 0.02, 'Step1', ha='center', fontsize=8, color='gray')
ax_mzm.text(5.5, 0.02, 'Step2', ha='center', fontsize=8, color='gray')

ax_mzm.set_xlabel(r'Time $t / \tau$', fontsize=12)
ax_mzm.set_ylabel('Overlap', fontsize=12)
ax_mzm.set_title(r'MZM ($E_1=0, t_1=0, \tau=' + f'{tau}$)', fontsize=13)
ax_mzm.legend(fontsize=8, loc='center right')
ax_mzm.set_ylim(-0.02, 1.05)
ax_mzm.grid(True, alpha=0.3)

# ── 右图: ABS braiding ──────────────────────────────────────────
ax_abs = axes[1]
for j in range(4):
    ls = '--' if j == 3 else '-'
    ax_abs.plot(all_t_abs/tau, all_o_abs[:,j], color=colors[j], ls=ls, lw=1.8,
                label=labels_abs[j])

for k in range(1, 6):
    ax_abs.axvline(x=k, color='gray', ls=':', lw=0.8, alpha=0.5)
ax_abs.text(1.5, 0.02, 'Step1', ha='center', fontsize=8, color='gray')
ax_abs.text(2.5, 0.02, 'Step2', ha='center', fontsize=8, color='gray')
ax_abs.text(3.5, 0.02, 'Step3', ha='center', fontsize=8, color='gray')
ax_abs.text(4.5, 0.02, 'Step1', ha='center', fontsize=8, color='gray')
ax_abs.text(5.5, 0.02, 'Step2', ha='center', fontsize=8, color='gray')

ax_abs.set_xlabel(r'Time $t / \tau$', fontsize=12)
ax_abs.set_ylabel('Overlap', fontsize=12)
ax_abs.set_title(r'ABS ($E_1=0.01, t_1=0, \tau=' + f'{tau}$)', fontsize=13)
ax_abs.legend(fontsize=8, loc='center right')
ax_abs.set_ylim(-0.02, 1.05)
ax_abs.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('abs_braid_evolution.png', dpi=150)
print("\nFigure saved: abs_braid_evolution.png")

# ── 也输出关键值 ──────────────────────────────────────────────────
print("\n=== Key results ===")
print(f"MZM  after step 3 (single swap):  o11={all_o_mzm[179,0]:.6f}, o12={all_o_mzm[179,1]:.6f}")
print(f"MZM  after step 6 (double swap):  o11={all_o_mzm[-1,0]:.6f}, o12={all_o_mzm[-1,1]:.6f}")
print(f"ABS  after step 3 (single swap):  o11={all_o_abs[179,0]:.6f}, o12={all_o_abs[179,1]:.6f}")
print(f"ABS  after step 6 (double swap):  o11={all_o_abs[-1,0]:.6f}, o12={all_o_abs[-1,1]:.6f}")
print(f"\nMZM:  non-Abelian swap complete? o12={all_o_mzm[-1,1]:.4f} (expect 1)")
print(f"ABS:  swap degraded by E1?      o11={all_o_abs[-1,0]:.4f}, o12={all_o_abs[-1,1]:.4f}")
