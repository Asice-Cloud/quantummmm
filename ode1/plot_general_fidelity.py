#!/usr/bin/env python3
"""
PRB111 一般情况 (E₁≠0, t₁≠0) 的保真度 vs τ 绘图。
单次 swap 和双次 swap，验证 γ₂↔γ₃ 是否互换。
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi
t_c = 0.3       # coupling peak
E0  = 0.3       # dot level

# ── 门控 profile ──────────────────────────────────────────────────
def fp(t, tau): return 0.5 * (1.0 + np.cos(pi * t / tau))
def fm(t, tau): return 0.5 * (1.0 - np.cos(pi * t / tau))

# ── 构建三段 A(t) ──────────────────────────────────────────────────
def build_A_step1(t, tau, E1v, t1v):
    A = np.zeros((5,5))
    A[0,1] =  2*E1v;  A[1,0] = -2*E1v
    t2 = t_c * fm(t, tau); A[1,3] = -2*t2; A[3,1] = 2*t2
    t1 = t1v * fm(t, tau); A[0,4] =  2*t1; A[4,0] = -2*t1
    Ed = E0 * fp(t, tau); A[3,4] =  2*Ed; A[4,3] = -2*Ed
    return A

def build_A_step2(t, tau, E1v, t1v):
    A = np.zeros((5,5))
    A[0,1] =  2*E1v;  A[1,0] = -2*E1v
    t2 = t_c * fp(t, tau); A[1,3] = -2*t2; A[3,1] = 2*t2
    t3 = t_c * fm(t, tau); A[2,3] =  2*t3; A[3,2] = -2*t3
    t1 = t1v * fp(t, tau); A[0,4] =  2*t1; A[4,0] = -2*t1
    return A

def build_A_step3(t, tau, E1v, t1v):
    A = np.zeros((5,5))
    A[0,1] =  2*E1v;  A[1,0] = -2*E1v
    t3 = t_c * fp(t, tau); A[2,3] =  2*t3; A[3,2] = -2*t3
    Ed = E0 * fm(t, tau); A[3,4] =  2*Ed; A[4,3] = -2*Ed
    t1 = t1v * np.ones_like(t); A[0,4] =  2*t1; A[4,0] = -2*t1
    return A

# ── 传播一步 ───────────────────────────────────────────────────────
def propagate_step(build_fn, tau, E1v, t1v):
    def rhs(t, y):
        A = build_fn(t, tau, E1v, t1v)
        R = y.reshape(5,5)
        return (A @ R).reshape(-1)
    y0 = np.eye(5).reshape(-1)
    sol = solve_ivp(rhs, (0, tau), y0, rtol=1e-10, atol=1e-14, method='RK45')
    return sol.y[:,-1].reshape(5,5)

# ── 三段全协议 ─────────────────────────────────────────────────────
def run_protocol(tau, E1v, t1v):
    R1 = propagate_step(build_A_step1, tau, E1v, t1v)
    R2 = propagate_step(build_A_step2, tau, E1v, t1v)
    R3 = propagate_step(build_A_step3, tau, E1v, t1v)
    R_one = R3 @ R2 @ R1        # 单次 swap
    R_two = R_one @ R_one        # 双次 swap
    return R_one, R_two

# ── 保真度计算 ─────────────────────────────────────────────────────
# 基于 Majorana 波函数：
# |ψ₁⁺⟩ = (γ₁ − iγ₂)/√2  初始态
# 单次 swap 期望： γ₂→γ₃，所以 |ψ_expected_single⟩ = (γ₁ − iγ₃)/√2
# 双次 swap 期望： γ₂→γ₃→−γ₂，γ₃→−γ₂→−γ₃，|ψ⟩ 回到自身

def fidelity_single_swap(R):
    """单次 swap 保真度：当前态 vs 期望 (γ₁−iγ₃)/√2"""
    # 初始态在 R 作用下： γ_i' = Σ_j R_{ij} γ_j
    # |ψ(t)⟩ ∝ γ₁' − i γ₂'
    # |⟨ψ_expected|ψ(t)⟩|² = |⟨(γ₁−iγ₃)|(γ₁'−iγ₂')⟩|² / 4
    # 用 Majorana 内积 ⟨γ_i|γ_j⟩ = δ_{ij}
    # γ₁' − iγ₂' = R[0,j]γ_j − i R[1,j]γ_j = (R[0,j] − iR[1,j])γ_j
    # ⟨γ₁−iγ₃| = γ₁† + iγ₃† → 投影到 (R[0,0]−iR[1,0]) 和 (R[0,2]−iR[1,2]) 上
    overlap = (R[0,0] - 1j*R[1,0]) + 1j*(R[0,2] - 1j*R[1,2])
    # 上面不太对，让我重新来
    # 初始 |ψ₀⟩ ∼ γ₁ − iγ₂，系数 c₁=1, c₂=−i, c₃=0
    # 末态 |ψ(t)⟩ ∼ γ₁' − iγ₂' = Σ_j (R[0,j] − iR[1,j]) γ_j
    # |ψ_expected⟩ ∼ γ₁ − iγ₃，系数 d₁=1, d₂=0, d₃=−i
    # overlap = Σ_j d_j* (R[0,j] − iR[1,j])
    d = np.array([1.0, 0.0, -1j, 0.0, 0.0])  # γ₁ − iγ₃
    c_evolved = R[0,:3] - 1j*R[1,:3]  # 只看 γ₁,γ₂,γ₃ 分量
    norm_evolved = np.sqrt(np.sum(np.abs(R[0,:3] - 1j*R[1,:3])**2))
    norm_expected = np.sqrt(2.0)
    overlap = np.dot(d[:3].conj(), c_evolved) / (norm_evolved * norm_expected)
    return np.abs(overlap)**2

def fidelity_double_swap(R):
    """双次 swap 保真度：两次 braid 后 γ₂→−γ₂，期望态 (γ₁+iγ₂)/√2"""
    d = np.array([1.0, 1j, 0.0, 0.0, 0.0])  # γ₁ + iγ₂ (两次 braid 的预期)
    c_evolved = R[0,:3] - 1j*R[1,:3]
    norm_evolved = np.sqrt(np.sum(np.abs(R[0,:3] - 1j*R[1,:3])**2))
    norm_expected = np.sqrt(2.0)
    overlap = np.dot(d[:3].conj(), c_evolved) / (norm_evolved * norm_expected)
    return np.abs(overlap)**2

# ── γ₂↔γ₃ 互换检查 ────────────────────────────────────────────────
def check_g2g3_swap(R):
    """检查 γ₂→γ₃, γ₃→−γ₂ 是否成立"""
    swap_quality = np.abs(R[1,2])  # γ₂ 中 γ₃ 的分量 (期望 ~1)
    sign_flip   = np.abs(R[2,1] + 1.0)  # γ₃ 中 −γ₂ 的分量 (期望 ~0)
    return swap_quality, sign_flip

# ── 主计算 ─────────────────────────────────────────────────────────
print("=== PRB111 一般情况保真度扫描 ===")

tau_vals = np.linspace(1, 100, 80)

# 三组参数
param_sets = [
    ("MZM: E₁=0, t₁=0",      0.0,   0.0,   'C0', '-'),
    ("ABS-small: E₁=0.01, t₁=0.005", 0.01,  0.005, 'C1', '--'),
    ("ABS-large: E₁=0.3, t₁=0.01",   0.3,   0.01,  'C2', '-.'),
]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax_fid_single = axes[0,0]
ax_fid_double = axes[0,1]
ax_swap       = axes[1,0]
ax_matrix     = axes[1,1]

results = {}

for label, E1v, t1v, color, ls in param_sets:
    print(f"\n--- {label} ---")
    fid_single = []
    fid_double = []
    swap_q = []
    sign_err = []
    tau_plot = []
    
    for tau in tau_vals:
        try:
            R1, R2 = run_protocol(tau, E1v, t1v)
            fs = fidelity_single_swap(R1)
            fd = fidelity_double_swap(R2)
            sq, se = check_g2g3_swap(R1)
            fid_single.append(fs)
            fid_double.append(fd)
            swap_q.append(sq)
            sign_err.append(se)
            tau_plot.append(tau)
        except Exception as e:
            pass
    
    results[label] = (np.array(tau_plot), np.array(fid_single), 
                       np.array(fid_double), np.array(swap_q), np.array(sign_err))
    
    # 保真度绘制
    ax_fid_single.plot(tau_plot, fid_single, color=color, ls=ls, lw=2, label=label)
    ax_fid_double.plot(tau_plot, fid_double, color=color, ls=ls, lw=2, label=label)
    ax_swap.plot(tau_plot, swap_q, color=color, ls=ls, lw=2, label=f'{label} (γ₂→γ₃)')
    ax_swap.plot(tau_plot, 1.0 - np.array(sign_err), color=color, ls=':', lw=1.5, 
                 label=f'{label} (γ₃→−γ₂)')
    
    # 打印几个关键 τ 值
    if len(tau_plot) > 0:
        print(f"  τ≈{tau_plot[0]:.1f}: fid_single={fid_single[0]:.4f}, fid_double={fid_double[0]:.4f}, γ₂→γ₃={swap_q[0]:.4f}")
    for t_idx in [20, 50, 100]:
        idx = np.argmin(np.abs(np.array(tau_plot) - t_idx))
        if idx < len(fid_single):
            print(f"  τ≈{tau_plot[idx]:.1f}: fid_single={fid_single[idx]:.4f}, "
                  f"fid_double={fid_double[idx]:.4f}, γ₂→γ₃={swap_q[idx]:.4f}")

# ── 子图 1: 单次 swap 保真度 ──────────────────────────────────────
ax_fid_single.axhline(y=1.0, color='gray', ls=':', lw=1, alpha=0.5)
ax_fid_single.annotate('ideal single-swap = 1', xy=(50, 0.97), fontsize=9, color='gray')
ax_fid_single.set_xlabel(r'$\tau$ (adiabatic parameter)', fontsize=12)
ax_fid_single.set_ylabel(r'Fidelity $|\langle\psi_{\rm exp}|\psi(\tau)\rangle|^2$', fontsize=11)
ax_fid_single.set_title('Single-swap fidelity', fontsize=13)
ax_fid_single.legend(fontsize=8, loc='lower right')
ax_fid_single.set_ylim(-0.02, 1.05)
ax_fid_single.grid(True, alpha=0.3)

# ── 子图 2: 双次 swap 保真度 ──────────────────────────────────────
ax_fid_double.axhline(y=1.0, color='gray', ls=':', lw=1, alpha=0.5)
ax_fid_double.annotate('ideal double-swap = 1', xy=(50, 0.97), fontsize=9, color='gray')
ax_fid_double.set_xlabel(r'$\tau$ (adiabatic parameter)', fontsize=12)
ax_fid_double.set_ylabel(r'Fidelity $|\langle\psi_{\rm exp}|\psi(2\tau)\rangle|^2$', fontsize=11)
ax_fid_double.set_title('Double-swap fidelity', fontsize=13)
ax_fid_double.legend(fontsize=8)
ax_fid_double.set_ylim(-0.02, 1.05)
ax_fid_double.grid(True, alpha=0.3)

# ── 子图 3: γ₂↔γ₃ 互换矩阵元 ──────────────────────────────────────
ax_swap.axhline(y=1.0, color='gray', ls=':', lw=1, alpha=0.5)
ax_swap.set_xlabel(r'$\tau$', fontsize=12)
ax_swap.set_ylabel('Matrix element magnitude', fontsize=11)
ax_swap.set_title(r'$\gamma_2\leftrightarrow\gamma_3$ swap check (single swap)', fontsize=13)
ax_swap.legend(fontsize=7, ncol=2, loc='lower right')
ax_swap.set_ylim(-0.02, 1.05)
ax_swap.grid(True, alpha=0.3)

# ── 子图 4: 特定 τ=50 时的旋转矩阵可视化 ───────────────────────────
tau_fixed = 50
E1_fixed, t1_fixed = 0.01, 0.005
R1_fixed, _ = run_protocol(tau_fixed, E1_fixed, t1_fixed)
im = ax_matrix.imshow(np.abs(R1_fixed[:3,:3]), cmap='Blues', vmin=0, vmax=1)
for i in range(3):
    for j in range(3):
        val = R1_fixed[i,j]
        color = 'white' if np.abs(val) > 0.5 else 'black'
        ax_matrix.text(j, i, f'{val:.3f}', ha='center', va='center', 
                       color=color, fontsize=9)
ax_matrix.set_xticks([0,1,2]); ax_matrix.set_yticks([0,1,2])
ax_matrix.set_xticklabels([r'$\gamma_1$', r'$\gamma_2$', r'$\gamma_3$'])
ax_matrix.set_yticklabels([r"$\gamma_1'$", r"$\gamma_2'$", r"$\gamma_3'$"])
ax_matrix.set_title(fr'MZM rotation matrix ($\tau$={tau_fixed}, $E_1$={E1_fixed}, $t_1$={t1_fixed})',
                    fontsize=12)
plt.colorbar(im, ax=ax_matrix, label='|matrix element|')

plt.tight_layout()
plt.savefig('general_fidelity.png', dpi=150)
print("\n\n图已保存为 general_fidelity.png")

# ── 也输出 γ₂↔γ₃ 是否成立 ──────────────────────────────────────
print("\n=== 双次 swap 时 γ₂↔γ₃ 是否互换？ ===")
for label, E1v, t1v, _, _ in param_sets:
    if label.startswith("MZM"):
        continue
    tau_test = [5, 20, 50, 100]
    print(f"\n{label}:")
    for tau in tau_test:
        try:
            R1, R2 = run_protocol(tau, E1v, t1v)
            print(f"  τ={tau}: R2[1,2]={R2[1,2]:.4f} (γ₂中γ₃分量), "
                  f"R2[2,1]={R2[2,1]:.4f} (γ₃中γ₂分量)")
        except:
            pass
