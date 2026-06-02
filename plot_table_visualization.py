#!/usr/bin/env python3
"""
绘制 solution.md §11 表格数据：旋转角 φ vs τ，四组参数对比。
解释 PRB105 "任意 Bloch 旋转"的物理含义。
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pi = np.pi
t_c = 0.3
E0 = 0.3

def fp(t, tau): return 0.5*(1.0+np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1.0-np.cos(pi*t/tau))

def build_A_step1(t, tau, E1v, t1v):
    A = np.zeros((5,5))
    A[0,1]=2*E1v; A[1,0]=-2*E1v
    t2=t_c*fm(t,tau); A[1,3]=-2*t2; A[3,1]=2*t2
    t1=t1v*fm(t,tau); A[0,4]=2*t1; A[4,0]=-2*t1
    Ed=E0*fp(t,tau); A[3,4]=2*Ed; A[4,3]=-2*Ed
    return A

def build_A_step2(t, tau, E1v, t1v):
    A = np.zeros((5,5))
    A[0,1]=2*E1v; A[1,0]=-2*E1v
    t2=t_c*fp(t,tau); A[1,3]=-2*t2; A[3,1]=2*t2
    t3=t_c*fm(t,tau); A[2,3]=2*t3; A[3,2]=-2*t3
    t1=t1v*fp(t,tau); A[0,4]=2*t1; A[4,0]=-2*t1
    return A

def build_A_step3(t, tau, E1v, t1v):
    A = np.zeros((5,5))
    A[0,1]=2*E1v; A[1,0]=-2*E1v
    t3=t_c*fp(t,tau); A[2,3]=2*t3; A[3,2]=-2*t3
    Ed=E0*fm(t,tau); A[3,4]=2*Ed; A[4,3]=-2*Ed
    return A

def propagate_step(build_fn, tau, E1v, t1v):
    def rhs(t,y):
        A=build_fn(t,tau,E1v,t1v); R=y.reshape(5,5); return (A@R).reshape(-1)
    y0=np.eye(5).reshape(-1)
    sol=solve_ivp(rhs,(0,tau),y0,rtol=1e-10,atol=1e-14,method='RK45')
    return sol.y[:,-1].reshape(5,5)

def run_protocol(tau, E1v, t1v):
    R1=propagate_step(build_A_step1,tau,E1v,t1v)
    R2=propagate_step(build_A_step2,tau,E1v,t1v)
    R3=propagate_step(build_A_step3,tau,E1v,t1v)
    return R3@R2@R1

def rotation_angle(R):
    tr = np.trace(R[:3,:3])
    phi = np.arccos(np.clip((tr-1)/2, -1, 1))
    return phi

# ── 精细扫描 ──────────────────────────────────────────────────────
tau_fine = np.linspace(1, 100, 150)

param_sets = [
    (r"$E_1=0,\ t_1=0$ (pure MZM)",     0.0,   0.0,   '#2c7bb6', '-', 1.4),
    (r"$E_1=0,\ t_1=0.01$ (ABS, $t_1$ only)", 0.0, 0.01, '#d7191c', '--', 1.4),
    (r"$E_1=0.01,\ t_1=0.005$ (ABS small)", 0.01, 0.005, '#fdae61', '-.', 1.4),
    (r"$E_1=0.3,\ t_1=0.01$ (ABS large)", 0.3, 0.01, '#5e3c99', ':', 1.8),
]

fig, axes = plt.subplots(2, 2, figsize=(15, 11))

# ── 子图 1: 四组 φ vs τ ───────────────────────────────────────────
ax1 = axes[0,0]
ax1.axhline(y=pi/2, color='gray', ls=':', lw=1, alpha=0.5)
ax1.annotate(r'$\pi/2$ (pure geometric braid)', xy=(55, pi/2+0.03), fontsize=9, color='gray')

for label, E1v, t1v, color, ls, lw in param_sets:
    phis = []
    for tau in tau_fine:
        R = run_protocol(tau, E1v, t1v)
        phis.append(rotation_angle(R))
    ax1.plot(tau_fine, phis, color=color, ls=ls, lw=lw, label=label)

ax1.set_xlabel(r'$\tau$ (adiabatic parameter)', fontsize=12)
ax1.set_ylabel(r'Rotation angle $\phi$ (rad)', fontsize=12)
ax1.set_title(r'SO(3) rotation angle $\phi$ vs $\tau$', fontsize=14)
ax1.legend(fontsize=8.5, loc='upper left')
ax1.set_ylim(0, pi+0.1)
ax1.set_yticks([0, pi/4, pi/2, 3*pi/4, pi])
ax1.set_yticklabels(['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])
ax1.grid(True, alpha=0.3)

# ── 子图 2: 四组 φ 的直方图 (展示唯一值分布) ──────────────────────
ax2 = axes[0,1]
tau_coarse = np.linspace(10, 100, 10)
bins = np.linspace(0, pi, 25)
colors_hist = ['#2c7bb6', '#d7191c', '#fdae61', '#5e3c99']

for idx, (label, E1v, t1v, color, ls, lw) in enumerate(param_sets):
    phis = []
    labels_tau = []
    for tau in tau_coarse:
        R = run_protocol(tau, E1v, t1v)
        phis.append(rotation_angle(R))
        labels_tau.append(f'{tau:.0f}')
    ax2.hist(phis, bins=bins, alpha=0.5, color=colors_hist[idx], 
             label=label.split('(')[0].strip(), edgecolor='white')
    
    # mark individual points below
    y_offset = -0.12 - idx*0.06
    ax2.scatter(phis, [y_offset]*len(phis), color=colors_hist[idx], 
                marker='|', s=100, alpha=0.8)

ax2.axvline(x=pi/2, color='gray', ls=':', lw=1, alpha=0.5)
ax2.set_xlabel(r'$\phi$ (rad)', fontsize=12)
ax2.set_ylabel('Count (10 tau points)', fontsize=12)
ax2.set_title('Distribution of rotation angles (tau=10:10:100)', fontsize=13)
ax2.legend(fontsize=8)
ax2.set_xticks([0, pi/4, pi/2, 3*pi/4, pi])
ax2.set_xticklabels(['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])
ax2.grid(True, alpha=0.3, axis='y')

# ── 子图 3: 旋转轴分量演化 (以 E₁=0, t₁=0.01 为例) ────────────────
ax3 = axes[1,0]
tau_axis = np.linspace(10, 100, 60)
nz_vals = []
phi_vals_axis = []
for tau in tau_axis:
    R = run_protocol(tau, 0.0, 0.01)
    # Extract rotation axis from SO(3) part
    R3 = R[:3,:3]
    # Rotation axis is eigenvector with eigenvalue 1
    # Or from skew-symmetric part: n ∝ (R - R^T)
    skew = R3 - R3.T
    nx = skew[2,1]; ny = skew[0,2]; nz = skew[1,0]
    norm = np.sqrt(nx**2 + ny**2 + nz**2)
    if norm > 1e-12:
        nz_vals.append(np.abs(nz/norm))
    else:
        nz_vals.append(0)
    phi_vals_axis.append(rotation_angle(R))

ax3.plot(tau_axis, phi_vals_axis, 'C0-', lw=2, label=r'$\phi$')
ax3.set_xlabel(r'$\tau$', fontsize=12)
ax3.set_ylabel(r'$\phi$ (rad)', fontsize=12, color='C0')
ax3.tick_params(axis='y', labelcolor='C0')

ax3b = ax3.twinx()
ax3b.plot(tau_axis, nz_vals, 'C3--', lw=2, label=r'$|\hat n_z|$')
ax3b.set_ylabel(r'$|\hat n_z|$ (axis tilt)', fontsize=12, color='C3')
ax3b.tick_params(axis='y', labelcolor='C3')
ax3b.set_ylim(-0.02, 1.02)

lines1, labels1 = ax3.get_legend_handles_labels()
lines2, labels2 = ax3b.get_legend_handles_labels()
ax3.legend(lines1+lines2, labels1+labels2, fontsize=9, loc='lower right')
ax3.set_title(r'$E_1=0,\ t_1=0.01$: rotation angle + axis tilt', fontsize=13)
ax3.grid(True, alpha=0.3)

# ── 子图 4: Bloch 球示意图 (物理含义解释) ─────────────────────────
ax4 = axes[1,1]
ax4.set_xlim(-1.3, 1.3); ax4.set_ylim(-1.3, 1.3)
ax4.set_aspect('equal')
ax4.axis('off')
ax4.set_title('Physical interpretation on Bloch sphere', fontsize=13, y=1.02)

# Draw circle
theta_c = np.linspace(0, 2*pi, 200)
ax4.plot(np.cos(theta_c), np.sin(theta_c), 'gray', lw=1, alpha=0.3)

# Draw axes
for angle, label, color in [(0, r'$x$ ($\sigma_x$)', 'C0'), 
                              (pi/2, r'$y$ ($\sigma_y$)', 'C2'),
                              (pi, r'$z$?', 'gray')]:
    pass
ax4.arrow(0, 0, 1.1, 0, head_width=0.05, head_length=0.05, fc='C0', ec='C0', lw=1.5)
ax4.arrow(0, 0, 0, 1.1, head_width=0.05, head_length=0.05, fc='C2', ec='C2', lw=1.5)
ax4.text(1.2, -0.08, r'$x$', fontsize=11, color='C0')
ax4.text(-0.08, 1.2, r'$y$', fontsize=11, color='C2')

# Pure MZM: fixed at x-axis
ax4.annotate('', xy=(1,0), xytext=(0,0),
            arrowprops=dict(arrowstyle='->', color='#2c7bb6', lw=3))
ax4.text(0.55, 0.18, r'Pure MZM: $\hat n\approx\hat x$', fontsize=9, color='#2c7bb6')
ax4.text(0.55, 0.05, r'$\phi\approx\pi/2$ fixed', fontsize=9, color='#2c7bb6')

# ABS: axis tilts, angle grows
for angle_idx, angle in enumerate([pi/8, pi/5, pi/3.5]):
    ax_n = angle + pi/2  # axis in xz plane → appears in upper half
    ax4.annotate('', xy=(0.85*np.cos(ax_n), 0.85*np.sin(ax_n)), xytext=(0,0),
                arrowprops=dict(arrowstyle='->', color='#d7191c', lw=1.5, alpha=0.5+0.15*angle_idx))

ax4.text(-1.0, 0.55, r'ABS ($t_1\neq 0$):', fontsize=10, color='#d7191c')
ax4.text(-1.0, 0.38, r'$\hat n$ tilts toward $z$', fontsize=9, color='#d7191c')
ax4.text(-1.0, 0.23, r'$\phi$ grows with $\tau$', fontsize=9, color='#d7191c')
ax4.text(-1.0, 0.08, r'$\to$ covers Bloch sphere', fontsize=9, color='#d7191c')

# Annotations for the table
ax4.text(-1.0, -0.3, 'Key insight:', fontsize=11, fontweight='bold')
ax4.text(-1.0, -0.55, r'$E_1=0,t_1=0$: $\phi$ pinned at $\pi/2$', fontsize=9)
ax4.text(-1.0, -0.75, r'$E_1=0,t_1\neq 0$: $\phi$ scans $[\pi/2,\infty)$', fontsize=9)
ax4.text(-1.0, -0.95, r'$E_1\neq 0,t_1\neq 0$: full SU(2) released', fontsize=9)

plt.tight_layout()
plt.savefig('prb105_table_visualization.png', dpi=150)
print("Figure saved: prb105_table_visualization.png")
