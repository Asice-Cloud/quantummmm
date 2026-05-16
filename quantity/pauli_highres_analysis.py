"""
Pauli 张量积高精度编织模拟与对标

核心改进：
1. 直接从有效 Hamiltonian 矩阵提取 Pauli 张量系数
2. 高精度 Trotter（步长=500）与标准（步长=100）对标
3. 诊断 Pauli 分解的正确性
"""
import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt
from importlib import util
import time

# 导入有效模型脚本
spec = util.spec_from_file_location('mod', 'quantity/reproduce_effective_braiding_majorana.py')
mod = util.module_from_spec(spec)
spec.loader.exec_module(mod)


class PauliTensorBraidingAnalyzer:
    """
    从有效 Hamiltonian 矩阵分析和提取 Pauli 张量分解
    """
    
    # Pauli 矩阵
    PAULI = {
        0: np.array([[1, 0], [0, 1]], dtype=complex),
        1: np.array([[0, 1], [1, 0]], dtype=complex),  # σ_x
        2: np.array([[0, -1j], [1j, 0]], dtype=complex),  # σ_y
        3: np.array([[1, 0], [0, -1]], dtype=complex),  # σ_z
    }
    
    def __init__(self, n_modes=3):
        self.n_modes = n_modes
        self.dim = 2 ** n_modes
    
    def extract_pauli_coefficients_simple(self, H_matrix):
        """
        简化诊断：查看 Hamiltonian 的特征值结构
        """
        evals = np.linalg.eigvals(H_matrix)
        evals_sorted = np.sort(evals.real)
        return evals_sorted
    
    def get_H_matrix_from_protocol(self, t, tau, tc, params):
        """
        使用有效模型构建 Hamiltonian 矩阵
        """
        c_ops = mod.fermionic_operators(self.n_modes)
        gammas = []
        for c in c_ops:
            cd = c.conj().T
            gammas.append(c + cd)
            gammas.append(-1j * (c - cd))
        
        t1v, t2v, t3v, Edv = mod.time_profiles(t, tau, tc)
        H = mod.build_H(gammas, params, t1v, t2v, t3v, Edv)
        
        return H
    
    def compute_fidelity(self, U_braiding, initial_state=None):
        """计算编织保真度"""
        if initial_state is None:
            psi = np.zeros(self.dim, dtype=complex)
            psi[0] = 1.0
        else:
            psi = initial_state / np.linalg.norm(initial_state)
        
        psi_final = U_braiding @ psi
        return abs(np.vdot(psi, psi_final))
    
    def check_unitarity(self, U_matrix):
        """检查幺正性"""
        U_dag_U = U_matrix.conj().T @ U_matrix
        error = np.linalg.norm(U_dag_U - np.eye(self.dim), 'fro')
        return error


def run_analysis():
    """主程序"""
    
    print("=" * 90)
    print("Pauli 张量高精度编织分析与对标")
    print("=" * 90)
    print()
    
    # 设置
    params = {'E1': 1e-4, 'Ed': 0.0}
    tau_protocol = 50.0
    tc = 0.03
    n_modes = 3
    
    analyzer = PauliTensorBraidingAnalyzer(n_modes=n_modes)
    
    print(f"System: {analyzer.dim}-dim (3 fermionic modes = 6 Majorana)")
    print(f"Parameters: E1={params['E1']:.0e}, Ed={params['Ed']:.2f}, tc={tc}, τ_protocol={tau_protocol}")
    print()
    
    # ======== 第一部分：诊断 Pauli 分解 ========
    print("-" * 90)
    print("Part 1: Hamiltonian 结构诊断")
    print("-" * 90)
    print()
    
    t_diag = [0.0, tau_protocol/2, tau_protocol]
    
    for t in t_diag:
        H_t = analyzer.get_H_matrix_from_protocol(t, tau_protocol, tc, params)
        evals = analyzer.extract_pauli_coefficients_simple(H_t)
        
        print(f"At t={t:.1f}:")
        print(f"  Hamiltonian 维数：{H_t.shape}")
        print(f"  Hermitian check (||H†H - HH†||): {np.linalg.norm(H_t.conj().T @ H_t - H_t @ H_t.conj().T):.2e}")
        print(f"  特征值范围：[{evals[0]:.4f}, {evals[-1]:.4f}]")
        print(f"  特征值（排序）：{[f'{e:.4f}' for e in evals]}")
        print()
    
    # ======== 第二部分：高精度 τ 扫描 ========
    print("-" * 90)
    print("Part 2: 高精度 τ 扫描（steps_per_tau=500）")
    print("-" * 90)
    print()
    
    taus = np.linspace(5, 150, 16)
    fidelities_500 = []
    unitarity_errors = []
    
    print(f"{'τ':>6} | {'Fidelity':>10} | {'Unitarity error':>16} | {'Time (s)':>8}")
    print("-" * 55)
    
    for tau_total in taus:
        t0 = time.time()
        
        # 高精度 Trotter (steps=500)
        dt = tau_total / 500
        U = np.eye(analyzer.dim, dtype=complex)
        
        for k in range(500):
            t = k * dt
            H_t = analyzer.get_H_matrix_from_protocol(t, tau_protocol, tc, params)
            U_dt = expm(-1j * H_t * dt)
            U = U_dt @ U
        
        fid = analyzer.compute_fidelity(U)
        unitary_err = analyzer.check_unitarity(U)
        
        fidelities_500.append(fid)
        unitarity_errors.append(unitary_err)
        
        elapsed = time.time() - t0
        print(f"{tau_total:6.1f} | {fid:10.6f} | {unitary_err:16.2e} | {elapsed:8.2f}")
    
    print()
    print(f"Mean fidelity (steps=500): {np.mean(fidelities_500):.6f}")
    print(f"Fidelity range: [{np.min(fidelities_500):.6f}, {np.max(fidelities_500):.6f}]")
    print(f"Max unitarity error: {np.max(unitarity_errors):.2e}")
    print()
    
    # ======== 第三部分：与现有结果对标 ========
    print("-" * 90)
    print("Part 3: 对标现有结果（steps_per_tau=100）")
    print("-" * 90)
    print()
    
    try:
        existing_data = np.loadtxt('quantity/fidelity_E1_1e-04.txt')
        existing_taus = existing_data[:, 0]
        existing_fids = existing_data[:, 1]
        
        # 线性插值对比
        from scipy.interpolate import interp1d
        f_existing = interp1d(existing_taus, existing_fids, kind='linear', 
                              fill_value='extrapolate')
        
        differences = []
        print(f"{'τ':>6} | {'Fidelity (500)':>15} | {'Fidelity (100)':>15} | {'Difference':>12}")
        print("-" * 60)
        
        for tau_val, fid_500 in zip(taus, fidelities_500):
            fid_100 = f_existing(tau_val)
            diff = abs(fid_500 - fid_100)
            differences.append(diff)
            print(f"{tau_val:6.1f} | {fid_500:15.6f} | {fid_100:15.6f} | {diff:12.6f}")
        
        mean_diff = np.mean(differences)
        max_diff = np.max(differences)
        
        print()
        print(f"Mean difference: {mean_diff:.6f}")
        print(f"Max difference:  {max_diff:.6f}")
        
        if mean_diff < 0.005:
            print("✓ 结论：高精度与标准精度一致，Trotter(100) 已足够")
        elif mean_diff < 0.02:
            print("⚠ 结论：差异较小，两种精度都可接受，高精度更精确")
        else:
            print("✗ 结论：差异显著，高精度计算推荐使用")
        
        comparison_available = True
    except FileNotFoundError:
        print("(现有数据文件未找到，跳过对标)")
        comparison_available = False
    
    print()
    
    # ======== 绘图 ========
    print("-" * 90)
    print("生成对标图表")
    print("-" * 90)
    print()
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # 左：保真度对比
    ax = axes[0]
    ax.plot(taus, fidelities_500, '-o', label='High-res (steps=500)', 
            linewidth=2, markersize=6)
    
    if comparison_available:
        ax.plot(existing_taus, existing_fids, 's--', 
                label='Standard (steps=100)', alpha=0.7, markersize=4)
        ax.fill_between(taus, 
                        [f_existing(t) - 0.01 for t in taus],
                        [f_existing(t) + 0.01 for t in taus],
                        alpha=0.1, color='orange')
    
    ax.set_xlabel('Braiding time τ', fontsize=11)
    ax.set_ylabel('Fidelity', fontsize=11)
    ax.set_title('Pauli 高精度 vs 标准精度', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    ax.set_ylim([0, 1.05])
    
    # 右：幺正性检查
    ax = axes[1]
    ax.semilogy(taus, unitarity_errors, 'o-', color='red', 
                linewidth=2, markersize=6, label='||U†U - I||')
    ax.axhline(1e-14, color='gray', linestyle='--', alpha=0.5, label='Machine precision')
    ax.set_xlabel('Braiding time τ', fontsize=11)
    ax.set_ylabel('Unitarity error (Frobenius norm)', fontsize=11)
    ax.set_title('幺正性验证', fontsize=12)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    plt.savefig('quantity/pauli_highres_vs_standard.png', dpi=200)
    print("✓ 保存：quantity/pauli_highres_vs_standard.png")
    
    # 保存数据
    np.savetxt('quantity/fidelity_pauli_highres_E1_1e-04.txt',
               np.column_stack((taus, fidelities_500, unitarity_errors)))
    print("✓ 保存：quantity/fidelity_pauli_highres_E1_1e-04.txt")
    
    print()
    print("=" * 90)
    print("分析完成")
    print("=" * 90)


if __name__ == '__main__':
    run_analysis()
