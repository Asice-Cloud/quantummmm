## verify 目录脚本索引与对应验证内容

### 一、R→Kitaev 映射（1D 有效模型）

- run:  kit1/
	- verify/majorana_braid_check.py
	- verify/build_micro_BdG.py
	- verify/run_micro_mapping_fit.py
	- verify/run_mapping_validation.py
- 验证内容：
	- 明确 R(a,b,c,d) 到 2-site Kitaev 链参数的解析映射：
		- t = b + c
		- Δ = b − c
		- μ = 4 d + μ₀
	- 用 BdG 显微链构造 build_micro_BdG 并抽取 (t_eff, Δ_eff, μ_eff)，数值拟合出与 (b+c, b−c, d) 的线性关系；
	- 在纯 2-site BdG 模型中再次扫描 (a,b,c,d)，验证解析映射与数值完全一致；
	- 对应文档位置：R_to_Kitaev.md、kit3-2.md、kit3-3.md、kit_all.md 第 2 章（2.2–2.3）。

### 二、1D Kitaev 链的拓扑与零模标度

- run:  kit1/
	- verify/run_spectrum_and_topology.py
	- verify/run_topology_and_scaling.py
- 验证内容：
	- 计算动量空间能隙 gap_kspace(t,Δ,μ)，并在开边界链上数零能态数目，确认 |μ| < 2|t| 区间为拓扑相、有端点零模；
	- 计算拓扑 winding 数与 Pfaffian 符号 vs μ，并在拓扑/平庸两点做零模能隙随链长 N 的指数标度拟合，提取关联长度 ξ；
	- 对应文档位置：kit3-4.md、kit4-1.md 的 warm-up、kit_all.md 第 2 章（2.3）关于 1D bulk–edge story 的图与数据。

### 三、R-矩阵的 YBE 约束与 R 空间几何

- run:  kit1/
	- verify/2.py  → 生成 verify/ybe_eqs.txt
	- verify/3.py
- 验证内容：
	- 用 Sympy 把 Pauli 形式的 R(a,b,c,d) 代入常数 YBE：R₁₂R₁₃R₂₃ = R₂₃R₁₃R₁₂，展开得到关于 (a,b,c,d) 的多项式方程组 ybe_eqs.txt；
	- 对这些多项式做因式分解/求解，辅助解析分类不同 YBE 解族（Ising-like、XY-like 等）；
	- 对应文档位置：kit3-5.md、kit_all.md 第 2.4 节（Hilbert 丛、平坦联络与 YBE 子流形）。

### 四、四-Majorana 玩具模型：half twist、Dehn twist 与 LQC 复杂度

- run:  kit2/
	- verify/run_instantaneous_braid_crosscheck.py
	- verify/run_dehn_twist_like_path.py
	- verify/run_dehn_twist_micro_vs_berry.py
	- verify/run_lqc_permutation_fit.py
	- verify/run_complexity_flatness_scan.py
	- verify/run_interacting_braid_check.py
- 验证内容：
	- 用 R→Kitaev 映射构造 4-Majorana Hamiltonian H，时间演化 U_time = exp(−i H τ)，与理想 braid 生成元 U_ideal = exp((π/4) γ₂γ₃) 精确对比（instantaneous braid cross-check）；
	- 在耦合空间构造闭合路径 H(φ)、H(θ₁,θ₂)，保证两维基态子空间始终存在，并通过相邻基向量重叠计算非阿贝尔 Berry holonomy（Dehn-twist-like path）；
	- 构造 half twist U_half = exp((π/4) γ_i γ_j)，将 U_half² 投影到几何基态子空间，并与几何 Berry holonomy 以及 Ising TQFT 的 (R^{σσ})² 做 SU(2) 共轭比较（micro vs Berry vs Ising Dehn twist）；
	- 把 4-Majorana 看成两比特，扫描浅层 LQC+SWAP 电路 (U₁ₐ⊗U₂ₐ)·SWAP·(U₁_b⊗U₂_b)，投影到基态子空间后与 U_Berry 比较，得到最佳 fidelity F，说明浅层电路可以实现该 Dehn-twist-like 操作（LQC+permutation fit）；
	- 在 H(φ) 上加入扰动 ε·(i/2)γ₁γ₂，得到 H(φ;ε)，对每个 ε 计算 U_Berry(ε)，用相同浅层 ansatz 拟合并得到 F(ε)，把 1−F(ε) 视为“曲率/复杂度” 指标，画出 complexity_flatness_F_vs_eps.png（complexity/flatness scan）；
	- 在 4-Majorana Fock 空间中加入真正的相互作用项 U n₁n₂ 和化学势 μ，计算 U_full 与 U_ideal 的差异以及在 (μ,U) 平面上的误差热图，考察相互作用对 braid 的稳健性（interacting braid check）；
	- 对应文档位置：kit3-5-2.md、kit3-5-3.md、kit3-5-4.md、kit_all.md 第 2.4–2.5 节和结论部分中关于 “从 R 到 half twist/Dehn twist 以及复杂度层”的讨论。

### 五、2D p+ip 模型：vortex Berry 与 Dehn twist 平台

- run:  kit3/
	- verify/run_pip_vortex_berry.py
	- verify/run_pip_vortex_scan.py
- 验证内容：
	- 在有限大小方格上构造 2D spinless p+ip BdG Hamiltonian，插入两个涡旋，沿参数 λ ∈ [0,1] 的路径让一个涡绕另一个涡一圈，提取零模子空间并计算 Berry holonomy U_Berry；
	- 在一系列 μ 上重复 vortex loop 实验，计算每个 μ 的 F_Dehn(μ)，即几何 U_Berry 与 Ising Dehn twist (R^{σσ})² 在 SU(2) 中的一致程度，并画出 F_Dehn(μ) 的平台结构；
	- 对应文档位置：kit4-1.md、kit4-3.md、kit_all.md 第 3 章（2D p+ip 涡旋 Dehn twist 与平台图）。

### 六、Honeycomb 模型：相图与 vison Berry

- run:  kit3/
	- verify/plot_fk_phase.py
	- verify/run_honeycomb_vison_berry.py
- 验证内容：
	- 计算 honeycomb 模型的 f(k) 并在 Jx+Jy+Jz=1 的简单形区域上扫描，标出 gapless/gapped 区域，画出 phase_diagram.png，同时在各向同性点 Jx=Jy=Jz=1/3 上画 |f(k)| 分布图（f_isotropic.png），作为 honeycomb 相图的基准；
	- 在 brick-wall 几何上构造 Majorana+Z₂ gauge 模型，枚举键与六边形 plaquette，给出一条简单的 Z₂ “vison loop” gauge 路径，并在该路径上计算低能子空间的 Berry holonomy，作为 honeycomb vison Berry 实验的脚手架；
	- 对应文档位置：kit4-2.md、kit4-3.md、kit4-4.md、kit_all.md 第 4 章（honeycomb 相图与 vison-loop Berry 比较）。

### 七、辅助/说明性文件

- verify/2.md
	- 说明性文本，用于早期 1D/Kitaev 练习与推导草稿，不参与主线数值扫描；
	- 可以作为阅读脚本时的注释补充。

