**XYZ / 8‑vertex：符号约束与 ABS 的严格链**

目标：对 XYZ/8‑vertex 型 R(u) 在 $u=0$ 处的线性化导数，求出哪些参数约束可使局域两格 Pauli 展开满足 free‑fermion（二次 BdG）条件，并说明在满足这些约束时如何得到非零配对项 $\Delta$ 并据此产生 ABS。

1) XYZ 线性化与 Pauli 系数（重述）
- 采用常见的 8‑vertex/XYZ 近似参数化（线性化）：写
$$
a(u)=A+a_1 u,\,\; b(u)=B+b_1 u,\; c(u)=C+c_1 u,\; d(u)=D+d_1 u,$$
并在基序 $|00\rangle,|01\rangle,|10\rangle,|11\rangle$ 下构造
$$
R(u)=\begin{pmatrix}a&0&0&d\\0&b&c&0\\0&c&b&0\\d&0&0&a\end{pmatrix}.
$$ 
对 $u$ 求导在 $u=0$ 得到 $\partial_u R|_{0}$ 的非零元为 $a_1,b_1,c_1,d_1$（按位置分别对应上文符号）。

- 将局域生成元取为常用规范
$$
h = P\frac{\partial_u R|_{0}}{\rho},
$$
并把 $h$ 在 Pauli 基上展开，得到（仅列非零分量）：
$$
    \begin{align*}
	c_{II} &= \frac{a_1+c_1}{2\rho},\\
	c_{XX} &= \frac{b_1+d_1}{2\rho},\\
	c_{YY} &= \frac{b_1-d_1}{2\rho},\\
	c_{ZZ} &= \frac{a_1-c_1}{2\rho}.
	\end{align*}
$$
对应的 BdG 映射（按本文第 5 节的线性关系）给出：
$$
	\begin{align*}
	t &= c_{XX}+c_{YY} = \frac{b_1}{\rho},\\
	\Delta &= c_{XX}-c_{YY} = \frac{d_1}{\rho},\\
	\mu &= -2(c_{Z0}+c_{0Z}) = 0.
	\end{align*}
$$

2) free‑fermion 约束（简化代数条件）
- 我们之前对一般两格 Pauli 系数做了 Groebner 消元（见 `tools/derive_constraints.py` 与 `results/free_fermion_constraints.txt`），得到若要把 $h$ 精确写成二次 BdG 必须满足的一组代数零条件。将这些通用约束代入本节的 XYZ 展开后，得到的必要条件可写为：

- 要求 $c_{ZZ}=0$（因为通用消元将 $c_{zz}$ 的实部列为必须为零的分量之一）。代入上式得到
$$
\frac{a_1-c_1}{2\rho}=0\quad\Longrightarrow\quad a_1 = c_1. \tag{C1}
$$
- 其余非零的允许分量（在消元允许的子空间内）可为 $c_{II}, c_{XX}, c_{YY}$（这些并不违背二次条件，且虚部应为零）。

因此对 XYZ 型 R，满足 free‑fermion（二次）條件的**代数约束**为 (C1)（并假定所有导数参数为实数，以满足虚部为零的条件）。

3) 由约束导出的 Δ 与 BdG 参量
- 在 (C1) 成立时，BdG 参量化为
$$
t=\frac{b_1}{\rho},\qquad \Delta=\frac{d_1}{\rho},\qquad \mu=0.$$ 
因此要产生非零配对 $\Delta\neq0$，只需保证
$$
d_1\neq0. \tag{C2}
$$

4) ABS 存在的判据（在此约化模型下的具体操作）
- 给定分段链（左段参数集 $(b_1^L,d_1^L)$，右段 $(b_1^R,d_1^R)$），在满足 (C1) 的前提下，可把左右段各自写成 p‑wave BdG：
$$
t_{L,R}=b_1^{L,R}/\rho,\quad \Delta_{L,R}=d_1^{L,R}/\rho,\quad \mu_{L,R}=0.
$$ 

- 按第 7 节的方法，寻找子隙内能量 $E$，解左右两段对 (*) 的 $X$ 根并选取 $X>1$ 的衰减根，再在界面上建立匹配矩阵 $M(E)$；若存在 $E$ 使 $\det M(E)=0$，则存在指数衰减的界面束缚态（ABS）。

- 直观条件：若左段具有显著配对 $\Delta_L\neq0$ 而右段配对小或为零（$\Delta_R\approx0$），且两侧体隙在某能量窗口重叠，则常常可找到满足上述匹配的 $E$，从而出现 ABS（这与我们用数值强制配对得到的示例一致）。

5) 示例（代数与数值提示）
- 代数示例：取 $\rho=1$，令 $a_1=c_1$（满足 (C1)），选 $b_1=1.0, d_1=0.5$，则
$$
t=1.0,\quad \Delta=0.5,\quad \mu=0.$$ 
用这些参数构造左右分段（例如左段 $\Delta_L=0.5$，右段 $\Delta_R=0$）并按本文第 7 节的匹配步骤，可以构造并求解 $\det M(E)=0$ 以定位 ABS。

- 数值建议：用 `tools/xyz_symbolic.py` 得到参数映射后，可把选定参数代入 `tools/pulse_abs.py` 或直接构造 BdG 并用 `numpy.linalg.eigh` 对角化来寻找束缚态；对候选本征态再做密度拟合得到衰减长 $\xi$ 以确认指数局域性。

6) 结论（归纳）
- XYZ/8‑vertex 家族因其更丰富的矩阵结构允许在 $\partial_u R|_0$ 的适当矩阵元（这里为 $d_1$）非零，从而在满足简单代数约束 $a_1=c_1$ 的情况下产生非零配对 $\Delta$ 并使局域两格 Hamiltonian 落入 free‑fermion 可表示的子空间。满足这些代数条件并在空间上构造配对差异（例如左段有配对、右段无配对）的情形，可严格导出 ABS 的判据并在数值上验证。

 

文件与脚本：已使用 `tools/xyz_symbolic.py` 做符号展开，`tools/derive_constraints.py` 给出了通用的 free‑fermion 消元约束；可以按需把上面的代数示例代入并运行 `tools/pulse_abs.py` / `tools/ed_extra_scan.py` 做数值验证并把结果追加到此文档。

**数值验证（示例演示）**

- 参数：`a1=c1=0.0`, `b1=1.0`, `d1=0.5`, `rho=1.0` → 对应单体映射为 `t=1.0`, 左段 `Delta_L=0.5`, 右段 `Delta_R=0.0`, `mu=0`。
- 运行 `tools/run_xyz_abs_demo.py`（构造 N=200 的 S–N 链并对角化）得到最低能级（按绝对值排序）显示：
	-2.621642e-16, 2.861543e-16, -3.116653e-02, 3.116653e-02, -6.227385e-02, 6.227385e-02
- 最低本征态密度 `results/xyz_demo_density0.txt` 的峰值位于格点 `i=199`，峰值约 `0.01982575`。用简单的右侧 1/e 判据无法在右侧找到 1/e 交叉点（因为峰在边界处），因此对该运行的简单衰减长度估计为 NaN（边界局域）。

- 说明：在这次示例里，最低能量态在系统右侧边界强烈局域化（不是严格位于中间界面）。这说明在实际数值检验时需要：
	- 确认界面位置（`tools/pulse_abs.build_chain_from_params` 中 interface_width 与分段索引），
	- 试验不同的 `interface_width` 或更中心化的界面参数，
	- 或使用较小系统更仔细地查看界面束缚态是否出现并定位其衰减长度 $\xi$。

	我们已把谱图保存为 `results/xyz_demo_spectrum.png`，密度及报告保存为 `results/xyz_demo_density0.txt` 与 `results/xyz_demo_report.txt`。

	**数值测试（右段带隙，μ_right=3.0）**

	- 目的：通过把右段设置为带隙平凡相（`Delta_right=0`, `mu_right=3.0`），检验是否在左配对/右带隙界面产生局域 ABS。
	- 参数（运行脚本）：`N=200`, `t_left=1.0, Delta_left=0.5, mu_left=0.0`, `t_right=1.0, Delta_right=0.0, mu_right=3.0`, `interface_width=4`。
	- 关键数值输出（已保存为 `results/try_mu3_vals.txt` 和 `results/try_mu3_density0.txt`）：
		- 最低绝对值本征能： ~1.53e-16, -4.33e-16（近零）
		- 接下来的若干对称能级： ±1.00096, ±1.00146, ±1.00149 …
		- 密度峰：`peak_idx=98`, `peak_val≈0.6317`，简单 1/e 估计 `xi_est=1`（非常局域）。

	- 结论：将右段变为带隙平凡相成功在界面附近产生了一个非常局域的近零态（峰值位于链中部附近 idx≈98，密度峰值 ≈0.63），说明当左右两段在拓扑性质上不同时，XYZ‑derived 配对确实能产生界面束缚态（ABS）。

	文件：`results/try_mu3_vals.txt`, `results/try_mu3_density0.txt`（谱与密度），结果已保存并可以用于后续拟合衰减长度或 Pauli‑Γ 展开分析。

	**MZM 证据汇总（当前最可信的零模）**

	我们用以下数值检验判断产生的近零态是否具有 Majorana 特征：

	- Majorana‑ness：`results/majorana_check.txt` 给出 `majorana_score=2.209885e-15`（≈0），以及 `abs_overlap_u_conjv=0.5`（规范相关），说明对最低模 `u ≈ conj(v)` 成立 —— 符合 Majorana 波函数的形态。 
	- 尺度不变性：`results/scale_summary.txt`（N=200,400,800）显示最低绝对能级分别为 ≈1.53e-16、≈5.6e-16、≈3.98e-17（均接近数值零，且随 N 不呈系统性漂移），这表明该零模随尺寸稳定。 
	- 参数鲁棒性：`results/mu_scan.txt` 中 μ_right 在 [2.5,3.5] 的扫描显示最低绝对能量保持在 ~1e-16–1e-15 量级，对 μ 的小幅变化稳健。 
	- 局域性（LDOS）：`results/try_mu3_ldos.png` 显示在界面处有高峰（`peak_val≈0.63, peak_idx≈98`），且对数图 `results/try_mu3_ldos_log.png` 表明密度在界面外指数小。 

	综合以上：在我们对 XYZ‑derived 参数的构造（左段配对 Δ>0，右段带隙平凡 μ_right large）下，产生的近零态同时满足 Majorana 形态、尺度稳定性与参数鲁棒性，故判定为 Majorana‑like 零模（MZM），而非脆弱的 ABS。

	**下一步：自动搜索 ABS 参数集（说明）**

	为了找到明显的 ABS（而非 MZM），我们实现了自动化搜索（脚本 `tools/search_abs.py`）去扫描“局域配对陷阱”类型的参数，即在导电（无配对）背景中嵌入短而强的配对区域，并判断候选态的“Majorana‑ness”（应较大）和对小扰动的敏感性（能量易漂移）。搜索结果保存在 `results/abs_search.json`（若找到候选点）。



