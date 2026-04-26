**从 R(u) 到 Andreev Bound States（ABS）的完整代数推导**

目标：给出如何从仓库中的两格 R(u)（谱参数的 YBE R 矩阵）严格构造出局域两体 Pauli 哈密顿、把它映射到费米 BdG 描述，并得到判断界面/短结产生 ABS 的代数判据与可执行数值流程。

一、概览与符号约定
- 局域两格哈密顿以 Pauli 基表示：对第 j 与 j+1 格定义
$$
	h_j \,=\,\sum_{\alpha,\beta\in\{0,x,y,z\}} c^{(j)}_{\alpha\beta}\;\sigma^\alpha_j\otimes\sigma^\beta_{j+1},
$$
其中 $\sigma^0=I$。我们的脚本（例如 [tools/xxz_R_and_H.py](tools/xxz_R_and_H.py)）按此 convention 展开 $4\times4$ 局域算符到 Pauli 系数 $c_{\alpha\beta}$。

二、从谱参数 R(u) 导出局域生成元
- 给定两格 R(u)，在 $u=0$ 处作导数并以 $R(0)$ 的尺度归一化，常用定义
$$
	h^{(R)}\;=\;P\cdot\frac{\partial_u R(u)\big|_{u=0}}{\rho},
$$

其中 $P$ 是两格置换算子，$\rho$ 是常数尺度（例如对三角 XXZ 族取 $\rho=\sin\eta$）。另一个常见生成元为
$$
	H_P(u) = -i R(u)^{-1}\partial_u R(u),
$$
但请注意数值上 $H_P(u)$ 在实谱参数上未必 Hermitian —— 需要按幺正门的需求选用适当规范（例如以 $R(0)$ 归一化或取 $u\to i u$）。

三、Pauli 系数到费米算符（Jordan–Wigner）
- 通过 Jordan–Wigner (JW) 变换，局域 Pauli 两体项可展开为费米算符的二次项与四次项：
	- 典型二次项（产生单粒子与配对）： $c_j^\dagger c_{j+1},\; c_j c_{j+1},\; c_j^\dagger c_j$。
	- 四费米项来自例如 $\sigma^z_j\sigma^z_{j+1}$ 等，会产生交互项 $n_j n_{j+1}$。
- 因此，写出 JW 展开并把系数按幂次分组：
$$
	H_{\text{JW}} = H^{(2)} + H^{(4)} + \cdots
$$
其中 $H^{(2)}$ 为所有二次项的线性组合（可写成 BdG 格式），$H^{(4)}$ 为四费米交互。

四、free‑fermion（BdG）可表示的代数条件
- 要能用精确的单粒子 BdG 描述，必须消去所有四费米项（或能通过变量替换/恒等式把它们并入常数）。这给出关于 $c_{\alpha\beta}$ 的一组多项式约束：
$$
	f_i(\{c_{\alpha\beta}\}) = 0,\quad i=1,\dots,m.
$$
在仓库中有 SymPy 脚本（见 [tools/derive_constraints.py](tools/derive_constraints.py)）用于类似消元工作。满足这些约束的局域 $h_j$ 才能精确对应某个 BdG 的二次形式。

五、构造 BdG 单粒子矩阵（若满足二次条件）
- 把所有二次项整理为一般 BdG 块表示（N 站点）：
$$
	H_{\mathrm{BdG}}=\begin{pmatrix} A & B\ -B^* & -A^T \end{pmatrix},
$$
其中 $A$ 来源于粒子型项（含跃迁 $t_j$ 与化学位 $\mu_j$），$B$ 来源于配对项 $\Delta_j$。具体映射可按脚本 [tools/verify_mzm.py](tools/verify_mzm.py) 中的公式执行：
$$
	t\sim c_{xx}+c_{yy}+i(c_{xy}-c_{yx}),\quad \Delta\sim c_{xx}-c_{yy}-i(c_{xy}+c_{yx}),\quad \mu\sim 4c_{zz}-2(c_{z0}+c_{0z}).
$$

六、均匀段的色散与子隙
- 在平移不变段代入 Bloch 解，得到模方程（near‑neighbor p‑wave 约定举例）
	$$
	E(k)^2 = \varepsilon(k)^2 + |\Delta(k)|^2,
	$$
	其中 $\varepsilon(k)=-\mu-2t\cos k,\;\Delta(k)=2i\Delta\sin k$。局域“体隙”定义为 $\min_k |E(k)|$。

七、束缚态的代数判据（界面问题）
- 对界面束缚态，寻找复波数 $k=k_r+i\kappa$ 使波函数在远端指数衰减。代入 $k=i\kappa$（或 $\pi+i\kappa$）得到关于 $X=\cosh\kappa$ 的代数方程（经代数整理得二次形式，p‑wave 近邻模型为例）：
$$
	(t^2+|\Delta|^2) X^2 + (\mu t) X + \frac{\mu^2-4|\Delta|^2-E^2}{4} = 0.\tag{*}
$$
- 给定候选能量 $E$（|E| 小于两侧能隙），若该二次方程对左右两段均有解且该解 $X>1$（即 $\kappa=\operatorname{arccosh}X>0$），则存在指数衰减的远端解。
- 接着在界面上写出有限维匹配方程（把左右衰减解在界面局部格点插入 BdG 差分方程），得到接口矩阵 $M(E)$；存在非平凡解等价于
$$
	\det M(E)=0.\tag{**}
$$
因此 ABS 的严格判据为：存在 $E$（在两个体隙内）满足 (*) 在左右两段都有 $X>1$ 的根，且 (**) 成立。

八、谱参数 R(u) 与含时脉冲的角色
- 若 R(u) 在小 u 处满足 $R(u)\approx R(0) e^{i u H_P}$（或规范化后等价），则可以把一段谱参数变化视为施加局域两体门 $e^{i\delta u H_P}$。通过在时间或空间上设计不同的 R(u) 脉冲序列，可以在左/右段产生不同的 Pauli 系数 $c^{(j)}_{\alpha\beta}(t)$，进而在平均或瞬时意义上生成非零 $\Delta_j(t)$ 与 $\mu_j(t)$。
- 注意数值/代数细节：直接用 $H_P(u)=-iR^{-1}\partial_u R$ 作为时间生成元时要保证 Hermiticity（我们的 XXZ 示例在实 u 下直接得出的 $H_P$ 可能非 Hermitian；可用尺度化 $h=P\partial_u R/\rho$ 或取虚谱参数 $u\mapsto i u$ 来得到 Hermitian 生成元）。

九、数值实现步骤（可直接运行）
1. 用 `tools/xxz_R_and_H.py` 或其它 R(u) 脚本得到 $h^{(R)}$，并用 `expand_on_pauli` 得到 $c_{\alpha\beta}$。
2. 检查 free‑fermion 约束（或人为施加 mean‑field 分解）。若不满足，可把 $H^{(4)}$ 做平均场近似或在脉冲设计中直接加入配对项。脚本参考：`tools/derive_constraints.py`。
3. 从 $c$ 计算 $t_j,\Delta_j,\mu_j$（见 `tools/verify_mzm.py:map_c_to_params`）。
4. 构造空间上分段的 BdG（左段有非零 $\Delta_L$、右段为常导 $\Delta_R\approx0$），用数值对角化得到本征态与能级（脚本：`tools/pulse_abs.py`）。
5. 对候选子隙能量用方程 (*) 解出 $X_{L,R}$ 并选 $X>1$ 的根，构造界面匹配矩阵 $M(E)$，计算 $\det M(E)$ 并求根（可用 SymPy 生成小宽度接口的显式多项式）。
6. 做数值验证：检查粒‑洞对称性、波函数指数拟合（见 `tools/abs_checks.py`），以及 Josephson 相位依赖拟合（短结公式）以确认 ABS 的物理行为。

十、与交互/拓扑零模的区分
- 若得到 $E\approx0$ 的模式，要区分 MZM（拓扑保护）与平凡 ABS：对系统尺寸 $N$ 做收敛测试（MZM 的 E 随 N 指数衰减到 0），对局域扰动测试（MZM 鲁棒，ABS 易漂移），并检查拓扑不变量（对可解模型）。

参考实现与示例
- 参见仓库脚本：
	- [tools/xxz_R_and_H.py](tools/xxz_R_and_H.py) — 从 XXZ R(u) 导出局域 h 并展开到 Pauli。 
	- [tools/pulse_abs.py](tools/pulse_abs.py) — 用 R 映射得到参数并构造 S–N 链与 Floquet 平均示例。
	- [tools/abs_checks.py](tools/abs_checks.py) — PH 对称、局域拟合、Josephson 扫描等自动检查工具。

结语：上述推导把“从 R(u) → Pauli h → JW → BdG → ABS 条件”连成了一个明确可执行的链。关键的严格性点在于是否能把局域 Pauli 展开精确写成二次费米形式（即 free‑fermion 约束）；若不满足则需施行 mean‑field 或在脉冲设计中外加显式配对以产生 ABS。若你要，我可以把第七步中接口矩阵 $M(E)$ 在“1‑site/2‑site 界面宽度”情形用 SymPy 符号化写出并求出 $\det M(E)$ 的显式多项式样例。

**附：从 Pauli 系数到 Majorana 模算符的逐步代数推导（可直接实现）**

下面给出完整、可直接编码的代数步骤，便于把数值得到的本征向量转成 Pauli 串表示的 Majorana 算符 $\Gamma$：

- 约定：链上有 $N$ 个物理格点，费米湮灭算符为 $c_j$ ($j=1\ldots N$)。Nambu 向量按 $\Psi=(c_1,\ldots,c_N,c_1^\dagger,\ldots,c_N^\dagger)^T$ 排列，BdG 本征问题为 $H_{\mathrm{BdG}}\Phi=E\Phi$。

- 步骤 1 — 从局域 Pauli 系数构建 BdG（若满足二次条件）：
	1. 用仓库中已有的 `expand_on_pauli` 得到每个键（bond）上的 $c^{(j)}_{\alpha\beta}$。
	2. 按照 JW 展开规则把二次成分识别出来并聚合为跃迁、配对与数项：得到离散化的 $t_j,\Delta_j,\mu_j$（参照 `tools/verify_mzm.py:map_c_to_params` 的具体代数映射）。
	3. 构造 $A$ 与 $B$ 块并拼成 $H_{\mathrm{BdG}}=\begin{pmatrix}A&B\\-B^*&-A^T\end{pmatrix}$。

- 步骤 2 — 求子隙本征态（数值）并构造费米模：
	1. 对 $H_{\mathrm{BdG}}$ 对角化，取某个低能量本征向量 $\Phi=(u_1,\dots,u_N,v_1,\dots,v_N)^T$ 对应能量 $E$（若为 Majorana 零模，$E\approx0$ 且可作规范化）。
	2. 该本征向量对应的费米算符（线性组合）为
	$$
	\gamma_\Phi = \sum_{j=1}^N \big( u_j c_j + v_j c_j^\dagger \big).
	$$
若要求 $\gamma_\Phi$ 为 Hermitian（Majorana），可把系数重构为实型组合：
$$
	\Gamma = \gamma_\Phi + \gamma_\Phi^\dagger = \sum_j\big( (u_j+v_j^*) c_j + (u_j^*+v_j) c_j^\dagger\big),
$$
然后对整体除以适当常数使 $\Gamma^\dagger=\Gamma$ 与 $\Gamma^2=1$（按需要近似归一）。

- 步骤 3 — Jordan–Wigner 展开成 Pauli 串（逐站点显式公式）：
	1. JW 基本映射（链的标准顺序）为
	$$
		 c_j = \Big(\prod_{k<j} Z_k\Big)\frac{X_j+iY_j}{2},\qquad
		 c_j^\dagger = \Big(\prod_{k<j} Z_k\Big)\frac{X_j-iY_j}{2}.
	$$
	2. 把 $\Gamma$ 中每一项替換上述表达式，得到一组带前因子的 Pauli 串（包含长的 Z 串）。展开并合并同一 Pauli 串的复系数得到
	$$
		 \Gamma = \sum_s a_s P_s,
	$$
其中 $P_s\in\{I,X,Y,Z\}^{\otimes N}$ 为 Pauli 串，$a_s\in\mathbb C$。由于原始 $\Gamma$ 要求 Hermitian，展开后系数满足 $a_s^*=a_s$。

- 步骤 4 — 截断与稀疏化（实用）
	1. 计算每条 Pauli 串的模 $|a_s|$，按大小排序并保留前 $M$ 项（例如 $M=50$ 或按累计贡献阈值 99%）。
	2. 统计每个 Pauli 串的支持（最大的站点间距），优先保留局域串以便物理实现。

- 逐项代数示例（单站点边界上的简单项）：
	假设在第 1 和第 2 格的局域哈密顿有一项 $h_{12}=\alpha\,X_1X_2+\beta\,Y_1Y_2$，则通过 JW：
	- $X_1X_2=(c_1+c_1^\dagger)(c_2+c_2^\dagger)\cdot(\prod_{k<1}Z_k)^2$，展开并收集得二次项中含 $c_1^\dagger c_2$ 与 $c_1 c_2$，其线性组合对应到 $t$ 與 $\Delta$。
	- 代数上可直接写出：
$$
		c_1^\dagger c_2 = \frac{1}{4}(X_1X_2+Y_1Y_2+i(X_1Y_2-Y_1X_2))\cdot(\prod_{k<1}Z_k)^0,
$$
$$
		c_1 c_2 = \frac{1}{4}(X_1X_2-Y_1Y_2+i(X_1Y_2+Y_1X_2)).
$$
从而直接读出 $t,\Delta$ 与 $\mu$ 的代数关系（与 `map_c_to_params` 中实现一致）。

- 实现建议（脚本化）
	1. 在 `tools/` 下新增 `pauli_to_majorana.py`：主流程读取 BdG 本征向量或直接读取 `pulse_abs` 产生的本征模向量，构造 $\Gamma$ 并把每个 $c_j, c_j^\dagger$ 用 JW 公式替換，最后把结果化简为 Pauli 串映射并按贡献排序输出 JSON/文本。
	2. 输出格式建议：每行 `coefficient, pauli-string`（例如 `0.312, ZZZIXI...`），并输出局域化度量（位址重心、归一化密度、支持长度）。
	3. 将前若干项（例如前 10 项）以 LaTeX 形式追加到 `steps.md` 或 `myabs.md` 作示例展示。

以上段落即为把数值本征模精确写为 Pauli 串所需的完整代数链与可执行实现细节。实现 `tools/pauli_to_majorana.py` 并对你希望的 pulse 案例运行后，我会把前 10 条 Pauli 展开与局域度量写回到 `myabs.md` 并提交结果文件到 `results/`。

**附录：pulse‑case 的 cosh 拟合示例**

我对 `results/vec_lowest.npy`（pulse 构造的 S–N 链最低本征向量）计算了费米模密度 $\rho_j=|u_j|^2+|v_j|^2$，并用 $\rho(j)\propto 1/\cosh((j-j_0)/\xi)^p$ 的形式做了单边 cosh 拟合，结果如下：

- 链长：N = 200。
- 峰值位置：`peak_idx = 98`，峰值密度 `max_rho = 1.027176e-02`。
- cosh 拟合（p=1）结果：
	- 左侧：`xi ≈ 200`（搜索上界，说明左侧衰减很慢或需扩大搜索区间），拟合残差 `res_rms ≈ 3.727e-03`。 
	- 右侧：`xi ≈ 33.0`，拟合残差 `res_rms ≈ 6.42e-04`。

文件与图像：
- 密度图与 cosh 拟合： [results/gamma_density_cosh.png](results/gamma_density_cosh.png)
- 拟合摘要： [results/gamma_density_cosh_summary.txt](results/gamma_density_cosh_summary.txt)

解读：右侧的 cosh 拟合良好并给出局域长度约 $\xi\approx33$ 个格点；左侧拟合达到搜索上界，表明那边的衰减更宽（或拟合空间需要扩大）。总体上，cosh 型拟合比单纯的单指数拟合在本例中更能同时描述界面附近的峰值与远端衰减。

如需，我可以把左侧 cosh 的 xi 搜索上界扩大并重新拟合以获得更准确的左侧 ξ 值，或把这些数值结果以 LaTeX 表格追加到本文的“数值示例”小节。


**关于其他 R 家族（XYZ / 8‑vertex）与 Δ 的符号结论**

- 要点：对任意两格通用矩阵 R(u) 的符号展开（见 `tools/generic_R_symbolic.py`），我们得到
	$$\Delta\propto d_{3}/\rho$$
	其中 $d_3$ 是 $\partial_u R(u)|_{u=0}$ 的特定矩阵元（按 00,01,10,11 的基序为第 4 个元素），$\rho$ 为规范化因子。
- 含义：只要该导数矩阵元 $d_3\neq0$，理论上即可产生非零配对项 $\Delta$。XXZ（6‑vertex）族在常用规范下恰好使该元为零，因此符号计算给出 $\Delta\equiv0$（这解释了前面的结论）。
- 因此要寻找“R(u 自然产生配对并可能导出 ABS”的候选，优先考虑参数更丰富的 R 家族，例如 XYZ 或 8‑vertex（Baxter）型，它们在某些参量处允许对应的导数矩阵元非零。
- 建议的下一步：对 XYZ / 8‑vertex R(u) 做符号展开（或对其导数矩阵元做代数条件分析），直接检验上述 $d_3$ 是否可非零，从而给出严格的 Δ≠0 的判据（这是下一项待办任务，已加入 TODO）。

我已把这一结论与建议加入本文，后续可把 XYZ 的符号推导结果并入同一节以完成严格证明链。


<!-- 
**补充：Groebner 消元得到的 free‑fermion 代数约束（实际结果）**

根据 `tools/derive_constraints.py` 对一般局域两体 Pauli 系数 $c_{\alpha\beta}$ 的 Groebner 消元，得到必须为零的系数分量（已写入 `results/free_fermion_constraints.txt`）。把消元结果用更物理的语句表述如下：

- 所有系数的虚部必须为零（即所有 $c_{\alpha\beta}$ 都是实数）：
	- `c_00_im`, `c_0x_im`, `c_0y_im`, `c_0z_im`, `c_x0_im`, `c_xx_im`, `c_xy_im`, `c_xz_im`, `c_y0_im`, `c_yx_im`, `c_yy_im`, `c_yz_im`, `c_z0_im`, `c_zx_im`, `c_zy_im`, `c_zz_im`。

- 实部被强制为零的分量（消元出了这些必须消失）：
	- `c_0x_re`, `c_0y_re`, `c_x0_re`, `c_xz_re`, `c_y0_re`, `c_yz_re`, `c_zx_re`, `c_zy_re`, `c_zz_re`。

因此，唯一允许非零的 Pauli 系数实部集合为：
`c_00_re`, `c_0z_re`, `c_xx_re`, `c_xy_re`, `c_yx_re`, `c_yy_re`, `c_z0_re`。

物理含义：上面结果表示若要把一般的两格 Pauli 哈密顿精确写成二次费米 (BdG) 形式，则该两格 Hamiltonian 的 Pauli 展开必须落在由上述允许分量张成的子空间内；否则必然存在四费米项或其它无法被线性二次 BdG 吸收的成分。

示例解集（演示性的具体参数，使哈密顿满足 free‑fermion 约束）

取以下实部值（其它实部与全部虚部均为零）：

- `c_00_re = 0.0`
- `c_0z_re = 0.10`
- `c_z0_re = 0.20`
- `c_xx_re = 1.00`
- `c_yy_re = -0.50`
- `c_xy_re = 0.00`
- `c_yx_re = 0.00`

用之前在文中给出的近似映射（见第五节）把这些 Pauli 系数转为 BdG 参数（此处假设所有虚部为零，使用简化的实数对应关系）：

- 跳跃项（实部）: $t \approx c_{xx}+c_{yy} = 1.00 + (-0.50) = 0.50$. 
- 配对幅度（实部）: $\Delta \approx c_{xx}-c_{yy} = 1.00 - (-0.50) = 1.50$. 
- 化学势（简化表达）: $\mu \approx -2(c_{z0}+c_{0z}) = -2(0.20+0.10) = -0.60$.

该示例给出一个显式的二次 BdG 解（所有虚部为零，且仅允许上面列出的实部非零），可直接用于构造 $H_{\mathrm{BdG}}$ 并检验谱与束缚态性质。

注：上面的数值映射是示范用的线性关系（基于第 5 节给出的成分识别公式）。在一般情形下应使用 `tools/verify_mzm.py:map_c_to_params` 中实现的精确代数表达式把 Pauli 系数映射到 `t, \Delta, \mu`，并对得到的 BdG 做数值对角化验证。
 -->
