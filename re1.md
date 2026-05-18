# re1: `ver5.md` 对应的代码、图像与文档对照

这份记录把 `ver5.md` 里的 braid / YBE / Majorana 路径推导，和实际代码、结果图对应起来。目标是回答三个问题：

1. 代码具体做了什么；
2. 图像具体说明了什么；
3. 这些东西如何回到 `ver5.md` 的数学推导与物理叙述。

---

## 1. 对应的核心代码

1.1 主脚本

- `tools/verify_path_braid.py`

这是 `ver5.md` 最直接对应的脚本。它把文档里的分段控制路径写成显式的 4×4 Hamiltonian，再按段做时间有序积分，最后投影到逻辑子空间并与理想 braid 门 / YBE 进行比较。



1.2 代码流程

脚本的主要步骤是：

1. 用 `path_segment(k, u)` 定义三段路径。
  - 第 1 段：`(u, 0, 1-u, 1)`
  - 第 2 段：`(1-u, u, 0, 1)`
  - 第 3 段：`(0, 1-u, u, 1)`
2. 用 `h4_from_path(...)` 把路径点映射成 4×4 Hamiltonian。
3. 用 `integrate_segment(...)` 对每一段做 midpoint product 近似，得到 `R1, R2, R3`。
4. 用 `compose_cycle(...)` 把三段组合成
  - `R_cycle = R3 @ R2 @ R1`
  - `R_total = R_cycle @ R_cycle`
5. 用 `logical_project(...)` 把 `R_total` 投影到逻辑子空间 `span{|01>, |10>}`，得到 `R_logical`。
6. 用 `logical_braid_target(axis)` 和 `best_logical_axis(...)` 和理想 braid 门比较。
7. 用 `braid_relation_residual(...)` 和 `ybe_residual(...)` 检查 braid relation 和 constant YBE。



1.3 代码和 `ver5.md` 的公式对应

`ver5.md` 里最关键的链条是：

`p(t) -> λ_k(t) -> h_{μν}(t) -> H_⊗(t) -> R(u) -> YBE`

这个脚本把这条链具体化了：

- `p(t)` 就是三段路径 `g_i(t)`；
- `λ_k(t)` 是路径决定的耦合幅度；
- `h_{μν}(t)` 是固定 Pauli 张量基下的系数；
- `H_⊗(t)` 是 4×4 Hamiltonian；
- `R(u)` 是每段积分得到的演化门；
- YBE / braid relation 则是最后的结构检查。

---

## 2. 对应的图像

本节已更新以反映最近的联合参数搜索（`cross_comp_2 × zeeman`）及其摘要图。主要可视化文件包括：

- 早期 summary figure：`results/path_braid/path_braid_steps80_T1_tc1_axisx_hc1.png`（低分辨率摘要），用于快速查看残差分量与逻辑子空间差异；
- 高分辨率对比：`results/path_braid/highres_steps2000_T2.png`（更细时间积分下的 summary），用于检验数值收敛对残差的影响；
- 本次联合搜索的最优摘要：`results/path_braid/scan_cross2_ze_best.png`，诊断数据保存在 `results/path_braid/scan_cross2_ze_best.npz`。

摘要图（`scan_cross2_ze_best.png`）的关键信息：

- 最优参数摘要：`tc=1.0`, `cross_comp_2=0.005`, `zeeman=0.00166667`, `T=3.175`, `steps=300`, `hc_factor=2.125`, `axis=z`。
- 运行指标：`score = 1.27574`，`target_res = 1.05790`，`braid_res = 0.197792`，`ybe_res = 0.020048`。

图像含义与对比要点：

- 左侧柱状残差（与早期 summary 的含义相同）显示在该候选点上 `target`（逻辑门逼近）和 `braid`/`YBE` 项都有不同程度的改善，说明小幅的 XX/YY 横向耦合（`cross_comp_2`）与微小 Zeeman 偏置在本模型中能部分补偿累积相位误差；
- 右侧逻辑子空间差异热图仍显示明显元素级别残差，表明矩阵元素的相位/幅值分布尚未完全对齐；总体 score 改善来自若干残差分量的综合下降，而非单一项的完全消除；
- 与此前的单变量网格相比（例如 `tc × cross_comp_2` 的最优 ≈ 1.29853），本次联合最优在数值上更优（score=1.27574），但仍需稳健性检验以排除对某些未建模静态偏差的敏感性。

数值与物理解读（简短）:

- `cross_comp_2`（段内 XX+YY）改变了两体交换的各向异性，有助于调节逻辑门的演化轴与累积相位；
- 微小 `zeeman`（ZI−IZ）项可作为逻辑能级不对称补偿段间相位差；
- 在该简化 4×4 模型里，两者联合作用产生了一个较浅的“最优盆地”，使得综合残差下降，但该盆地的厚度（鲁棒性）仍需通过交叉验证和 gap/泄漏监控确认。

结论（关于图像的局限性）:

- summary figure 是诊断工具，不应单独作为“复现”证明：需要把该候选点放入稳健性测试（步数变化、静态偏移、噪声平均、gap_min 监控）才可判断其可用性；
- 高分辨率图（如 `highres_steps2000_T2.png`）可以帮助排除数值伪影，但不能替代对物理扰动（例如微小 Zeeman 或 μ 漂移）的检验；
- 因此在把图像写进论文或报告前，应把摘要图的结论与交叉验证结果一起列出。

---

## 3. 和 `ver5.md` 的对应关系

3.1 `ver5.md` 的主线:

`ver5.md` 主要是：

1. 从门控路径 `g_i(t)` 出发；
2. 得到时变微观 Hamiltonian `H_M(t)` / `H_⊗(t)`；
3. 将其写成 Pauli 张量形式；
4. 把整条路径积分成门 `R_total`；
5. 投影到逻辑子空间后，与理想 braid 门比较；
6. 再检查 braid relation / YBE。

这和 `verify_path_braid.py` 的实现是一一对应的。


3.2 文档中的关键数学对象和代码映射:

- `ver5.md` 的 `H_⊗(t)`  -> `verify_path_braid.py` 的 `h4_from_path(...)`
- `ver5.md` 的分段路径 `g^(1), g^(2), g^(3)` -> `path_segment(k, u)`
- `ver5.md` 的总门 `R_{⊗,F}` -> `compose_cycle(...)` 里的 `R_total`
- `ver5.md` 的逻辑子空间 `P_L` -> `logical_project(...)`
- `ver5.md` 的理想 braid 门 `U_braid` -> `logical_braid_target(axis)`
- `ver5.md` 的 braid/YBE 检验 -> `braid_relation_residual(...)` 和 `ybe_residual(...)`



3.3 文档里真正被代码支持的结论:

这套代码支持的结论是：

- 路径可以生成非平凡门；
- 逻辑子空间投影后可以比较 braid-like 行为；
- 但并不是任意路径都会自动给出理想 braid 或严格 YBE 解；
- 更大的积分分辨率可以改善数值稳定性，但不会自动把一个不合适的路径变成理想 braid。

## 4. 图像分析

结论先说清楚：**这是论文思路的结构性复现，不是当前参数下的精确数值复现。**



4.1 图像:

`results/path_braid/path_braid_steps80_T1_tc1_axisx_hc1.png` 和 `results/path_braid/highres_steps2000_T2.png` 都是在看同一件事：

- 逻辑子空间投影后的门，和理想 braid 门有多像；
- 三体嵌入后是否满足 braid relation / YBE；
- 路径离论文里“理想非阿贝尔 braid”有多远。

从图上看：

- 低分辨率图里，逻辑门和理想 braid 门的残差明显不小；
- 高分辨率图里，YBE 残差显著变小，但 braid target 残差仍然偏大；
- 也就是说，这条路径更像“已经有 YBE 结构，但还没有完全贴到论文理想 braid 门”的结果。



4.2 和论文最后那个 $U$ 是否一致:

论文里真正想得到的是：

- 对真 MZM，braiding 结果应当主要由拓扑决定，理想情况下和 braiding 时间无关；
- 对 ABS，结果会带上明显的动力学相，并且依赖路径的具体分段和时间长度。

而我们现在这套数值结果显示：

- 门是非平凡的，确实有 braid-like 结构；
- 但当前路径和参数并没有把逻辑门压到论文里那个“理想、时间无关、接近标准 braid 算符”的极限。

因此，**它和论文的形式是一致的，但数值上还不是精确复现**。



4.3 复现部分:

1. 复现了论文的理论链条：路径 -> Hamiltonian -> 演化门 -> braid/YBE 检验；
2. 复现了论文使用的两体张量结构和逻辑子空间投影方法；
3. 复现了“braid 结果可以通过残差图来判断”的分析框架。


图里会出现“YBE 更接近，但 braid target 仍不够小”的现象：

- 说明张量结构和三体嵌入方向正确；
- 但门控路径和参数还需要进一步调优，才能逼近论文里最终的理想 $U$## 4. 图像分析（基于联合搜索结果）

下面的分析针对最近的联合网格最优结果（参见 `results/path_braid/scan_cross2_ze_best.png` 与 `scan_cross2_ze_best.npz`），把图像的诊断量与物理解读结合起来给出可操作结论。

4.1 关键数值概览

- 最优参数：`tc=1.0`, `cross_comp_2=0.005`, `zeeman=0.00166667`, `T=3.175`, `steps=300`。
- 诊断指标：`score = 1.27574`，`target_res = 1.05790`，`braid_res = 0.197792`，`ybe_res = 0.020048`。

4.2 图像要点（what the figure shows）

- 残差柱状图：`target` 分量仍占主导但比基线降低，`braid` 与 `YBE` 也出现下降，说明该候选在逻辑门逼近与代数一致性上都有所受益；
- 逻辑子空间差异热图：矩阵元素层面仍有显著局部差异（相位与幅值不均），表明改善不是把每个矩阵元精确对齐，而是若干残差项的综合减少；
- 局部盆地迹象：在该模型空间内，`cross_comp_2` 与 `zeeman` 的小组合产生了一个浅槽（score 较低），意味着存在可被利用的补偿自由度，但盆地较浅，鲁棒性待验证。

4.3 物理解读（why it helps）

- `cross_comp_2`（XX+YY 段内耦合）改变了交换项的各向异性，等效于微调逻辑门的演化轴，从而能部分抵消由于路径分段引入的角度/相位偏差；
- 小的 `zeeman` 项相当于在逻辑子空间引入能级不对称，能局部调节 |01> 与 |10> 的相位累积，协同 `cross_comp_2` 达到更好的全局对齐；
- 两者配合能在不显著破坏能隙的前提下减小 target/braid 分量，但若这类补偿过度依赖精确幅值则在实验中难以复现。

4.4 局限性与鲁棒性风险（what to watch for）

- Gap 与泄漏：若改善来自于缩小最小瞬时能隙（gap_min），则该“改进”可能是伪像，会增加非绝热泄漏；必须实时记录 `gap_min` 与投影到非逻辑子空间的泄漏概率来排除这一情形；
- 静态偏差敏感性：早期交叉验证显示某些局部最优对 ±0.01 的 Zeeman 偏移非常敏感，因此需要用更宽的静态偏移集合（例如 ±0.005 或 ±0.01）做平均检验；
- 数值伪影排查：用不同 `steps`（例如 80 与 300、甚至更高）复算并对比，确认不是积分步长导致的虚假最优。

4.5 行动建议（基于图像的可执行步骤）

1. 立即对 `scan_cross2_ze_best` 做三组交叉验证：`steps={80,300}` × `zeeman={0, ±0.005}`，输出均值与方差，并加入 `gap_min` 与泄漏投影统计（若代码中无泄漏投影，先加入计算 `1 - Tr(P_L ρ)` 的监测）；
2. 在候选点附近做局部细化网格（cross_comp_2 ∈ [0.0,0.01], zeeman ∈ [ -0.005, 0.005 ]，分辨率 11×11）绘制 score 热图，以判断盆地的厚度与形状；
3. 若盆地鲁棒，继续把该策略（小幅 cross_comp_2 + small zeeman）纳入带噪声平均的鲁棒优化目标（用平均 score 替代单点 score，并约束 `gap_min > g_thresh`）；
4. 把最终的摘要图（`scan_cross2_ze_best.png`）与交叉验证表格一起放入报告，明确标注“在 X 范围内的均值与标准差”，避免单点宣称。

4.6 简要结论

基于当前摘要图的图像分析：联合微小的横向交叉耦合与轻微 Zeeman 偏置在本 4×4 简化模型内能带来数值上可观的综合残差下降（score 从 ~1.348 降到 1.2757），但该改进目前仍表现为浅盆地，存在对静态偏差敏感和可能的 gap 缩小风险。下一步应优先进行交叉验证与 gap/泄漏监控以确认该候选的物理可行性。。



## 5. Fig.3/Fig.4/Fig.5 趋势图（paper-style 坐标）

为了复现趋势，新增了一个独立脚本：

- [tools/reproduce_trend_figs.py](tools/reproduce_trend_figs.py)

它使用当前简化模型，统一输出三张 paper-style 趋势图，坐标按论文常见写法处理：

- 能谱面板使用 `E/Δ` 对 `t/T`；
- 重叠度面板使用 `T(100/Δ)` 作为横轴；
- 结论定位为“趋势复现”，不宣称逐点定量一致。

本次运行生成文件如下：

- [results/paper_trends/paper_trend_fig3.png](results/paper_trends/paper_trend_fig3.png)
- [results/paper_trends/paper_trend_fig4.png](results/paper_trends/paper_trend_fig4.png)
- [results/paper_trends/paper_trend_fig5.png](results/paper_trends/paper_trend_fig5.png)
- [results/paper_trends/paper_trend_fig3.npz](results/paper_trends/paper_trend_fig3.npz)
- [results/paper_trends/paper_trend_fig4.npz](results/paper_trends/paper_trend_fig4.npz)
- [results/paper_trends/paper_trend_fig5.npz](results/paper_trends/paper_trend_fig5.npz)

这三张图对应的趋势解释是：

1. Fig.3-like：MZM-like 分支对时间更稳，ABS-like 分支随 braiding time 更明显变化。
2. Fig.4-like：引入正弦调制后，最终重叠度对时间的依赖显著减弱，体现动力学相抵消趋势。
3. Fig.5-like：在较绝热区间，不同振幅曲线彼此接近，体现对振幅扰动的低敏感性。

因此，当前仓库已经具备“按论文风格坐展示 Fig.3/4/5 关键趋势”的可复现实例；其定位仍是机理与趋势复现，而非严格的逐点定量复现。

- MZM / ABS 趋势构建流程（实现要点）：
  1. 用 `tools/embed_kitaev.build_bdg(...)` 构造简化的 Kitaev-chain BdG 快照并提取低能谱分支（脚本中函数 `compute_bdg_trace` 封装了时间网格与分支提取）。
  2. 用 `tools/tetron_path_sim.theta_from_time(...)` 与 `H_eff_from_theta(...)` 构造 2×2 有效哈密顿，并在短时间步上用 `expm(-1j * H * dt)` 累积得到总演化算符（脚本中函数 `compute_overlap_curve` 实现）。
  3. 把 BdG 低能谱（随 $t/T$）与最终 overlap($T$) 曲线并列绘制，以比较 MZM（`delta=0`）与 ABS（`delta>0`）在时间上的稳定性差异。MZM 情形一般表现为对 $T$ 更不敏感，ABS 更易随 $T$ 波动。
  4. 自动参数搜索由 `tools/reproduce_trend_figs.py` 的 `run_auto_scan` 实现，评分函数 `trend_score_fig3/4/5` 用来挑选最“paper-like”的参数组合并导出对比图与 `.npz` 数据。

- 复现命令（在仓库根目录，确保 `PYTHONPATH=.`）：

```bash
PYTHONPATH=. .venv/bin/python tools/reproduce_trend_figs.py --outdir results/paper_trends
# 如果需要自动搜索最优参数：
PYTHONPATH=. .venv/bin/python tools/reproduce_trend_figs.py --outdir results/paper_trends --auto-scan --scan-n-per-step 180
```

- 主要脚本与接口：
  - `tools/reproduce_trend_figs.py`：生成 Fig.3/4/5 及 auto-scan（评分并导出 `paper_trend_auto_scan_best.npz`）。
  - `tools/tetron_path_sim.py`：时间网格、`theta_from_time`、`H_eff_from_theta`、短步 propagator（2×2 有效模型）。
  - `tools/embed_kitaev.py`：toy BdG 构造与零能 LDOS / 本征值计算（用于谱线面板）。
  - `quantity/pauli_tensor_nonabelian_reproduction_script.py`：Pauli 张量、三体嵌入与 YBE / 非阿贝尔度量示例。

- 说明：
  - 这些图是使用仓库内的简化 “toy” 模型生成的，目的在于趋势与机理的复现；若要做到论文中严格的 BdG 数值复现，需要替换/扩展 BdG 实现、提高时间积分精度并加入 gap/leakage 约束与噪声平均等鲁棒性检验。

## 具体逻辑（代码 → 数据 → 绘图）

下面把“我们是如何把两个简化模型并行起来生成 Fig.3/4/5”的具体逻辑按代码、产生的数据和绘图步骤逐条写清：

- 入口函数：`tools/reproduce_trend_figs.py` 的 `plot_fig3(outdir)`, `plot_fig4(outdir)`, `plot_fig5(outdir)`（可选 `run_auto_scan` 用于参数搜索）。

- 共用底层接口：
  - 时间网格与分段路径：`tetron.make_time_grid(T_step, n_per_step)`、`tetron.gates_at(step,s)`、`tetron.theta_from_time(step,s)`。
  - 有效两能级模型（2×2）：`tetron.H_eff_from_theta(theta, delta)` 返回 $H=\cos\theta\,\sigma_x+\sin\theta\,\sigma_y+\delta\,\sigma_z$；短步 propagator 使用 `expm(-1j*H*dt)` 累积得到总演化算符（`compute_overlap_curve` 实现）。
  - toy BdG (Kitaev 链)：`embed_kitaev.build_bdg(mu,t_links,Delta_links)` 构造 BdG 矩阵，`np.linalg.eigvalsh(H)` / `embed_kitaev.compute_zero_ldos` 提取低能谱（`compute_bdg_trace` 封装）。

- Fig.3 具体流程（能谱 + MZM/ABS overlap）：
  1. 调用 `compute_bdg_trace(T_step=400.0, n_per_step=300, delta_mod=None)`：
     - 内部用 `tetron.make_time_grid` 得到时间点；对每个时间点用 `tetron.gates_at` 得到 `g1..g4`；`map_gates_to_links` 把门控映射为链上 `mu,t_links,Delta_mod`；`embed_kitaev.build_bdg` 构造 BdG 并用 `eigvalsh` 取四个最低本征值 → 返回 `t_over_T` 与 `branches`（shape=(N,4)）。
  2. 在显示窗口 `Tscan = linspace(FIG3_DISPLAY_TMIN, FIG3_DISPLAY_TMAX, FIG3_DISPLAY_NPTS)` 对每个 `TT`：
     - 调用 `compute_overlap_curve(T_step=TT, delta=0.0, n_per_step=...)`（MZM）和 `compute_overlap_curve(..., delta=FIG3_DISPLAY_ABS_DELTA)`（ABS）；函数返回随时间的 `overlaps` 序列，取 `abs(overlaps[-1])` 作为该 `TT` 的 final-overlap。
  3. 绘图：上面板画 `branches` vs `t_over_T`，下面板画 `Tscan/P.T_UNIT` vs 两条 final-overlap 曲线；保存 `paper_trend_fig3.png` 并用 `np.savez` 导出 `paper_trend_fig3.npz`（键：`t_over_T, branches, Tscan, mzm_overlaps, abs_overlaps`）。

- Fig.4 具体流程（调制下的动力学相抵消展示）：
  1. BdG 谱（case a/b）：分别用 `compute_bdg_trace(..., delta_mod=0.59*P.DELTA, amp=0.02*P.DELTA)` 和 `compute_bdg_trace(..., delta_mod=0.57*P.DELTA, amp=0.02*P.DELTA)`，在 `compute_bdg_trace` 中 `VD_here = P.VD*(1.0 + delta_mod + amp*cos(pi*t/T_step))` 被用作 QD 深度代理，从而在 BdG 中体现正弦调制。
  2. Overlap（case a/b）：对 `Tscan = linspace(300,700,17)`，分别调用 `compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=250, modulation=(d0,amp))`（其中 `compute_overlap_curve` 内部把 `delta_eff = d0 + amp * cos(pi * t / T_step)`，直接影响 2×2 瞬时能量与动态相累积）。
  3. 绘图：上面板显示 case a 的 BdG 分支，下面板并列绘制 case a 与 case b 的 final-overlap 曲线；保存 `paper_trend_fig4.png`/`.npz`。

- Fig.5 具体流程（绝热区下振幅不敏感性）：
  1. 设 `Tscan = linspace(300,1000,20)`；对两个振幅 `amp in {0.01,0.03}` 调用 `compute_overlap_curve(T_step=TT, delta=0.2, n_per_step=250, modulation=(0.59, amp))`，收集每个 `TT` 的 final-overlap。
  2. 绘图：在同一坐标绘制两条 amp 曲线并对比尾部（大 T）是否靠拢；保存 `paper_trend_fig5.png`/`.npz`。

- 自动参数搜索与评分：`run_auto_scan(outdir, n_per_step=...)` 枚举参数格点，使用 `trend_score_fig3/4/5` 对应评分函数评估每组参数的“paper-like”程度，保存 `paper_trend_auto_scan_best.npz` 与 `paper_trend_auto_scan_summary.txt`。

- 复现与检查：运行主脚本后，可用 `np.load('results/paper_trends/paper_trend_fig3.npz')` 检查数组与绘图输入；若要严格验证数值稳定性，应增加 `n_per_step/steps`、计算 `gap_min` 与 leakage 投影并对参数做平均/交叉验证。

## 6. 结果图解读

下面对三张趋势图给出更细的解读与对照逻辑：


**6.1 Fig.3-like**（能谱演化 + 时间依赖对比）:

对应图像：

- [results/paper_trends/paper_trend_fig3.png](results/paper_trends/paper_trend_fig3.png)

如果你要看“MZM 和 ABS 到底差在哪里”，默认的 Fig.3-like 图已经切到一个更清晰的显示窗口了；如果想把差异看得更强，再对照这张自动挑选参数后的对比图：

- [results/paper_trends/paper_trend_auto_compare.png](results/paper_trends/paper_trend_auto_compare.png)

默认版现在已经能看出 MZM-like 更平、ABS-like 更容易弯折的趋势；自动扫描版则把这种差异进一步放大，更适合拿来做细读和对照。

**(a) 能谱演化面板**

这张面板展示的是简化 BdG 模型在门控路径下的低能谱分支随 $t/T$ 的变化。可以读出三点：

1. 谱线在分段节点附近有明显的“结构性拐折”，说明门控路径对低能谱有分段驱动作用；这对应论文里“分段交换路径导致的结构性演化”。
2. 低能谱分支整体仍围绕 $E=0$ 对称分布，体现的是“准零能态在路径驱动下的演化”，对应论文讨论的 MZM/ABS 区分背景。
3. 相比纯 MZM 极限，该谱仍有可见的能量幅度，说明这是“ABS-like 的可动能谱”，而不是严格零能扁平。趋势层面上与论文 Fig.3 的 ABS 情况一致。

“ABS 的拐点”：不是一个真正的相变点，也不是能隙突然闭合的点，而是一个**crossover（交叉过渡）**。在这个简化模型里，拐点意味着两件事发生了切换：

1. 低能态不再只是“跟着路径平滑地走”，而开始明显感受到有限能量分支的弯折。
2. 绝热跟随与动力学相积累之间的平衡被打破，导致谱线斜率改变、曲线出现折弯。

换句话说，这个拐点说明 ABS 不是像 MZM 那样被拓扑钉住的零能模；它会随着控制路径和时间尺度改变自己的瞬时能量，所以它对 braiding time 更敏感。

**(b) 重叠度 vs T 面板**

这一面板体现的核心趋势是：

1. MZM-like 曲线在 $T$ 扫描中变化更小，说明结果更接近“时间无关”的拓扑主导行为。
2. ABS-like 曲线随 $T$ 变化更明显，体现“动力学相主导”的时间依赖性。
3. 两条曲线平均值的分离，反映了“拓扑保护 vs 动力学相污染”的趋势差异。


在我们的简化模型里：

1. MZM-like 曲线更平，表示最终重叠度对 $T$ 的依赖弱，交换结果更接近“路径决定、时间不敏感”。
2. ABS-like 曲线出现折弯/拐点，表示随着 $T$ 改变，动力学相不再能被完全抵消，结果开始明显偏离 MZM-like 的稳定区。
3. 这个折弯就是“ABS 不稳定”的可视化信号：它告诉你，ABS 不是一个完全拓扑保护的 braiding 对象，时间尺度一变，交换结果就会跟着变。

所以，如果把这张图压缩成一句话，它说明的是：**MZM 的交换更像拓扑控制，ABS 的交换更像动力学控制。**

这与论文 Fig.3 的主旨一致：MZM 情况更稳定，ABS 情况更依赖 braiding time。

Fig.3 这张图在区分“交换结果主要由拓扑决定”与“交换结果主要由动力学相决定”。

更具体一点：

1. MZM-like 曲线更平，说明它对 braiding time 不敏感，交换后的结果更像是由路径拓扑决定的。
2. ABS-like 曲线出现拐点/折弯，说明它的低能态不是被钉死在零能附近，而是会随着时间尺度变化发生可见的能量和相位重排。
3. 因此ABS 也能做交换，但它的结果不是天然稳定的；当时间尺度改变时，它的交换输出会跟着变。


**6.2 Fig.4-like**（调制导致动力学相抵消）:

对应图像：

- [results/paper_trends/paper_trend_fig4.png](results/paper_trends/paper_trend_fig4.png)

如果把论文原意说得更直接一点，Fig.4 不是在展示“又一组好看的谱线”，而是在演示一个控制技巧：**用正弦调制去反转 ABS 的能级符号，从而把动态相尽量抵消掉**。论文分了两层意思：先只消掉动态相，再把剩余的方位相也一起消掉。

**(a) 调制下的能谱演化**

论文里这一幅对应的是 $V_x = [0.59 + 0.02\cos(t\pi/T)]\,\Delta$。它想表达的不是“能谱本身长得多复杂”，而是：调制以后，ABS 的能量在 braiding 过程中会跨过零点，正负能量段在积分时可以互相抵消，所以动态相 $\theta_2$ 近似消失。

更具体地说，论文这里关心的是“积分后的净相位”，不是瞬时能量的某一个点。也就是说：

1. 能谱可以一度上升或下降，但只要正负两段在总积分里相互抵掉，动态相就被消除了。
2. 这一步消掉后，braiding 的结果已经不再像 Fig.3 那样强烈依赖过程时长。
3. 但这还没到完全恢复 MZM 稳定性的程度，因为方位相 $\theta_1-\theta_3$ 还可能残留。

我们的简化图把这件事压缩成了“调制后能谱/重叠曲线更平、更接近彼此”。它没有把论文里“动态相 $\theta_2$ 先被消掉”的相位结构单独画出来，但趋势上是在表达同一件事：**调制让动力学项失去主导地位**。

**(b) 调制后重叠度 vs T**

论文里这一部分有两个层次：

1. 在 $V_x=[0.59+0.02\cos(t\pi/T)]\,\Delta$ 时，$\psi_{1}^{-}(6T)$ 的**幅度**已经基本不依赖 $T$，但**相位**仍会随 $T$ 振荡。
2. 在 $V_x=[0.57+0.02\cos(t\pi/T)]\,\Delta$ 且只作用于 $t\in[4T,6T]$ 时，连 $\theta_1-\theta_3$ 也被消掉，于是幅度和相位都对 $T$ 不敏感，接近真正的 MZM 式稳定非阿贝尔交换。

我们的图目前更接近论文的第一层：它能看出调制后最终重叠曲线变平、对 $T$ 的敏感性减弱，也就是说“动态相被压低”这件事是对上的。但它还没有像论文那样把“幅度稳定但相位还在振荡”与“幅度和相位都稳定”这两种情形分成两个清晰子面板，所以需要把它理解成**趋势对应，而不是完整相位结构的逐项复刻**。

换句话说，Fig.4 对我们的图的要求其实很明确：

1. 先让结果看起来不再强烈依赖 $T$，这表示动态相正在被取消。
2. 再看调制是否足够强，能不能把剩下的方位相也压掉。
3. 只要这两个方向的趋势成立，就说明我们的简化模型已经抓住论文 Fig.4 的核心物理。

这与论文 Fig.4 的要点一致：合理的调制会让最终的非阿贝尔结果更加接近拓扑决定的极限。



**6.3 Fig.5-like**（绝热区振幅不敏感）:

对应图像：

- [results/paper_trends/paper_trend_fig5.png](results/paper_trends/paper_trend_fig5.png)

如果说 Fig.4 关注的是“调制能不能把动态相抵消掉”，那么 Fig.5 关注的就是更实际的一点：**同样的调制机制，换一个振幅，结果还稳不稳。** 论文的回答是：只要还在绝热区，结果是稳的。

这一图的关键趋势是：

1. 在大 $T$（绝热区）里，不同振幅曲线彼此靠近，说明对振幅扰动的响应减弱。
2. 曲线尾部趋于平坦，体现绝热极限下的“参数不敏感性”。

论文原文给出的具体比较是 $V_x=[0.59+0.01\cos(t\pi/T)]\,\Delta$ 和 $V_x=[0.59+0.03\cos(t\pi/T)]\,\Delta$。它要说明的不是这两个振幅本身有多特别，而是：**只要 braiding 足够慢，振幅从 0.01 变到 0.03，最终交换统计几乎不变。** 这意味着去掉动态相的机制不是一个“精细调参才成立”的技巧，而是一个在绝热条件下相当鲁棒的机制。

我们的简化图现在表达的就是这个趋势：两条振幅不同的曲线在大 $T$ 端逐渐贴近，说明随着绝热性增强，调制振幅对最终结果的影响被压低了。需要诚实地说的是，我们这里画的是最终重叠度的趋势，不是论文中完整的复幅相位/矩阵元素分析，所以它更适合被读成“绝热鲁棒性”的验证，而不是“完整非阿贝尔矩阵元”的复现。

所以 Fig.5 这张图对我们的图最准确的要求就是：

1. 不同振幅的曲线要在大 $T$ 区域越来越接近。
2. 这种接近不是偶然的数值巧合，而是绝热条件下动力学相被压制的表现。
3. 只要这一点成立，就说明我们的简化模型已经把论文 Fig.5 的核心结论抓住了。

这正是论文 Fig.5 的核心趋势：足够绝热时，braiding 结果对振幅变化不敏感。



6.4 总结

从上述三张图的趋势可以明确给出结论：

1. 简化模型已经能复现论文 Fig.3/4/5 的**定性趋势**。
2. 具体数值幅度、曲线细节和参数映射仍与论文完整 BdG 模型存在差距。
3. 因此最准确的表述是：**趋势一致，但定量不一致**。

## 7. 总结

这次针对论文 Fig.3/Fig.4/Fig.5 的复现，最终可以总结为:已经在简化模型中复现了论文想表达的趋势和机理，但还没有做到完整 BdG 模型下的定量一致。

更具体地说：

1. Fig.3-like 主要复现了“ABS 随 braiding time 更敏感、MZM 更稳定”的对比趋势；默认版已经能看出差异，自动扫描版把差异放得更明显。
2. Fig.4-like 复现了“通过正弦调制抵消动力学相，使 braiding 结果更接近时间无关”的趋势。
3. Fig.5-like 复现了“在绝热区，不同调制振幅下的结果趋于一致”的鲁棒性趋势。

对应脚本：

- [tools/reproduce_trend_figs.py](tools/reproduce_trend_figs.py)

对应结果路径：

- [results/paper_trends/paper_trend_fig3.png](results/paper_trends/paper_trend_fig3.png)
- [results/paper_trends/paper_trend_fig4.png](results/paper_trends/paper_trend_fig4.png)
- [results/paper_trends/paper_trend_fig5.png](results/paper_trends/paper_trend_fig5.png)
- [results/paper_trends/paper_trend_auto_compare.png](results/paper_trends/paper_trend_auto_compare.png)
- [results/paper_trends/paper_trend_auto_scan_summary.txt](results/paper_trends/paper_trend_auto_scan_summary.txt)
- [results/paper_trends/paper_trend_auto_scan_best.npz](results/paper_trends/paper_trend_auto_scan_best.npz)

