# Fig.8-10 Comparison Summary

## 论文图片详细解读

本文件先对论文中的 Fig.8、Fig.9 和 Fig.10 做逐图解读，明确每块子图的物理含义、对比点和关键趋势。之后再对照我们复现图，逐项评估相似性与差异。

- Fig.8 主要用两种势配置（uniform / inhomogeneous）对比能级谱、zero-bias conductance 和 braiding fidelity，重点在于展示不均匀势增强非阿贝尔 braiding 的稳定性。
- Fig.9 通过 6 种不同 nanowire / 边界条件配置，考察 zero-bias conductance 在不同场景下随磁场变化的峰值是否稳定、是否分裂，以及是否能区分 Majorana-like 和 ABS-like 行为。
- Fig.10 以 disorder 配置为背景，将谱图、局域估计器与零偏 conductance 结合起来，进一步把重点放在选取 3 个典型 B 点后，braiding fidelity 随时间成本 τ 的变化，以及量子点能级 Ed 对谱图的影响。

每个图中的子面板都有明确对应关系：
- 能级谱面板用于判断是否存在近零能态或分裂态；
- conductance 面板用于观察 zero-bias peak 的出现条件和稳定性；
- fidelity 面板用于检验条件下是否真正具备可靠 braiding 性能；
- Ed 扫描面板用于区分内插态与量子点诱导的 ABS 结构。

在我们的复现中，差异往往出现在“是否真的计算了 braiding fidelity”以及“是否使用了与论文一致的横轴变量（如 Ed）”这两点上。

## 论文原图核心含义

### Fig.8
- (a)/(d): uniform / inhomogeneous 潜势下的能级谱随磁场 B 变化。
- (b)/(e): 对应的 tunneling differential conductance，重点是 zero-bias peak (ZBP) 的出现与稳定性。
- (c)/(f): 在固定 B 点下，braiding fidelity 随 time cost τ 的变化；inhomogeneous 情况下 fidelity 更高、更加稳定。

### Fig.9
- 6 种不同 nanowire/边界结构配置下的 tunneling conductance。
- 重点比较不同情形下 zero-bias conductance 峰是否随 B 变化保持稳固或发生分裂。
- 展示不同失配 / quantum-dot / disorder 情况下 ZBP 行为的差异。

### Fig.10
- (a) disorder configuration 示意。
- (b) 能谱与局部估计器。
- (c) differential conductance vs (E,B)，强调 3 个 B 点处的 ZBP。
- (d)-(f) 3 个 B 点在不同 τ 下的 braiding fidelity。
- (g)-(i) 3 个 B 点对应的 dot-nanowire 能谱 vs QD 能级 Ed。

## 我们的复现图对比

### Fig.8 对比

论文原图侧重：
- uniform / inhomogeneous 的能级谱差异；
- 真实 differential conductance 的 ZBP；
- braiding fidelity 在 inhomogeneous 下更高更稳定。

我们复现：
- (a)/(d) 也画了 uniform-like / inhomogeneous-like 能级谱；
- (b)/(e) 现在改成了带探针权重的 zero-bias conductance 近似图，形式更接近论文；
- (c)/(f) 仍是 kernel `\mathcal{N}` 随 `u` 变化，而不是 fidelity 曲线。

判断：
- 结构上类似，但核心物理含义不同；
- 趋势上部分一致：出现了 uniform vs inhomogeneous 的对比；
- 但并不能说“完全一致”，因为我们没有直接复现论文里“braiding fidelity 更高”的结论。

### Fig.9 对比

论文原图侧重：
- 6 种情形下 ZBP 稳定性与分裂行为的差异。

我们复现：
- 也有 6 个子图；
- 现在的 G 计算已更接近零偏 conductance；
- 但图之间差异比论文弱，表现为“同类情形的 conductance 变化”，而不是论文那样的强烈区分。

判断：
- 形式上很像；
- 趋势上也有“不同情形下 conductance 变化”的方向一致性；
- 但论文中每个场景的 ZBP 稳定性差异在我们图中表现得不够明显。

### Fig.10 对比

论文原图侧重：
- disorder 配置 + 能谱 + zero-bias conductance；
- 3 个 B 点的 braiding fidelity vs τ；
- dot-nanowire 能谱 vs Ed。

我们复现：
- (a) disorder configuration 及 (b) spectrum vs B 与 3 个 B 标记类似；
- (c) 也有 zero-bias conductance 近似图；
- 但 (d)-(f) 不是 fidelity，而是 kernel `\mathcal{N}(τ)`；
- (g)-(i) 不是 Ed 扫描，而是 spectrum vs `δ`。

判断：
- 前半部分结构相似；
- 后半部分物理量和横轴不同，因此不能说趋势一致。

## 我们复现图现象与分析

### Fig.8-like 现象与分析
- 能级谱面板体现了两种势条件下低能态分布的差异：uniform-like 情形低能态更对称、分裂较小，而 inhomogeneous-like 情形出现了更明显的近零能态和非平衡能级重排。
- conductance 面板显示在 zero-bias 附近存在峰值增强，但峰形更宽、更平缓，缺少论文中明确的尖峰与稳定平顶特征。
- kernel `\mathcal{N}` 面板呈现出随路径参数 `u` 的变化趋势，说明我们当前的“braiding指示器”在参数空间中有变化，但它不是直接的操作保真度，因此不能直接用来判定论文所说的 braiding 成功性。
- 结论：我们的 Fig.8-like 图展示了类似的“能级对比 + zero-bias 响应”结构，但在“是否真正对应高保真 braiding”这一点上仍然不够精确。

### Fig.9-like 现象与分析
- 六个子图总体保持了不同条件下 conductance 随 B 变化的趋势差异，反映了模型对参数变化的敏感性。
- 但从现象上看，ZBP 的出现/消失与分裂主要表现为连续的谱强度改变，而不是论文中那种明显的分裂线条和稳定峰值区分。
- 由于我们对 G 的定义基于 eigenvector 权重和零偏谱密度，图像更像“谱强度映射”而不是“差分电导”; 因此它可以给出趋势判断，但不能严格复刻论文里的实验 zero-bias conductance 细节。
- 结论：Fig.9-like 图能够捕捉到“场强与 conductance 关联”的现象，但缺失了论文中区分不同 ABS/Majorana 情形的典型特征。

### Fig.10-like 现象与分析
- (a)-(c) 部分：我们的 disorder 配置图和 spectrum/vs-B 图呈现出类似的背景结构，并且 zero-bias conductance 近似图能找到若干候选 B 点。
- (d)-(f) 部分：当前使用 kernel `\mathcal{N}(τ)` 的变化来替代 fidelity 曲线，使得这组面板更多表现为时间成本参数对非阿贝尔估计器的影响，而不是论文所要求的 braiding fidelity 动力学。
- (g)-(i) 部分：扫描变量是模型参数 `δ`，而不是论文强调的 quantum dot 能级 `Ed`，因此谱图仍然表达了参数调整下的能级演化，但缺少特定 dot-纳米线耦合下的 ABS vs Majorana 区分。
- 结论：Fig.10-like 图在结构上接近论文，但在关键物理变量和“braiding fidelity vs Ed”判断上存在本质差别。

## 综合结论

- 现在的复现图在“图像结构”和“zero-bias conductance 含义”上已经拉近了与论文的距离。
- 但严格来说，仍然是“同类现象的 Pauli-tensor 模型类比复现”，而不是论文原始 Majorana/ABS fidelity 图的逐点复刻。
- 如果要判断“最起码趋势是否一样”，结果是：
  - Fig.8: 部分趋势一致，但 fidelity 结论未直接复现；
  - Fig.9: 形式一致、趋势弱一致；
  - Fig.10: 只有前半部分趋势一致，后半部分不一致。

## 进一步建议

- 若要增强一致性，下一步应把 `(c)/(f)` 的 kernel `\mathcal{N}` 改成更接近论文的 braiding fidelity；
- 并把 Fig.10 的 `(g)-(i)` 改成 dot-nanowire 的 `Ed` 扫描。

## 脚本实现思路和方法

### 核心模型：Pauli 张量有效哈密顿

脚本使用四维 Pauli 张量空间中的有效哈密顿量：

$$H_4(u,\delta) = \cos(u) \sigma_x \otimes \sigma_x + \frac{1}{2}\sin(u)(\sigma_y \otimes \sigma_x - \sigma_x \otimes \sigma_y) + \frac{\delta}{2}(\sigma_z \otimes I - I \otimes \sigma_z)$$

其中 $u$ 是路径参数（通过 $u = 2\pi B / B_{max}$ 与磁场 $B$ 关联），$\delta$ 是不均匀势强度。

**公式来源说明**：这个哈密顿 **不是** 直接从某个标准 Yang-Baxter R 解通过 $H(u) = i\partial_u R(u) R(u)^{-1}$ 推导而来的。相反，它是基于以下原则的参数化构造：

1. **Pauli 张量结构**：采用 $\sigma_\alpha \otimes \sigma_\beta$ 作为基础，确保系统具有两体 Pauli 相互作用。
2. **非阿贝尔性**：通过 $\sigma_y \otimes \sigma_x - \sigma_x \otimes \sigma_y$ 这类反对称组合引入非交换性，这是 Yang-Baxter 关系所要求的。
3. **路径参数化**：用 $u$ 参数化演化，使得 $\cos(u)$ 和 $\sin(u)$ 项产生周期性与旋转。
4. **物理启发**：$\delta$ 项类似于论文中不均匀势的作用，通过打破对称性来调控能谱与零偏 conductance。

因此，这个公式更准确地说是一个 **受 Yang-Baxter 思想启发、基于 Pauli 张量代数构造的模型参数化**，而非严格的 R 矩阵投影。

**与 8-vertex R 矩阵的关系**：理论上，若从标准 8-vertex R 矩阵出发，可以通过对数导数 $H(u) = i\partial_u R(u) R(u)^{-1}$ 直接推导出对应的生成元。例如，标准 8-vertex 模型：

$$R_{\rm 8v}(u)=\begin{pmatrix} \cos u & 0 & 0 & i\delta\sin u\\
0 & \sin u & \cos u & 0\\
0 & \cos u & \sin u & 0\\
i\delta\sin u & 0 & 0 & \cos u\end{pmatrix}$$

通过投影和对数导数可得到：

$$H_{\rm eff}(u) \propto \begin{pmatrix} \delta & e^{-iu}\\ e^{iu} & -\delta\end{pmatrix}$$

对应 Pauli 系数 $(\cos u, \sin u, \delta)$。这与我们当前 H 的结构部分类似，但我们采用的是更一般化、更易于参数调整的参数化形式，以便于在数值中快速探索不同物理情景。


### 路径有序演化与 Yang-Baxter 关系

脚本通过离散化步骤计算路径有序幺正演化：

1. **路径有序 $R$ 矩阵**：将 $u$ 从 0 扫到终点值，每步计算 $R(u_i) = \exp(-i H_4(u_i) \Delta u) @ R(u_{i-1})$，最终得到 $R_u$。
2. **嵌入到 6 维空间**：通过 $R^{(12)}_u = R_u \otimes I_2$ 和 $R^{(23)}_u = I_2 \otimes R_u$ 将 4 维操作复制到相邻对。
3. **Yang-Baxter 偏差计算**：评估 $\Delta = R^{(12)}_u R^{(23)}_{uv} R^{(12)}_v - R^{(23)}_v R^{(12)}_{uv} R^{(23)}_u$，反映非阿贝尔结构。

### Kernel 非阿贝尔度量 $\mathcal{N}_{kernel}$

脚本定义了一个基于第三阶 kernel 积分的非阿贝尔估计器：

$$\mathcal{N}_{kernel} = \sqrt{\sum_{abmn} \left| \int_0^u ds \int_0^{u+v} dt \, K(s,t;u,v) \, h_a^b(s) \, h_m^n(t) \right|^2}$$

其中：
- $h_a^b(u)$ 是哈密顿系数（Pauli 基展开）
- $K(s,t;u,v)$ 是分段线性 kernel，在 $s,t$ 平面内编码时序与路径结构
- 求和遍历 Pauli 矩阵索引 $\{I,x,y,z\}$

这个度量反映了"路径依赖性"和"非交换程度"，但 **不是** 论文中的 braiding fidelity；它更多是一个结构性指示器。

### Zero-bias 导电谱的计算

脚本实现了一种基于本征态权重的零偏导电近似：

1. **本征态获取**：对每个 B 点，计算 $H_4(u)$ 的本征值 $E_n$ 和本征向量 $v_n$。
2. **权重分配**：根据选定的 lead（例如第一个 Pauli site）上的概率密度计算权重 $w_n = |v_n[lead]|^2$，归一化使 $\sum_n w_n = 1$。
3. **Lorentzian 求和**：对每个能量 $E$ 计算 conductance：
   $$G(E,B) = \sum_n w_n \frac{\gamma^2}{(E-E_n)^2 + \gamma^2}$$
   其中 $\gamma$ 是 broadening 宽度（通常 0.009–0.012 eV）。

这个方法模拟了实验中通过 tunneling spectroscopy 观测到的导电峰，但由于缺少真实 tunneling matrix element，本质上是"本征态项"的加权求和，而不是精确的差分电导。

### 三张图的生成流程

**Fig.8-like**：
- 参数对比：uniform-like ($\delta = 0.10$) vs inhomogeneous-like ($\delta = 0.04$)
- 左列 (a)/(d)：能级谱 vs B
- 中列 (b)/(e)：zero-bias conductance 热图 vs (E,B)
- 右列 (c)/(f)：kernel $\mathcal{N}$ vs 路径参数 $u$

**Fig.9-like**：
- 6 个子图分别对应 6 组参数组合（虚拟的"不同纳米线配置"）
- 每个子图展示该参数下的 zero-bias conductance 热图
- 通过改变 $\delta$ 值来模拟"不同条件"的差异

**Fig.10-like**：
- (a)：模拟 disorder 配置示意（虽然我们的模型中 disorder 不显式存在）
- (b)：能级谱与标记的 3 个 B 点
- (c)：整体 zero-bias conductance 热图
- (d)-(f)：3 个选定 B 点上，kernel $\mathcal{N}$ 随 $u$ 的变化（替代 fidelity vs τ）
- (g)-(i)：3 个 B 点的能级谱 vs 模型参数 $\delta$（替代 vs Ed）

### 关键参数与精度

- **网格分辨率**：B 扫描通常 160 点，E 扫描 220 点，路径参数扫描 80–160 点，Kernel 积分采样 120×160
- **Broadening**：$\gamma = 0.009$–$0.012$ eV，选择平衡谱分辨率与图像清晰度
- **路径步数**：通常 140–220 步，确保路径有序演化的收敛
- **Kernel 网格**：$s$ 和 $t$ 独立采样，重积分通过矩形法则近似

### 与论文的差异根源

本质上，脚本复现了"Pauli tensor 模型下的能谱与非阿贝尔特征"现象，而非直接复现论文的"Majorana/ABS 结构与 braiding fidelity"。差异主要来自：

1. **模型差异**：Pauli tensor vs Majorana fermion / ABS
2. **物理量差异**：kernel $\mathcal{N}$ vs braiding fidelity；spectrum vs (Ed,τ) vs (δ,τ)
3. **操作差异**：路径有序演化 + kernel 积分 vs 论文的完整 braiding 过程
4. **度量差异**：基于 tunneling DOS 的 zero-bias conductance vs 基于 tunneling probability 的实验 $dI/dV$

