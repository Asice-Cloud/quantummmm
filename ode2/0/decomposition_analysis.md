# SO(5)/Sp(2) 几何-动力学分解条件与抵消策略

> 核心问题：对于 $\dot U = KU$，$U(T) \in Sp(2) \cong Spin(5)$，
> 何时 $U(T) = G \cdot D$（几何 $\perp$ 动力学），不可分解时如何抵消 $D$？

---

## 一、预备：问题的形式化

### 1.1 演化方程

$$
\dot U(t) = K(t) U(t),\quad K \in \mathfrak{sp}(2),\quad U \in Sp(2)
$$

$$
K = \begin{pmatrix} A & B \\ C & D \end{pmatrix},\quad
\begin{aligned}
A &= \frac{E_1+|t_3|}{2}\mathbf i + \frac{|t_2|}{2}\mathbf j \\[4pt]
D &= \frac{-E_1+|t_3|}{2}\mathbf i + \frac{|t_2|}{2}\mathbf j \\[4pt]
B &= \frac{|t_1|}{2} + \frac{E_d}{2}\mathbf k \\[4pt]
C &= -\frac{|t_1|}{2} + \frac{E_d}{2}\mathbf k
\end{aligned}
$$

### 1.2 参数的物理分类

| 参数 | 类别 | 作用 |
|------|------|------|
| $t_2, t_3, E_d$ | **控制参数**（几何驱动） | 时变脉冲，闭合回路 → holonomy |
| $E_1$ | **寄生参数**（动力学驱动） | MZM 杂化能，几乎不可控 |
| $t_1$ | **半可控参数**（动力学驱动） | ancilla-MZM 耦合，部分可控 |

### 1.3 什么是"分解"

**严格分解**：$U(T) = U_{\text{geo}} \cdot U_{\text{dyn}}$，满足

1. $U_{\text{geo}}$ 仅依赖于控制参数在参数空间中的**路径几何**（路径的像，非参数化速度）
2. $U_{\text{dyn}}$ 仅依赖于哈密顿量的**时间积分** $\int H(t)dt$
3. $U_{\text{geo}}$ 和 $U_{\text{dyn}}$ 各自属于 $Sp(2)$ 的可区分子群

**弱分解**：$U(T)$ 的某些分量（如 MZM 子空间投影 $R_{123}$）可写成纯几何 × 纯动力学的形式。

---

## 二、分解判据：三条等价条件

### 2.1 判据 A：哈密顿量的对易子闭包

> **充要条件**：存在 Lie 代数分解 $\mathfrak{g} = \mathfrak{g}_{\text{geo}} \oplus \mathfrak{g}_{\text{dyn}}$（直和），
> 且对任意 $t, s \in [0,T]$，$K_{\text{geo}}(t) \in \mathfrak{g}_{\text{geo}}$，$K_{\text{dyn}}(s) \in \mathfrak{g}_{\text{dyn}}$，
> 满足 $[\mathfrak{g}_{\text{geo}}, \mathfrak{g}_{\text{dyn}}] = 0$。

在此条件下：
$$U(T) = \mathcal{T}\exp\!\left(\int_0^T K_{\text{geo}}\,dt\right) \cdot
         \mathcal{T}\exp\!\left(\int_0^T K_{\text{dyn}}\,dt\right)$$

**这是最强的分解形式。** 仅有 $[\mathfrak{g}_{\text{geo}}, \mathfrak{g}_{\text{dyn}}] = 0$ 的这个条件，
才能保证时间编序指数可因子化。

证明：
$$\mathcal{T}\exp\!\left(\int (K_{\text{geo}}+K_{\text{dyn}})dt\right)
= \mathcal{T}\exp\!\left(\int K_{\text{geo}}\,dt\right)
  \cdot \mathcal{T}\exp\!\left(\int K_{\text{dyn}}\,dt\right)$$

当且仅当 $[K_{\text{geo}}(t_1), K_{\text{dyn}}(t_2)] = 0$ 对所有 $t_1, t_2$ 成立。
这等价于它们生成的子代数对易。

### 2.2 判据 B：Maurer-Cartan 形式的平坦性（纤维丛视角）

在纤维丛 $Sp(2) \to \mathbb{HP}^1$（纤维 $U(2)$）中：

$$U = \sigma(q) \cdot \begin{pmatrix} u & 0 \\ 0 & v \end{pmatrix}$$

Maurer-Cartan 形式 $\omega_u = \operatorname{Im}(A + Bq)$ 是纤维联络。

> **判据**：若对某参数 $\lambda$ 有 $\partial \omega_u / \partial \lambda = 0$ 沿整条路径，
> 则 $\lambda$ 不贡献动力学相位（纤维联络在 $\lambda$ 方向**平坦**）。

**$E_1=0$ 时**：$\partial\omega_u/\partial t_1 = 0$（已证，见 report §7.5 和 fiber_orbit_framework.md §4）。
→ 改变 $t_1$ 不改变纤维 holonomy。但 $t_1$ 仍改变基空间轨迹 $q(t)$。

### 2.3 判据 C：李代数的理想结构

> 紧半单李代数 $\mathfrak{g}$ 可分解为非平凡理想的直和 $\Longleftrightarrow$ $\mathfrak{g}$ 不是**单李代数**。

| 李代数 | 单性 | 非平凡理想 | 可直和分解？ |
|--------|------|-----------|------------|
| $\mathfrak{so}(3) \cong \mathfrak{su}(2)$ | 单 | 无 | ✗ |
| $\mathfrak{so}(4) \cong \mathfrak{su}(2) \oplus \mathfrak{su}(2)$ | 非单 | 两个 $\mathfrak{su}(2)$ | ✓（抽象层面） |
| $\mathfrak{so}(5) \cong \mathfrak{sp}(2)$ | **单** | **无** | ✗ |
| $\mathfrak{u}(1) \oplus \mathfrak{su}(2)$ | 非单 | $\mathfrak{u}(1)$ 和 $\mathfrak{su}(2)$ | ✓ |

---

## 三、三段协议的逐段分析

### 3.1 Step 3：可分解 ✓（最优情况）

**活跃生成元**：$t_3$（$\Sigma_{34}$）、$E_d$（$\Sigma_{45}$）、$E_1$（$\Sigma_{12}$ 仅在 $E_1\neq 0$）

**李代数闭包**：$\mathfrak{u}(1) \oplus \mathfrak{su}(2)$

**分解**：
$$\underbrace{K_{\text{step3}}}_{\mathfrak{u}(1)\oplus\mathfrak{su}(2)} =
\underbrace{K_{E_d}}_{ \mathfrak{u}(1)} + \underbrace{K_{t_3}, K_{E_1}}_{\mathfrak{su}(2)}$$

$[\mathfrak{u}(1), \mathfrak{su}(2)] = 0$ 恒成立。

$$U_{\text{step3}} = \underbrace{e^{\int E_d(t)\,\Sigma_{45}\,dt}}_{\text{动力学相位 }U_{\text{dyn}}}
                    \cdot
                    \underbrace{\mathcal{T}\exp\!\left(\int (K_{t_3}+K_{E_1})dt\right)}_{\text{几何旋转 }U_{\text{geo}}}$$

**关键性质**：
- $U_{\text{dyn}}$ 是 ancilla 子空间的纯 U(1) 相位 → 不影响 MZM 子空间的编织结果
- $U_{\text{geo}}$ 携带完整的非阿贝尔几何信息
- **这是三段协议中唯一严格可分解的步骤**

### 3.2 Step 2：不可分解 ✗（最差情况）

**活跃生成元**：全部 5 项 → 李代数闭包 = $\mathfrak{so}(5)$（满 10 维）

**障碍**：$\mathfrak{so}(5)$ 是 **B₂ 型单李代数**，秩 2，**无任何非平凡理想**。

> 没有非平凡理想 → 无法将任何生成元分离到对易子闭子代数中 →
> **所有 5 项在时间编序指数中不可逆地混合。**

**物理含义**：
- $E_1$ 和 $t_1$ 的动力学效应与 $t_2, t_3, E_d$ 的几何编织**不可分离地纠缠**
- 这就是 Magnus 展开中高阶项不收敛、相互作用绘景中"微扰"不闭合的根本原因
- **不是计算方法的问题，是群论结构不允许**

**残余希望**（仅限 $E_1=0$ 时的弱分解）：

Riccati ODE $\dot q = C + [A,q] - qBq$（$A=D$）在 $E_1=0$ 时：
- $[A,q]$ 提供纯旋转（来自 $t_2, t_3$ 的几何驱动）
- $C$ 和 $qBq$ 提供 ancilla 耦合（来自 $t_1, E_d$）
- 三步积分后 $\int \omega_u^{(k)} dt = 0$ → 纤维 holonomy 不依赖 $t_1$

**这不是严格的 $U = G\cdot D$，而是**：在 SO(3) 投影（MZM 子空间）上，
$t_1$ 的影响等效于改变旋转轴和旋转角，但 $\sigma_z$ 方向净效应为零。
这是 report §7.5 中 $\hat n_z = 0$ 的深层代数原因。

### 3.3 Step 1：名义可分解但物理不对齐 ⚠️

**活跃生成元**：$E_1, t_2, t_1, E_d$ → 李代数闭包 = $\mathfrak{so}(4)$（6 维）

**表面上的好消息**：$\mathfrak{so}(4) \cong \mathfrak{su}(2) \oplus \mathfrak{su}(2)$ **是直和**。

**实际的问题**：物理生成元在 $\mathfrak{su}(2) \oplus \mathfrak{su}(2)$ 分解中的位置：

$$\mathfrak{so}(4) = \underbrace{\mathfrak{su}(2)_L}_{\text{左等倾旋转}} \oplus \underbrace{\mathfrak{su}(2)_R}_{\text{右等倾旋转}}$$

但 $K_{\text{geo}}$（含 $t_2, E_d$）和 $K_{\text{dyn}}$（含 $E_1, t_1$）的生成元
**在两个 $\mathfrak{su}(2)$ 因子上都有投影**——它们并不分别落在 $\mathfrak{su}(2)_L$ 和 $\mathfrak{su}(2)_R$ 中。

具体检验（在 Sp(2) 表示中）：

$$\begin{aligned}
\text{几何生成元：}& \quad \Sigma_{24}(t_2) = \begin{pmatrix}\mathbf j/2&0\\0&\mathbf j/2\end{pmatrix},\quad
\Sigma_{45}(E_d) = \begin{pmatrix}0&\mathbf k/2\\ \mathbf k/2&0\end{pmatrix} \\[4pt]
\text{动力学生成元：}& \quad \Sigma_{12}(E_1) = \begin{pmatrix}\mathbf i/2&0\\0&-\mathbf i/2\end{pmatrix},\quad
\Sigma_{15}(t_1) = \begin{pmatrix}0&-1/2\\1/2&0\end{pmatrix}
\end{aligned}$$

计算对易子：
$$[\Sigma_{24}, \Sigma_{12}] = -\Sigma_{14} \neq 0,\quad
[\Sigma_{45}, \Sigma_{15}] = \Sigma_{14} \neq 0$$

**几何生成元和动力学生成元不对易** → 时间编序指数**不可因子化**。

> 尽管 $\mathfrak{so}(4)$ 抽象层面可分解为直和，但物理上几何项和动力学项
> 跨两个 $\mathfrak{su}(2)$ 因子分布 → **无实用分解**。

---

## 四、分解条件总结

| 条件 | Step 1 | Step 2 | Step 3 |
|------|--------|--------|--------|
| 李代数闭包 | $\mathfrak{so}(4)$ | $\mathfrak{so}(5)$ | $\mathfrak{u}(1)\oplus\mathfrak{su}(2)$ |
| 单李代数？ | 否（直和） | **是** | 否 |
| $[\mathfrak{g}_{\text{geo}}, \mathfrak{g}_{\text{dyn}}] = 0$？ | ✗（不对齐） | ✗（无理想） | **✓** |
| 严格分解 $U = G\cdot D$ | ✗ | ✗ | **✓** |
| 弱分解（MZM 子空间） | $E_1=0$：$\hat n_z=0$ | $E_1=0$：$\hat n_z=0$ | 自动成立 |
| 解析解 | 仅 $E_1=0$ | 仅 $E_1=0$ | 可写闭式 |

**核心结论**：

1. **精确的 $U = G \cdot D$ 分解仅在 Step 3 成立**，条件是 $[\mathfrak{g}_{\text{geo}}, \mathfrak{g}_{\text{dyn}}] = 0$

2. **最根本的障碍是 $\mathfrak{so}(5)$ 的单性**——Step 2 中整个 10 维李代数被激活，
   没有任何非平凡理想可供分解。这是群论层面的不可约性，非任何计算技巧可绕过。

3. **$E_1=0$ 提供了弱分解**：虽然全空间 $U$ 不能因子化，但 MZM 子空间投影
   $R_{123} \in SO(3)$ 中动力学影响被限制在 $xy$ 平面（$\hat n_z=0$），
   根源是纤维联络 $\omega_u$ 在 $t_1$ 方向平坦。

---

## 五、不可分解时的动力学抵消策略

当 $\mathfrak{g}$ 是单李代数（Step 2）或物理生成元不对齐（Step 1）时，
无法通过代数分解消除动力学效应。以下策略从**群论/控制论**层面设计抵消方案。

### 5.1 策略一：复合脉冲补偿（已实现 ✓）

**原理**：在编织协议后追加补偿协议 $\bar\Gamma$，使 $U(\bar\Gamma) \cdot U(\Gamma) \sim U_{\text{ideal}}$。

**自由度**：$\bar\Gamma$ 含 6 个参数 $\{\bar\tau_1, \bar\tau_2, \bar\tau_3, \bar t_1^{(1)}, \bar t_1^{(2)}, \bar t_1^{(3)}\}$。

**匹配条件**：$U(\bar\Gamma) \cdot U(\Gamma)$ 和 $U_{\text{ideal}}$ 共轭等价 → 2 个本征值条件。

$6 \gg 2$ → 解空间是 4 维流形，问题**良定**。

**李代数解释**：$\bar\Gamma$ 的三个独立步在 $Sp(2)$ 中生成 6 维子流形。
当 $E_1 \neq 0$（$A \neq D$），完整 $\mathfrak{sp}(2)$ 被激活，补偿协议的
切空间满秩 → 可达所有需要的共轭类。

**数值验证**：`compensation_clean.py` 和 `compensation_range.py` 证实：
- $E_1/t_1 = 0.1 \sim 10$ 范围内，补偿后本征值偏差 $< 0.05$
- 改善倍数 $10\times \sim 50\times$

**群论视角**：
令 $U_0 = U(\Gamma)$，$\bar U = U(\bar\Gamma)$。线性化：
$$\bar U \cdot U_0 \approx (I + \sum_i \delta_i M_i) U_0
= U_0 (I + \sum_i \delta_i U_0^{-1} M_i U_0)$$

$M_i$ 是 $\partial \bar U / \partial p_i$ 在恒等元处的值。$U_0^{-1} M_i U_0$ 通过
伴随作用 $Ad_{U_0}$ 旋转到 $U_0$ 的"局部坐标系"。对 Cartan 子代数投影
即得 6×2 矩阵 $J$，$J\delta = b$ 的解 $\delta$ 给出补偿参数修正量。
$\operatorname{rank}(J) = 2$ 由 $E_1 \neq 0$ 保证。

### 5.2 策略二：零面积脉冲（E₁=0 时的自然抵消）

**原理**：利用门控函数 $f_\pm$ 在时间上的对称性。

$$\int_0^{3\tau} f_-(t)\,dt = \int_0^{3\tau} f_+(t)\,dt = \frac{3\tau}{2}$$

当 $E_1=0$ 时，$A=D$ 保证 Riccati 结构中正反向 ancilla 通道完全对称。
$\omega_u^{(k)} = \frac{|t_1|}{2}q_3 + \frac{E_d}{2}q_0$ 在三步脉冲下的时间积分
自动为零——**无需主动设计，是系统对称性的自然结果**。

**推广**：对于一般的 $t_1(t)$ 包络，若设计 $\int_0^{3\tau} t_1(t) f(t) dt = 0$
（其中 $f(t)$ 是合适的权重函数），动力学效应可被**部分抵消**。

### 5.3 策略三：动力学去耦（Dynamical Decoupling, DD）

**原理**：在编织过程中插入快速 $\pi$ 脉冲，使寄生哈密顿量 $H_{\text{dyn}}$ 在
脉冲间隔内被平均为零。

**在 Sp(2) 表示中的实现**：

寄生项 $H_{\text{dyn}} = iE_1\gamma_1\gamma_2 - i|t_1|\gamma_b\gamma_1$ 对应的生成元：
- $\Sigma_{12}$（$E_1$）：对角块，MZM 子空间的 $\sigma_z$
- $\Sigma_{15}$（$t_1$）：非对角块，ancilla-MZM 耦合

DD 序列设计：
1. 选择去耦群 $\mathcal{G}_{DD} \subset Sp(2)$，使 $\frac{1}{|\mathcal{G}_{DD}|}\sum_{g} g H_{\text{dyn}} g^{-1} = 0$
2. 以周期 $\Delta t \ll 1/||H_{\text{dyn}}||$ 施加群元素脉冲
3. 在脉冲间隔内 $H_{\text{geo}}$（编织驱动）正常演化

**候选去耦群**：
$$\mathcal{G}_{DD} = \{I, e^{\pi \Sigma_{45}}\} = \{I, e^{\pi E_d \Sigma_{45}}\}$$

$e^{\pi \Sigma_{45}}$ 将 $\Sigma_{12} \mapsto -\Sigma_{12}$，$\Sigma_{15} \mapsto -\Sigma_{15}$，
实现对寄生项的一阶平均消除。

**局限**：DD 脉冲与编织脉冲 $t_2, t_3$ 不对易 → 需要精心设计时间顺序以避免
干扰几何编织。DD 频率必须 $\gg$ 编织速度 → 需要 $\tau \gg 1/||H_{\text{dyn}}||$。

### 5.4 策略四：对称回波（Echo）

**原理**：正向协议 + 反向协议，利用时间反演对称性。

若 $H_{\text{dyn}}$ 在时间反演下奇对称（或可通过规范变换实现），则：
$$U_{\text{forward}} \cdot U_{\text{backward}} = U_{\text{geo}}^2$$

动力学相位在往返过程中抵消。

**在编织协议中的实现**：
- 正向：$\gamma_2 \to \gamma_3$（标准三步协议）
- 反向：$\gamma_3 \to \gamma_2$（交换 $t_2$ 和 $t_3$ 的角色）

两次编织后的净效应：
$$U_{\text{echo}} = U(\gamma_3\leftarrow\gamma_2) \cdot U(\gamma_2\leftarrow\gamma_3)$$

若 $\gamma_1$ 的耦合 $E_1$ 在两次编织间不变号，则动力学相位在 $U_{\text{echo}}$ 中
以对易子形式残留（而非直接相加），可显著抑制其影响。

### 5.5 策略五：最优化控制（Quantum Optimal Control）

**原理**：将抵消问题表述为泛函优化。

$$\min_{t_2(t), t_3(t), E_d(t), t_1(t)} \; \|U(T) - U_{\text{target}}\|^2 + \lambda \mathcal{R}[controls]$$

使用 GRAPE / CRAB / Krotov 等方法在 $\mathfrak{sp}(2)$ 流形上搜索最优控制脉冲。

**优势**：
- 不依赖于代数分解假设
- 可同时优化保真度和鲁棒性（对 $E_1$ 涨落不敏感）

**代价**：
- 数值求解，无解析理解
- 需要梯度的伴随方法来处理 10 维流形

### 5.6 策略六：余伴随轨道跃迁操控

**原理**（来自 `fiber_orbit_framework.md` §5 和 `math_tools_for_SO5_evolution.md` §3）：

改变 $E_1/|t_2|$ 的**比值**而非绝对值，控制 $U(3\tau)$ 的共轭类（余伴随轨道类型）：

| 比值范围 | 共轭类 | 本征值简并度 | 稳定性 |
|---------|--------|------------|--------|
| $E_1 \ll |t_2|,|t_3|$ | 特殊轨道 | 2（二重简并） | 对 $t_1$ 涨落不敏感 |
| $E_1 \sim |t_2|,|t_3|$ | **一般轨道** | 1（全非简并） | 敏感 |
| $E_1 \gg |t_2|,|t_3|$ | 特殊轨道（恢复） | 2 | RWA 对称恢复 |

**策略**：操作 $t_2, t_3$（它们是可控的！）使 $E_1/|t_{2,3}| \ll 1$，
保持系统在**特殊轨道**上。此时即使 $E_1 \neq 0$，只要 $E_1$ 相对于编织耦合足够小，
$t_1$ 的影响仍然被几何抑制（flat direction）。

**这意味着**：不需要 $E_1$ 严格为零，只需 $E_1 \ll t_c$（编织耦合强度）。
当 $E_1/t_c \lesssim 0.1$ 时，动力学相位被有效抑制。

---

## 六、策略对比与推荐

| 策略 | 代数前提 | 实验难度 | $\tau$ 要求 | 对 $E_1\neq 0$ 有效？ | 解析可控？ |
|------|---------|---------|------------|---------------------|-----------|
| **复合脉冲补偿** | 无 | 中（需独立控制三步 $t_1$） | 任意 | **✓** | 半解析 |
| **零面积脉冲** | $E_1=0$（$A=D$） | 低（自动满足） | 任意 | ✗ | **✓** |
| **动力学去耦** | 去耦群存在 | 高（需快速脉冲） | $\tau \gg$ | **✓** | 微扰 |
| **对称回波** | 时间反演对称 | 低（只需反转协议） | 任意 | 部分 | **✓** |
| **最优化控制** | 无 | 高（需梯度优化） | 任意 | **✓** | ✗ |
| **轨道操控** | 无 | 低（增大 $t_c$） | 任意 | **✓**（$E_1/t_c \ll 1$） | **✓** |

### 推荐方案（按优先级）

1. **增大 $t_c$**（轨道操控）：保持 $E_1/t_c \lesssim 0.1$，最简单且最有效
2. **复合脉冲补偿**：当 $E_1/t_c$ 无法降低时，用 6 参数补偿协议
3. **对称回波**：如果只需验证非阿贝尔统计（而非实现特定量子门），
   双次编织天然提供了回波结构

---

## 七、补充：B₂ 根系统的深层解释

$\mathfrak{so}(5) \cong \mathfrak{sp}(2)$ 是 B₂ = C₂ 型，秩 2，8 个根：

```
         α₁+2α₂  (短根，最高根 ↔ σ_z)
              ↑
 α₁ ←──●──→ α₁+α₂ ←──●──→ 2α₁+α₂
              ↑
             α₂  (长根 ↔ σ_x, σ_y 平面)
```

- **2 个短根**：$\pm\alpha_1, \pm(\alpha_1+2\alpha_2)$ → $i\gamma_1\gamma_2$（$\sigma_z$）
- **6 个长根**：$\pm\alpha_2, \pm(\alpha_1+\alpha_2), \pm(2\alpha_1+\alpha_2)$ → $i\gamma_2\gamma_3, i\gamma_3\gamma_1, \dots$

**$E_1=0$ 意味着短根方向的生成元 $\Sigma_{12}$ 未被激活** →
Weyl 群作用仅在长根平面内 → 旋转轴锁定 $xy$ 平面。

**$E_1 \neq 0$ 激活短根** → 满 B₂ 根系统 → 所有 3 个 Pauli 方向释放。

这就是 **Cartan 分解**视角下的根本区别：

$$\mathfrak{so}(5) = \underbrace{\mathfrak{so}(4)}_{\mathfrak{k} \text{ (极大紧子代数)}} \oplus \underbrace{\mathfrak{p}}_{\text{4 个"真正非阿贝尔"方向}}$$

$E_1=0$：$K(t) \in \mathfrak{k} = \mathfrak{so}(4)$ → 动力学被约束在极大紧子代数内
$E_1\neq 0$：$K(t) \in \mathfrak{k} \oplus \mathfrak{p} = \mathfrak{so}(5)$ → 满代数，无约束

---

## 八、结论

### 可直接回答的三个问题

**Q1: SO(5)/Sp(2) 能分解出几何项 + 动力学项吗？**

**A**: **仅在特定子代数中可以。**
- Step 3（$\mathfrak{u}(1)\oplus\mathfrak{su}(2)$）：**可以严格分解**，因为 $[\mathfrak{g}_{\text{geo}}, \mathfrak{g}_{\text{dyn}}] = 0$
- Step 1（$\mathfrak{so}(4)$）：**名义上可以**（$\mathfrak{so}(4)$ 是直和），但**物理生成元不对齐**，几何项和动力学项跨两个 $\mathfrak{su}(2)$ 因子分布，不可因子化
- Step 2（$\mathfrak{so}(5)$）：**绝对不可以**——$\mathfrak{so}(5)$ 是单李代数，无任何非平凡理想

**Q2: 分解的充要条件是什么？**

**A**: 充要条件是 $[\mathfrak{g}_{\text{geo}}, \mathfrak{g}_{\text{dyn}}] = 0$。
- 充分性：对易子为零 → 时间编序指数可因子化
- 必要性：若不对易，Magnus 展开中 Baker-Campbell-Hausdorff 级数产生无穷多交叉项

等价地，在纤维丛视角下：纤维联络 $\omega_u$ 在动力学参数方向**平坦**（$\partial\omega_u/\partial\lambda = 0$）。

**Q3: 不能分解时如何抵消动力学项？**

**A**: 六种策略，推荐优先级：
1. **轨道操控**（增大 $t_c$ 使 $E_1/t_c \ll 1$）——最简单
2. **复合脉冲补偿**（6 参数 $\bar\Gamma$ 协议）——最通用
3. **对称回波**（往返编织）——天然抵消
4. **动力学去耦**（快速 $\pi$ 脉冲）——需要高速控制
5. **最优化控制**（GRAPE/CRAB）——数值方法
6. **零面积脉冲设计**——仅 $E_1=0$ 时自动满足

### 最核心的群论事实

> **$\mathfrak{so}(5)$ 是 B₂ 型单李代数。它没有非平凡理想 → 
> 代数层面的几何-动力学分离在 Step 2 中原则上不可能。**
> 
> 这是群论的基本事实，不是近似、截断或计算方法的问题。
> 所有"抵消"策略本质上都是在更大的参数空间中寻找路径，
> 使 $U_{\text{total}}$ 落在所需共轭类上，而非真的从代数中"去除"动力学生成元。
