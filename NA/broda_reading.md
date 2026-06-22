# Broda 非 Abel Stokes 定理 —— 详细解读

> 基于 `NA/na0.txt`（hep-th/9511150，10页原始版本）
> 和 `NA/na1.txt`（math-ph/0012035，46页综述版本）
> 的完整内容。

---

## 一、定理要解决什么问题？

### 1.1 Abel 版本回顾

标准 Stokes 定理（Abel 版本）：

$$\int_{\partial M} \omega = \int_M d\omega$$

对 $U(1)$ 规范场：$F = dA$，所以 $\oint_C A = \int_S F$。

全局版本（用 holonomy）：

$$\exp\left(i\oint_C A\right) = \exp\left(i\int_S F\right)$$

### 1.2 非 Abel 版本的核心困难

对 Yang-Mills 场：$F = dA - iA \wedge A$，场强不是恰当形式。且积分值取在 Lie 代数 $\mathfrak{g}$ 中——不同点的值**不对易**。

Broda 的目标：将路径编序指数

$$W[C] = \mathcal{P}\exp\left(i\oint_C A\right) \in G$$

表达为曲面 $S$（$\partial S = C$）上的某种积分。

---

## 二、三种等价的表述形式

Broda 在综述（na1）中明确区分了三种方案：

| 方案 | 核心思想 | 优点 | 缺点 |
|------|---------|------|------|
| **算子形式** | 直接定义曲面编序 $\mathcal{P}_S$ | 最直观，与 Abel 形式最接近 | 曲面编序难以严格定义 |
| **相干态路径积分** | 用群相干态改写 Wilson loop | 几何直觉清晰，适合半经典近似 | 依赖参考态选择 |
| **全纯路径积分** | 通过 Fock 空间实现 Lie 代数 | 最系统，给出严格的泛函积分 | 涉及 BRS 量子化，较抽象 |

三种方案**等价**——它们给出同一个 Wilson loop 的不同表示。

---

## 三、算子形式

### 3.1 表述

$$\boxed{\mathcal{P} \exp\left(i \oint_{\partial S} A_i dx^i\right) = \mathcal{P}_S \exp\left(\frac{i}{2} \int_S \tilde{F}_{ij} dx^i \wedge dx^j\right)}$$

其中 $\mathcal{P}_S$ 是**曲面编序**算子，$\tilde{F}_{ij}(x)$ 是**路径依赖的曲率**：

$$\tilde{F}_{ij}(x) = U^{-1}(x, O) F_{ij}(x) U(x, O)$$

这里 $U(x, O)$ 是从基点 $O$ 沿曲面 $S$ 内的路径 $\mathcal{L}$ 到点 $x$ 的平行移动算符：

$$U(x, O) = \mathcal{P} \exp\left(i \int_{\mathcal{L}} A_i dy^i\right)$$

### 3.2 几何含义

- $F_{ij}$ 是局部曲率，但 $\tilde{F}_{ij}$ 通过平行移动把所有点的曲率"拉回"到基点 $O$ 的 Lie 代数中
- 曲面编序 $\mathcal{P}_S$ 对分格后的面元按偏序排列——这是对路径编序 $\mathcal{P}$ 的二维推广
- 规范不变性：$\tilde{F}$ 的定义保证了曲面上不同点的值可以放在同一个 Lie 代数中比较和排序

### 3.3 关键困难

$\mathcal{P}_S$ 需要定义曲面上的偏序关系。对于平面回路这相对简单，但对于打结回路（knotted loops），Broda 与合作者 Duniec 在 `math-ph/0109028` 中专门处理了这个问题——需要将曲面分解为 Seifert 曲面并谨慎定义面元排序。

---

## 四、相干态路径积分形式

### 4.1 核心思想

用规范群 $G$ 的相干态将 Wilson loop 重写为路径积分，然后利用 Abelian Stokes 定理。

### 4.2 相干态构造（五步标准流程）

**第1步**：取 $\mathfrak{g}$ 的 Cartan 基 $\{H_i, E_\alpha, E_{-\alpha}\}$

**第2步**：选一个幺正不可约表示 $R$，取最高权态 $|R\rangle$（满足 $E_\alpha |R\rangle = 0$）

**第3步**：最大稳定子群 $H$：$h|R\rangle = |R\rangle e^{i\phi(h)}$，$h \in H$

**第4步**：陪集分解 $g = \xi h$，$\xi \in G/H$

**第5步**：相干态 $|\xi, R\rangle = \xi |R\rangle$

相干态的核心性质：
- 归一化：$\langle g,R|g,R\rangle = 1$
- **单位分解**：$\int |g,R\rangle d\mu(g) \langle g,R| = I$（这是路径积分的基础）

### 4.3 路径积分推导

相干态之间的 transition amplitude：

$$\langle g'',R | \mathcal{P}\exp\left(i \int_{t'}^{t''} A(t) dt\right) | g',R \rangle$$

利用单位分解 $N$ 次插入，取 $N \to \infty$ 极限，得到泛函积分。

单个无穷小步的振幅：

$$\langle g_n,R|(1 + i\epsilon A_n)|g_{n-1},R\rangle = \exp\left(\langle R| -g^\dagger \dot g + i g^\dagger A g |R\rangle \epsilon + O(\epsilon^2)\right)$$

### 4.4 最终公式

Wilson loop 的相干态路径积分形式：

$$W[C] = \int \mathcal{D}\mu(g) \exp\left(i \oint_C \mathcal{L}\right)$$

其中 Lagrangian 为：

$$\boxed{\mathcal{L} = \frac{1}{\kappa} \text{Tr}\left[m \cdot H \; A_g(t)\right]}$$

而 $A_g(t) = i g^\dagger \dot g + g^\dagger A g$ 是**规范变换后的联络**。

**关键**：$\mathcal{L} = B_i dx^i$ 是一个 **Abelian 微分形式**（标量值！）。因此可以直接应用标准（Abel）Stokes 定理：

$$\oint_C B = \int_S dB$$

于是：

$$\boxed{W[C] = \int \mathcal{D}\mu(g) \exp\left(i \int_S dB\right)}$$

### 4.5 物理含义

- $g(t) \in G$ 是辅助的"量子力学"自由度，定义在回路 $C$ 上
- $\mathcal{D}\mu(g)$ 是 $G/H$ 上的 Haar 测度
- $m \cdot H$ 编码了表示 $R$ 的最高权——不同的表示给出不同的 $m_i$
- $dB$ 涉及曲率 $F$ 和联络 $A$——**这就是为什么"路径依赖曲率"出现在算子形式中**

---

## 五、全纯路径积分形式（最系统化）

这是 Broda 最精心发展的方案，也是与我们的 Sp(2) Riccati ODE 联系最紧密的方案。

### 5.1 核心思路

1. 在辅助 **Fock 空间**中实现 Lie 代数 $\mathfrak{g}$ 的表示
2. 将 parallel-transport operator 表达为全纯路径积分
3. 将积分从回路 $C$ "膨胀"到曲面 $S$
4. 用 2D 拓扑量子场论 + BRS 量子化处理曲面上的冗余自由度
5. 用 **Abelian Stokes 定理**完成证明

### 5.2 辅助 Fock 空间实现

产生/湮灭算符：$[a_k, a_l^\dagger]_\mp = \delta_{kl}$

Lie 代数生成元的 Fock 实现：
$$\hat{T}^a = T^a_{kl} a_k^\dagger a_l$$

验证：$[\hat{T}^a, \hat{T}^b] = i f^{abc} \hat{T}^c$（忠实表示）

**关键**：粒子数算符 $\hat{N} = \delta_{kl} a_k^\dagger a_l$ 与 Hamiltonian 对易：
$$[\hat{N}, \hat{H}] = 0$$

这意味着**如果初始态在单粒子子空间，演化永远留在单粒子子空间**。

### 5.3 辅助力学问题

经典 Lagrangian：

$$L(z, \bar{z}) = i \bar{z} D_t z$$

其中 $D_t = \frac{d}{dt} - i \dot{x}^i A_i^a T^a$。

对应的 Hamiltonian（在 Fock 空间中）：

$$\hat{H}(t) = -\dot{x}^i(t) A_i^a[x(t)] T^a_{kl} a_k^\dagger a_l$$

### 5.4 全纯路径积分表示的 Parallel-Transport Operator

单粒子态之间的 transition amplitude：

$$\boxed{U_{rs}(t'', t') = \int \mathcal{D}\bar{z} \mathcal{D}z \; \exp\left(-\bar{z}(t')z(t') + i\int_{t'}^{t''} L[\bar{z}, z] dt\right) \; z_r(t'')\bar{z}_s(t')}$$

这里 $z_k(t)$ 是全纯坐标（Bargmann-Fock 表示中的相干态参数）。

### 5.5 曲面上的 2D 拓扑场论

在曲面 $S$（参数化 $x^i(\sigma^1, \sigma^2)$，$t' \le \sigma^1 = t \le t''$，$0 \le \sigma^2 = s \le 1$）上定义作用量：

$$\boxed{S_{\text{cl}} = \int_S \epsilon^{AB} \left(i D_A \bar{\psi} D_B \psi + \frac{1}{2} \bar{\psi} F_{AB} \psi\right) d^2\sigma}$$

其中 $D_A = \partial_A x^i D_i$，$F_{AB} = \partial_A x^i \partial_B x^j F_{ij}$。

这个理论具有**拓扑规范对称性**：

$$\delta\psi(x) = \theta(x), \quad \delta\bar{\psi}(x) = \xi(x)$$

其中 $\theta, \xi$ 在内部任意但在边界 $\partial S = C$ 上为零。

**关键性质**：利用分部积分和 **Abelian Stokes 定理**，这个 2D 作用量退化为边界上的 1D 作用量：

$$S_{\text{cl}} = i \int_{\partial S} \bar{\psi} D_t \psi dt$$

即 $S_{\text{cl}}$ 在边界上与原始 Lagrangian $L$ 一致！

### 5.6 BRS 量子化

由于理论有规范对称性（$\delta\psi = \theta$，$\delta\bar\psi = \xi$），需要规范固定。

引入 BRS 算子 $s$（$s^2 = 0$）和鬼场 $\phi, \chi, \bar\phi, \bar\chi$，规范固定项：

$$S' = s\left(\int_S (\bar\phi \Delta\psi \pm \bar\psi \Delta\chi) d^2\sigma\right)$$

这给出约束 $\Delta\psi = 0$，$\Delta\bar\psi = 0$——即 $\psi$ 和 $\bar\psi$ 在曲面内部是调和的，完全由边界值决定。**这消除了曲面内部的冗余泛函积分。**

### 5.7 最终的全纯非 Abel Stokes 定理

$$\boxed{\int \mathcal{D}^2 z \; z_k z_\ell \exp\left(-\bar{z}z + i\oint_C i\bar{z}Dz\right) = \int \mathcal{D}^2 z \; z_k z_\ell \exp\left(-\bar{z}z + i S_{\text{cl}}\right)}$$

或用参数化形式：

$$\boxed{\begin{aligned}
&\int \mathcal{D}^2 z \; z_k(t'')\bar{z}_\ell(t') \exp\left(-\bar{z}(t')z(t') + i\int_{t'}^{t''} L[\bar{z}, z] dt\right) \\
= &\int \mathcal{D}^2 z \; z_k(t'', 1)\bar{z}_\ell(t', 1) \exp\left(-\bar{z}(t', 1)z(t', 1) + i\int_0^1 \int_{t'}^{t''} L_{\text{cl}} dt ds\right)
\end{aligned}}$$

**等号两边**：
- **左边**：回路 $C$ 上的全纯路径积分 = parallel-transport operator（Wilson loop 的矩阵元）
- **右边**：曲面 $S$ 上的全纯路径积分，被积函数是 2D 经典作用量

**积分测度集中在边界上**，两边的测度完全相同。

### 5.8 证明逻辑链

```
辅助 Fock 空间中的 Schrödinger 问题
    ↓
Parallel-transport operator 的全纯路径积分表示 (1D，在 C 上)
    ↓
2D 拓扑场论在曲面 S 上的构造 (S_cl 在边界 = L)
    ↓
BRS 量子化 → 消除曲面内部的冗余自由度
    ↓
Abelian Stokes 定理 → 证明 LHS (1D) = RHS (2D)
    ↓
非 Abel Stokes 定理成立
```

**核心技巧**：整个证明中，非 Abel 结构被编码在 Fock 空间的算符实现中。曲面 $S$ 上的作用量**本身不涉及非 Abel 编序**——编序自动被全纯路径积分中的 $z_k(t'')\bar{z}_\ell(t')$ 因子处理。Abelian Stokes 定理只在**最后一步**被使用——将边界上的全纯作用量等价于曲面上的 2D 作用量。

---

## 六、三种方案的比较与联系

### 6.1 对应关系

| | 算子形式 | 相干态 PI | 全纯 PI |
|---|---|---|---|
| **核心对象** | $\mathcal{P}_S \exp(\int \tilde{F})$ | $\int \mathcal{D}\mu(g) e^{i\oint B}$ | $\int \mathcal{D}^2z \, z_k z_\ell e^{iS_{\text{cl}}}$ |
| **辅助自由度** | 无（直接处理） | $g \in G/H$（陪集空间） | $z_k \in \mathbb{C}^{\dim R}$（全纯坐标） |
| **"Abel化"方式** | 路径依赖曲率 $\tilde{F}$ | $\mathcal{L} = \text{Tr}[m\cdot H A_g]$ 是标量 | 2D 拓扑 QFT + BRS |
| **Stokes 用在哪里** | 在"曲面编序"的定义中 | $\oint B = \int dB$ | $S_{\text{cl}}|_{\partial S} = i\oint \bar\psi D\psi$ |

### 6.2 统一的数学本质

三种方案的共同点：
1. 都用一个**扩张的构造**（路径依赖曲率 / 相干态 / Fock 空间）来处理非交换性
2. 都最终依赖于 **Abelian Stokes 定理**——非 Abel → Abel 的桥接
3. 都给出 Wilson loop 的等价表示

区别在于"非 Abel 信息存储在哪里"：
- 算子形式：存储在 $\tilde{F}$ 的路径依赖性中
- 相干态 PI：存储在 $g(t)$ 的路径积分和 $m \cdot H$ 中
- 全纯 PI：存储在 $z_k(t'')\bar{z}_\ell(t')$ 的编序因子中

---

## 七、推广：打结回路

Broda 与 Duniec（math-ph/0109028）将定理推广到打结回路。核心困难：

- 打结回路的 Seifert 曲面不是唯一的
- 曲面自相交时面元排序不明确

解决方案：使用 Seifert 曲面的分层分解，对每层定义独立的偏序。

---

## 八、应用（综述中介绍的部分）

### 8.1 2D Yang-Mills 理论

非 Abel Stokes 定理给出 Wilson loop 期望值的严格表达式，可用于验证 Migdal-Makeenko 圈方程。

### 8.2 Chern-Simons 拓扑场论

这是**最重要的应用**。Chern-Simons 理论的 Wilson loop 期望值给出扭结不变量（Jones 多项式等）。非 Abel Stokes 定理的非微扰性质对于推导 skein 关系至关重要。

Broda 自己在这个方向有大量工作（1990-1994）。路径积分形式允许将 Wilson loop 表达为曲面上的泛函积分，进而用微扰/非微扰方法计算。

### 8.3 高维推广

Broda 提出了一般配方：
1. 在 $\partial M$（$d$ 维流形 $M$ 的边界）上构造外场 $A$ 中的拓扑场论
2. 量子化该理论 → 得到路径积分形式的配分函数
3. 应用 Abelian Stokes 定理将作用量从 $\partial M$ 推到 $M$ 内部
4. 得到高维非 Abel Stokes 定理

---

## 九、与我们的 Sp(2) Riccati ODE 的联系（预告）

Broda 的全纯路径积分中的核心量是：

$$U_{rs} = \int \mathcal{D}^2z \; z_r \bar{z}_s \; \exp(iS_{\text{cl}})$$

对于我们的 Sp(2) 系统：
- $r,s = 1,2$（2 维四元数表示 → 4 维复表示）
- $a_k, a_k^\dagger$ 是 4 对产生/湮灭算符
- $z = (z_1, z_2, z_3, z_4) \in \mathbb{C}^4 \cong \mathbb{H}^2$
- 单粒子子空间 = 4 维复 = 我们的旋量空间

Riccati 变量 $q \in \mathbb{H}$ 自然出现在全纯坐标的比值中：

$$q = z_{\text{ancilla}} / z_{\text{MZM}}$$

而 Broda 的 2D 经典作用量 $S_{\text{cl}}$ 在我们的设定中简化为：

$$S_{\text{cl}} \to \int_0^{3\tau} \text{Tr}[K(t) \cdot \Omega(q(t))] dt$$

其中 $\Omega(q) = |q\rangle\langle q|$ 是 Sp(2) 相干态的投影算符。**鞍点条件 $\delta S_{\text{cl}}/\delta \bar{q} = 0$ 正是我们的 Riccati ODE。**

（后续文档将详细展开这个对应。）
