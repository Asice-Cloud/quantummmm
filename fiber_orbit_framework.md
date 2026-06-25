# 5-Majorana 编织的纤维丛-余伴随轨道描述

> 从 `ode1/report.md` 的物理系统和 Sp(2) 表示出发，
> 用纤维丛和余伴随轨道重新表述多体演化。

---

## 一、物理系统

### 1.1 5 个 Majorana 模式

$\gamma_1, \gamma_2, \gamma_3$：MZM，$\gamma_a, \gamma_b$：ancilla。编织目标：交换 $\gamma_2 \leftrightarrow \gamma_3$。

$$H_{EM} = iE_d\gamma_a\gamma_b + iE_1\gamma_1\gamma_2 + i|t_2|\gamma_a\gamma_2 - i|t_1|\gamma_b\gamma_1 - i|t_3|\gamma_a\gamma_3$$

| 参数 | 物理含义 | 可控？ |
|------|---------|-------|
| $E_1$ | $\gamma_1$–$\gamma_2$ 杂化能 | ❌ 材料常数 |
| $t_1$ | $\gamma_1$–ancilla 耦合 | ⚠️ 部分可控 |
| $t_2, t_3$ | 门控编织耦合 | ✅ 时变脉冲 |
| $E_d$ | 量子点能级 | ✅ 时变脉冲 |

### 1.2 三段编织协议

门控函数 $f_\pm(t) = (1 \pm \cos(\pi t/\tau))/2$，三步各持续 $\tau$：

- Step 1：$t_2$ 打开，$E_d \to 0$ → $\gamma_2$ 移入量子点
- Step 2：$t_3$ 打开，$t_2$ 关闭 → $\gamma_3$ 接管
- Step 3：$t_3$ 关闭，$E_d$ 恢复 → 回到初态，$\gamma_2 \leftrightarrow \gamma_3$ 交换

### 1.3 李代数结构

双线性生成元 $X_{ij} = i\gamma_i\gamma_j$（$C_2^5 = 10$ 个）满足 $\mathfrak{so}(5)$ 对易关系。

$$\mathfrak{so}(5) \cong \mathfrak{sp}(2) \quad\text{(10 维实简单李代数同构)}$$

$Sp(2) \cong Spin(5)$ 是 $SO(5)$ 的双重覆盖。

---

## 二、Sp(2) 旋量表示

### 2.1 四元数分块

演化算符 $U(t) \in Sp(2)$ 满足 $\dot U = KU$，$K \in \mathfrak{sp}(2)$：

$$K = \begin{pmatrix} A & B \\ C & D \end{pmatrix}, \qquad U = \begin{pmatrix} X & Y \\ Z & W \end{pmatrix}$$

$A, D$ 为纯虚四元数（对应对角块，MZM 和 ancilla 各自的 $\mathfrak{su}(2)$ 旋转）。
$B, C$ 为一般四元数（对应非对角块，ancilla–MZM 耦合），满足 $\mathfrak{sp}(2)$ 条件 $C = -\bar B$。

### 2.2 物理参数对应的分块

$$\boxed{\begin{aligned}
A(t) &= \frac{E_1 + |t_3|}{2}\,\mathbf i + \frac{|t_2|}{2}\,\mathbf j \\[4pt]
D(t) &= \frac{-E_1 + |t_3|}{2}\,\mathbf i + \frac{|t_2|}{2}\,\mathbf j \\[4pt]
B(t) &= \frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf k \\[4pt]
C(t) &= -\frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf k
\end{aligned}}$$

**关键**：$A - D = E_1\mathbf{i}$。$E_1 = 0 \Leftrightarrow A = D$。

---

## 三、纤维丛分解：$Sp(2) \to \mathbb{HP}^1$

### 3.1 主丛结构

$Sp(2)$ 是总空间。基空间 $\mathbb{HP}^1$（四元数射影线，4 维）由 $q = ZX^{-1} \in \mathbb{H}$ 参数化。
纤维 $U(2)$（6 维）是 ancilla 和 MZM 的内部相位空间。

任意 $U \in Sp(2)$ 唯一分解为：

$$U = \sigma(q) \cdot h, \qquad \sigma(q) = \frac{1}{\sqrt{1+|q|^2}}\begin{pmatrix} I & -q^\dagger \\ q & I \end{pmatrix}, \qquad h = \begin{pmatrix} u & 0 \\ 0 & v \end{pmatrix}$$

- $\sigma(q)$：截面（基空间提升）
- $u, v \in Sp(1)$：单位四元数，纤维上的运动

### 3.2 MZM 纤维的 Maurer-Cartan 形式

**定理**：$\omega_u \equiv \dot u u^{-1} = \operatorname{Im}(A + Bq)$。

**证明**：
$u = X/|X|$。$\dot u u^{-1} = K_{\text{eff}} - \frac{|\dot X|}{|X|}I$。
$u$ 是单位四元数 $\Rightarrow$ 标量部分恒零 $\Rightarrow$ $|\dot X|/|X| = \operatorname{Re}(K_{\text{eff}})$。
此等式由 $\mathfrak{sp}(2)$ 条件 $C = -\bar B$ 保证（sympy 验证）。
$\therefore \omega_u = \operatorname{Im}(K_{\text{eff}})$。$\square$

### 3.3 分量展开

$$\boxed{\begin{aligned}
\omega_u^{(i)} &= \frac{E_1 + |t_3|}{2} + \frac{|t_1|}{2}q_1 - \frac{E_d}{2}q_2 \quad (\sigma_x)\\[4pt]
\omega_u^{(j)} &= \frac{|t_2|}{2} + \frac{|t_1|}{2}q_2 + \frac{E_d}{2}q_1 \quad (\sigma_y)\\[4pt]
\omega_u^{(k)} &= \frac{|t_1|}{2}q_3 + \frac{E_d}{2}q_0 \quad (\sigma_z)
\end{aligned}}$$

**$\omega_u^{(k)}$ 不显含 $E_1$**。$E_1$ 仅出现在 $\omega_u^{(i)}$ 中。

ancilla 纤维同理：$\omega_v = \operatorname{Im}(D + Cp)$，其中 $p = YW^{-1}$。

---

## 四、$E_1 = 0$：平坦联络

### 4.1 平坦的表述

$E_1 = 0 \Rightarrow A = D$。$\omega_u^{(k)} = \frac{|t_1|}{2}q_3 + \frac{E_d}{2}q_0$。

Riccati 方程 $\dot q = C + [A,q] - qBq$ 保证 $\dot\alpha/|t_2| = \dot\beta/|t_3|$，
从而 $|t_3|\beta - |t_2|\alpha \equiv 0$。三步脉冲对称性（$\int f_m = \int f_p$）导致：

$$\int_0^{3\tau} \omega_u^{(k)}(t)\,dt = 0, \quad \forall\, t_1$$

$$\boxed{\left.\frac{\partial\theta_u}{\partial t_1}\right|_{E_1=0} = 0}$$

**纤维联络在 $t_1$ 方向平坦**——改变 $t_1$ 不改变纤维 holonomy。这是 $\sigma_z$ 保护的深层几何根源。

### 4.2 $E_1 > 0$：联络弯曲

$A \neq D \Rightarrow \dot\alpha/|t_2| \neq \dot\beta/|t_3|$。$\int \omega_u^{(k)} dt \neq 0$ 且依赖 $t_1$。

弯曲的结构形式：
$$\theta_u(E_1, t_1) = \theta_u^0 + E_1 \cdot f(t_1) + O(E_1^2)$$

$f(0) = 0$，$f(t_1)$ 无解析闭式（非可积系统），可从数值提取。

---

## 五、余伴随轨道：共轭类分类

### 5.1 SO(5) 的轨道分层

$\xi \in \mathfrak{so}(5)$，轨道 $O_\xi = \{g\xi g^{-1} : g \in SO(5)\}$。由两个特征值 $(\lambda_1, \lambda_2)$ 分类：

| 特征值 | 稳定子 | 轨道维数 |
|--------|--------|---------|
| $(0, 0)$ | $SO(5)$ (10D) | 0 |
| $(\lambda, 0)$ | $U(2)$ (4D) | 6 |
| $(\lambda_1, \lambda_2)$, $\lambda_1 \neq \lambda_2$ | $T^2$ (2D) | **8**（一般轨道） |

### 5.2 与 $U(3\tau)$ 本征值的对应

$U(3\tau) \in Sp(2)$ 的 $4 \times 4$ 复矩阵本征值成对出现 $(\lambda, \bar\lambda, \lambda^{-1}, \bar\lambda^{-1})$。

本征值的简并度决定 $U(3\tau)$ 的共轭类：

| 独立本征值个数 | 稳定子 | 共轭类维数 | 对应轨道类型 |
|-------------|--------|----------|------------|
| 2（各二重简并） | $U(2)$ (4D) | 6 | 特殊 |
| 4（全部不简并） | $T^2$ (2D) | **8** | **一般** |
| 3（一个二重简并） | 介于两者之间 | 7 | 过渡 |

### 5.3 轨道跃迁

| $E_1$ 相对 $|t_2|,|t_3|$ | 独立本征值 | 共轭类 | 含义 |
|---|---|---|---|---|---|
| $E_1 = 0$ | 2（二重简并） | 特殊 | $t_1$ 方向平坦 |
| $E_1 \ll |t_2|,|t_3|$ | 2 | 特殊 | 仍在保护区内 |
| $E_1 \sim |t_2|,|t_3|$ | **3–4（简并破缺）** | **一般** | **轨道跃迁** |
| $E_1 \gg |t_2|,|t_3|$ | 2（恢复） | 特殊 | RWA 对称恢复 |

---

## 六、统一图景

```
                     Sp(2) 纤维丛
                  Sp(2) → ℍℙ¹ (U(2) 纤维)
                 /                    \
         基空间 q                   纤维 u, v
     (MZM/ancilla 比值)        (内部相位 holonomy)
                 \                    /
              ω_u = Im(A + Bq)      ω_v = Im(D + Cp)
                       \              /
                    E₁=0 时 t₁ 方向平坦
                    E₁≠0 时联络弯曲
                           |
                    余伴随轨道分层
                    U(3τ) 的本征值
                    简并度 = 稳定子大小
                           |
              E₁ ≪ |t₂|,|t₃|: 特殊轨道（可规范变换）
              E₁ ∼ |t₂|,|t₃|: 一般轨道（需复合脉冲）
              E₁ ≫ |t₂|,|t₃|: 特殊轨道（编织已破坏）
```

### 核心结论

1. **$E_1=0$ 的保护是几何的**：纤维联络在 $t_1$ 方向平坦——不是分量巧合，是丛结构保证
2. **$E_1 \neq 0$ 的效应是拓扑的**：改变了 $U(3\tau)$ 的共轭类——不是"多了一个分量"，是轨道类型变了
3. **判据是整体的**：$U(3\tau)$ 的本征值简并度直接告诉你系统在哪种轨道上——不需要拆分任何分量

---

## 七、复合脉冲补偿理论

### 7.1 问题定义

给定 $E_1 \neq 0$ 的主协议 $\Gamma$（$\tau, t_1$ 固定），$U(\Gamma)$ 偏离理想共轭类。

**目标**：找补偿协议 $\bar\Gamma$（3 步，每步独立 $\bar\tau_n, \bar t_1^{(n)}$），使

$$\boxed{U(\bar\Gamma) \cdot U(\Gamma) \sim U_{\text{ideal}}}$$

其中 $\sim$ 表示共轭等价（本征值相等）。

### 7.2 自由度分析

$\bar\Gamma$ 参数：$2 \times 3 = 6$ 个。$Sp(2)$ 的共轭类由 2 个独立本征值确定。匹配条件：2 个方程。$6 > 2$，解空间为 4 维流形——问题良定。

### 7.3 存在性

$E_1 \neq 0$ 时 $A \neq D$，完整 $\mathfrak{sp}(2)$ 被激活。$\bar\Gamma$ 三个独立步的像在 $Sp(2)$ 中是 6 维子流形，在恒等元附近的切空间维数与 $Sp(2)$ 匹配，覆盖所需共轭类。补偿解存在。

### 7.4 与纤维联络的关系

纤维联络 $\omega_u = \operatorname{Im}(A+Bq)$ 解释**为什么** $E_1 \neq 0$ 时补偿可行（联络弯曲创造了额外的可达方向），但补偿条件本身不需要引用 $\omega_u$——只需匹配 $U_{\text{total}}$ 的本征值。

### 7.5 数值实现

固定 $E_1, \Gamma$，在 $\bar\Gamma$ 的 6 维参数空间中搜索使 $\operatorname{evals}(U(\bar\Gamma) \cdot U(\Gamma))$ 与 $\operatorname{evals}(U_{\text{ideal}})$ 距离最小化的点。这是标准优化问题。
