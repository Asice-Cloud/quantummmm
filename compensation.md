# 引入新生成元能否实现抵消条件

> 出发点,在不破坏四元数 Riccati 框架的前提下，引入新耦合能否使 $\epsilon_y = 0$ 或 $\epsilon_z = 0$？

---

## 1. Riccati 框架的保持条件

new3.md §2 的推导依赖 $K \in \mathfrak{sp}(2)$ 的 $2\times2$ 四元数分块结构：

$$
K = \begin{pmatrix}A & B \\ C & D\end{pmatrix},\quad A,D \in \operatorname{Im}\mathbb{H},\quad B,C \in \mathbb{H}.
$$

任何新耦合必须也是 $\mathfrak{sp}(2)$ 的元素，即对应某个 $\Sigma_{ij}$。否则 $H_{\rm eff} = i(A+Bq)$ 的推导整体失效。

---

## 2. SO(5) 生成元的分块分类

report.md §3.3 给出了全部 10 个生成元 $\Sigma_{ij}$ 的分块形式，按对角/非对角分类：

**对角块（进入 $A$、$D$）**——两指标均在 MZM 侧或均在 ancilla 侧：

| 生成元 | 对应耦合 | $A$ 的四元数分量 | $D$ 的四元数分量 | 当前模型 |
|---|---|---|---|---|
| $\Sigma_{12}$，$i\gamma_1\gamma_2$ | $E_1$ | $+\mathbf{i}/2$ | $-\mathbf{i}/2$ | ✓ |
| $\Sigma_{24}$，$i\gamma_a\gamma_2$ | $\|t_2\|$ | $+\mathbf{j}/2$ | $+\mathbf{j}/2$ | ✓ |
| $\Sigma_{34}$，$i\gamma_a\gamma_3$ | $\|t_3\|$ | $-\mathbf{i}/2$ | $-\mathbf{i}/2$ | ✓ |
| $\Sigma_{13}$，$i\gamma_1\gamma_3$ | $\Lambda_y$（新） | $+\mathbf{j}/2$ | $-\mathbf{j}/2$ | ✗ |
| $\Sigma_{23}$，$i\gamma_2\gamma_3$ | $\Lambda_z$（新） | $-\mathbf{k}/2$ | $+\mathbf{k}/2$ | ✗ |
| $\Sigma_{44}$，$i\gamma_a\gamma_a$（平凡） | — | — | — | — |

**非对角块（进入 $B$、$C$）**——跨越 MZM–ancilla 边界：

| 生成元 | 对应耦合 | 当前模型 |
|---|---|---|
| $\Sigma_{15}$，$i\gamma_b\gamma_1$ | $\|t_1\|$ | ✓ |
| $\Sigma_{45}$，$i\gamma_a\gamma_b$ | $E_d$ | ✓ |
| $\Sigma_{25}$，$i\gamma_b\gamma_2$ | — | ✗ |
| $\Sigma_{35}$，$i\gamma_b\gamma_3$ | — | ✗ |
| $\Sigma_{14}$，$i\gamma_a\gamma_1$ | — | ✗ |

引入非对角块的新生成元会改变 $B$、$C$ 的结构，从而改变 $H_{\rm eff}$ 的整个表达式，不在本文讨论范围。

---

## 3. 新生成元的显式分块

由 report.md §3.2 的 Gamma 矩阵直接计算：

**$\Sigma_{13} = \frac{1}{4}[\Gamma_1,\Gamma_3]$：**

$$
\Gamma_1\Gamma_3 = \begin{pmatrix}\mathbf{j}&0\\0&-\mathbf{j}\end{pmatrix},\quad
\Gamma_3\Gamma_1 = \begin{pmatrix}-\mathbf{j}&0\\0&\mathbf{j}\end{pmatrix}
\implies
\Sigma_{13} = \begin{pmatrix}\mathbf{j}/2&0\\0&-\mathbf{j}/2\end{pmatrix}
$$

**$\Sigma_{23} = \frac{1}{4}[\Gamma_2,\Gamma_3]$：**

$$
\Gamma_2\Gamma_3 = \begin{pmatrix}-\mathbf{k}&0\\0&\mathbf{k}\end{pmatrix},\quad
\Gamma_3\Gamma_2 = \begin{pmatrix}\mathbf{k}&0\\0&-\mathbf{k}\end{pmatrix}
\implies
\Sigma_{23} = \begin{pmatrix}-\mathbf{k}/2&0\\0&\mathbf{k}/2\end{pmatrix}
$$

两者均为纯对角块，$B$、$C$ 不变，Riccati 框架保持。

---

## 4. 引入新耦合后 $A$、$D$ 的修改

引入 $\Lambda_y(t)\,i\gamma_1\gamma_3$ 和 $\Lambda_z(t)\,i\gamma_2\gamma_3$，对应 $K$ 中加入
$\Lambda_y\,\Sigma_{13} + \Lambda_z\,\Sigma_{23}$，则：

$$
A \;\to\; A + \frac{\Lambda_y}{2}\,\mathbf{j} - \frac{\Lambda_z}{2}\,\mathbf{k}
$$

$$
D \;\to\; D - \frac{\Lambda_y}{2}\,\mathbf{j} + \frac{\Lambda_z}{2}\,\mathbf{k}
$$

$B$、$C$ 不变。

---

## 5. 对 $\epsilon_y$、$\epsilon_z$ 的影响

修改后，new3.md §4.2 的对应表变为：

$$
\epsilon_y = \underbrace{\frac{|t_2|}{2} + \frac{\Lambda_y}{2}}_{\text{几何项（来自 }A\text{）}} + [Bq]_\mathbf{j}
$$

$$
\epsilon_z = \underbrace{-\frac{\Lambda_z}{2}}_{\text{几何项（来自 }A\text{，原为零）}} + [Bq]_\mathbf{k}
$$

**关键变化：**

- $\Lambda_z$ 给 $\epsilon_z$ 增加了一个**直接可控的几何偏置** $-\Lambda_z/2$。这是当前模型里没有的自由度——$\epsilon_z$ 此前完全由 ancilla 反馈 $[Bq]_\mathbf{k}$ 决定，无法通过任何现有参数单独调控。

- $\Lambda_y$ 给 $\epsilon_y$ 的几何项增加了偏置 $\Lambda_y/2$，与 $|t_2|/2$ 的作用类似，但对 $D$ 的符号相反（$\Sigma_{13}$ 的 $D$ 块符号与 $\Sigma_{24}$ 不同），会改变 ancilla 侧的动力学，从而间接改变 $q(t)$。

---

## 6. 抵消条件的形式

在绝热极限下，$q \to q^*(t)$ 由修改后的不动点方程决定（$A$、$D$ 已包含 $\Lambda_y$、$\Lambda_z$）：

$$
0 = C + D'q^* - q^*A' - q^*Bq^*
$$

其中 $A' = A + \frac{\Lambda_y}{2}\mathbf{j} - \frac{\Lambda_z}{2}\mathbf{k}$，$D' = D - \frac{\Lambda_y}{2}\mathbf{j} + \frac{\Lambda_z}{2}\mathbf{k}$。

**$\epsilon_z = 0$ 的抵消条件：**

$$
\frac{\Lambda_z}{2} = [Bq^*]_\mathbf{k} = \frac{|t_1|}{2}q_z^* + \frac{E_d}{2}q_0^*
$$

其中 $q^*$ 依赖 $\Lambda_z$（通过修改后的不动点方程），这是一个关于 $\Lambda_z(t)$ 的**隐式自洽方程**。

**$\epsilon_y = 0$ 的抵消条件：**

$$
\frac{|t_2|}{2} + \frac{\Lambda_y}{2} = -[Bq^*]_\mathbf{j}
$$

同样是隐式方程，$q^*$ 依赖 $\Lambda_y$（通过 $A'$、$D'$）。

---

## 7. 可行性分析

两个抵消条件均为隐式自洽方程，能否有解取决于不动点方程的解对 $\Lambda_y$、$\Lambda_z$ 的依赖性。

**$\epsilon_z$ 方向（$\Lambda_z$）：**

目前 $A$、$D$ 均无 $\mathbf{k}$ 分量，$\Lambda_z$ 是引入 $\mathbf{k}$ 方向的**唯一来源**。不动点方程中，$\mathbf{k}$ 方向的耦合来自 $B$、$C$（含 $E_d\,\mathbf{k}$）和新的 $A'$、$D'$（含 $\Lambda_z\,\mathbf{k}$）。两者的 $\mathbf{k}$ 方向结构不同（$B$、$C$ 是非对角块，$\Lambda_z$ 在对角块），因此 $\Lambda_z$ 对 $q_z^*$、$q_0^*$ 的影响需要数值分析。

**$\epsilon_y$ 方向（$\Lambda_y$）：**

$\Sigma_{13}$ 的 $\mathbf{j}$ 方向与 $\Sigma_{24}$（即 $|t_2|$）在 $A$ 中相同，但在 $D$ 中符号相反。这意味着 $\Lambda_y$ 会破坏 $A$ 和 $D$ 的 $\mathbf{j}$ 分量之间的对称性，对 ancilla 动力学有额外影响，$q_y^*$、$q_x^*$ 的变化需要数值验证。

---

## 8. 结论

在保持四元数 Riccati 框架的前提下，**可以引入的 MZM 内部新耦合只有两个**：

- $\Lambda_z(t)\,i\gamma_2\gamma_3$：给 $\epsilon_z$ 提供直接可控的几何偏置，是当前模型中**缺失的唯一针对 $\epsilon_z$ 的自由度**
- $\Lambda_y(t)\,i\gamma_1\gamma_3$：给 $\epsilon_y$ 的几何项增加偏置，但与 $|t_2|$ 的作用部分重叠

这两个耦合对应器件中 $\gamma_2\gamma_3$ 和 $\gamma_1\gamma_3$ 的直接耦合，在实验上是否可调，取决于具体器件结构（见 image.png）。

抵消条件本身是隐式自洽方程，解的存在性和唯一性需要对修改后的不动点方程做数值分析。
