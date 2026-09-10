# Model

---

## 一、出发点：PRB111 有效模型

### 1.1 物理场景

PRB111 (Zhang et al., Phys. Rev. B **111**, 205411, 2025) 研究 ABS 存在时
Majorana 的 braiding 性质。系统包含 5 个 Majorana 模式：

- $\gamma_1,\gamma_2,\gamma_3$：MZM 模式（$\gamma_2$ 和 $\gamma_3$ 将被交换）
- $\gamma_a,\gamma_b$：量子点 ancilla 模式 

有效哈密顿量：

$$
H_{EM}(t) = iE_d\gamma_a\gamma_b + iE_1\gamma_1\gamma_2 + i|t_2|\gamma_a\gamma_2 - i|t_1|\gamma_b\gamma_1 - i|t_3|\gamma_a\gamma_3
$$

| 参数 | 含义 | 理想 MZM | ABS |
|---|---|---|---|
| $E_1$ | $\gamma_1$–$\gamma_2$ 杂化能 | 0 | $\neq 0$ |
| $t_1$ | $\gamma_1$–ancilla 耦合 | 0 | $\neq 0$ |
| $t_2,t_3$ | 门控编织耦合 | 时变 | 时变 |
| $E_d$ | 量子点能级 | 时变 | 时变 |

### 1.2 三段 braiding 协议

每段 $\tau$ 时长，门控函数 $f_\pm(t)=\frac{1\pm\cos(\pi t/\tau)}{2}$：

```
Step 1 (0→τ):   G1 关 → t₂ 打开, E_d → 0         γ₂ 移入量子点
Step 2 (τ→2τ):  G2 关, G1 开 → t₃ 打开, t₂ 关闭  γ₃ 接管, γ₂ 退出
Step 3 (2τ→3τ): G2 开 → t₃ 关闭, E_d 恢复        γ₃ 退出, 回到初态
                                    ─────────────
                                    结果: γ₂↔γ₃ 交换
```

---

## 二、为什么是 so(5)？

### 2.1 双线性生成元的李代数闭包

定义 $X_{ij}=i\gamma_i\gamma_j$（共 $C_2^5=10$ 个）。它们满足标准 $so(5)$ 对易关系：

$$[X_{ij}, X_{kl}] = 2i(\delta_{jk}X_{il} - \delta_{ik}X_{jl} - \delta_{jl}X_{ik} + \delta_{il}X_{jk})$$

### 2.2 论文哈密顿量在 so(5) 基底中的投影

$H_{EM}$ 恰好是 5 个生成元的线性组合：

$$H_{EM} = E_1 X_1 + |t_2| X_5 - |t_1| X_7 - |t_3| X_6 + E_d X_{10}$$

**三段协议本质上是在 10 维 so(5) 李代数里的分段时变演化，不预设二能级投影。**

### 2.3 

### 2.4 为什么研究李代数闭包？

演化算符 $U(t)=\mathcal T\exp(\int H(t)dt)$ 由哈密顿量 $H(t)$ 决定。
核心事实：

> **$H(t)$ 在李代数 $\mathfrak g$ 中 ⟹ $U(t)$ 被限制在李群 $G$ 中。**

李代数闭包 $\mathfrak g$ 的维数就是演化所需的最小参数个数。

| 段 | 闭包 $\mathfrak g$ | 李群 $G$ | 维数 | 解析含义 |
|---|---|---|---|---|
| Step 1 | $so(4)$ | $SO(4)$ | 6 | $SU(2)\times SU(2)$，两个 $su(2)$ 一般不对易 → 不可因子化 |
| Step 2 | $so(5)$ | $SO(5)$ | 10 | 维数满，时变轴 → commutant 平凡 → 无捷径 |
| Step 3 | $u(1)\oplus su(2)$ | $U(1)\times SU(2)$ | 4 | $u(1)$ 与 $su(2)$ 对易 → 可分解为 $e^{u(1)}\cdot e^{su(2)}$ → **可写闭式** |

**这就是为什么前面试过的 interaction picture、Magnus 展开都无法让 Step 2 降维——
不是因为方法不对，而是 Step 2 的代数本身就是满的 10 维 $so(5)$，里面没有更小的
不变子代数可用。** 李代数闭包分析让我们在动手算之前就知道：Step 3 可以简化，
Step 1/2 不可能。

---

## 三、so(5)表示

### 3.1 李代数同构与李群关系

李代数层面：$\mathfrak{so}(5) \cong \mathfrak{sp}(2)$（10 维实简单李代数同构）。

李群层面：$\text{Sp}(2) \cong \text{Spin}(5)$ 是 $\text{SO}(5)$ 的双重覆盖：

$$\text{SO}(5) \cong \text{Sp}(2)/\{\pm I\}$$

即 $U \in \text{Sp}(2)$ 和 $-U$ 映射到同一个 $\text{SO}(5)$ 矩阵。所有物理可观测量
（fidelity、Bloch 矢量、SO(5) 旋转矩阵）对 $\pm I$ 商自动不变，因此两套代码
（直接 SO(5) 和 Sp(2) 四元数）给出完全相同的结果（偏差 $<10^{-9}$）。

**$\mathfrak{sp}(2)$ 的定义**（四元数形式）：

$$\mathfrak{sp}(2) = \left\{ \begin{pmatrix} u & q \\ -\bar q & v \end{pmatrix}
\;\Big|\; u,v\in\operatorname{Im}\mathbb H,\; q\in\mathbb H \right\}$$

- $u,v$：纯虚四元数，各有 3 个实自由度 → $3+3=6$
- $q$：一般四元数，4 个实自由度
- 总维数：$6+4 = 10$，与 $so(5)$ 一致

等价地，$\mathfrak{sp}(2)$ 是所有 $2\times2$ 四元数反厄米矩阵：
$X^\dagger = -X$（$\dagger$ 为共轭转置）。

### 3.2 Cl(5) Gamma 矩阵

5 个 $2\times2$ 四元数矩阵：

$$\Gamma_1=\begin{pmatrix}0&1\\1&0\end{pmatrix},\;
\Gamma_2=\begin{pmatrix}0&-\mathbf i\\ \mathbf i&0\end{pmatrix},\;
\Gamma_3=\begin{pmatrix}0&-\mathbf j\\ \mathbf j&0\end{pmatrix},\;
\Gamma_4=\begin{pmatrix}0&-\mathbf k\\ \mathbf k&0\end{pmatrix},\;
\Gamma_5=\begin{pmatrix}1&0\\0&-1\end{pmatrix}$$

满足 $\{\Gamma_i,\Gamma_j\}=2\delta_{ij}$。

### 3.3 旋量生成元：$\Sigma_{ij} = \frac14[\Gamma_i, \Gamma_j]$

以 $\Sigma_{12}$ 为例显式计算：

$$\Gamma_1\Gamma_2 = \begin{pmatrix}0&1\\1&0\end{pmatrix}\begin{pmatrix}0&-\mathbf i\\\mathbf i&0\end{pmatrix}
= \begin{pmatrix}\mathbf i&0\\0&-\mathbf i\end{pmatrix}$$

$$\Gamma_2\Gamma_1 = \begin{pmatrix}0&-\mathbf i\\\mathbf i&0\end{pmatrix}\begin{pmatrix}0&1\\1&0\end{pmatrix}
= \begin{pmatrix}-\mathbf i&0\\0&\mathbf i\end{pmatrix}$$

$$[\Gamma_1,\Gamma_2] = \Gamma_1\Gamma_2 - \Gamma_2\Gamma_1 = \begin{pmatrix}2\mathbf i&0\\0&-2\mathbf i\end{pmatrix}$$

$$\Sigma_{12} = \frac14[\Gamma_1,\Gamma_2] = \begin{pmatrix}\mathbf i/2&0\\0&-\mathbf i/2\end{pmatrix}$$

同理可算得全部 5 个活跃生成元的分块形式：

| $\Sigma_{12}(E_1)$ | $\Sigma_{24}(|t_2|)$ | $\Sigma_{15}(-|t_1|)$ | $\Sigma_{34}(-|t_3|)$ | $\Sigma_{45}(E_d)$ |
|---|---|---|---|---|---|
| $\begin{pmatrix}\mathbf i/2&0\\0&-\mathbf i/2\end{pmatrix}$ | $\begin{pmatrix}\mathbf j/2&0\\0&\mathbf j/2\end{pmatrix}$ | $\begin{pmatrix}0&-1/2\\1/2&0\end{pmatrix}$ | $\begin{pmatrix}-\mathbf i/2&0\\0&-\mathbf i/2\end{pmatrix}$ | $\begin{pmatrix}0&\mathbf k/2\\\mathbf k/2&0\end{pmatrix}$ |

### 3.4 从 Majorana 哈密顿量到旋量演化 $\dot U = KU$

**第一步**：物理 Schrödinger 方程。
态矢量 $|\psi(t)\rangle = U(t)|\psi(0)\rangle$ 满足 $i\partial_t|\psi\rangle = H|\psi\rangle$，因此
$$\dot U(t) = -iH(t)\,U(t).$$

**第二步**：哈密顿量在旋量表示中的形式。
物理哈密顿量 $H_{EM} = \sum h_{ij}(t)\,(i\gamma_i\gamma_j)$。
在旋量表示中，将 Majorana 算符替换为 Gamma 矩阵 $\gamma_i \to \Gamma_i$：

$$H_{\text{spinor}}(t) = \sum h_{ij}(t)\,(i\Gamma_i\Gamma_j).$$

对 $i\neq j$，利用反对易关系 $[\Gamma_i,\Gamma_j] = \Gamma_i\Gamma_j - \Gamma_j\Gamma_i
= 2\Gamma_i\Gamma_j$，有 $\Gamma_i\Gamma_j = \frac12[\Gamma_i,\Gamma_j] = 2\Sigma_{ij}$。

代入得 $H_{\text{spinor}} = \sum h_{ij}\,(i\cdot 2\Sigma_{ij}) = 2i\sum h_{ij}\Sigma_{ij}$。

**第三步**：定义 $K(t) := \sum h_{ij}(t)\,\Sigma_{ij}$，则
$$\dot U = -i H_{\text{spinor}} U = -i(2iK)U = 2KU.$$

**归一化：**其中因子 2 可被吸收进生成元的归一化约定中。在我们的约定下，直接取
$\boxed{\dot U(t) = K(t)\,U(t)}$，$K\in\mathfrak{sp}(2)$。

**第四步**：将 $K$ 和 $U$ 按 $2\times2$ 四元数分块：

$$K = \begin{pmatrix}A&B\\ C&D\end{pmatrix},\qquad
U = \begin{pmatrix}X&Y\\ Z&W\end{pmatrix},\qquad
A,B,C,D,X,Y,Z,W\in\mathbb H.$$

$A,D$ 为纯虚四元数，对应于 $\mathfrak{sp}(2)$ 的对角块；
$B,C$ 为一般四元数，对应于非对角色块。

将 生成元分量代入 $K = \sum h_{ij}\Sigma_{ij}$，读出 $A,B,C,D$ 的显式：

$$\boxed{\begin{aligned}
A(t) &= \frac{E_1 + |t_3|}{2}\,\mathbf i + \frac{|t_2|}{2}\,\mathbf j, &
D(t) &= \frac{-E_1 + |t_3|}{2}\,\mathbf i + \frac{|t_2|}{2}\,\mathbf j,\\[4pt]
B(t) &= \frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf k, &
C(t) &= -\frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf k.
\end{aligned}}$$



---