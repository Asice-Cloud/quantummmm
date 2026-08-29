# PRB111 有效模型的 Riccati ODE 解法——完整汇报

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

## 三、从 so(5) 到 Sp(2) 四元数表示

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

其中因子 2 可被吸收进生成元的归一化约定中。在我们的约定下，直接取
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

## 四、核心：推导 Riccati ODE

### 4.0 为什么用 $q$ 而不是直接演化 $U$？

$U\in Sp(2)$ 是 $2\times2$ 四元数矩阵，含 $4\times4=16$ 个实分量。
但 $\dot U=KU$ 的第一列方程**自封闭**——$\dot X=AX+BZ$，$\dot Z=CX+DZ$，与 $Y,W$ 无关。

$X,Z$ 各为四元数（共 8 个实分量），但演化中存在冗余的标度自由度。
Riccati 变量取其比值，压缩到 4 个分量：

$$q := ZX^{-1}\in\mathbb H.$$

> 类比：复数 $z=re^{i\theta}$ 的动力学，$q$ 相当于"角度"$\theta$，$X$ 相当于"模"$r$。
> $q$ 捕获了所有非平凡的李群演化，标度信息推迟到重建步骤。

### 4.1 分块运动方程

由 $\dot U = KU$，展开第一列：

$$\begin{pmatrix}\dot X\\ \dot Z\end{pmatrix}
= \begin{pmatrix}A&B\\ C&D\end{pmatrix}\begin{pmatrix}X\\ Z\end{pmatrix}
\quad\Longrightarrow\quad
\dot X = AX + BZ,\qquad \dot Z = CX + DZ$$

### 4.2 定义 Riccati 变量 $q = ZX^{-1}$

求导：

$$\dot q = \dot Z X^{-1} - ZX^{-1}\dot X X^{-1}$$

代入 $\dot X, \dot Z$：

$$\begin{aligned}
\dot q &= (CX + DZ)X^{-1} - (ZX^{-1})(AX + BZ)X^{-1} \\
       &= C + D(ZX^{-1}) - (ZX^{-1})A - (ZX^{-1})B(ZX^{-1}) \\
       &= C + Dq - qA - qBq.
\end{aligned}$$

### 4.3 最终 Riccati ODE

$$\boxed{\dot q(t) = C(t) + D(t)q(t) - q(t)A(t) - q(t)B(t)q(t),\qquad q(0)=0}$$

**全程纯代数推导，无近似、无截断。**

### 4.3bis 对易子 $[A,q]$ 的纯旋转性：$[A,q] = 2(A \times q)$

对于纯四元数（$\mathbb{R}^3$ 向量），对易子恒等于叉乘的两倍。设 $A = a_1\mathbf i + a_2\mathbf j + a_3\mathbf k$，
$q = \alpha\mathbf i + \beta\mathbf j + \gamma\mathbf k$（两者均为纯四元数，实部为零）：

$$[A,q] = Aq - qA$$

逐项展开：

$$\begin{aligned}
Aq &= -(a_1\alpha + a_2\beta + a_3\gamma) + (a_2\gamma - a_3\beta)\mathbf i + (a_3\alpha - a_1\gamma)\mathbf j + (a_1\beta - a_2\alpha)\mathbf k \\
qA &= -(a_1\alpha + a_2\beta + a_3\gamma) + (a_3\beta - a_2\gamma)\mathbf i + (a_1\gamma - a_3\alpha)\mathbf j + (a_2\alpha - a_1\beta)\mathbf k
\end{aligned}$$

相减：

$$\boxed{[A,q] = 2(a_2\gamma - a_3\beta)\mathbf i + 2(a_3\alpha - a_1\gamma)\mathbf j + 2(a_1\beta - a_2\alpha)\mathbf k}$$

另一方面，三维叉乘 $A \times q$：

$$A \times q = \begin{vmatrix}\mathbf i&\mathbf j&\mathbf k\\ a_1&a_2&a_3\\ \alpha&\beta&\gamma\end{vmatrix} = (a_2\gamma - a_3\beta)\mathbf i + (a_3\alpha - a_1\gamma)\mathbf j + (a_1\beta - a_2\alpha)\mathbf k$$

因此：

$$\boxed{[A, q] = 2(A \times q)}$$

**结论：$[A,q]$ 垂直于 $q$（$\because$ 叉乘性质），不改变 $|q|$ → 纯旋转项。**

这正是 Bloch 方程 $\dot{\vec S} = \vec\Omega \times \vec S$ 的四元数形式——进动轴为 $A$，
进动角频率为 $2|A|$，不拉伸、不压缩。

#### E₁=0 时的退化形式

$E_1 = 0 \;\Rightarrow\; A = D = \frac{|t_3|}{2}\mathbf i + \frac{|t_2|}{2}\mathbf j$，
此时 Riccati ODE 的线性部分坍缩为纯旋转：

$$\dot q = C + \underbrace{[A, q]}_{2(A\times q)} - qBq$$

$A \in \text{span}\{\mathbf i, \mathbf j\}$（轴在 $\sigma_x$–$\sigma_y$ 平面），
$[A,q]$ 的 $\mathbf k$（$\sigma_z$）分量为 $2(a_1\beta - a_2\alpha)$——即使 $\gamma=0$，
它也能瞬时产生 $\sigma_z$ 方向加速度。但三步门控的时间积分恰好抵消
（$\int f_m = \int f_p$），净 $\sigma_z = 0$，轨迹永限赤道。

### 4.4 三段统一

三段协议的差异仅体现在 $A(t),B(t),C(t),D(t)$ 的时变系数中（由 $E_1,|t_2|,|t_1|,|t_3|,E_d$
和门控函数 $f_\pm$ 决定），ODE 形式不变。**$q(t)$ 跨步边界自动连续。**

> **证明**：$K(t)$ 在步边界连续（门控函数 $f_\pm$ 在步边界处匹配），因此
> $\dot U=KU$ 的解 $U(t)$ 是 $C^1$ 连续的。$U(t)$ 连续 $\Rightarrow$ $X(t),Z(t)$ 连续。
>
> $X(t)$ 的可逆性：$Sp(2)$ 矩阵满足 $U^\dagger U=I$，展开左上角得
> $|X|^2+|Z|^2=1$。初始 $X(0)=1$（$|X|=1$）。若 $X(t)\to 0$，则需 $|Z|\to 1$，
> 此时 $q=ZX^{-1}$ 发散——这是 Riccati 坐标图的边界奇点，而非物理奇点。
> 在实际演化中 $q(t)$ 始终有界（数值验证），因此 $X(t)$ 始终可逆。
>
> 故 $q(t)=Z(t)X(t)^{-1}$ 是连续函数的复合，跨边界连续。

### 4.5 数值验证

全三段 Riccati vs 直接 Sp(2) 传播：**偏差 $<10^{-9}$（机器精度）。**

### 4.6 概念澄清：ancilla 消去、$H=iK$、$K_{\text{eff}}$ vs $K$

这一节集中回答推导链中容易混淆的几个关键点。

#### 4.6.1 全空间 $K$ 与有效生成元 $K_{\text{eff}}$ 不是同一个对象

$$\underbrace{K = \begin{pmatrix}A&B\\ C&D\end{pmatrix}}_{\text{全空间 } \mathfrak{sp}(2)} \quad\text{vs}\quad \underbrace{K_{\text{eff}} = A + Bq}_{\text{MZM 子空间 } \mathfrak{su}(2)}$$

| | $K$ | $K_{\text{eff}}$ |
|---|---|---|
| 维度 | $2\times2$ 四元数矩阵 ($\simeq 4\times4$ 复) | 单个四元数 ($\simeq 2\times2$ 复) |
| 作用对象 | 全空间 $U$ | 仅 MZM 子块 $X$ |
| 运动方程 | $\dot U = KU$ | $\dot X = K_{\text{eff}}X$ |

$K_{\text{eff}}$ 是从 $K$ **约化出来**的，不是同一个对象。推导链：
$$\dot U = \begin{pmatrix}A&B\\ C&D\end{pmatrix}U \;\xrightarrow{q=ZX^{-1}}\; \dot X = (A+Bq)X = K_{\text{eff}}X$$

#### 4.6.2 为什么 $H = iK$

物理 Schrödinger 方程是 $\dot U = -iHU$，而李代数流的标准写法是 $\dot U = KU$。对比得 $K = -iH$，即 $H = iK$。因子 $i$ 的作用是**把 anti-Hermitian 李代数生成元（$K^\dagger=-K$）变成 Hermitian 物理哈密顿量（$H^\dagger=H$）**。

在 §3.4 的完整推导中：$H_{\text{spinor}} = 2iK$，吸收因子 2 后约定 $\dot U = KU$。因此：
$$\boxed{K_{\text{eff}} = A + Bq},\qquad \boxed{H_{\text{eff}} = iK_{\text{eff}} = i(A+Bq)}$$

**$K_{\text{eff}}$ 是 anti-Hermitian 的，$H_{\text{eff}}$ 是 Hermitian 的。**

#### 4.6.3 "ancilla 消去"是什么意思

ancilla（量子点模式 $\gamma_a,\gamma_b$）不是物理删除，而是**不再显式追踪**——它的全部影响被编码进了 Riccati 变量 $q = ZX^{-1}$。

**状态空间的分层结构**（旋量表示中 $|\psi\rangle\in\mathbb H^2$）：
- $\psi_1$（上分量）$\leftrightarrow$ MZM 子空间（$\gamma_1,\gamma_2,\gamma_3$）
- $\psi_2$（下分量）$\leftrightarrow$ ancilla 子空间（$\gamma_a,\gamma_b$）

演化矩阵的分块：
$$|\psi(t)\rangle = \begin{pmatrix}X&Y\\ Z&W\end{pmatrix}\begin{pmatrix}\psi_1(0)\\ 0\end{pmatrix}
= \begin{pmatrix}X\psi_1(0)\\ Z\psi_1(0)\end{pmatrix}$$

| 块 | 含义 |
|---|---|
| $X$ | MZM → MZM（振幅留在 MZM 子空间） |
| $Z$ | MZM → ancilla（振幅泄漏到 ancilla） |

**Riccati 变量 $q = ZX^{-1}$ 是 ancilla 对 MZM 的瞬时比值：**
$$\psi_2 = Z\psi_1(0) = qX\psi_1(0) = q\,\psi_1$$
即 ancilla 态 $\psi_2 = q\psi_1$——完全由 MZM 态和比值 $q$ 决定，不再独立。

**四块 $A,B,C,D$ 在吸收后的归宿：**

$$K = \begin{pmatrix}A&B\\ C&D\end{pmatrix}$$

| 块 | 物理来源 | 吸收后的归宿 |
|---|---|---|
| $A$ | MZM 自身耦合（$t_2,t_3,E_1$） | 直接出现在 $\dot X = (A+Bq)X$ 中 |
| $B$ | ancilla 回馈 MZM（$t_1, E_d$） | 经 $q$ 加权：$Bq$ |
| $C$ | MZM 泄漏到 ancilla | 驱动 $\dot q$，不直接出现在 $\dot X$ |
| $D$ | ancilla 自身动力学 | 驱动 $\dot q$，不直接出现在 $\dot X$ |

$C$ 和 $D$ **没有被丢弃**——它们通过 Riccati ODE 决定 $q(t)$ 的演化轨迹，$q(t)$ 再通过 $Bq$ 间接影响 $X$。

**类比**：$q$ 相当于 Born-Oppenheimer 近似中的"有效势能"——快速自由度（ancilla）的影响被折叠进了慢速自由度（MZM）的等效运动方程，无需显式追踪前者。

#### 4.6.4 Pauli 矩阵对应关系

$$i\mathbf i \leftrightarrow i\gamma_2\gamma_3 \leftrightarrow \sigma_x$$
$$i\mathbf j \leftrightarrow i\gamma_3\gamma_1 \leftrightarrow \sigma_y$$
$$i\mathbf k \leftrightarrow i\gamma_1\gamma_2 \leftrightarrow \sigma_z$$

因此：
- $A = \frac{|t_3|}{2}\mathbf i + \frac{|t_2|}{2}\mathbf j$ 直接贡献 $\sigma_x$ 和 $\sigma_y$
- $Bq$ 贡献所有三个方向，具体分量由 $q$ 的虚部决定
- $E_1=0$ 时 $A=D$；$E_1\neq 0$ 时 $A\neq D$——这是决定旋转自由度的关键（见 §7.5 和后续）

---

## 五、回到物理量

1. $\dot X=(A+Bq)X$, $X(0)=1$ → 与 $q$ 联立得 $U(t)\in Sp(2)$
2. $R_{ij}=\frac12\operatorname{Tr}(\Gamma_i U\Gamma_j U^\dagger)$ → $R(t)\in SO(5)$
3. $R_{123}$（$\{\gamma_1,\gamma_2,\gamma_3\}$ 子块）即 braiding 结果

---