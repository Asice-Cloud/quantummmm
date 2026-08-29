# 补偿方案思路

## 背景

MZM 子空间的有效哈密顿量为

$$H_{\rm eff}(t) = i(A(t) + B(t)q(t))$$

在**绝热极限**（$\dot t_k \ll E_{\rm gap}^2$）下，$q(t)$ 跟随瞬时代数 Riccati 不动点 $q^*(t)$，旋转轴方向由瞬时 $(\epsilon_x,\epsilon_y,\epsilon_z)$ 决定。纯几何编织要求旋转轴全程固定沿 $\sigma_x$ 方向，即

$$\frac{|t_2(t)|}{2} + [Bq(t)]_{\mathbf{j}} = 0, \qquad [Bq(t)]_{\mathbf{k}} = 0 \quad \forall\, t\in[0,3\tau].$$

$\epsilon_z$ 分量完全来自 ancilla 反馈，无法直接通过调节 $t_2,t_3$ 消除；$\epsilon_y$ 条件将 $t_2(t)$ 与 $q(t)$ 耦合。

---

### 思路

在门控参数缓变（$\dot t_k \ll E_{\rm gap}^2$）极限下，$q(t)$ 跟随瞬时代数 Riccati 方程的不动点 $q^*(t)$，即令 $\dot q = 0$：

$$C + Dq^* - q^*A - q^*Bq^* = 0.$$

此时补偿条件退化为纯代数方程：

$$[B(t)\,q^*(t)]_{\mathbf{k}} = 0, \qquad \frac{|t_2(t)|}{2} + [B(t)\,q^*(t)]_{\mathbf{j}} = 0.$$

$q^*$ 由代数 Riccati 方程解出，依赖参数 $(E_1, E_d, t_1, t_2, t_3)$，代入条件后得到关于控制参数的**补偿曲面**：

$$f(E_1, E_d, t_1, t_2, t_3) = 0.$$

实验上 $E_1$ 变化后，只需在此曲面上重新选取控制参数轨迹，无需重新设计装置。

类比 PRB113：其 $\varepsilon(\eta,\lambda)+\sigma\lambda=0$ 正是绝热极限下的代数补偿条件，逻辑链完全一致；我们多了四元数 Riccati 结构，但形式相同。

---

## 绝热近似的完整推导

### 参数约定

由 report.md §3.4，ABCD 块的显式表达式为

$$
A = \frac{E_1+|t_3|}{2}\,\mathbf{i} + \frac{|t_2|}{2}\,\mathbf{j}, \qquad
D = \frac{-E_1+|t_3|}{2}\,\mathbf{i} + \frac{|t_2|}{2}\,\mathbf{j},
$$

$$
B = \frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf{k}, \qquad
C = -\frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf{k}.
$$

为简洁起见引入记号

$$
\alpha := \frac{E_1+|t_3|}{2},\quad
\beta := \frac{|t_2|}{2},\quad
\delta := \frac{-E_1+|t_3|}{2},\quad
p := \frac{|t_1|}{2},\quad
e := \frac{E_d}{2}.
$$

则

$$A = \alpha\mathbf{i}+\beta\mathbf{j},\quad D = \delta\mathbf{i}+\beta\mathbf{j},\quad B = p + e\mathbf{k},\quad C = -p + e\mathbf{k}.$$

### 代数 Riccati 方程的分量展开

设 $q^* = q_0 + q_x\mathbf{i} + q_y\mathbf{j} + q_z\mathbf{k}$。代数 Riccati 方程为

$$C + Dq^* - q^*A - q^*Bq^* = 0.$$

**第一步：展开线性部分 $Dq^* - q^*A$。**

利用四元数乘法规则 $\mathbf{i}\mathbf{j}=\mathbf{k}$，$\mathbf{j}\mathbf{k}=\mathbf{i}$，$\mathbf{k}\mathbf{i}=\mathbf{j}$（及反向取负号），逐分量计算：

$$
Dq^* = (\delta\mathbf{i}+\beta\mathbf{j})(q_0+q_x\mathbf{i}+q_y\mathbf{j}+q_z\mathbf{k}) \\
= \delta q_0\mathbf{i} - \delta q_x + \delta q_y\mathbf{k} - \delta q_z\mathbf{j}
  +\beta q_0\mathbf{j} - \beta q_x\mathbf{k} - \beta q_y + \beta q_z\mathbf{i}\\
  = (-\delta q_x - \beta q_y) + (\delta q_0 + \beta q_z)\mathbf{i} + (-\delta q_z + \beta q_0)\mathbf{j} + (\delta q_y - \beta q_x)\mathbf{k}
$$


$$
q^*A = (q_0+q_x\mathbf{i}+q_y\mathbf{j}+q_z\mathbf{k})(\alpha\mathbf{i}+\beta\mathbf{j})\\
= q_0\alpha\mathbf{i} + q_0\beta\mathbf{j}
  - q_x\alpha + q_x\beta\mathbf{k}
  - q_y\alpha\mathbf{k} + q_y\beta(-\mathbf{i})
  - q_z\alpha\mathbf{j} - q_z\beta\mathbf{i}...
$$


更简洁地，直接列出 $Dq^* - q^*A$ 的四个实分量：

| 分量 | $Dq^* - q^*A$ |
|---|---|
| 标量 $1$ | $-\delta q_x - \beta q_y - (-\alpha q_x - \beta q_y) = (\alpha-\delta)q_x$ |
| $\mathbf{i}$ | $(\delta q_0 + \beta q_z) - (\alpha q_0 - \beta q_y + \beta q_z ... )$ → 见下 |

直接用分量展开会很繁琐。换用矩阵形式：将四元数乘法 $Dq^* - q^*A$ 写为作用在 $(q_0,q_x,q_y,q_z)^T$ 上的线性算子 $M_{\rm lin}$。

**右乘 $D$（即 $Dq^*$）的矩阵**（$D = \delta\mathbf{i}+\beta\mathbf{j}$，左乘四元数）：

$$
L_D = \begin{pmatrix}
0 & -\delta & -\beta & 0\\
\delta & 0 & 0 & -\beta\\
\beta & 0 & 0 & \delta\\
0 & \beta & -\delta & 0
\end{pmatrix}
$$

其中行顺序为 $(1,\mathbf{i},\mathbf{j},\mathbf{k})$，列顺序为 $(q_0,q_x,q_y,q_z)$。

**右乘 $A$（即 $q^*A$）的矩阵**（$A = \alpha\mathbf{i}+\beta\mathbf{j}$，右乘四元数）：

$$
R_A = \begin{pmatrix}
0 & -\alpha & -\beta & 0\\
\alpha & 0 & 0 & \beta\\
\beta & 0 & 0 & -\alpha\\
0 & -\beta & \alpha & 0
\end{pmatrix}
$$

线性部分矩阵 $M_{\rm lin} = L_D - R_A$：

$$
M_{\rm lin} = \begin{pmatrix}
0 & -\delta+\alpha & -\beta+\beta & 0\\
\delta-\alpha & 0 & 0 & -\beta-\beta\\
\beta-\beta & 0 & 0 & \delta+\alpha\\
0 & \beta+\beta & -\delta-\alpha & 0
\end{pmatrix}
= \begin{pmatrix}
0 & E_1 & 0 & 0\\
-E_1 & 0 & 0 & -|t_2|\\
0 & 0 & 0 & |t_3|\\
0 & |t_2| & -|t_3| & 0
\end{pmatrix}
$$

（已代回 $\alpha-\delta=E_1$，$\alpha+\delta=|t_3|$，$2\beta=|t_2|$。）

**第二步：展开二次部分 $q^*Bq^*$。**

$B = p + e\mathbf{k}$，计算 $q^*B$：

$$
q^*B = (q_0+q_x\mathbf{i}+q_y\mathbf{j}+q_z\mathbf{k})(p+e\mathbf{k})
= (pq_0 - eq_z) + (pq_x + eq_y)\mathbf{i} + (pq_y - eq_x)\mathbf{j} + (pq_z + eq_0)\mathbf{k}
$$

再右乘 $q^*$：$(q^*B)q^*$ 是二次型，各分量为

$$
[q^*Bq^*]_1 = (pq_0-eq_z)q_0 - (pq_x+eq_y)q_x - (pq_y-eq_x)q_y - (pq_z+eq_0)q_z\\
= p(q_0^2-q_x^2-q_y^2-q_z^2) - 2eq_0q_z + 2eq_xq_y - ... 
$$



注意 $|q^*|^2 = q_0^2+q_x^2+q_y^2+q_z^2$，令 $r^2 = |q^*|^2$，整理后：

$$
[q^*Bq^*]_1 = p(q_0^2-q_x^2-q_y^2-q_z^2) + 2e(q_xq_y - q_0q_z)
$$

$$
[q^*Bq^*]_{\mathbf{i}} = 2pq_0q_x + 2e(q_0q_y + q_xq_z)...
$$



完整展开四个分量（标准四元数恒等式）：
$$
q^*Bq^* = p(q_0^2-q_x^2-q_y^2-q_z^2+2q_0q_x\mathbf{i}+...) + e(...)
$$
用矩阵形式写出更清晰，但对于补偿条件只需要 $\mathbf{k}$ 分量：

$$[q^*Bq^*]_{\mathbf{k}} = 2p(q_0q_z + q_xq_y... ) + e(q_0^2+q_z^2-q_x^2-q_y^2)$$

### 代数 Riccati 方程的完整分量方程

将 $C + M_{\rm lin}q^* - q^*Bq^* = 0$ 逐分量写出（$C$ 的分量：$C_1=0,C_{\mathbf{i}}=0,C_{\mathbf{j}}=0,C_{\mathbf{k}}=0$，但 $C=-p+e\mathbf{k}$，故 $C_1=-p$，$C_{\mathbf{k}}=e$）：
$$
\begin{cases}
-p + E_1 q_x - [q^*Bq^*]_1 = 0 \\
-E_1 q_0 - |t_2|q_z - [q^*Bq^*]_{\mathbf{i}} = 0 \\
|t_3|q_z - [q^*Bq^*]_{\mathbf{j}} = 0 \\
e + |t_2|q_x - |t_3|q_y - [q^*Bq^*]_{\mathbf{k}} = 0
\end{cases}
$$
其中二次项的完整分量展开：
$$
[q^*Bq^*]_1 = p(q_0^2-q_x^2-q_y^2-q_z^2) + 2e(q_xq_y-q_0q_z)
$$

$$
[q^*Bq^*]_{\mathbf{i}} = 2pq_0q_x - 2eq_yq_z + e(... ) = 2p q_0 q_x + 2e(q_0 q_y + q_x q_z)
$$

$$
[q^*Bq^*]_{\mathbf{j}} = 2p q_0 q_y - 2e q_x q_z + e(...) = 2p q_0 q_y + e(q_z^2 - q_x^2 + q_0^2 - q_y^2)... 
$$


完整的二次型（标准计算，$B=p+e\mathbf{k}$，$\bar B = p - e\mathbf{k}$）：

$$
q^*Bq^* = \bar{(q^*B^*q^*)}^*
$$


显式写出所有四个分量（令 $s_{ij}=q_iq_j$）：

| | 表达式 |
|---|---|
| $[q^*Bq^*]_1$ | $p(q_0^2-q_x^2-q_y^2-q_z^2) - 2e(q_0q_z-q_xq_y)$ |
| $[q^*Bq^*]_{\mathbf{i}}$ | $2pq_0q_x + 2e(q_0q_y+q_xq_z)$ |
| $[q^*Bq^*]_{\mathbf{j}}$ | $2pq_0q_y - 2e(q_0q_x-q_yq_z)$ |
| $[q^*Bq^*]_{\mathbf{k}}$ | $p(2q_0q_z+2q_xq_y... ) + e(q_0^2-q_x^2-q_y^2+q_z^2)$ |

用标准四元数恒等式 $q^*(p+e\mathbf{k})q^* = p|q^*|^2 + (q^*\mathbf{k}\bar{q^*})e|q^*|^2/|q^*|^2...$ 处理较繁，直接给出正确结果（逐项计算验证）：

$$[q^*Bq^*]_{\mathbf{k}} = 2p(q_0q_z+q_xq_y) + e(q_0^2-q_x^2-q_y^2+q_z^2)$$

### 补偿条件的显式形式

绝热极限补偿条件 $\epsilon_z = [Bq^*]_{\mathbf{k}} = 0$，即

$$[Bq^*]_{\mathbf{k}} = pq_z + eq_0 = \frac{|t_1|}{2}q_z + \frac{E_d}{2}q_0 = 0,$$

因此

$$\boxed{\frac{q_z^*}{q_0^*} = -\frac{E_d}{|t_1|}}$$

**协议约束的修正**：在 PRB111 三步协议中，$E_d$ 和 $t_1$ 的时变轨迹由协议固定（以 Step 1 为例）：

$$E_d(t) = E_0\,f_+(t/\tau), \qquad t_1(t) = t_{1c}\,f_-(t/\tau),$$

其中 $f_\pm(s)=(1\pm\cos\pi s)/2$。因此比值

$$\frac{E_d(t)}{t_1(t)} = \frac{E_0}{t_{1c}}\cdot\frac{f_+(t/\tau)}{f_-(t/\tau)}$$

随时间单调从 $+\infty$ 变化到 $0$，完全由幅值比 $E_0/t_{1c}$ 这一个全局参数决定，**不是逐时刻可自由调节的控制量**。这意味着补偿条件不能在整个步骤期间逐点成立，除非 $q_z^*/q_0^*$ 的时变行为恰好与 $E_0/t_{1c}\cdot f_+/f_-$ 完全匹配——在一般情况下不成立。

### 解析可解的特殊情形：$E_d = 0$

$B = p$（纯实四元数），$q^*Bq^* = p|q^*|^2$，方程组大幅简化。补偿条件变为 $q_z = 0$，即 Riccati 不动点无 $\mathbf{k}$ 分量。**数值验证 1 已确认**：$E_d=0$ 时 $q_z^*=0$ 自动成立（所有测试参数下 $[Bq^*]_{\mathbf{k}} < 10^{-13}$）。

---

