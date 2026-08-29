# $\gamma_1,\gamma_2,\gamma_3$ 有效子空间的哈密顿量

> 出发点：report.md 的全部推导，不引入任何额外近似。
> 目标：把 ancilla $\gamma_a,\gamma_b$ 的效应精确地折叠进一个作用在
> MZM 子空间的有效哈密顿量，给出其显式结构，并指出动力学相位的代数来源。

---

## 1. 全空间演化方程

物理哈密顿量（report §1.1）：

$$
H_{EM}(t) = iE_d\,\gamma_a\gamma_b + iE_1\,\gamma_1\gamma_2
+ i|t_2|\,\gamma_a\gamma_2 - i|t_1|\,\gamma_b\gamma_1 - i|t_3|\,\gamma_a\gamma_3
$$

在旋量表示（report §3.4）下，演化算符 $U(t)\in\mathrm{Sp}(2)$ 满足 (其中系数2吸收到K里面)

$$
\dot U(t) = K(t)\,U(t), \qquad U(0)=\mathbb{1},
$$

其中 $K(t)\in\mathfrak{sp}(2)$ 是 $H_{EM}$ 在 $2\times2$ 四元数矩阵表示下的像：

$$
K(t) = \begin{pmatrix} A(t) & B(t) \\ C(t) & D(t) \end{pmatrix}, \qquad
A,B,C,D \in \mathbb{H},
$$

系数（report §3.4）：

$$
\boxed{
\begin{aligned}
A(t) &= \frac{E_1 + |t_3|}{2}\,\mathbf{i} + \frac{|t_2|}{2}\,\mathbf{j}, \\
D(t) &= \frac{-E_1 + |t_3|}{2}\,\mathbf{i} + \frac{|t_2|}{2}\,\mathbf{j}, \\
B(t) &= \frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf{k}, \\
C(t) &= -\frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf{k}.
\end{aligned}
}
$$

**分块的物理含义**（report §4.6.3）：

| 块 | 耦合的 Majorana 对 | 作用 |
|---|---|---|
| $A$ | $t_2\!\to\!\gamma_a\gamma_2$，$t_3\!\to\!\gamma_a\gamma_3$，$E_1\!\to\!\gamma_1\gamma_2$ | MZM 自身 + 编织门控 |
| $D$ | 同上但符号不同（ancilla 自旋） | ancilla 自身动力学 |
| $B$ | $t_1\!\to\!\gamma_b\gamma_1$，$E_d\!\to\!\gamma_a\gamma_b$ | ancilla → MZM 回馈 |
| $C$ | 同上反向 | MZM → ancilla 泄漏 |

旋量上分量 $\psi_1$ 对应 MZM 子空间，下分量 $\psi_2$ 对应 ancilla 子空间。
$U$ 的第一列满足自封闭方程（report §4.1）：
$$
\dot X = A X + B Z, \qquad \dot Z = C X + D Z.
$$

---

## 2. 证明：$i(A+Bq)$ 是 MZM 子空间的有效哈密顿量

### 2.1 旋量空间的物理分解

旋量表示将 $H_{EM}$ 映射为 $2\times2$ 四元数矩阵方程
$\dot U = K U$，$U\in\mathrm{Sp}(2)$，旋量空间为 $\mathbb{H}^2$（两分量四元数列向量）。

**上下分量的物理对应来自 Gamma 矩阵的明确分块结构**（report §3.2）：
$$
\Gamma_5 = \begin{pmatrix}1 & 0 \\ 0 & -1\end{pmatrix}.
$$

$\Gamma_5$ 的本征值 $+1$（上分量）和 $-1$（下分量）对应两个不同的 Majorana 奇偶扇区。
具体地，对角生成元 $\Sigma_{12}$ 和 $\Sigma_{34}$ 分别由 $E_1$（耦合 $\gamma_1\gamma_2$）和 $|t_3|$（耦合 $\gamma_a\gamma_3$）驱动，它们的矩阵块为：
$$
\Sigma_{12} = \begin{pmatrix}\mathbf{i}/2 & 0 \\ 0 & -\mathbf{i}/2\end{pmatrix}, \quad
\Sigma_{34} = \begin{pmatrix}-\mathbf{i}/2 & 0 \\ 0 & -\mathbf{i}/2\end{pmatrix},
$$

而非对角生成元 $\Sigma_{15}$（耦合 $\gamma_b\gamma_1$，即 $t_1$）和 $\Sigma_{45}$（耦合 $\gamma_a\gamma_b$，即 $E_d$）为：

$$
\Sigma_{15} = \begin{pmatrix}0 & -1/2 \\ 1/2 & 0\end{pmatrix}, \quad
\Sigma_{45} = \begin{pmatrix}0 & \mathbf{k}/2 \\ \mathbf{k}/2 & 0\end{pmatrix}.
$$

由此可读出：**对角块 $A,D$ 由纯 MZM 内部耦合（$t_2,t_3,E_1$）驱动，非对角块 $B,C$ 由跨越 MZM–ancilla 边界的耦合（$t_1,E_d$）驱动**。上分量即为初态在 MZM 扇区的分量。

### 2.2 初始条件与态的演化

物理上，初始态制备在 MZM 扇区：

$$
|\psi(0)\rangle = \begin{pmatrix}\psi_1(0) \\ 0\end{pmatrix}, \quad \psi_1(0) \in \mathbb{H}\setminus\{0\}.
$$

演化算符 $U(t)\in\mathrm{Sp}(2)$ 按列分块写为

$$
U(t) = \begin{pmatrix}X(t) & Y(t) \\ Z(t) & W(t)\end{pmatrix},
$$

则 $t$ 时刻的态为：

$$
|\psi(t)\rangle = U(t)|\psi(0)\rangle = \begin{pmatrix}X(t) \\ Z(t)\end{pmatrix}\psi_1(0).
$$

其中 $X(t)\psi_1(0)$ 是留在 MZM 扇区的分量，$Z(t)\psi_1(0)$ 是泄漏到 ancilla 扇区的分量。$Y,W$ 列对应 ancilla 初态，与本协议无关。

### 2.3 $X(t)$ 满足的方程

由 $\dot U = KU$ 展开第一列（$Y,W$ 列不出现）：

$$
\begin{pmatrix}\dot X \\ \dot Z\end{pmatrix}
= \begin{pmatrix}A & B \\ C & D\end{pmatrix}
\begin{pmatrix}X \\ Z\end{pmatrix}
\implies
\dot X = AX + BZ, \quad \dot Z = CX + DZ. \tag{1}
$$

这是精确方程，无任何截断或近似。方程 (1) 的两行是**耦合的**：$X$ 的演化依赖 $Z$，$Z$ 的演化依赖 $X$。

### 2.4 Riccati 变量消去 $Z$

**命题**：设 $X(t)$ 可逆，令 $q(t):=Z(t)X(t)^{-1}$，则 $q$ 满足
$$
\dot q = C + Dq - qA - qBq, \quad q(0)=0, \tag{2}
$$

且 $X(t)$ 满足

$$
\dot X = (A + Bq)\,X. \tag{3}
$$

**证明**：对 $q = ZX^{-1}$ 求导：
$$
\dot q = \dot Z X^{-1} + Z\frac{d}{dt}(X^{-1}) = \dot Z X^{-1} - Z X^{-1} \dot X X^{-1}.
$$

将 (1) 代入 $\dot Z = CX+DZ$ 和 $\dot X = AX+BZ$：

$$
\dot q = (CX+DZ)X^{-1} - ZX^{-1}(AX+BZ)X^{-1}
= C + DZX^{-1} - ZX^{-1}A - ZX^{-1}BZX^{-1}
= C + Dq - qA - qBq. \quad\square
$$

将 $Z = qX$ 代入 (1) 的上行：

$$
\dot X = AX + BqX = (A+Bq)X. \quad\square
$$

两个推导均为纯四元数代数恒等变换

**$X(t)$ 的可逆性**：$U\in\mathrm{Sp}(2)$ 满足 $U^\dagger U=\mathbb{1}$，展开左上角得 $X^\dagger X + Z^\dagger Z = \mathbb{1}$，故 $|X|^2 + |Z|^2 = 1$。初始 $X(0)=\mathbb{1}$，$|X(0)|=1$。在实际协议中 $q(t)$ 有界（数值验证 report §4.5，偏差 $<10^{-9}$），因此 $X(t)$ 全程可逆。

### 2.5 $i(A+Bq)$ 是有效哈密顿量

方程 (3) 完整描述了 MZM 扇区分量 $X(t)$ 的演化，且是精确的——ancilla 扇区 $Z(t)$ 的信息被**完整地**折叠进了 $q(t)$，无信息丢失。

$$
\boxed{H_{\rm eff}(t) = i\,K_{\rm eff}(t) = i\bigl(A(t) + B(t)\,q(t)\bigr).}
$$

**$H_{\rm eff}$ 是 Hermitian 的**：$K_{\rm eff} = A+Bq$ 作为单个四元数，其对应的 $2\times2$ 复矩阵是 $\mathfrak{su}(2)$ 元素（纯虚四元数对应 $\mathfrak{su}(2)$ 的反 Hermitian 生成元），故 $iK_{\rm eff}\in\mathfrak{su}(2)$ 是 Hermitian 矩阵。

---

## 3. MZM 有效哈密顿量的显式结构

$K_{\rm eff}(t)\in\mathbb{H}$（单个四元数），对应的**有效哈密顿量**（report §4.6.2）：

$$
H_{\rm eff}(t) = i\,K_{\rm eff}(t) = i\bigl(A(t) + B(t)\,q(t)\bigr).
$$

代入 $A,B$ 的显式：

$$
H_{\rm eff}(t) = i\left[
\underbrace{\frac{E_1+|t_3|}{2}\,\mathbf{i} + \frac{|t_2|}{2}\,\mathbf{j}}_{A(t)}
+ \underbrace{\left(\frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf{k}\right)}_{B(t)} q(t)
\right].
$$

$Bq$ 各分量展开（$q = q_0 + q_x\mathbf{i}+q_y\mathbf{j}+q_z\mathbf{k}$）：

$$
[Bq]_{\mathbf{i}} = \frac{|t_1|}{2}q_x + \frac{E_d}{2}q_y, \quad
[Bq]_{\mathbf{j}} = \frac{|t_1|}{2}q_y - \frac{E_d}{2}q_x, \quad
[Bq]_{\mathbf{k}} = \frac{|t_1|}{2}q_z + \frac{E_d}{2}q_0.
$$


---

## 4. 三 MZM 子空间的独立有效哈密顿量

### 4.1 子空间结构

三个 Majorana 零模 $\gamma_1,\gamma_2,\gamma_3$ 满足

$$
\{\gamma_i,\gamma_j\} = 2\delta_{ij}, \quad \gamma_i^\dagger = \gamma_i.
$$

它们的二次型算符 $i\gamma_i\gamma_j$（$i<j$）构成 $\mathfrak{so}(3)\cong\mathfrak{su}(2)$ 的三个生成元，与 Pauli 矩阵的对应为

$$
\sigma_x \leftrightarrow i\gamma_2\gamma_3, \quad
\sigma_y \leftrightarrow i\gamma_3\gamma_1, \quad
\sigma_z \leftrightarrow i\gamma_1\gamma_2.
$$

在费米子奇偶守恒的约束下，MZM 子空间的哈密顿量只能是这三个生成元的实线性组合：

$$
\boxed{H_3 = \epsilon_x\,i\gamma_2\gamma_3 + \epsilon_y\,i\gamma_3\gamma_1 + \epsilon_z\,i\gamma_1\gamma_2}
$$

其中 $\epsilon_x,\epsilon_y,\epsilon_z\in\mathbb{R}$ 

### 4.2 显式对应

将 $H_{\rm eff}(t)=i(A(t)+B(t)q(t))$ 展开为 Pauli 分量，对照 $H_3$ 的结构，三个耦合强度由下表给出：

| 分量 | $H_3$ 参数 | 我们模型中的显式表达式 | 来源 |
|---|---|---|---|
| $i\gamma_2\gamma_3$（$\sigma_x$） | $\epsilon_x$ | $\dfrac{E_1+\|t_3\|}{2} + [Bq]_{\mathbf{i}}$ | $A$ 的 $\mathbf{i}$ 分量 + ancilla 反馈 |
| $i\gamma_3\gamma_1$（$\sigma_y$） | $\epsilon_y$ | $\dfrac{\|t_2\|}{2} + [Bq]_{\mathbf{j}}$ | $A$ 的 $\mathbf{j}$ 分量 + ancilla 反馈 |
| $i\gamma_1\gamma_2$（$\sigma_z$） | $\epsilon_z$ | $[Bq]_{\mathbf{k}}$ | 纯 ancilla 反馈 |

其中 $Bq$ 的各分量展开为（$q=q_0+q_x\mathbf{i}+q_y\mathbf{j}+q_z\mathbf{k}$）：

$$
[Bq]_{\mathbf{i}} = \frac{|t_1|}{2}q_x + \frac{E_d}{2}q_y, \quad
[Bq]_{\mathbf{j}} = \frac{|t_1|}{2}q_y - \frac{E_d}{2}q_x, \quad
[Bq]_{\mathbf{k}} = \frac{|t_1|}{2}q_z + \frac{E_d}{2}q_0.
$$

