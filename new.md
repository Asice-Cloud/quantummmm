# 从 Sp(2) Riccati 到纯 MZM 七维封闭系统

## 目标

消去 ancilla 自由度（$\gamma_a, \gamma_b$），得到仅含三个 MZM（$\gamma_1, \gamma_2, \gamma_3$）的封闭动力学系统。

---

## 一、起点：Sp(2) Riccati 框架（回顾）

### 1.1 有效哈密顿量

$$H_{EM}(t) = iE_d\gamma_a\gamma_b + iE_1\gamma_1\gamma_2 + i|t_2|\gamma_a\gamma_2 - i|t_1|\gamma_b\gamma_1 - i|t_3|\gamma_a\gamma_3$$

### 1.2 Sp(2) 旋量表示

通过 $so(5) \cong sp(2)$ 同构，将 Majorana 算符替换为 $2\times2$ 四元数 Gamma 矩阵。演化矩阵 $U \in Sp(2)$ 满足：

$$\dot U = KU, \quad K = \begin{pmatrix} A & B \\ C & D \end{pmatrix} \in \mathfrak{sp}(2)$$

其中 $A,B,C,D \in \mathbb{H}$ 的显式（来自生成元 $\Sigma_{ij}$ 的投影）：

$$\boxed{\begin{aligned}
A(t) &= \frac{E_1+|t_3|}{2}\,\mathbf i + \frac{|t_2|}{2}\,\mathbf j \\[4pt]
D(t) &= \frac{-E_1+|t_3|}{2}\,\mathbf i + \frac{|t_2|}{2}\,\mathbf j \\[4pt]
B(t) &= \frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf k \\[4pt]
C(t) &= -\frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf k
\end{aligned}}$$

**关键性质**：
- $A, D \in \operatorname{Im}\mathbb{H}$（纯虚四元数，各 3 实分量，但仅 $\mathbf i, \mathbf j$ 分量非零）
- $B, C \in \operatorname{span}\{1, \mathbf k\}$（标量 + $\mathbf k$，各 2 实分量）
- $C = -\bar B$（即 $C_0 = -B_0$, $C_{\mathbf k} = B_{\mathbf k}$）
- $A = D \iff E_1 = 0$

### 1.3 Riccati 变量与 ODE

将 $U = \begin{pmatrix} X & Y \\ Z & W \end{pmatrix}$ 的第一列提取出来，定义：

$$q(t) := Z(t)X(t)^{-1} \in \mathbb{H}, \quad q(0) = 0$$

$q$ 满足 Riccati ODE：

$$\boxed{\dot q(t) = C(t) + D(t)q(t) - q(t)A(t) - q(t)B(t)q(t), \qquad q(0) = 0}$$

$q$ 有 4 个实分量，编码了 ancilla 对 MZM 的全部瞬时影响。

### 1.4 MZM 子空间演化

MZM 子空间的演化由 $X(t)$ 描述：

$$\dot X(t) = K_{\text{eff}}(t)\,X(t), \qquad K_{\text{eff}}(t) := A(t) + B(t)q(t)$$

$X(0) = 1$（初始化时 MZM 子空间未受扰动）。

$X \in \mathbb{H}$ 是 MZM → MZM 的振幅。归一化来自 $U^\dagger U = I$：

$$|X|^2 + |Z|^2 = 1 \;\Rightarrow\; |X|^2 = \frac{1}{1+|q|^2}$$

---

## 二、消去 ancilla：定义 $K := K_{\text{eff}} = A + Bq$

### 2.1 变量替换

定义：

$$\boxed{K(t) := A(t) + B(t)q(t) \in \mathbb{H}}$$

则 $q = B^{-1}(K - A)$（当 $B \neq 0$ 时）。注意到 $B(t)$ 仅在 $|t_1(t)| = E_d(t) = 0$ 时奇异——这在三步协议中最多发生在孤立的步边界瞬刻，可处理为极限。

$K$ 有 4 个实分量。将其分解为标量部分与矢量部分：

$$\boxed{K = k_0 + \Omega}, \quad k_0 := \operatorname{Re}(K), \quad \Omega := \operatorname{Im}(K) \in \operatorname{span}\{\mathbf i, \mathbf j, \mathbf k\} \cong \mathbb{R}^3$$

$\Omega$ 即 MZM 的有效磁场（旋转轴 × 角频率）。

### 2.2 $K$ 的自洽 ODE

对 $K = A + Bq$ 求导：

$$\dot K = \dot A + \dot B q + B\dot q$$

代入 Riccati ODE $\dot q = C + Dq - qA - qBq$，并用 $q = B^{-1}(K - A)$：

$$\begin{aligned}
\dot K &= \dot A + \dot B B^{-1}(K-A) \\
       &\quad + B\Big[C + D B^{-1}(K-A) - B^{-1}(K-A)A - B^{-1}(K-A)B B^{-1}(K-A)\Big]
\end{aligned}$$

利用 $B B^{-1} = 1$ 化简最后一项：$B^{-1}(K-A)B B^{-1}(K-A) = B^{-1}(K-A)^2$。

两边左乘 $B$ 后整理，得到 $K$ 的封闭 ODE：

$$\boxed{\dot K = \dot A + \dot B B^{-1}(K-A) + BC + BD B^{-1}(K-A) - (K-A)A - (K-A)^2}$$

**每一项都只含 $K$ 和已知的时间函数 $A,B,C,D,\dot A,\dot B$。ancilla 变量 $q$ 已被完全消去。**

初始条件：由 $q(0) = 0$ 得

$$\boxed{K(0) = A(0)}$$

（$A(0)$ 由 $t=0$ 时的门控参数给定。在三步协议中，Step 1 开始时 $t_2(0)=t_3(0)=E_1$（若 $E_1\neq 0$），故 $A(0)$ 可能非零。）

### 2.3 分拆为标量与矢量方程

令 $A = \mathbf a$, $D = \mathbf d$（纯矢量），$B = b_0 + b_3\mathbf k$, $C = -b_0 + b_3\mathbf k$。

四元数 $B$ 与 $B^{-1}$ 均在 $\operatorname{span}\{1,\mathbf k\}$ 内，记：

$$B^{-1} = \frac{b_0 - b_3\mathbf k}{b_0^2 + b_3^2} =: \beta_0 + \beta_3\mathbf k, \quad \beta_0 := \frac{b_0}{|B|^2},\; \beta_3 := -\frac{b_3}{|B|^2}$$

**乘积 $BD B^{-1}$**：$D$ 是纯矢量 $\mathbf d = d_1\mathbf i + d_2\mathbf j$。$B$ 对其的共轭作用是一个 $\mathbf i$–$\mathbf j$ 平面内的旋转：

$$BD B^{-1} = \tilde d_1 \mathbf i + \tilde d_2 \mathbf j$$

其中
$$\tilde d_1 = \frac{(b_0^2-b_3^2)d_1 - 2b_0b_3 d_2}{|B|^2}, \qquad \tilde d_2 = \frac{2b_0b_3 d_1 + (b_0^2-b_3^2)d_2}{|B|^2}$$

**乘积 $BC$**：$BC = (b_0 + b_3\mathbf k)(-b_0 + b_3\mathbf k) = -b_0^2 + b_3^2$（标量）。

**乘积 $\dot B B^{-1}$**：在 $\operatorname{span}\{1,\mathbf k\}$ 内，记作 $\gamma_0 + \gamma_3\mathbf k$。

$$\dot B B^{-1} = \frac{\dot b_0 b_0 + \dot b_3 b_3}{|B|^2} + \frac{\dot b_3 b_0 - \dot b_0 b_3}{|B|^2}\,\mathbf k$$

**项 $-(K-A)A$**：$K-A = k_0 + \Omega - \mathbf a$，$A = \mathbf a$（纯矢量）。

$$-(K-A)A = -(k_0 + \Omega - \mathbf a)\mathbf a = -k_0\mathbf a + \mathbf a\cdot(\Omega - \mathbf a) - \mathbf a \times (\Omega - \mathbf a)$$

（此处 $\mathbf a \times (\Omega - \mathbf a) = \mathbf a \times \Omega$，因为 $\mathbf a \times \mathbf a = 0$。）

标量部分：$\mathbf a \cdot (\Omega - \mathbf a)$，矢量部分：$-k_0\mathbf a - \mathbf a \times \Omega$。

**项 $-(K-A)^2$**：对于一般四元数，

$$(K-A)^2 = \big(k_0^2 - |\Omega - \mathbf a|^2\big) + 2k_0(\Omega - \mathbf a)$$

因此 $-(K-A)^2$ 的标量部分为 $-(k_0^2 - |\Omega - \mathbf a|^2)$，矢量部分为 $-2k_0(\Omega - \mathbf a)$。

将它们写成完整形式：分解 $\dot K = \dot k_0 + \dot\Omega$：

$$\boxed{\begin{aligned}
\dot k_0 &= \operatorname{Re}\!\big[\dot B B^{-1}(K-A)\big] + (b_3^2 - b_0^2) + \mathbf a\cdot\Omega - |\mathbf a|^2 - k_0^2 + |\Omega - \mathbf a|^2 \\[6pt]
\dot\Omega &= \dot{\mathbf a} + \operatorname{Im}\!\big[\dot B B^{-1}(K-A)\big] + \tilde{\mathbf d} - k_0\mathbf a - \mathbf a\times\Omega - 2k_0(\Omega - \mathbf a)
\end{aligned}}$$

其中 $\tilde{\mathbf d} = BD B^{-1}$ 的矢量部分，$\operatorname{Re}/\operatorname{Im}$ 分别取四元数的标量/矢量部分。

---

## 三、MZM 纯态的封闭演化

### 3.1 纯 MZM 四元数

定义单位四元数 $r(t)$ 为 $X(t)$ 的方向部分：

$$r(t) := \frac{X(t)}{|X(t)|} \in \mathbb{H}, \qquad |r(t)| = 1$$

$r \in SU(2)$，仅 3 个实自由度，完全编码了三个 MZM 的量子态。其与 SO(3) 编织矩阵 $R_{123}$ 的关系是标准旋量覆盖：

$$(R_{123})_{ij} = \frac{1}{2}\operatorname{Tr}(\sigma_i\, r\, \sigma_j\, r^\dagger)$$

其中 $\sigma_1, \sigma_2, \sigma_3$ 对应 $i\gamma_2\gamma_3, i\gamma_3\gamma_1, i\gamma_1\gamma_2$。

### 3.2 $r$ 的运动方程

由 $\dot X = KX$ 及 $X = \rho r$（$\rho = |X|$）：

$$\dot\rho\, r + \rho\, \dot r = K \rho r$$

两边右乘 $r^{-1}$，分别取实部与虚部：

- 实部：$\dot\rho / \rho = k_0$（标量部分仅缩放模长，已被 $|X|^2 + |Z|^2 = 1$ 自洽约束）
- 虚部：

$$\boxed{\dot r(t) = \Omega(t)\, r(t)}$$

**这是 Bloch 方程的四元数形式：$\dot r = \frac{1}{2}(\vec\omega \cdot \vec\sigma)\, r$，有效磁场为 $\Omega(t)$。**

### 3.3 完整七维封闭系统

$$\boxed{\begin{aligned}
\dot k_0 &= f_0(k_0, \Omega; t) \\[4pt]
\dot\Omega &= \mathbf{f}(k_0, \Omega; t) \\[4pt]
\dot r &= \Omega\, r
\end{aligned}}$$

| 层次 | 变量 | 维数 | 自治？ | 含义 |
|------|------|------|--------|------|
| 有效驱动层 | $(k_0, \Omega)$ | 4 | ✅ | ancilla 历史的痕迹，仅由门控参数驱动 |
| MZM 响应层 | $r$ | 3 | ❌（被 $\Omega$ 驱动） | 纯 MZM 旋转 |

**ancilla 从未出现在任何方程中。其唯一痕迹是初始条件 $K(0) = A(0)$，出现一次后消失。**

初始条件汇总：
$$K(0) = A(0), \quad r(0) = 1$$

---

## 四、退化极限

### 4.1 $E_1 = 0$：$A = D$

$$\dot\Omega = \dot{\mathbf a} + (\text{旋转项}) - \mathbf a \times \Omega - 2k_0(\Omega - \mathbf a)$$

此时 $A = D$ 使得系统具有额外对称性，$\Omega(t)$ 的轨迹约束在 $\mathbf i$–$\mathbf j$ 平面。$\Omega(t)$ 可以预先积分为门控参数的显式泛函（report §7 的解析解）。

### 4.2 $B = 0$（$t_1 = E_d = 0$）

$K = A$，$\dot r = \mathbf a(t)\, r$——平凡。

### 4.3 瞬时 $B = 0$（步边界）

$\dot B B^{-1}$ 会发散。此时应回归 Riccati 表述（$q$ 的 ODE）跨越该瞬刻。数值实现可检测 $|B|$ 并切换。

---

## 五、变量对比

| 方案 | 变量数 | ancilla 显式出现？ | 方程自治？ |
|------|--------|-------------------|------------|
| SO(5) 直接 | 25 | ✅ | ✅ |
| Sp(2) 直接 | 10 | ✅ | ✅ |
| Riccati ($q$) | 4 | ✅（$q$ 即 ancilla 比值） | ✅ |
| **$K$ 有效生成元** | **4** | **❌** | **✅** |
| **$(K, r)$ 联合** | **7** | **❌** | **✅** |

---

---

## 七、PRB113 路线：SO(5) 正交旋转尝试解耦

PRB113 (Nitsch et al., 2026) 通过**连续正交旋转 Majorana 基**，将对角化后的哈密顿量中解耦的 Majorana 对分离出来，得到低能有效扇区。本节将此方法应用于我们的 5-Majorana 模型。

### 7.1 SO(5) 反对称矩阵形式

将 $H_{EM}$ 写为 $H = \frac{i}{2} \sum_{ij} M_{ij} \gamma_i \gamma_j$，其中基底排序为 $\{\gamma_1, \gamma_2, \gamma_3, \gamma_a, \gamma_b\}$：

$$M(t) = \begin{pmatrix}
0 & E_1 & 0 & 0 & |t_1| \\
-E_1 & 0 & 0 & -|t_2| & 0 \\
0 & 0 & 0 & |t_3| & 0 \\
0 & |t_2| & -|t_3| & 0 & E_d \\
-|t_1| & 0 & 0 & -E_d & 0
\end{pmatrix}$$

（已逐项验证符号正确性。）

$M(t)$ 是 $5 \times 5$ 实反对称矩阵，总是可以通过一个 $SO(5)$ 旋转 $O(t)$ 对角化为：

$$O^T(t) M(t) O(t) = \begin{pmatrix}
0 & & & & \\
& 0 & \varepsilon_1 & & \\
& -\varepsilon_1 & 0 & & \\
& & & 0 & \varepsilon_2 \\
& & & -\varepsilon_2 & 0
\end{pmatrix}$$

其中 $\pm i\varepsilon_1, \pm i\varepsilon_2$ 是 $M$ 的特征值（一个特征值恒为零，对应瞬时零模）。

### 7.2 与 PRB113 的关键区别

| | PRB113 | 本模型 |
|---|---|---|
| Majorana 数 | 6（3 对） | 5（2 对 + 1 个） |
| 瑕疵参数 | $\zeta$（静态常数） | $E_1, t_1$（时变） |
| 编织耦合 | $\epsilon_1, t_{12}, t_{13}$（时变） | $t_2, t_3, E_d$（时变） |
| 修正项 | $\Lambda$（调谐至简并） | 无独立修正参数 |

**PRB113 能成功对角化并解耦的核心原因**：瑕疵 $\zeta$ 是常数 → 在协议第二步（$\theta = \pi/2$）期间，对角化角度 $\theta_\alpha, \theta_\mu$ 为常数 → 被积函数与自身对易 → 时间排序 $\mathcal{T}$ 可忽略。

**本模型的困难**：$E_1, t_1$ 随门控函数 $f_\pm(t)$ 时变 → 对角化旋转在所有时刻都在变化 → 没有"常数角度"区间。

### 7.3 第一步：$E_1 = t_1 = 0$ 极限下的 SO(4) 结构

此时 $M$ 退化为：

$$M_0(t) = \begin{pmatrix}
0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & -|t_2| & 0 \\
0 & 0 & 0 & |t_3| & 0 \\
0 & |t_2| & -|t_3| & 0 & E_d \\
0 & 0 & 0 & -E_d & 0
\end{pmatrix}$$

$\gamma_1$ **完全解耦**。剩余 $\{\gamma_2, \gamma_3, \gamma_a, \gamma_b\}$ 构成一个 SO(4) 子系统：

$$M_0^{(4)} = \begin{pmatrix}
0 & 0 & -|t_2| & 0 \\
0 & 0 & |t_3| & 0 \\
|t_2| & -|t_3| & 0 & E_d \\
0 & 0 & -E_d & 0
\end{pmatrix}$$

这恰好是纯净编织的极限——report 中已验证此时编织保真度 $\to 1$（$\tau \ge 50$），$\gamma_2 \leftrightarrow \gamma_3$ 交换精确成立。

### 7.4 第二步：PRB113 风格的方向提取

定义 ancilla $\gamma_a$ 所耦合的"有效方向"：

$$\boxed{\gamma_\Delta(t) := \frac{|t_2|\gamma_2 - |t_3|\gamma_3 + E_d\gamma_b}{N(t)},\quad N(t) := \sqrt{|t_2|^2 + |t_3|^2 + E_d^2}}$$

则 ancilla 相关的项简并为：

$$i|t_2|\gamma_a\gamma_2 - i|t_3|\gamma_a\gamma_3 + iE_d\gamma_a\gamma_b = iN(t)\,\gamma_a\gamma_\Delta$$

这正是 PRB113 的 $i\gamma_1\gamma_\Delta$ 项的对应物（我们的 $\gamma_a$ 对应他们的 $\gamma_1$，我们的 $\gamma_\Delta$ 对应他们的编织方向）。

在 $\{\gamma_2, \gamma_3, \gamma_b\}$ 的三维空间中，$\gamma_\Delta$ 是一个时变方向。补全正交基：

$$\boxed{\begin{aligned}
\gamma_{\theta'} &:= \frac{-|t_3|\gamma_2 - |t_2|\gamma_3}{\sqrt{|t_2|^2 + |t_3|^2}} \\[4pt]
\gamma_{\phi'} &:= \frac{E_d|t_2|\gamma_2 - E_d|t_3|\gamma_3 - (|t_2|^2+|t_3|^2)\gamma_b}{N(t)\sqrt{|t_2|^2 + |t_3|^2}}
\end{aligned}}$$

三者构成标准正交基：$\{\gamma_\Delta, \gamma_{\theta'}, \gamma_{\phi'}\}$（均为合法 Majorana 算符）。

在此基下，SO(4) 子块的哈密顿量为：

$$H_0 = iN(t)\,\gamma_a\gamma_\Delta$$

$\gamma_{\theta'}, \gamma_{\phi'}$ 与 $\gamma_a$ 无耦合——它们仅在基旋转变化时通过 Berry 联络产生效应。

### 7.5 第三步：引入 $E_1, t_1 \neq 0$

现在将 $\gamma_1$ 重新纳入。原始耦合项：

$$iE_1\gamma_1\gamma_2 - i|t_1|\gamma_b\gamma_1$$

用逆变换表达 $\gamma_2, \gamma_b$ 为 $\{\gamma_\Delta, \gamma_{\theta'}, \gamma_{\phi'}\}$：

$$\begin{aligned}
\gamma_2 &= \frac{|t_2|}{N}\gamma_\Delta - \frac{|t_3|}{\sqrt{t_2^2+t_3^2}}\gamma_{\theta'} + \frac{E_d|t_2|}{N\sqrt{t_2^2+t_3^2}}\gamma_{\phi'} \\[4pt]
\gamma_b &= \frac{E_d}{N}\gamma_\Delta - \frac{\sqrt{t_2^2+t_3^2}}{N}\gamma_{\phi'}
\end{aligned}$$

代入后，$\gamma_1$ 的耦合项变为：

$$iE_1\gamma_1\gamma_2 - i|t_1|\gamma_b\gamma_1 = i\gamma_1\!\left[\,\alpha_\Delta\gamma_\Delta + \alpha_{\theta'}\gamma_{\theta'} + \alpha_{\phi'}\gamma_{\phi'}\,\right]$$

其中：

$$\boxed{\begin{aligned}
\alpha_\Delta &= E_1\frac{|t_2|}{N} + |t_1|\frac{E_d}{N} \\[4pt]
\alpha_{\theta'} &= -E_1\frac{|t_3|}{\sqrt{t_2^2+t_3^2}} \\[4pt]
\alpha_{\phi'} &= E_1\frac{E_d|t_2|}{N\sqrt{t_2^2+t_3^2}} + |t_1|\frac{\sqrt{t_2^2+t_3^2}}{N}
\end{aligned}}$$

### 7.6 全哈密顿量在新基下

$$\boxed{H = iN(t)\,\gamma_a\gamma_\Delta + i\gamma_1\!\left[\alpha_\Delta\gamma_\Delta + \alpha_{\theta'}\gamma_{\theta'} + \alpha_{\phi'}\gamma_{\phi'}\right]}$$

**结构解读**：

- **第一项** $iN\gamma_a\gamma_\Delta$：$\gamma_a$ 仅与 $\gamma_\Delta$ 耦合——与 PRB113 结构一致
- **第二项**：$\gamma_1$ 同时与 $\gamma_\Delta, \gamma_{\theta'}, \gamma_{\phi'}$ 三个方向耦合

**这与 PRB113 的关键差异**：在 PRB113 中，$\gamma_1$ 仅与 $\gamma_\Delta$ 耦合（$i\gamma_1\gamma_\Delta$），因为他们的瑕疵 $\tilde{H}$ 涉及的 $\tilde{\gamma}_1$ 已经在 $\gamma_\Delta$ 方向中。我们的模型中，$\gamma_1$ 同时耦合到三个正交方向——这是因为 $E_1$ 和 $t_1$ 作用在 $\gamma_2$ 和 $\gamma_b$ 上，而这两者在旋转基中都有 $\gamma_{\theta'}, \gamma_{\phi'}$ 的分量。

### 7.7 能否进一步解耦？

现在哈密顿量涉及 5 个 Majorana：$\{\gamma_1, \gamma_a, \gamma_\Delta, \gamma_{\theta'}, \gamma_{\phi'}\}$。耦合结构为：

$$\begin{aligned}
\gamma_a &\leftrightarrow \gamma_\Delta \quad (\text{强度 } N) \\
\gamma_1 &\leftrightarrow \{\gamma_\Delta, \gamma_{\theta'}, \gamma_{\phi'}\} \quad (\text{强度 } \alpha_{\Delta}, \alpha_{\theta'}, \alpha_{\phi'})
\end{aligned}$$

这等价于一个 $4 \times 4$ 反对称子矩阵（$\gamma_1, \gamma_\Delta, \gamma_{\theta'}, \gamma_{\phi'}$）与 $\gamma_a$ 通过 $N$ 耦合：

$$M' = \begin{pmatrix}
0 & \alpha_\Delta & \alpha_{\theta'} & \alpha_{\phi'} & 0 \\
-\alpha_\Delta & 0 & 0 & 0 & N \\
-\alpha_{\theta'} & 0 & 0 & 0 & 0 \\
-\alpha_{\phi'} & 0 & 0 & 0 & 0 \\
0 & -N & 0 & 0 & 0
\end{pmatrix}$$

这里排序为 $\{\gamma_1, \gamma_\Delta, \gamma_{\theta'}, \gamma_{\phi'}, \gamma_a\}$。

**$\gamma_{\theta'}$ 和 $\gamma_{\phi'}$ 仅通过 $\gamma_1$ 与系统其余部分耦合**——它们不与 $\gamma_a$ 或 $\gamma_\Delta$ 直接耦合。

这暗示可以进一步旋转 $\{\gamma_{\theta'}, \gamma_{\phi'}\}$ 平面，使 $\gamma_1$ 只与其中一个耦合：

$$\gamma_\perp := \frac{\alpha_{\theta'}\gamma_{\theta'} + \alpha_{\phi'}\gamma_{\phi'}}{\sqrt{\alpha_{\theta'}^2 + \alpha_{\phi'}^2}},\qquad
\gamma_{\perp'} := \frac{\alpha_{\phi'}\gamma_{\theta'} - \alpha_{\theta'}\gamma_{\phi'}}{\sqrt{\alpha_{\theta'}^2 + \alpha_{\phi'}^2}}$$

则：

$$H = iN\gamma_a\gamma_\Delta + i\gamma_1\!\left[\alpha_\Delta\gamma_\Delta + \alpha_\perp\gamma_\perp\right]$$

其中 $\alpha_\perp := \sqrt{\alpha_{\theta'}^2 + \alpha_{\phi'}^2}$。**$\gamma_{\perp'}$ 完全解耦！**

最终有效哈密顿量仅涉及 4 个 Majorana：

$$\boxed{H_{\text{eff}} = iN(t)\,\gamma_a\gamma_\Delta + i\alpha_\Delta(t)\,\gamma_1\gamma_\Delta + i\alpha_\perp(t)\,\gamma_1\gamma_\perp}$$

### 7.8 结构的物理含义

$$\underbrace{iN\gamma_a\gamma_\Delta}_{\text{编织驱动}} + \underbrace{i\alpha_\Delta\gamma_1\gamma_\Delta}_{\text{E₁ 沿编织方向的泄漏}} + \underbrace{i\alpha_\perp\gamma_1\gamma_\perp}_{\text{E₁+t₁ 正交方向的泄漏}}$$

| 项 | 系数 | 物理来源 |
|---|---|---|
| $iN\gamma_a\gamma_\Delta$ | $N = \sqrt{t_2^2+t_3^2+E_d^2}$ | 编织门控 |
| $i\alpha_\Delta\gamma_1\gamma_\Delta$ | $\alpha_\Delta = \frac{E_1|t_2| + |t_1|E_d}{N}$ | 沿编织方向的瑕疵耦合 |
| $i\alpha_\perp\gamma_1\gamma_\perp$ | $\alpha_\perp = \sqrt{\alpha_{\theta'}^2+\alpha_{\phi'}^2}$ | 正交于编织方向的瑕疵耦合 |

**至此实现了从 5 Majorana → 4 Majorana 的降维**（$\gamma_{\perp'}$ 解耦）。剩余 4 个 Majorana $\{\gamma_1, \gamma_a, \gamma_\Delta, \gamma_\perp\}$ 构成一个 SO(4) 子系统——恰好是 PRB113 处理的那种结构。

### 7.9 与 PRB113 的对齐

PRB113 在 $\theta = \pi/2$（第二步）时，将剩余耦合写为 $i\gamma_1^D\gamma_\Delta^D$ 和 $i\gamma_\eta^D\gamma_{\theta'}^D$ 的形式。我们现在的 $\{\gamma_1, \gamma_a, \gamma_\Delta, \gamma_\perp\}$ 恰好对应他们的低能四 Majorana 子系统。可以进一步做类似于他们 Eq.(S17) 的旋转 ansatz：

$$\begin{aligned}
\gamma_1^D &= \cos\vartheta_1\,\gamma_1 + \sin\vartheta_1\,\gamma_a \\
\gamma_a^D &= \cos\vartheta_1\,\gamma_a - \sin\vartheta_1\,\gamma_1 \\
\gamma_\Delta^D &= \cos\vartheta_2\,\gamma_\Delta + \sin\vartheta_2\,\gamma_\perp \\
\gamma_\perp^D &= \cos\vartheta_2\,\gamma_\perp - \sin\vartheta_2\,\gamma_\Delta
\end{aligned}$$

选择合适的 $\vartheta_1(t), \vartheta_2(t)$ 可以使哈密顿量对角化。$\vartheta_1, \vartheta_2$ 的时变由 $N, \alpha_\Delta, \alpha_\perp$ 决定。

**关键**：在协议第二步（$t_2$ 关闭，$t_3$ 主导），$\alpha_\Delta \to E_1|t_3|/N$，$\alpha_\perp$ 有确定的时间依赖。如果 $\alpha_\Delta$ 和 $\alpha_\perp$ 的比值在第二步保持恒定，则 $\vartheta_1, \vartheta_2$ 为常数 → 可忽略 $\mathcal{T}$ → 得到类似 PRB113 的解析 Berry 相位。

### 7.10 结论：正交旋转路线的成果与局限

| 成果 | 说明 |
|------|------|
| ✅ 5→4 降维 | $\gamma_{\perp'}$ 完全解耦 |
| ✅ 结构对齐 | 剩余 4 Majorana 结构与 PRB113 一致 |
| ✅ 哈密顿量对角化 | 可通过角度 $\vartheta_1, \vartheta_2$ 实现 |
| ⚠️ 不同 | 我们的"瑕疵"全部时变——$\vartheta_1, \vartheta_2$ 未必为常数 |
| ⚠️ 不同 | 无独立调谐参数 $\Lambda$ 来恢复简并 |

**与 Riccati 路线的对比**：

| | Riccati（Sp(2)） | 正交旋转（SO(5)） |
|---|---|---|
| 降维方式 | $q = ZX^{-1}$（ancilla 比值） | 连续正交旋转 + 解耦 |
| 最终变量数 | 4（$q$） | 4（$\vartheta_1, \vartheta_2$ 参数化的 SO(4) 旋转） |
| ancilla 显式 | ✅（$q$ 是 ancilla 比值） | ❌（完全消去） |
| 代数结构 | $\mathfrak{sp}(2)$ 四元数 | $\mathfrak{so}(5)$ 实矩阵 |
| 优势 | 方程紧凑，数值高效 | MZM 演化在"裸"旋转框架中更直观 |
| 劣势 | $q$ 物理含义间接 | 旋转角 $\vartheta_{1,2}$ 的 ODE 比 Riccati 更复杂 |

**核心发现**：两条路线是完全对偶的。Riccati 用 4 维四元数 $q$ 编码 ancilla 记忆，正交旋转用 4 维角度 $(\vartheta_1, \vartheta_2)$ + Berry 相位编码同样的信息。区别在于：Riccati 的方程更简单（矩阵 Riccati 是规范形），但物理图像不直观；正交旋转的物理图像清晰（直接看到哪个 Majorana 对解耦），但角度 ODE 更复杂。

$E_1 = t_1 = 0$ 时两条路线都退化到平凡情况：$\gamma_1$ 完全自由，编织完美。$E_1 \neq 0$ 时两条路线都不可进一步简化——这是物理本身决定的，不是方法的局限。

---

### 7.11 绝热消除：$H_{\text{eff}} = h_0 I + \mathbf{h} \cdot \boldsymbol{\sigma}$

现在最关键的一步：从 4-Majorana 子系统绝热消除快变量，得到 MZM 的 $2\times2$ 有效哈密顿量。

#### 7.11.1 4-Majorana 的矩阵表示

取 $\{\gamma_1, \gamma_\perp, \gamma_\Delta, \gamma_a\}$ 的 Pauli 矩阵表示（4 维 Hilbert 空间 = 2 量子比特）：

$$\begin{aligned}
\gamma_1 &= \sigma_x \otimes I, & \gamma_\perp &= \sigma_y \otimes I \\
\gamma_\Delta &= \sigma_z \otimes \sigma_x, & \gamma_a &= \sigma_z \otimes \sigma_y
\end{aligned}$$

验证：全部 4 个两两反对易，$\gamma_i^2 = I$。✓

哈密顿量 $H = iN\gamma_a\gamma_\Delta + i\alpha_\Delta\gamma_1\gamma_\Delta + i\alpha_\perp\gamma_1\gamma_\perp$ 在此表示下：

$$\boxed{H = N\,(I \otimes \sigma_z) + \alpha_\Delta\,(\sigma_y \otimes \sigma_x) - \alpha_\perp\,(\sigma_z \otimes I)}$$

在基 $\{|00\rangle, |01\rangle, |10\rangle, |11\rangle\}$（第一量子比特 = MZM，第二 = ancilla 快模）下：

$$H = \begin{pmatrix}
N - \alpha_\perp & 0 & 0 & -i\alpha_\Delta \\
0 & -(N + \alpha_\perp) & -i\alpha_\Delta & 0 \\
0 & i\alpha_\Delta & N + \alpha_\perp & 0 \\
i\alpha_\Delta & 0 & 0 & -(N - \alpha_\perp)
\end{pmatrix}$$

#### 7.11.2 分块结构

$H$ 自然分裂为两个总宇称扇区：

- **偶宇称** $\{|00\rangle, |11\rangle\}$：$H_{\text{even}} = \begin{pmatrix} N-\alpha_\perp & -i\alpha_\Delta \\ i\alpha_\Delta & -(N-\alpha_\perp) \end{pmatrix}$
- **奇宇称** $\{|01\rangle, |10\rangle\}$：$H_{\text{odd}} = \begin{pmatrix} -(N+\alpha_\perp) & -i\alpha_\Delta \\ i\alpha_\Delta & N+\alpha_\perp \end{pmatrix}$

#### 7.11.3 绝热极限 $N \gg \alpha_\Delta, \alpha_\perp$

$N(t)$ 在编织过程中 $\sim 0.3$ meV（门控振幅），而 $E_1, t_1 \ll 0.3$（瑕疵参数）。因此 ancilla 模（$i\gamma_a\gamma_\Delta = I\otimes\sigma_z$）是快变量。

偶宇称扇区低能态 $|E_-\rangle$ 对应 $I\otimes\sigma_z = -1$（即第二个量子比特为 $|1\rangle$）：

$$H_{\text{even}}|E_-\rangle \approx -(N - \alpha_\perp) - \frac{\alpha_\Delta^2}{2(N-\alpha_\perp)}$$

奇宇称扇区低能态 $|O_-\rangle$ 同样对应 $I\otimes\sigma_z = -1$：

$$H_{\text{odd}}|O_-\rangle \approx -(N + \alpha_\perp) - \frac{\alpha_\Delta^2}{2(N+\alpha_\perp)}$$

两者之间的能隙：

$$\Delta E = E_{\text{even}} - E_{\text{odd}} \approx 2\alpha_\perp$$

#### 7.11.4 有效 $2\times2$ 哈密顿量

在以 $\{|E_-\rangle, |O_-\rangle\}$ 为基的 2 维低能流形中：

$$\boxed{H_{\text{eff}} = -\left(N + \frac{\alpha_\Delta^2}{2N}\right) I + \alpha_\perp\,\sigma_z + O\!\left(\frac{1}{N^2}\right)}$$

**这就是你想要的 $H_{\text{eff}} = h_0 I + \mathbf{h} \cdot \boldsymbol{\sigma}$ 形式！**

| 项 | 系数 | 物理来源 |
|---|---|---|
| $h_0 I$ | $-(N + \alpha_\Delta^2/2N)$ | 整体能量平移（含 $\alpha_\Delta$ 的二阶微扰修正） |
| $h_z \sigma_z$ | $\alpha_\perp$ | $E_1, t_1$ 在垂直于编织方向的投影 |
| $h_x \sigma_x$ | 0 | 瞬时哈密顿量不产生此分量 |
| $h_y \sigma_y$ | 0 | 瞬时哈密顿量不产生此分量 |

#### 7.11.5 编织从哪来？Berry 相位

瞬时 $H_{\text{eff}}$ 只有 $\sigma_z$ 分量——但这**不是**编织的结果。编织（$\gamma_2 \leftrightarrow \gamma_3$ 交换 = $\sigma_x$ 旋转）来自 **Berry 联络**：

$$A_\mu = \langle\psi_-|\partial_\mu|\psi_-\rangle$$

低能本征态 $|\psi_-(t)\rangle$ 自身随参数 $\{N, \alpha_\Delta, \alpha_\perp\}$ 旋转。这个旋转在参数空间的闭环积分给出非对易的 Berry 相位（Wilczek-Zee 联络）：

$$U_{\text{braid}} = \mathcal{P} \exp\!\left(-\oint A_\mu\, d\lambda^\mu\right)$$

展开后，$U_{\text{braid}}$ 一般包含 $\sigma_x, \sigma_y, \sigma_z$ 全部三个分量。

| | 瞬时 $H_{\text{eff}}$ | 完整演化 $U$ |
|---|---|---|
| $\sigma_z$ | ✅（来自 $\alpha_\perp$） | ✅ |
| $\sigma_x$ | ❌ | ✅（来自 Berry 相位） |
| $\sigma_y$ | ❌ | ✅（来自 Berry 相位） |

#### 7.11.6 $E_1 = 0$ 时的验证

$E_1 = 0, t_1 \neq 0$ 时：

$$\begin{aligned}
\alpha_\Delta &= \frac{|t_1|E_d}{N} \quad (\text{仅 } t_1 \text{ 贡献}) \\[4pt]
\alpha_\perp &= \frac{|t_1|\sqrt{t_2^2+t_3^2}}{N}
\end{aligned}$$

瞬时 $H_{\text{eff}} = h_0 I + \alpha_\perp \sigma_z$。$\alpha_\perp \propto |t_1|$ 在小 $t_1$ 下是小量。Berry 相位给出 $\sigma_x$ 分量（编织），$\sigma_z$ 分量来自 $t_1$ 的动态相位——与 report §7 的解析解 $R_{123} = \exp(-i\phi\,\hat n\cdot\vec\sigma)$ 完全一致，其中 $\hat n = (\cos\alpha, \sin\alpha, 0)$ 在 $xy$ 平面。

#### 7.11.7 $E_1 \neq 0$ 时

$$\alpha_\Delta = \frac{E_1|t_2| + |t_1|E_d}{N}, \quad \alpha_\perp = \sqrt{\alpha_{\theta'}^2 + \alpha_{\phi'}^2}$$

其中 $\alpha_{\theta'} = -E_1|t_3|/\sqrt{t_2^2+t_3^2},\;\alpha_{\phi'} = (E_1 E_d|t_2| + |t_1|(t_2^2+t_3^2))/(N\sqrt{t_2^2+t_3^2})$.

瞬时 $H_{\text{eff}}$ 仍是 $h_0 I + \alpha_\perp \sigma_z$，但 $\alpha_\perp$ 同时含 $E_1$ 和 $t_1$。两者在时变门控下不对易 → $\alpha_\perp(t)$ 的积分产生复杂干涉 → 这正是 report 中 Fig 1(d) 的物理。

---

### 7.12 与 Riccati 路线的最终对比

| | Riccati（Sp(2) 四元数） | 正交旋转（SO(5)） |
|---|---|---|
| 从 5 Majorana 到 | $q \in \mathbb{H}$（4 维） | $\{\gamma_1, \gamma_a, \gamma_\Delta, \gamma_\perp\}$（4 Majorana） |
| ancilla 的归宿 | 编码在 $q = ZX^{-1}$ 中 | $\gamma_a$ 被 $N$ 拉开，绝热消除 |
| MZM 有效理论 | $\dot r = \Omega\, r$（Bloch 球） | $H_{\text{eff}} = h_0 I + h_z \sigma_z$ + Berry 联络 |
| 数值复杂度 | 4 维 Riccati ODE（极简） | 4 维角度 ODE（较繁） |
| 物理透明度 | $q$ 含义间接 | ✅ 直观：能看到哪个 Majorana 对解耦 |
| **能否写 $H_{\text{eff}} = \mathbf{h}\cdot\boldsymbol{\sigma}$** | 间接（需从 $\Omega$ 反推） | ✅ **直接**：$h_z = \alpha_\perp$，$h_x, h_y$ 来自 Berry |

**答案**：是的，通过 SO(5) 正交旋转 + 绝热消除，我们明确得到了 $H_{\text{eff}} = h_0 I + \alpha_\perp \sigma_z$。瞬时哈密顿量只有 $\sigma_z$ 分量；编织的 $\sigma_x, \sigma_y$ 分量来自 Berry 相位——这恰好是 PRB113 的核心机制。

---

### 路线一：Sp(2) Riccati → K_eff → MZM

```
物理: H_EM (5 Majorana: γ₁,γ₂,γ₃,γ_a,γ_b)
  ↓ Sp(2) 表示
四元数: q = ZX⁻¹ (ancilla 记忆, 4 维)
  ↓ 变量替换 K = A + Bq
四元数: K = k₀ + Ω (有效驱动, 4 维, 无 ancilla)
  ↓
Ω 驱动 r: ṙ = Ω r (纯 MZM, 3 维)
  ↓
通往: 编织矩阵 R₁₂₃
```

### 路线二：SO(5) 正交旋转 → 解耦 → 有效 4-Majorana

```
物理: H_EM (5 Majorana)
  ↓ SO(5) 反对称矩阵 M(t)
旋转基: {γ_Δ, γ_θ', γ_ϕ'} (跟随编织方向)
  ↓ 进一步旋转解耦 γ_{⊥'}
4 Majorana: {γ₁, γ_a, γ_Δ, γ_⊥} (SO(4) 子系统)
  ↓ ansatz 对角化 (ϑ₁, ϑ₂)
Berry 相位 → 编织矩阵 R₁₂₃
```

### 两条路线的统一

| | Riccati | 正交旋转 |
|---|---|---|
| 数学结构 | $\mathfrak{sp}(2) \subset \mathbb{H}^{2\times2}$ | $\mathfrak{so}(5) \subset \mathbb{R}^{5\times5}$ |
| 降维机制 | 四元数比值 $q=ZX^{-1}$ | Gram-Schmidt 正交化 + 解耦 |
| 最终自由度 | $(k_0, \Omega, r)$ = 7 | $(\vartheta_1, \vartheta_2, r)$ = 7 |
| $E_1=0$ 极限 | $A=D$ → 解析解 | $\alpha_\perp=0$, $\gamma_1$ 半解耦 |
| $E_1\neq 0$ | 4 维 Riccati ODE（最简） | 4 维角度 ODE |

**共同结论**：

1. ancilla 可以完全消去，MZM 存在 7 维封闭动力学
2. $E_1=0$ 时可解析（两条路线均给出闭式解）
3. $E_1\neq 0$ 时需解 4 维 ODE（ancilla 记忆不可消除）
4. 4 维是物理本身要求的极小表示——Riccati 和正交旋转是对偶描述

---

## 八点五、Sp(2) 描述的进一步开发空间

SO(5) 正交旋转在「理解」上胜出，但 Sp(2) 在「深度」上有独特的开发潜力。

### 8.5.1 $K_{\text{eff}}$ 自身的 Riccati 结构

将 $\dot K$ 的方程（§2.2）与 Riccati ODE 对比：

$$\dot q = C + Dq - qA - qBq \qquad \text{(q 的 Riccati)}$$

$$\dot K = \dot A + \dot B B^{-1}(K-A) + BC + BD B^{-1}(K-A) - (K-A)A - (K-A)^2 \qquad \text{(K 的方程)}$$

整理 $\dot K$ 为 $K$ 的二次型：

$$\boxed{\dot K = \tilde C + \tilde D K - K \tilde A - K \tilde B K}$$

其中：

$$\begin{aligned}
\tilde C &= \dot A - \dot B B^{-1}A + BC - BD B^{-1}A + A^2 \\
\tilde D &= \dot B B^{-1} + BD B^{-1} + A \\
\tilde A &= \dot B B^{-1} + A \\
\tilde B &= 1 \quad \text{(!)}
\end{aligned}$$

**$K$ 也满足 Riccati ODE，且 $\tilde B = 1$ 意味着非线性项极简：$-K^2$。**

这是一个重大简化——$K$ 的 Riccati 方程中，非线性项没有任何参数依赖，就是朴素的 $-K^2$。对于四元数 $K = k_0 + \Omega$：

$$-K^2 = -(k_0^2 - |\Omega|^2) - 2k_0\Omega$$

**可开发方向**：$K$-Riccati 的非线性项与参数无关，这是 $q$-Riccati 不具备的性质。它意味着 $K$ 的动力学在所有参数区域有统一的非线性结构——特别适合做不动点分析、稳定性分析和相图。

### 8.5.2 几何解释：$q$ 作为齐性空间坐标

$q = ZX^{-1}$ 的定义在数学上有精确的几何含义。$U = \begin{pmatrix} X & Y \\ Z & W \end{pmatrix} \in Sp(2)$，而 $q$ 参数化了商空间：

$$q \in Sp(2) / (Sp(1) \times Sp(1)) \cong S^4$$

具体来说，$Sp(2)$ 对第一列 $(X, Z)^T$ 的作用通过右下角的 $Sp(1) \times Sp(1)$ 稳定子群产生商空间——这正是四元数射影空间 $\mathbb{H}P^1 \cong S^4$。

| 对象 | 几何含义 |
|------|---------|
| $Sp(2)$ | 全对称群（10 维） |
| $q \in \mathbb{H}$ | $S^4$ 上的仿射坐标（4 维） |
| $q = 0$ | ancilla 未占据的「原点」 |
| $|q| \to \infty$ | ancilla 完全占据（$X \to 0$） |
| $|X|^2 = 1/(1+|q|^2)$ | $S^4$ 上的径向函数 |

**可开发方向**：$S^4$ 上的 Riccati 流具有 $Sp(2)$ 等距群的对称性。这意味着存在守恒量（Casimir 不变量），可以降维。$S^4$ 上的测地线对应最简单的演化——$E_1 = t_1 = 0$ 时 $q$ 的轨迹是否沿测地线？

### 8.5.3 $E_1 = 0$ 对称性及其破缺的系统分类

$E_1 = 0 \iff A = D$。此时 Riccati ODE 对称群的李代数从一般的 $\mathfrak{sp}(2)$ 约化为一个子代数：

$$\dot q = C + [A, q] - qBq$$

线性部分退化为纯对易子 $[A,q]$——这是 $\mathfrak{sp}(1) \cong \mathfrak{su}(2)$ 的伴随作用，仅旋转 $q$ 的矢量部分。

**对称性层级**：

| 条件 | Riccati 形式 | 对称性 | 守恒量 |
|------|-------------|--------|--------|
| $E_1=t_1=0$ | $\dot q = [A,q] - qBq$（$B,C$ 纯 $\mathbf k$） | 最大 | $|q|$ 的演化受限 |
| $E_1=0, t_1\neq 0$ | $\dot q = C + [A,q] - qBq$ | 中等 | 轴锁定 $xy$ 平面 |
| **$E_1\neq 0, t_1\neq 0$** | $\dot q = C + Dq - qA - qBq$（$A\neq D$） | **最小** | 无 |

**可开发方向**：用对称性破缺参数 $\delta = A-D = E_1\mathbf i$ 做微扰展开。$E_1 = 0$ 的解析解作为零阶，$\delta$ 展开给出 $E_1 \neq 0$ 的近似解析公式——不依赖数值 ODE。

### 8.5.4 控制论视角：$\Omega(t)$ 作为 Bloch 球的操控函数

关键观察：MZM 的运动方程 $\dot r = \Omega(t) r$ 中，$\Omega(t)$ 就是控制函数。问题是：通过选择门控参数 $(E_1, t_1, \tau, \text{门控形状})$，$\Omega(t)$ 能覆盖 $\mathbb{R}^3$ 中的哪些曲线？

这正是控制论中的**可达集（reachable set）** 问题：

| 控制参数 | 可控的 $\Omega$ 分量 |
|----------|---------------------|
| $t_2(t)$ | $\Omega_x, \Omega_y$（编织驱动） |
| $t_3(t)$ | $\Omega_x, \Omega_y$（编织驱动） |
| $E_d(t)$ | $\Omega_z$（通过 $Bq$ 的 $\mathbf k$ 分量） |
| $E_1$ | $\Omega$ 全部分量（打破 $A=D$） |
| $t_1$ | $\Omega$ 全部分量（驱动 $q$ 的幅值） |

**可开发方向**：
1. 给定目标 Bloch 矢量 $\mathbf{n}_{\text{target}}$，是否存在门控协议驱动 $r(0) \to r_{\text{target}}$？
2. $E_1 = 0$ 时可达集是 $S^2$ 上的 1 维曲线（已确认）；$E_1 \neq 0$ 时是否覆盖整个球面？（report 数值表明是）
3. 最优控制：最小 $\tau$ 或最小 $t_1$ 下实现目标旋转

### 8.5.5 从 Sp(2) 回到 SO(5)：双重覆盖的显式利用

$\text{Sp}(2) \cong \text{Spin}(5)$ 是 $\text{SO}(5)$ 的双重覆盖。这意味着：
- SO(5) 中的正交旋转 $O(t)$ 对应 Sp(2) 中的两个矩阵 $\pm U(t)$
- SO(5) 中 $\gamma_{\perp'}$ 的解耦对应 Sp(2) 中某个子矩阵的块对角化

**可开发方向**：将 SO(5) 正交旋转的结果「翻译」回 Sp(2) 语言。特别地，$\gamma_{\perp'}$ 的解耦应对应 $q$ 的某个分量为零或守恒。这将把两条路线统一在一个框架内。

### 8.5.6 总结：Sp(2) 未开发的潜力

| 方向 | 难度 | 潜在收益 |
|------|------|---------|
| $K$ 的简化 Riccati（$\tilde B=1$） | ⭐ 低 | 统一非线性结构，便于不动点分析 |
| $S^4$ 几何与守恒量 | ⭐⭐ 中 | 可能发现新的运动积分，降维 |
| $E_1$ 微扰展开 | ⭐⭐ 中 | $E_1\neq 0$ 的半解析公式 |
| 控制论可达集 | ⭐⭐ 中 | 指导门控协议设计 |
| Sp(2)↔SO(5) 统一 | ⭐⭐⭐ 高 | 两条路线完美统一 |

**最值得优先推进的**：$E_1$ 微扰展开（实用性最高——直接给近似解析公式）和 $K$ 的简化 Riccati（结构最优雅——$\tilde B=1$ 是意外之喜）。

---

## 九、综合评判：Riccati vs SO(5) 正交旋转

### 9.1 信息量对比

| 维度 | Riccati | SO(5) 正交旋转 |
|------|---------|---------------|
| **ancilla 的归宿** | 编码在 $q$ 的 4 个分量中 | $\gamma_{\perp'}$ 解耦，$\gamma_a$ 被 $N$ 拉开 |
| **MZM 有效理论** | $\dot r = \Omega\, r$（Bloch 方程） | $H_{\text{eff}} = h_0 I + \alpha_\perp \sigma_z$ + Berry 联络 |
| **解耦的 Majorana** | 不可见 | ✅ $\gamma_{\perp'}$ 显式解耦 |
| **快/慢模分离** | 不可见 | ✅ $i\gamma_a\gamma_\Delta$（快，能隙 $2N$）vs $\{\gamma_1,\gamma_\perp\}$（慢） |
| **Berry 联络结构** | 隐含在 $\Omega(t)$ 的时变中 | ✅ 可直接从基矢旋转计算 $A_\mu$ |
| **参数的角色** | $A,B,C,D$ 四块混在一起 | ✅ $\alpha_\Delta$（沿编织方向）vs $\alpha_\perp$（正交方向）清晰分离 |

**结论**：SO(5) 正交旋转提供的信息远超 Riccati。Riccati 擅长「算」，SO(5) 擅长「理解」。

### 9.2 消除动力学相位、恢复纯几何驱动的条件

SO(5) 框架下，低能有效哈密顿量：

$$H_{\text{eff}} = h_0 I + \alpha_\perp \sigma_z$$

动力学相位来自 $\int \alpha_\perp(t)\,dt$。纯几何驱动要求此积分为零（或为 $2\pi$ 的整数倍）。

$\alpha_\perp$ 的显式结构：

$$\alpha_\perp = \sqrt{\alpha_{\theta'}^2 + \alpha_{\phi'}^2}$$

$$\alpha_{\theta'} = -\frac{E_1|t_3|}{\sqrt{t_2^2+t_3^2}}, \quad
\alpha_{\phi'} = \frac{E_1 E_d|t_2| + |t_1|(t_2^2+t_3^2)}{N\sqrt{t_2^2+t_3^2}}$$

由此直接读出消除动力学相位的条件：

| 条件 | 物理含义 |
|------|---------|
| $E_1 = 0$ **且** $t_1 = 0$ | 纯净 MZM 极限，$\alpha_\perp \equiv 0$——完美几何编织 |
| $\int_0^{3\tau} \alpha_\perp(t)\,dt = 2\pi n$ | 动力学相位恰好绕整数圈——编织保真度恢复 |
| $E_1 = 0$，$\int_0^{3\tau} \alpha_\perp(t)\,dt = 0$ | $t_1$ 的时间积分抵消——需要特殊的门控形状 |

与 PRB113 的对比：

| | PRB113 | 本模型 |
|---|---|---|
| 瑕疵参数 | $\zeta$ 静态常数 | $E_1, t_1$ 时变 |
| 恢复简并的手段 | 调谐 $\Lambda$（独立自由参数） | 调谐 $t_1, E_1$（但它们是「瑕疵」本身） |
| 操控自由度 | 瑕疵 + 修正 = 两个独立旋钮 | 所有参数都参与编织协议 |

**关键洞察**：PRB113 有独立修正参数 $\Lambda$ 来抵消瑕疵 $\zeta$ 的动力学效应——这是他们能"恢复简并"的核心。我们的模型没有这个自由度——$t_1$ 和 $E_1$ 本身就是编织门控的一部分，无法独立调谐。**SO(5) 框架清楚地揭示了这个结构性的不对称。**

### 9.3 Yang-Baxter 方程的适用性

编织算符 $B_{23}$（交换 $\gamma_2 \leftrightarrow \gamma_3$）应满足辫群关系：

$$B_{12} B_{23} B_{12} = B_{23} B_{12} B_{23}$$

在 SO(5) 框架中分析：

1. **纯净 MZM 极限（$E_1=t_1=0$）**：$\gamma_1$ 完全解耦，编织算符 $B_{23} \in SO(3)_{\gamma_2,\gamma_3}$ 是对 $\{\gamma_2,\gamma_3\}$ 的标准 $\pi/2$ 旋转。此时 Yang-Baxter 平凡成立（因为 $B_{12}$ 和 $B_{23}$ 作用在不同的 Majorana 对上，类似直积结构）。

2. **ABS 存在（$E_1\neq 0$ 或 $t_1\neq 0$）**：编织算符不再只是 $\gamma_2,\gamma_3$ 的旋转——$\gamma_1$ 通过 $\alpha_\perp$ 参与演化。$B_{23}$ 在低能 2 维流形上的表示 $U \in SU(2)$ 一般不是辫群的标准表示。

   SO(5) 框架可以直接检验：Berry 联络 $A_\mu$ 在参数空间中的路径积分给出 $U = \mathcal{P}\exp(-\oint A_\mu d\lambda^\mu)$。Yang-Baxter 要求 $U_{12} U_{23} U_{12} = U_{23} U_{12} U_{23}$，这等价于 $U$ 属于辫群的某个表示。

3. **SO(5) 框架的优势**：Berry 联络 $A_\mu$ 的曲率 $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu]$ 在低能流形上是一个 $2\times2$ 矩阵值 2-形式。Yang-Baxter 等价于曲率在参数空间某些面上的积分为零（平坦联络对应辫群表示）。SO(5) 框架可以直接计算这些量。

### 9.4 最终评判：两条路线的分工

```
         SO(5) 正交旋转                  Riccati (Sp(2))
              │                              │
     ┌────────┼────────┐              ┌──────┼──────┐
     │        │        │              │      │      │
   物理     解析     理论            数值    快速    验证
   洞察     条件     连接            计算    扫描    对照
     │        │        │              │      │      │
     ▼        ▼        ▼              ▼      ▼      ▼
  解耦的    消除动   Yang-          4变量   参数    论文
  Majorana  力学相   Baxter         ODE    空间    复现
  对        位的     方程                   热图
            条件
```

**推荐工作流**：

1. **用 SO(5) 做分析**：找出解耦的 Majorana 对，确定 $\alpha_\perp$ 的显式，分析消除动力学相位的条件，检查 Yang-Baxter 相容性
2. **用 Riccati 做计算**：4 维 ODE 高效数值积分，产出保真度热图、Bloch 球轨迹、参数扫描
3. **交叉验证**：两个框架在数值上必须一致（偏差 $<10^{-9}$）

**一句话总结**：SO(5) 正交旋转是「理论家的工具」——告诉你系统为什么长这样、哪些参数控制哪些效应、还能往哪个方向改进。Riccati 是「计算家的工具」——用最小的变量数得出最精确的数值结果。两者不竞争，互补。
