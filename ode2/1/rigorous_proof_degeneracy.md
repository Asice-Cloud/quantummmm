# PRB111 模型中动力学项抵消：显式条件与严格证明

> 对标 PRB113 (Nitsch et al., 2026) 的 $\varepsilon+\sigma\lambda=0$ 抵消方案

---

## 一、PRB113 的抵消逻辑（对照基准）

PRB113 核心方程（SI Eq. S31）：

$$\varepsilon(\eta, \lambda) + \sigma\lambda = 0$$

**物理含义**：通过调节 $\lambda(t)$，使低能两个本征态在**每一时刻**严格简并。

$$\lambda_1(t) = \lambda_2(t) = \lambda(t) \quad \forall t$$

因此动力学相位 $\exp(-i\int_0^T \lambda_i dt)$ 对两个态完全相同 → 相对动力学相位恒为 $1$。

---

## 二、PRB111 能否做到同类瞬时抵消？不能。

### 2.1 $K$ 的本征值显式

$K \in \mathfrak{sp}(2)$，在 $4\times4$ 复表示下本征值为 $\pm i\lambda_1, \pm i\lambda_2$。

$\lambda_{1,2}^2$ 是二次方程 $\mu^2 - a_2\mu + a_0 = 0$ 的根：

$$\boxed{\lambda_{1,2}^2 = \frac{a_2 \pm \sqrt{\Delta}}{2}}$$

其中：

$$\boxed{a_2 = \frac{E_1^2 + E_d^2 + t_1^2 + t_2^2 + t_3^2}{2}}$$

$$\boxed{\begin{aligned}
a_0 = &\ \frac{E_1^4}{16} - \frac{E_1^2 E_d^2}{8} + \frac{E_1^2 t_1^2}{8} + \frac{E_1^2 t_2^2}{8} - \frac{E_1^2 t_3^2}{8} + \frac{E_1 E_d t_2 |t_1|}{2} \\
      &+ \frac{E_d^4}{16} + \frac{E_d^2 t_1^2}{8} + \frac{E_d^2 t_2^2}{8} + \frac{E_d^2 t_3^2}{8} \\
      &+ \frac{t_1^4}{16} - \frac{t_1^2 t_2^2}{8} - \frac{t_1^2 t_3^2}{8} + \frac{t_2^4}{16} + \frac{t_2^2 t_3^2}{8} + \frac{t_3^4}{16}
\end{aligned}}$$

$$\boxed{\Delta = a_2^2 - 4a_0 = E_1^2 E_d^2 + E_1^2 t_3^2 - 2E_1 E_d t_2 |t_1| + t_1^2 t_2^2 + t_1^2 t_3^2}$$

### 2.2 瞬时简并条件

瞬时简并 $\Longleftrightarrow \lambda_1 = \lambda_2 \Longleftrightarrow \Delta = 0$：

$$\boxed{E_1^2 E_d^2 + E_1^2 t_3^2 - 2E_1 E_d t_2 |t_1| + t_1^2 t_2^2 + t_1^2 t_3^2 = 0}$$

当 $E_1 \neq 0$ 时，$E_1^2(E_d^2 + t_3^2)$ 为正。要使 $\Delta = 0$，必须有 $t_3 = 0$（否则 $E_1^2 t_3^2$ 无法被 $-2E_1 E_d t_2 |t_1|$ 完全补偿，因为 $t_1^2(t_2^2+t_3^2) \ge 0$）。

> **定理一**：$E_1 \neq 0 \;\land\; |t_3| \neq 0 \;\Longrightarrow\; \Delta > 0 \;\Longrightarrow\; \lambda_1 \neq \lambda_2$。
>
> 在编织 Step 2（$|t_3| > 0$）中，瞬时能级简并**原则上不可行**。

对比：PRB113 的 $\Delta_{\text{PRB113}}$ 含一个可独立调节的 $\lambda$，使其可不强制编织参数为零。

---

## 三、全局抵消的显式条件

既然不能逐点简并，只能要求**终点积累**正确。$U(T) \in Sp(2)$，$4\times4$ 复本征值为 $e^{\pm i\theta_1}, e^{\pm i\theta_2}$。

**抵消条件**：两个独立本征角等于理想值（模 $2\pi$）：

$$\boxed{\theta_1(T) \equiv \theta_1^{\text{ideal}} \pmod{2\pi}}$$

$$\boxed{\theta_2(T) \equiv \theta_2^{\text{ideal}} \pmod{2\pi}}$$

等价地用迹表示：

$$\boxed{\operatorname{Tr}[U(T)] = \operatorname{Tr}[U_{\text{ideal}}]}$$

$$\boxed{\operatorname{Tr}[U(T)^2] = \operatorname{Tr}[U_{\text{ideal}}^2]}$$

### 3.1 为什么这是「显式条件」

迹是**可直接计算**的量。在 Riccati 框架中：

$$U = \begin{pmatrix} X & Y \\ Z & W \end{pmatrix},\quad Z = qX$$

$$\operatorname{Tr}[U] = 2\,\big(\operatorname{Re}(X) + \operatorname{Re}(W)\big)$$

$$\operatorname{Tr}[U^2] = 2\,\big(\operatorname{Re}(X^2 + YZ) + \operatorname{Re}(ZY + W^2)\big)$$

$q(t)$ 由 Riccati ODE $\dot q = C + Dq - qA - qBq$ 决定。$X(t), W(t)$ 由 $\dot X = (A+Bq)X$ 决定。

**两个迹方程 = 控制参数必须满足的 2 个实方程。** 这与 PRB113 的 1 个方程 $\varepsilon+\sigma\lambda=0$ 完全类似——区别在于 PRB111 的方程涉及完整的时间积分（因为 $\mathfrak{so}(5)$ 不可约），而非瞬时值。

### 3.2 理想参考值

对于 $E_1=t_1=0, t_c=E_0=0.3, \tau=50$（绝热），数值计算得：

$$\theta_{1,2}^{\text{ideal}} = \{0.0304,\; 1.5404\} \;\text{rad}$$

对应 $\operatorname{Tr}[U_{\text{ideal}}] = 2(\cos 0.0304 + \cos 1.5404) \approx 2.0654$。

---

## 四、PRB113 式 $\varepsilon+\sigma\lambda=0$ 与 PRB111 式 $\theta_i = \theta_i^{\text{ideal}}$ 对比

| | PRB113 | PRB111 |
|---|---|---|
| **方程数** | 1 | 2 |
| **方程类型** | 瞬时代数方程 | 终态泛函方程 |
| **变量** | $\lambda(t)$ 可逐点解 | $t_2(t), t_3(t), E_d(t), t_1(t)$ 需全程联立 |
| **能写成闭式？** | ✅（二分图哈密顿量可解析对角化） | ❌（$\mathfrak{so}(5)$ 是单李代数，无不变子空间） |
| **数值求解** | 不需要 | Riccati ODE（4 分量，极快） |
| **物理含义** | 每时刻 $\Delta E = 0$ | 终点 $\int\Delta E\,dt \equiv 0 \pmod{2\pi}$ |
| **本质** | 无动力学相位积累 | 动力学相位恰为 $2\pi$ 的整数倍 |

---

## 五、绝热情形下的简化

大 $\tau$ 时，动力学相位近似为 $\int \lambda_i(t) dt$。抵消条件近似为：

$$\boxed{\int_0^{3\tau} \lambda_1(t)\,dt \approx \theta_1^{\text{ideal}} + 2\pi n_1}$$

$$\boxed{\int_0^{3\tau} \lambda_2(t)\,dt \approx \theta_2^{\text{ideal}} + 2\pi n_2}$$

其中 $\lambda_{1,2}(t)$ 是上面 §2.1 给出的**显式代数函数**。$n_1, n_2 \in \mathbb{Z}$。

**实操**：对三步协议，每步参数化后（放开 $t_c$、$E_0$、$t_1$、$\tau$ 的每步约束），有 $\ge 6$ 个自由参数。2 个方程 $\ll 6$ 个参数 → 解空间非空。

---

## 六、定理二：全局抵消的存在性（Chow 可控性）

**定理**：对于 $E_1 \neq 0$ 的 PRB111 系统，存在控制波形 $\{t_2(t), t_3(t), E_d(t), t_1(t)\}$ 使 $U(T)$ 与 $U_{\text{ideal}}$ 共轭等价。

**证明**：

1. 系统为 $\dot U = (K_{\text{drift}} + \sum u_i K_i)U$，$K_{\text{drift}} = \frac{E_1}{2}\Sigma_{12}$
2. $\{K_{\text{drift}}, \Sigma_{24}, \Sigma_{34}, \Sigma_{45}, \Sigma_{15}\}$ 的李括号闭包 $= \mathfrak{sp}(2)$（满 10 维，report §2.4 已验证）
3. $Sp(2)$ 是紧连通李群，LARC 成立 → 由 Chow-Rashevsky 定理，从恒等元出发的可达集 $=$ 整个 $Sp(2)$
4. $U_{\text{ideal}} \in Sp(2)$ → 存在控制波形和时间 $T$ 使 $U(T) = U_{\text{ideal}}$ ∎

**注**：存在性证明不提供显式波形，但保证了数值搜索（如 `compensation_clean.py`）不会搜索空集。

---

## 七、总结：PRB111 的抵消条件

$$\boxed{\begin{aligned}
&\text{条件 I:}\quad \operatorname{Tr}\big[\mathcal{T}e^{\int_0^T K(t)dt}\big] = \operatorname{Tr}[U_{\text{ideal}}] \\[8pt]
&\text{条件 II:}\quad \operatorname{Tr}\big[\big(\mathcal{T}e^{\int_0^T K(t)dt}\big)^2\big] = \operatorname{Tr}[U_{\text{ideal}}^2]
\end{aligned}}$$

其中 $K(t) = \begin{pmatrix}A&B\\ C&D\end{pmatrix}$，$A,B,C,D$ 是 $E_1, t_1(t), t_2(t), t_3(t), E_d(t)$ 的已知四元数函数（report §3.4）。

**与 PRB113 的本质差异**：PRB113 可将其条件化为 $\varepsilon+\sigma\lambda=0$（瞬时代数方程），因为其哈密顿量有二分图结构允许解析对角化。PRB111 不能——因为 $E_1|t_3| \neq 0$ 期间 $\mathfrak{so}(5)$ 被完全激活，没有不变子空间可资分块。条件必须是终态的，不是瞬时的。这不是推导技巧问题，是 $\mathfrak{so}(5)$ 为单李代数的直接推论。
