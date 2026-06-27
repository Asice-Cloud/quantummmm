# Riccati 框架下的动力学抵消条件

> 给出用 $q(T)$ 和 $X(T)$ 直接表达的抵消方程，可直接用于参数设计

---

## 一、Riccati 变量和重建

演化：
$$\dot q = C + Dq - qA - qBq,\quad q(0)=0$$
$$\dot X = (A+Bq)X,\quad X(0)=1$$
$$\dot W = (D+Cp)W,\quad p = -\bar q,\quad W(0)=1$$

其中 $A,B,C,D$ 是 $E_1, t_1(t), t_2(t), t_3(t), E_d(t)$ 的已知四元数函数（report §3.4）。

全矩阵重建：
$$Z = qX,\qquad Y = -\bar q W$$

$$U(T) = \begin{pmatrix} X(T) & Y(T) \\ Z(T) & W(T) \end{pmatrix}$$

---

## 二、本征角与 $X, W, q$ 的关系

$U(T) \in Sp(2)$，$4\times4$ 复本征值：$e^{\pm i\theta_1}, e^{\pm i\theta_2}$。

$$\boxed{\cos\theta_1 + \cos\theta_2 = X_0 + W_0}$$

$$\boxed{\cos 2\theta_1 + \cos 2\theta_2 = (X^2 - \bar q W q X)_0 + (W^2 - q X \bar q W)_0}$$

其中 $(\cdot)_0$ 表示四元数的实部，所有的量取 $t=T$ 时刻的值。

**推导**：
$$\operatorname{Tr}[U] = 2(X_0+W_0) = 2(\cos\theta_1+\cos\theta_2)$$
$$\operatorname{Tr}[U^2] = \operatorname{Tr}\begin{pmatrix} X^2+YZ & * \\ * & ZY+W^2 \end{pmatrix} = 2\big((X^2+YZ)_0 + (ZY+W^2)_0\big)$$
代入 $Z=qX$，$Y=-\bar qW$ 即得。$\operatorname{Tr}[U^2] = 2(\cos 2\theta_1 + \cos 2\theta_2)$。

---

## 三、抵消条件（设计方程）

设理想编织 $(E_1=0, t_1=0)$ 的终态本征角为 $\theta_1^{\text{ideal}}, \theta_2^{\text{ideal}}$。

定义目标常数：
$$S_1^{\text{ideal}} := \cos\theta_1^{\text{ideal}} + \cos\theta_2^{\text{ideal}}$$
$$S_2^{\text{ideal}} := \cos 2\theta_1^{\text{ideal}} + \cos 2\theta_2^{\text{ideal}}$$

**抵消条件**：

$$\boxed{X_0(T) + W_0(T) = S_1^{\text{ideal}}}$$

$$\boxed{(X^2)_0(T) + (W^2)_0(T) - \big(\bar q W q X\big)_0(T) - \big(q X \bar q W\big)_0(T) = S_2^{\text{ideal}}}$$

对 $E_1=t_1=0, t_c=E_0=0.3, \tau=50$（绝热），数值参考值：
$$S_1^{\text{ideal}} \approx 1.0327,\qquad S_2^{\text{ideal}} \approx 0$$

（注意：$S_2^{\text{ideal}}$ 在绝热极限下精确为 0，因为理想编织的两个本征角满足 $\theta_1 + \theta_2 = \pi/2$。）

---

## 四、参数设计方程

将控制参数向量记为 $\mathbf{p}$（每步独立的 $\tau, t_c, E_0, t_1$ 等，$\ge 6$ 维）。

定义残差函数：

$$\boxed{R_1(\mathbf{p}) := X_0(T;\mathbf{p}) + W_0(T;\mathbf{p}) - S_1^{\text{ideal}}}$$
$$\boxed{R_2(\mathbf{p}) := (X^2)_0(T;\mathbf{p}) + (W^2)_0(T;\mathbf{p}) - (\bar qWqX)_0(T;\mathbf{p}) - (qX\bar qW)_0(T;\mathbf{p}) - S_2^{\text{ideal}}}$$

其中 $X(T;\mathbf{p}), W(T;\mathbf{p}), q(T;\mathbf{p})$ 由 Riccati ODE 对参数 $\mathbf{p}$ 传播得到。

**设计任务**：求解

$$\boxed{\mathbf{R}(\mathbf{p}) = \mathbf{0}}$$

这是 2 个方程，$\ge 6$ 个未知数。解空间是 $\ge 4$ 维流形。可用 Newton 法或 Levenberg-Marquardt 求解。

---

## 五、梯度（Newton 法所需）

$\partial R_i / \partial p_j$ 可通过 Riccati ODE 的切线变分方程并行传播：

$$\frac{\partial}{\partial p_j}\begin{pmatrix} q(T) \\ X(T) \\ W(T) \end{pmatrix} = \text{变分方程的解}$$

联合传播状态 + 变分状态（$4 + 4 + 4$ 分量 $\times$ $n$ 个参数）即可得 Jacobian。用 Newton 法迭代：

$$\mathbf{p}^{(k+1)} = \mathbf{p}^{(k)} - J^+(\mathbf{p}^{(k)})\,\mathbf{R}(\mathbf{p}^{(k)})$$

其中 $J^+$ 是 $2 \times n$ Jacobian 的 Moore-Penrose 伪逆（选最小范数修正）。

---

## 六、与 PRB113 的对应

| | PRB113 | PRB111（本方案） |
|---|---|---|
| 设计方程 | $\varepsilon(\eta,\lambda)+\sigma\lambda=0$ | $\mathbf{R}(\mathbf{p})=\mathbf{0}$ |
| 方程数 | 1 | 2 |
| 未知数 | 1（$\lambda$，每时刻） | $\ge 6$（全程参数） |
| 逐点可解？ | ✅ | ❌（$\mathfrak{so}(5)$ 单性） |
| 如何解 | 直接代入闭式 | Newton 迭代 + Riccati 传播 |
| 每次迭代代价 | O(1) | O(N)（Riccati ODE，4+8 分量，极轻量） |

---

## 七、绝热近似下的简化

大 $\tau$ 时，动力学贡献主导，$\theta_i \approx \int \lambda_i(t) dt$。条件近似为：

$$\int_0^T \lambda_1(t;\mathbf{p})\,dt \approx \theta_1^{\text{ideal}} + 2\pi n_1$$
$$\int_0^T \lambda_2(t;\mathbf{p})\,dt \approx \theta_2^{\text{ideal}} + 2\pi n_2$$

其中 $\lambda_{1,2}(t)$ 是 §2.1 给出的**显式代数函数**。这两条方程不需求解 Riccati ODE——被积函数直接由控制参数表达。可用作 Newton 法的初始猜测。
