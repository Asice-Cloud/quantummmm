# 补偿方案：从 Riccati ODE 出发

> 仅使用 ODE 结构，不引用纤维丛理论。

---

## 一、ODE 体系

演化方程 $\dot U = KU$，$U \in Sp(2)$，$K \in \mathfrak{sp}(2)$。

$$K = \begin{pmatrix} A & B \\ C & D \end{pmatrix}, \qquad U = \begin{pmatrix} X & Y \\ Z & W \end{pmatrix}$$

$$A = \frac{E_1+|t_3|}{2}\mathbf{i} + \frac{|t_2|}{2}\mathbf{j},\quad D = \frac{-E_1+|t_3|}{2}\mathbf{i} + \frac{|t_2|}{2}\mathbf{j}$$
$$B = \frac{|t_1|}{2} + \frac{E_d}{2}\mathbf{k},\quad C = -\frac{|t_1|}{2} + \frac{E_d}{2}\mathbf{k}$$

$|t_2|,|t_3|,E_d$ 按余弦脉冲 $f_\pm(t) = (1 \pm \cos(\pi t/\tau))/2$ 时变。三步协议各持续 $\tau$。

Riccati 约化：$q = ZX^{-1}$，$\dot q = C + Dq - qA - qBq$。

---

## 二、$E_1=0$ 时的结构

**引理 1**：$E_1 = 0 \Rightarrow A = D$。Riccati 方程 $\dot q = C + [A,q] - qBq$ 给出

$$\frac{\dot\alpha}{|t_2|} = \frac{\dot\beta}{|t_3|}$$

从而 $|t_3|\beta - |t_2|\alpha \equiv 0$。$\alpha, \beta$ 是 $q$ 的 $\mathbf{i}, \mathbf{j}$ 分量。

**推论**：$\operatorname{Im}(A+Bq)$ 的 $\mathbf{k}$ 分量 $(Bq)_k = \frac{|t_1|}{2}q_3 + \frac{E_d}{2}q_0$ 在三步脉冲下的时间积分为零，与 $t_1$ 无关。

$$E_1 = 0 \;\Longrightarrow\; \frac{\partial U(3\tau)}{\partial t_1} \text{ 不改变共轭类}$$

---

## 三、$E_1 \neq 0$ 时的破缺

$A - D = E_1\mathbf{i} \neq 0$。引理 1 不成立。$\int_0^{3\tau} (Bq)_k dt \neq 0$，且依赖 $t_1$。

**引理 2**：$E_1 \neq 0$ 时，$U(3\tau)$ 对 $(\tau, t_1)$ 的偏导生成完整的 $\mathfrak{sp}(2)$。

---

## 四、补偿协议的 ODE 表述

主协议 $\Gamma$：$(\tau, t_1) = (\tau_0, t_1^0)$。$U(\Gamma) = \mathcal{T}\exp(\int_0^{3\tau_0} K dt)$。

目标 $U_{\text{ideal}}$：$E_1 = 0$ 时的演化矩阵。

补偿协议 $\bar\Gamma$：三步，参数 $(\bar\tau_1, \bar\tau_2, \bar\tau_3, \bar t_1^{(1)}, \bar t_1^{(2)}, \bar t_1^{(3)})$，共 6 个。

条件：

$$\boxed{U(\bar\Gamma) \cdot U(\Gamma) \sim U_{\text{ideal}}}$$

---

## 五、线性化补偿方程

**线性 ODE 变分公式**：对 $\dot U = K(t)U$，参数 $\lambda$ 的变分

$$\boxed{\frac{\partial U(T)}{\partial \lambda} = U(T) \int_0^T U(t)^{-1} \frac{\partial K(t)}{\partial \lambda} U(t)\,dt}$$

令 $U_0 = U(\Gamma)$。在恒等元附近 $U(\bar\Gamma) \approx I + \sum \delta_i M_i$，$M_i$ 由上式取参考协议的值计算。

乘积线性化：

$$U(\bar\Gamma) \cdot U_0 \approx U_0 + \sum \delta_i M_i U_0 = U_0\!\left(I + \sum \delta_i \, U_0^{-1} M_i U_0\right)$$

设 $\pi : \mathfrak{sp}(2) \to \mathfrak{h}$ 为到 Cartan 子代数（对角纯虚四元数，2 维）的投影。

$$\boxed{J \delta = b}$$

$$J_{ij} = \pi(U_0^{-1} M_i U_0)_j,\qquad b_j = \pi(U_0^{-1} U_{\text{ideal}} - I)_j$$

$J$ 为 $6 \times 2$ 矩阵。引理 2 保证 $\operatorname{rank}(J) = 2$。解 $\delta = J^+ b + \ker(J)$。

---

## 六、$M_i$ 的显式

$\partial K/\partial t_1$：$t_1$ 只出现在 $B, C$ 中：

$$\frac{\partial B}{\partial t_1} = \frac{1}{2} \cdot f(t), \qquad \frac{\partial C}{\partial t_1} = -\frac{1}{2} \cdot f(t)$$

$f(t)$ 是步内的脉冲包络（$f_r$ 或 $f_f$）。

$\partial K/\partial \bar\tau_n$：通过链式法则从 $f_\pm(t/\bar\tau_n)$ 对 $\bar\tau_n$ 求导。这是已知的初等函数。

代入变分积分即得 6 个 $M_i$。

---

## 七、算法

1. 传播 $U(\Gamma)$（一次 RK4）
2. 选参考补偿参数，传播参考 $U_{\text{ref}}(t)$（一次 RK4）
3. 计算 6 个 $M_i$（6 次数值积分，被积函数为已知函数）
4. 构造 $J, b$，解 $J\delta = b$ → 得 $\delta$（补偿参数修正量）

---

## 八、与纤维丛理论的关系

$J\delta = b$ 中的 $M_i$ 张成的方向恰为 $\mathfrak{g}_{\text{ctrl}}$（控制李代数）。$\pi$ 投影对应余伴随轨道的 Cartan 坐标。$\operatorname{rank}(J)=2$ 对应 $6 > 2$ 的自由度计数。两条路线同一结论。
