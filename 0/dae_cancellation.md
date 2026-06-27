# DAE 抵消方案：基于 Riccati ODE 的逐点路径设计

> 不分解群，不全局补偿，直接在 Riccati ODE 层面要求 $K_{\text{eff}} = K_{\text{eff}}^0$

---

## 一、核心想法

MZM 子空间演化完全由 $K_{\text{eff}} = A + Bq$ 决定：

$$\dot X = K_{\text{eff}} X,\quad X(0)=1$$

如果能在**每一时刻**使 $K_{\text{eff}}(t) = K_{\text{eff}}^0(t)$（$E_1=0, t_1=0$ 的理想值），
则 $X(t) = X^0(t)$，MZM 编织严格理想——无论 $E_1$ 多大。

这给出 4 个代数方程，联立 Riccati ODE 的 4 个 ODE，构成一个 **DAE 系统**。

---

## 二、$K_{\text{eff}}$ 的显式展开

$$A = \frac{E_1+|t_3|}{2}\mathbf i + \frac{|t_2|}{2}\mathbf j,\qquad B = \frac{|t_1|}{2} + \frac{E_d}{2}\mathbf k$$

$$q = q_0 + q_1\mathbf i + q_2\mathbf j + q_3\mathbf k$$

$$Bq = \left(\frac{|t_1|q_0}{2} - \frac{E_d q_3}{2}\right)
      + \left(\frac{|t_1|q_1}{2} - \frac{E_d q_2}{2}\right)\mathbf i
      + \left(\frac{|t_1|q_2}{2} + \frac{E_d q_1}{2}\right)\mathbf j
      + \left(\frac{|t_1|q_3}{2} + \frac{E_d q_0}{2}\right)\mathbf k$$

$$\boxed{K_{\text{eff}} = A + Bq}$$

$$\boxed{\begin{aligned}
S(\mathbf{u},q) &= \frac{|t_1|}{2}q_0 - \frac{E_d}{2}q_3 \\[4pt]
I(\mathbf{u},q) &= \frac{E_1+|t_3|}{2} + \frac{|t_1|}{2}q_1 - \frac{E_d}{2}q_2 \\[4pt]
J(\mathbf{u},q) &= \frac{|t_2|}{2} + \frac{|t_1|}{2}q_2 + \frac{E_d}{2}q_1 \\[4pt]
K(\mathbf{u},q) &= \frac{|t_1|}{2}q_3 + \frac{E_d}{2}q_0
\end{aligned}}$$

$\mathbf{u} = (|t_1|, |t_2|, |t_3|, E_d)^T$。

---

## 三、理想参考 $K_{\text{eff}}^0$

在 $E_1=0, t_1=0$ 的标准协议下，预计算 $q^0(t)$（Riccati 轨迹）和 $K_{\text{eff}}^0(t)$：

$$K_{\text{eff}}^0(t) = \underbrace{-\frac{E_d^0(t) q_3^0(t)}{2}}_{S^0}
                     + \underbrace{\left(\frac{|t_3^0(t)|}{2} - \frac{E_d^0(t) q_2^0(t)}{2}\right)}_{I^0}\mathbf i
                     + \underbrace{\left(\frac{|t_2^0(t)|}{2} + \frac{E_d^0(t) q_1^0(t)}{2}\right)}_{J^0}\mathbf j
                     + \underbrace{\frac{E_d^0(t) q_0^0(t)}{2}}_{K^0}\mathbf k$$

$t_2^0(t), t_3^0(t), E_d^0(t)$ 按 $f_\pm(t)$ 余弦脉冲时变。$q^0(t)$ 由 $\dot q^0 = C^0 + D^0 q^0 - q^0 A^0 - q^0 B^0 q^0$ 决定。

---

## 四、匹配方程 $K_{\text{eff}} = K_{\text{eff}}^0$

$$\boxed{\begin{aligned}
\frac{|t_1|}{2}q_0 - \frac{E_d}{2}q_3 &= -\frac{E_d^0}{2}q_3^0 \qquad\text{(S)}\\[4pt]
\frac{E_1+|t_3|}{2} + \frac{|t_1|}{2}q_1 - \frac{E_d}{2}q_2 &= \frac{|t_3^0|}{2} - \frac{E_d^0}{2}q_2^0 \qquad\text{(I)}\\[4pt]
\frac{|t_2|}{2} + \frac{|t_1|}{2}q_2 + \frac{E_d}{2}q_1 &= \frac{|t_2^0|}{2} + \frac{E_d^0}{2}q_1^0 \qquad\text{(J)}\\[4pt]
\frac{|t_1|}{2}q_3 + \frac{E_d}{2}q_0 &= \frac{E_d^0}{2}q_0^0 \qquad\text{(K)}
\end{aligned}}$$

---

## 五、控制函数的显式解

方程 (S) 和 (K) 是 $|t_1|, E_d$ 的 $2\times2$ 线性系统：

$$\begin{pmatrix} q_0 & -q_3 \\ q_3 & q_0 \end{pmatrix}
\begin{pmatrix} |t_1|/2 \\ E_d/2 \end{pmatrix}
= \frac{E_d^0}{2} \begin{pmatrix} -q_3^0 \\ q_0^0 \end{pmatrix}$$

行列式 $\det = q_0^2 + q_3^2$。对 $|q| > 0$ 可逆，解得：

$$\boxed{\begin{pmatrix} |t_1| \\ E_d \end{pmatrix}
= \frac{E_d^0}{q_0^2+q_3^2}
\begin{pmatrix} -q_0 q_3^0 + q_3 q_0^0 \\[4pt] q_0 q_0^0 + q_3 q_3^0 \end{pmatrix}}$$

代回 (I) 和 (J) 得 $|t_2|, |t_3|$：

$$\boxed{\begin{aligned}
|t_3| &= |t_3^0| - E_1 - |t_1|q_1 + E_d q_2 - E_d^0 q_2^0 \\[4pt]
|t_2| &= |t_2^0| - |t_1|q_2 - E_d q_1 + E_d^0 q_1^0
\end{aligned}}$$

---

## 六、$q$ 的自洽 Riccati ODE

将 $\mathbf{u} = \mathbf{u}(q)$ 代入 Riccati ODE：

$$\dot q = C(\mathbf{u}(q)) + D(\mathbf{u}(q)) q - q A(\mathbf{u}(q)) - q B(\mathbf{u}(q)) q$$

得到 $\dot q = F(q; K_{\text{eff}}^0)$，右端是 $q$ 和已知函数 $K_{\text{eff}}^0(t)$ 的显式函数。

**$F(q)$ 的具体构造**：
1. 由上述公式用 $q, q^0$ 算出 $|t_1|, E_d, |t_2|, |t_3|$
2. 构造 $A, B, C, D$：
   $$A = \frac{E_1+|t_3|}{2}\mathbf i + \frac{|t_2|}{2}\mathbf j,\quad D = \frac{-E_1+|t_3|}{2}\mathbf i + \frac{|t_2|}{2}\mathbf j$$
   $$B = \frac{|t_1|}{2} + \frac{E_d}{2}\mathbf k,\quad C = -\frac{|t_1|}{2} + \frac{E_d}{2}\mathbf k$$
3. 代入 Riccati 右端：$F(q) = C + Dq - qA - qBq$

---

## 七、可解性分析

### 7.1 $t=0$ 的奇点

$q(0)=0$ 时矩阵退化（$q_0=q_3=0$）。由 Riccati ODE：

$$\dot q(0) = C(0) = -\frac{|t_1(0)|}{2} + \frac{E_d(0)}{2}\mathbf k = \frac{E_0}{2}\mathbf k$$

（$t_1(0)=0$ 因 Step 1 开始 G1 关断）。故 $q(\epsilon) \approx \epsilon \cdot \frac{E_0}{2}\mathbf k$。
$q_0^2+q_3^2 \propto \epsilon^2$，分母二次零 → 控制解在 $t=0$ 处**一阶极点**。

**处理方案**：在 $[0, \delta]$ 的初始段使用原始控制协议（$t_1=0$ 或使用理想值），
$\delta$ 后切换到 DAE 控制。$\delta$ 可选为 $q_0^2+q_3^2 \ge \varepsilon$ 的首达时刻。

### 7.2 全局可解性

自治 ODE $\dot q = F(q)$ 的右端在 $\{q: q_0^2+q_3^2 > 0\}$ 上**实解析**（有理函数 + 平方根最多来自 $|t_1|$ 等绝对值，但 $t_c, E_0, t_1$ 取绝对值后是连续可微远离零的区域）。

由标准 ODE 理论，若 $q(t)$ 在 $[0,T]$ 上有界，则解存在。核心问题是 $F(q)$ 的二次非线性是否导致有限时间 blow-up。这需要通过数值探索。

### 7.3 与标准 Riccati 方程的差异

标准 Riccati $\dot q = C + Dq - qA - qBq$ 中 $A,B,C,D$ 是**外生控制函数**。
自洽 Riccati $\dot q = F(q)$ 中 $A,B,C,D$ 是 **$q$ 的函数**——形成反馈回路。

反馈可能稳定（$|t_1|, E_d$ 自动调节抑制 $q$ 增长），也可能不稳定（正反馈放大）。
这类似于控制系统中的"模型参考自适应控制"——$K_{\text{eff}}^0$ 是参考模型，
$\mathbf{u}(q)$ 是自适应律。

---

## 八、与 PRB113 对比

| | PRB113 | PRB111 DAE 方案 |
|---|---|---|
| 抵消方式 | 恢复瞬时简并 | 匹配 $K_{\text{eff}}$ |
| 条件 | $\varepsilon + \sigma\lambda = 0$ | $K_{\text{eff}} = K_{\text{eff}}^0$ |
| 条件类型 | 1 个代数方程 | 4 个代数方程 + 4 个 ODE |
| 涉及 $q$？ | 不需要（二分图可解析对角化） | 需要（$q$ 通过 Riccati ODE 决定） |
| 控制显式 | $\lambda = f(\varepsilon,\sigma)$ | $\mathbf{u} = \mathbf{u}(q, q^0)$ |
| 额外自由度 | 富余 MBS 对 | 无——Riccati 框架自身闭合 |
| 奇点 | 无 | $q=0$ 处 $|t_1|, E_d$ 奇异 |

---

## 九、实施路线

1. **预计算**：传播理想 Riccati ODE，输出 $q^0(t)$ 和 $K_{\text{eff}}^0(t)$ 的时间序列
2. **积分自洽 ODE**：从 $q(\delta) = q^0(\delta)$ 开始（初始段用理想控制），
   用自适应 RK45 积分 $\dot q = F(q; K_{\text{eff}}^0)$ 到 $T$
3. **反算控制**：从 $q(t)$ 用 §5 公式计算 $|t_1(t)|, |t_2(t)|, |t_3(t)|, E_d(t)$
4. **验证**：用计算出的控制波形跑完整 Sp(2) 传播，检查 $R_{123} = R_{123}^{\text{ideal}}$
