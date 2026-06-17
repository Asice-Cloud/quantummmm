## $A \neq D$ 对称性破缺释放第三分量：Riccati 严格证明

1: 设

$$A = a_1\mathbf{i} + a_2\mathbf{j}, \quad D = d_1\mathbf{i} + d_2\mathbf{j}, \quad C = c_1\mathbf{i} + c_2\mathbf{j}, \quad B = b_1\mathbf{i} + b_2\mathbf{j}$$

$$a_1 = \frac{E_1+|t_3|}{2},\; a_2 = \frac{|t_2|}{2},\; d_1 = \frac{-E_1+|t_3|}{2},\; d_2 = \frac{|t_2|}{2}$$

$$c_1 = \frac{|t_2|}{2},\; c_2 = \frac{-E_1+|t_3|}{2},\; b_1 = -\frac{|t_2|}{2},\; b_2 = -\frac{E_1+|t_3|}{2}$$

$$q = \alpha\mathbf{i} + \beta\mathbf{j} + \gamma\mathbf{k}$$

关键恒等式：$A-D = E_1\mathbf{i}$。$E_1=0 \Leftrightarrow A=D$。



2 约定:

将 $q$ 展开为标准基：

$$q = \alpha\mathbf{i} + \beta\mathbf{j} + \gamma\mathbf{k},\qquad \alpha,\beta,\gamma \in \mathbb{R}$$

基乘法：$\mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = -1$，$\mathbf{i}\mathbf{j} = \mathbf{k} = -\mathbf{j}\mathbf{i}$，$\mathbf{j}\mathbf{k} = \mathbf{i} = -\mathbf{k}\mathbf{j}$，$\mathbf{k}\mathbf{i} = \mathbf{j} = -\mathbf{i}\mathbf{k}$。

对任意纯四元数 $p = p_1\mathbf{i} + p_2\mathbf{j} + p_3\mathbf{k}$，定义 $(p)_k = p_3$ 为 $\mathbf{k}$ 基的系数。

### 3  $qBq$ 的 $\mathbf{k}$ 分量

$B = b_1\mathbf{i} + b_2\mathbf{j}$，$b_1 = -|t_2|/2$，$b_2 = -(E_1+|t_3|)/2$。

先算 $Bq$：

$$Bq = (b_1\mathbf{i} + b_2\mathbf{j})(\alpha\mathbf{i} + \beta\mathbf{j} + \gamma\mathbf{k})$$

逐项：
$$\begin{aligned}
b_1\mathbf{i} \cdot \alpha\mathbf{i} &= -b_1\alpha \\
b_1\mathbf{i} \cdot \beta\mathbf{j} &= b_1\beta\mathbf{k} \\
b_1\mathbf{i} \cdot \gamma\mathbf{k} &= -b_1\gamma\mathbf{j} \\
b_2\mathbf{j} \cdot \alpha\mathbf{i} &= -b_2\alpha\mathbf{k} \\
b_2\mathbf{j} \cdot \beta\mathbf{j} &= -b_2\beta \\
b_2\mathbf{j} \cdot \gamma\mathbf{k} &= b_2\gamma\mathbf{i}
\end{aligned}$$

$$Bq = -(b_1\alpha + b_2\beta) + b_2\gamma\,\mathbf{i} - b_1\gamma\,\mathbf{j} + (b_1\beta - b_2\alpha)\mathbf{k}$$

记 $Bq = s + v_1\mathbf{i} + v_2\mathbf{j} + v_3\mathbf{k}$，其中 $s = -(b_1\alpha+b_2\beta)$，$v_1 = b_2\gamma$，$v_2 = -b_1\gamma$，$v_3 = b_1\beta-b_2\alpha$。

现在算 $q(Bq) = (\alpha\mathbf{i} + \beta\mathbf{j} + \gamma\mathbf{k})(s + v_1\mathbf{i} + v_2\mathbf{j} + v_3\mathbf{k})$，只提取 $\mathbf{k}$ 分量：

$$(qBq)_k = \gamma s + \alpha v_2 - \beta v_1$$

代入 $s, v_1, v_2$：

$$\begin{aligned}
(qBq)_k &= \gamma(-b_1\alpha - b_2\beta) + \alpha(-b_1\gamma) - \beta(b_2\gamma) \\
&= -b_1\alpha\gamma - b_2\beta\gamma - b_1\alpha\gamma - b_2\beta\gamma \\
&= -2\gamma(b_1\alpha + b_2\beta)
\end{aligned}$$

代入 $b_1 = -|t_2|/2$, $b_2 = -(E_1+|t_3|)/2$：

$$\boxed{(qBq)_k = \gamma\bigl(|t_2|\alpha + (E_1+|t_3|)\beta\bigr)}$$

**关键性质**：$(qBq)_k \propto \gamma$。$\gamma = 0 \;\Rightarrow\; (qBq)_k = 0$——非线性项不产生初始 $\mathbf{k}$ 分量，只能放大已有的。

### 4  $C, Dq, qA$ 的 $\mathbf{k}$ 分量

$$C = c_1\mathbf{i} + c_2\mathbf{j} = \frac{|t_2|}{2}\mathbf{i} + \frac{-E_1+|t_3|}{2}\mathbf{j}$$

$(C)_k = 0$。

$Dq - qA$ 的 $\mathbf{k}$ 分量已在 §4.3bis 算出：$(Dq-qA)_k = |t_3|\beta - |t_2|\alpha$。

### 5 完整分量 ODE

$$\boxed{\dot\alpha = \frac{|t_2|}{2} + |t_2|\gamma + \frac{|t_2|}{2}(\alpha^2 + \beta^2)} \tag{1}$$

$$\boxed{\dot\beta = \frac{-E_1 + |t_3|}{2} - |t_3|\gamma + \frac{E_1 + |t_3|}{2}(\alpha^2 + \beta^2)} \tag{2}$$

$$\boxed{\dot\gamma = |t_3|\beta - |t_2|\alpha - \gamma\bigl(|t_2|\alpha + (E_1 + |t_3|)\beta\bigr)} \tag{3}$$

初始：$\alpha(0)=\beta(0)=\gamma(0)=0$。

### 6 $E_1=0$ 时：$A=D$，对称性保护 $\gamma \equiv 0$

令 $E_1=0$，$\gamma=0$。(1)(2)(3) 简化为：

$$\dot\alpha = \frac{|t_2|}{2}(1 + \alpha^2 + \beta^2)$$

$$\dot\beta = \frac{|t_3|}{2}(1 + \alpha^2 + \beta^2)$$

$$\dot\gamma = |t_3|\beta - |t_2|\alpha$$

由前两式：

$$\frac{\dot\alpha}{|t_2|} = \frac{\dot\beta}{|t_3|} \tag{4}$$

(4) 是 $A=D$ 的直接后果——对称的 $A$ 和 $D$ 对 $\alpha$ 和 $\beta$ 施加了等比例的驱动力。

(4) 意味着 $\alpha,\beta$ 的演化保持比例 $\alpha/|t_2| = \beta/|t_3|$（初始均为0，故保持比例恒等）。由此 $|t_3|\beta - |t_2|\alpha \equiv 0$，$\dot\gamma \equiv 0$，$\gamma \equiv 0$ 自洽。

**$E_1=0$ 时 $\gamma \equiv 0$ 是 ODE 的不动解。第三分量永不激活。**

### 7 $E_1 \neq 0$ 时：$A \neq D$，对称性破缺

$E_1 \neq 0$ 时，(2) 变为：

$$\dot\beta = \frac{-E_1 + |t_3|}{2} + \frac{E_1 + |t_3|}{2}(\alpha^2 + \beta^2) \quad (\gamma=0)$$

(1) 不变。此时：

$$\frac{\dot\beta}{|t_3|} = \frac{1}{2}(1 + \alpha^2 + \beta^2) + \frac{E_1}{2|t_3|}\bigl(-1 + (\alpha^2 + \beta^2)\bigr)$$

$$\frac{\dot\alpha}{|t_2|} = \frac{1}{2}(1 + \alpha^2 + \beta^2)$$

$$\frac{\dot\alpha}{|t_2|} \neq \frac{\dot\beta}{|t_3|}$$

**(4) 被 $E_1$ 破坏。** $\alpha,\beta$ 不再保持等比例演化 → $|t_3|\beta - |t_2|\alpha$ 不再恒零 → $\dot\gamma \neq 0$ → $\gamma$ 开始生长。

### 8 $qBq$ 非线性放大

一旦 $\gamma \neq 0$，由 §12.3 的 $(qBq)_k$ 公式：

$$(qBq)_k = \gamma(|t_2|\alpha + (E_1 + |t_3|)\beta)$$

$$\dot\gamma = \underbrace{|t_3|\beta - |t_2|\alpha}_{\text{线性种子（被 }E_1\text{ 点燃）}} \;-\; \underbrace{\gamma(|t_2|\alpha + (E_1 + |t_3|)\beta)}_{\text{非线性放大（}E_1\text{ 增强增益）}}$$

$\gamma \uparrow \;\Rightarrow\; |(qBq)_k| \uparrow \;\Rightarrow\; |\dot\gamma| \uparrow$ —— **正反馈**。$E_1$ 既在 $b_2 = -(E_1+|t_3|)/2$ 中增强放大系数，又在 (2) 中维持不对称。

### 9 结论

$$\boxed{E_1=0 \;\Leftrightarrow\; A=D \;\Rightarrow\; \frac{\dot\alpha}{|t_2|} = \frac{\dot\beta}{|t_3|} \;\Rightarrow\; \gamma \equiv 0}$$

$$\boxed{E_1 \neq 0 \;\Leftrightarrow\; A \neq D \;\Rightarrow\; \frac{\dot\alpha}{|t_2|} \neq \frac{\dot\beta}{|t_3|} \;\xrightarrow{\;\dot\gamma\neq 0\;}\; \gamma > 0 \;\xrightarrow{\;qBq\;}\; \text{正反馈放大}}$$

**$E_1 \neq 0$ 打破 $A=D$ 对称性 → $\alpha,\beta$ 演化比例失调 → $\gamma$ 获得种子 → $qBq$ 非线性正反馈 → 第三分量完全释放。**