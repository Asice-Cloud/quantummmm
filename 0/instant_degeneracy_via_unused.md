# 添加未用双线性实现瞬时简并：$i\gamma_1\gamma_3$ 方案

> 对标 PRB113 的 $\varepsilon + \sigma\lambda = 0$，在 PRB111 中用未用自由度实现同类抵消

---

## 一、方案

在原始哈密顿量上加未用双线性 $i\gamma_1\gamma_3$（系数 $\delta$）：

$$H_{\text{total}} = H_{EM} + \delta \cdot i\gamma_1\gamma_3$$

在 Sp(2) 表示中，$\Sigma_{13} = \frac14[\Gamma_1, \Gamma_3] = \begin{pmatrix}\mathbf j/2 & 0 \\ 0 & -\mathbf j/2\end{pmatrix}$。

这使得 $A, D$ 的 $\mathbf j$ 分量不再相等：

$$A = \frac{E_1+|t_3|}{2}\mathbf i + \frac{|t_2|+\delta}{2}\mathbf j$$
$$D = \frac{-E_1+|t_3|}{2}\mathbf i + \frac{|t_2|-\delta}{2}\mathbf j$$

$$A - D = E_1\mathbf i + \delta\mathbf j$$

---

## 二、瞬时简并条件

sympy 计算给出 $K^2 \propto I$ 的条件：

$$\boxed{E_1|t_3| + \delta|t_2| = 0\quad\text{（对角色块差为零）}}$$

$$\boxed{E_d\delta + |t_3||t_1| = 0,\quad E_1E_d = |t_2||t_1| \quad\text{（非对角块为零）}}$$

由此解得：

$$\boxed{\delta = -\frac{E_1|t_3|}{|t_2|},\qquad |t_1| = \frac{E_1 E_d}{|t_2|}}$$

在 $|t_2| \neq 0$ 且所有参数为正的物理条件下，这使 $K^2 \propto I$ 严格成立。

---

## 三、与 PRB113 的对应

| | PRB113 | 本方案 |
|---|---|---|
| 修正项 | $H_{\text{corr}} = \Lambda(i\gamma_2\tilde\gamma_2 + i\gamma_3\tilde\gamma_3)$ | $\delta \cdot i\gamma_1\gamma_3$ |
| 简并条件 | $\varepsilon + \sigma\lambda = 0$ | $\delta = -E_1|t_3|/|t_2|$，$|t_1| = E_1E_d/|t_2|$ |
| 方程数 | 1 | 2（耦合的） |
| 额外自由度来源 | 富余 MBS 对 | PRB111 未用的双线性 |
| 涉及的物理耦合 | MBS 对内部杂化（可静电调） | $\gamma_1$–$\gamma_3$ 耦合 + $t_1$ 自适应 |

---

## 四、物理可行性

### 4.1 $i\gamma_1\gamma_3$ 能否物理实现？

$i\gamma_1\gamma_3$ 是 $\gamma_1$ 和 $\gamma_3$ 之间的直接耦合。在纳米线几何中，如果 $\gamma_1$ 和 $\gamma_3$ 的波函数有空间交叠，或通过一个额外的静电门可调，则原则上可实现。这与 $t_1$（$\gamma_b$–$\gamma_1$ 耦合）和 $t_2$（$\gamma_a$–$\gamma_2$ 耦合）的机制类似——需要 $\gamma_1$ 和 $\gamma_3$ 之间的可调隧道耦合。

**实验上这不一定容易**：$\gamma_1$ 和 $\gamma_3$ 是不同 MZM，可能空间分离较远。但不排除通过第三个量子点中介。

### 4.2 $t_1$ 被锁定

条件 $|t_1| = E_1E_d/|t_2|$ 意味着 $t_1$ 不再是自由参数——它必须自适应跟踪 $t_2$ 和 $E_d$。这实际上是**将 $t_1$ 升级为主动控制通道**，而非被动误差源。

### 4.3 $t_2 \to 0$ 的奇点

当 $|t_2| \to 0$（Step 1 初始和 Step 3），$\delta$ 和 $|t_1|$ 都会发散。这是方案的**主要物理限制**：

- **Step 1 起始**：$t_2$ 从 0 开始爬升 → 初始阶段 $\delta$ 很大，$t_1$ 很大。$t_1$ 有最大物理值（门耦合上限），所以方案在 $t_2$ 极小阶段不可行。
- **Step 3**：$t_2 = 0$ 全程 → 方案完全失效。Step 3 的动力学相位无法通过此机制消除。

### 4.4 可行区间

方案有效的区间：$t_2$ 足够大，使 $|t_1| = E_1E_d/t_2$ 和 $|\delta| = E_1|t_3|/t_2$ 都在物理可及范围内。

设 $t_1^{\max}$ 为 $t_1$ 的最大物理值（门耦合强度上限），则要求：

$$t_2 \ge \frac{E_1 E_d}{t_1^{\max}},\quad |t_3| \le \frac{|\delta|^{\max} t_2}{E_1}$$

---

## 五、与理想编织的差异

理想编织要求 $E_1 = 0, t_1 = 0$。本方案中：

- $E_1 \neq 0$（接受缺陷）
- $t_1 \neq 0$（用作控制通道）
- 增加 $\delta \cdot i\gamma_1\gamma_3 \neq 0$（修正项）

**这不是"让 $E_1$ 看起来为零"的魔法。这是"用 $t_1$ 和 $\delta$ 构建一个等价动力学"——一个不同的路径，但到达相同的终点。**

代价：
1. 需要物理上实现 $i\gamma_1\gamma_3$ 耦合
2. $t_1$ 失去参数自由度，被简并条件锁定
3. $t_2 \approx 0$ 区域不可行

---

## 六、残余问题：Step 3

Step 3 中 $t_2 = 0$，$i\gamma_1\gamma_3$ 方案失效。可能补充方案：

1. **缩短 Step 3**：使 $\tau_3$ 极小，Step 3 累积的动力学相位可忽略
2. **前两步预补偿**：在 Step 1/2 中刻意累积相反的动力学相位，使三步总和为零
3. **在 Step 3 换用不同未用双线性**：如 $i\gamma_3\gamma_b$（影响非对角块），需要单独分析

---

## 七、总结

| 项目 | 状态 |
|------|:---:|
| 能否瞬时简并（$t_2 \neq 0$）？ | ✅ 严格证明 |
| 简并条件是否有显式？ | ✅ $\delta = -E_1|t_3|/|t_2|$，$|t_1| = E_1E_d/|t_2|$ |
| 是否需要额外 Majorana？ | ❌ 只用 PRB111 已有的 5 个 |
| $t_2=0$ 区间 | ❌ 方案失效 |
| 物理可实现性 | ⚠️ 需 $i\gamma_1\gamma_3$ 耦合 + $t_1$ 自适应跟踪 |
