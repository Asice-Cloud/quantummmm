# 绝热极限下的不动点与纯几何编织条件

> 出发点：new3.md 的严格推导结果，不引入额外近似。
> 目标：(1) 证明绝热极限下 Riccati 方程存在稳定不动点；(2) 给出纯几何编织的代数含义；(3) 分析其可实现性。

---

## 1. 出发点：Riccati 方程

由 new3.md §2.4，$q(t) = Z(t)X(t)^{-1}$ 满足精确的四元数 Riccati 方程

$$
\dot q = C + Dq - qA - qBq, \tag{R}
$$

系数为（$t_3=0$ 情形暂不设，保持一般性）

$$
A = \frac{E_1+|t_3|}{2}\,\mathbf{i} + \frac{|t_2|}{2}\,\mathbf{j},\quad
D = \frac{-E_1+|t_3|}{2}\,\mathbf{i} + \frac{|t_2|}{2}\,\mathbf{j},
$$

$$
B = \frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf{k},\quad
C = -\frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf{k}.
$$

---

## 2. 绝热极限与不动点方程

**绝热条件**：参数 $t_k(t), E_d(t)$ 变化缓慢，满足

$$
|\dot A|, |\dot B|, |\dot C|, |\dot D| \ll E_{\rm gap}^2,
$$

其中 $E_{\rm gap}$ 是 ancilla 的能隙（由 $B, C, E_d$ 决定）。

在此条件下，$q(t)$ 能跟随参数的瞬时不动点 $q^*(t)$，即满足

$$
\boxed{0 = C + Dq^* - q^*A - q^*Bq^*.} \tag{FP}
$$

这是一个四元数代数 Riccati 方程（QARE），给定参数 $(A,B,C,D)$ 后，$q^*$ 由此方程决定。

### 2.1 不动点的存在性

将 $q^* = q_0 + q_x\mathbf{i} + q_y\mathbf{j} + q_z\mathbf{k}$ 代入 (FP)，展开为实数方程组。设参数为实数 $\alpha = (E_1+|t_3|)/2$、$\beta = |t_2|/2$、$p = |t_1|/2$、$e = E_d/2$，则

$$
A = \alpha\mathbf{i}+\beta\mathbf{j},\quad D = \delta\mathbf{i}+\beta\mathbf{j}\;(\delta=(-E_1+|t_3|)/2),\quad B = p + e\mathbf{k},\quad C = -p + e\mathbf{k}.
$$

展开 (FP) 的四个实分量方程（$\mathbf{1},\mathbf{i},\mathbf{j},\mathbf{k}$ 各一个），得到关于 $(q_0,q_x,q_y,q_z)$ 的耦合非线性方程组。

**数值观察**（来自 verify_adiabatic.py 的系统扫描）：对协议参数的全范围扫描，代数 Riccati 方程均存在唯一稳定实四元数解 $q^*$，且 $|q^*| < 1$（保证 $X(t)$ 可逆）。

**关键数值结论**：对所有扫描参数，不动点均满足

$$
q_z^* = 0, \qquad q_0^* = 0.
$$

这是模型结构（$B, C$ 的特殊形式）的内禀性质，不依赖参数取值。

---

## 3. 纯几何编织的代数含义

### 3.1 $H_{\rm eff}$ 的分量结构

由 new3.md §4.2，绝热极限下 $H_{\rm eff}$ 的三个 MZM 耦合强度为

$$
\epsilon_x = \frac{E_1+|t_3|}{2} + [Bq^*]_\mathbf{i},\quad
\epsilon_y = \frac{|t_2|}{2} + [Bq^*]_\mathbf{j},\quad
\epsilon_z = [Bq^*]_\mathbf{k}.
$$

其中

$$
[Bq^*]_\mathbf{k} = \frac{|t_1|}{2}q_z^* + \frac{E_d}{2}q_0^*.
$$

由数值结论 $q_z^*=q_0^*=0$，立即得到

$$
\epsilon_z^* = 0 \quad \text{（自动成立，无需额外调参）}.
$$

### 3.2 纯几何编织的定义与必要条件

编织操作的目标幺正变换（对 $\gamma_2,\gamma_3$ 交换两次）为

$$
U_{\rm ideal} = e^{i\frac{\pi}{4}i\gamma_2\gamma_3} = e^{i\frac{\pi}{4}\sigma_x}.
$$

**定义**：若整个编织过程中

$$
\epsilon_y(t) = 0 \quad \text{且} \quad \epsilon_z(t) = 0 \quad \forall\, t\in[0,6\tau],
$$

则 $H_{\rm eff}(t) = \epsilon_x(t)\,i\gamma_2\gamma_3$，演化算符退化为

$$
U(6\tau) = \exp\!\left(-i\,i\gamma_2\gamma_3\int_0^{6\tau}\epsilon_x(t)\,dt\right).
$$

此时旋转轴固定为 $\sigma_x$，编织结果只依赖积分 $\int\epsilon_x\,dt$，与路径的时间参数化无关（在绝热极限下 $\epsilon_x$ 由协议几何决定）——这是"纯几何"的含义。

**注意**：$\epsilon_y = \epsilon_z = 0$ 是旋转轴固定为 $\sigma_x$ 的**充分条件**，但它是否构成理想编织的必要条件，需要进一步分析 $U(6\tau)$ 对 $\epsilon_y, \epsilon_z$ 的灵敏度。此处将其作为充分条件使用。

---

## 4. 两个条件的性质差异

$$
\epsilon_z^* = [Bq^*]_\mathbf{k} = \frac{|t_1|}{2}q_z^* + \frac{E_d}{2}q_0^*
$$

$$
\epsilon_y^* = \frac{|t_2|}{2} + [Bq^*]_\mathbf{j} = \frac{|t_2|}{2} + \frac{|t_1|}{2}q_y^* - \frac{E_d}{2}q_x^*
$$

| 条件 | 来源 | 在不动点处的状态 | 可控性 |
|---|---|---|---|
| $\epsilon_z^* = 0$ | 纯 ancilla 反馈 $[Bq^*]_\mathbf{k}$ | **自动满足**（$q_z^*=q_0^*=0$） | 无需调参 |
| $\epsilon_y^* = 0$ | 几何项 $|t_2|/2$ + ancilla 反馈 $[Bq^*]_\mathbf{j}$ | **不自动满足**，依赖 $t_2, q_y^*, q_x^*$ | 需要约束 $t_2(t)$ |

因此两个条件的物理地位截然不同：

- $\epsilon_z^* = 0$ 是模型的**内禀性质**，由 Riccati 方程的结构保证
- $\epsilon_y^* = 0$ 是对协议参数 $t_2(t)$ 的**自洽约束**，需要在设计协议时显式满足

---

## 5. $\epsilon_y^* = 0$ 的自洽约束

将 $\epsilon_y^* = 0$ 代入表达式：

$$
\frac{|t_2|}{2} = -[Bq^*]_\mathbf{j} = -\frac{|t_1|}{2}q_y^* + \frac{E_d}{2}q_x^*.
$$

但 $q^*(t)$ 本身由不动点方程 (FP) 决定，其中包含 $t_2$（通过 $A$ 和 $D$ 的 $\mathbf{j}$ 分量）。因此这是一个**隐式自洽方程**：

$$
\frac{|t_2|}{2} = F\!\left(E_1, |t_1|, |t_2|, |t_3|, E_d\right),
$$

其中右端通过 $q^*(t_2)$ 隐式依赖 $t_2$。

这个方程是否有解、解是否唯一，需要数值分析。

---

## 6. 补偿方案的框架（借鉴 PRB113）

PRB113 的补偿思路：在哈密顿量里引入额外项 $H_{\rm corr} = \Lambda(t)\,i\gamma_2\tilde\gamma_2 + \Lambda(t)\,i\gamma_3\tilde\gamma_3$，使得总哈密顿量的不动点自动满足简并条件。额外参数 $\Lambda(t)$ 是补偿的自由度。

在我们的模型里，对应的补偿思路是：**引入一个直接作用在 MZM 子空间、针对 $\epsilon_y$ 的补偿耦合**，使得修改后的 $\epsilon_y^*$ 在绝热不动点处自动为零，而不需要对 $t_2(t)$ 施加隐式自洽约束。

具体可行方案和补偿参数的选取，见 plan_1.md（待补充）。

---

## 附：不动点方程的分量展开

设 $q^* = q_0 + q_x\mathbf{i} + q_y\mathbf{j} + q_z\mathbf{k}$，将 $C + Dq^* - q^*A - q^*Bq^* = 0$ 展开：

四元数乘法规则：$\mathbf{i}^2=\mathbf{j}^2=\mathbf{k}^2=-1$，$\mathbf{ij}=\mathbf{k}$，$\mathbf{jk}=\mathbf{i}$，$\mathbf{ki}=\mathbf{j}$（及反交换）。

$Dq^*$ 的各分量：

$$
Dq^* = (\delta\mathbf{i}+\beta\mathbf{j})(q_0+q_x\mathbf{i}+q_y\mathbf{j}+q_z\mathbf{k})
$$

$$
= \delta q_0\mathbf{i} - \delta q_x + \delta q_y\mathbf{k} - \delta q_z\mathbf{j}
+ \beta q_0\mathbf{j} - \beta q_x\mathbf{k} - \beta q_y + \beta q_z\mathbf{i}
$$

$q^*A$ 的各分量：

$$
q^*A = (q_0+q_x\mathbf{i}+q_y\mathbf{j}+q_z\mathbf{k})(\alpha\mathbf{i}+\beta\mathbf{j})
$$

$$
= \alpha q_0\mathbf{i} - \alpha q_x + \alpha q_y\mathbf{k} - \alpha q_z\mathbf{j}
+ \beta q_0\mathbf{j} - \beta q_x\mathbf{k} - \beta q_y + \beta q_z\mathbf{i}
$$

$q^*Bq^*$ 需展开 $B = p + e\mathbf{k}$：

$$
Bq^* = (p+e\mathbf{k})(q_0+q_x\mathbf{i}+q_y\mathbf{j}+q_z\mathbf{k})
= p q^* + e\mathbf{k}q^*
$$

$$
e\mathbf{k}q^* = e(-q_z + q_y\mathbf{i} - q_x\mathbf{j} + q_0\mathbf{k})
$$

故

$$
Bq^* = (pq_0 - eq_z) + (pq_x+eq_y)\mathbf{i} + (pq_y-eq_x)\mathbf{j} + (pq_z+eq_0)\mathbf{k}
$$

再左乘 $q^*$：

$$
q^*(Bq^*) = q^* \cdot [(pq_0-eq_z) + (pq_x+eq_y)\mathbf{i} + (pq_y-eq_x)\mathbf{j} + (pq_z+eq_0)\mathbf{k}]
$$

记 $Bq^* = r_0 + r_x\mathbf{i}+r_y\mathbf{j}+r_z\mathbf{k}$，则

$$
[q^*Bq^*]_{\mathbf{1}} = q_0 r_0 - q_x r_x - q_y r_y - q_z r_z,
$$
$$
[q^*Bq^*]_{\mathbf{i}} = q_0 r_x + q_x r_0 + q_y r_z - q_z r_y,
$$
$$
[q^*Bq^*]_{\mathbf{j}} = q_0 r_y - q_x r_z + q_y r_0 + q_z r_x,
$$
$$
[q^*Bq^*]_{\mathbf{k}} = q_0 r_z + q_x r_y - q_y r_x + q_z r_0.
$$

将以上代入 $C + Dq^* - q^*A - q^*Bq^* = 0$ 的四个分量，得到关于 $(q_0,q_x,q_y,q_z)$ 的四个耦合实方程，可数值求解。
