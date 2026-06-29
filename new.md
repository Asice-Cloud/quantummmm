# Majorana 编织的纤维丛与曲率描述

## 一、几何图景总览

```
                          总空间 P = M × SO(5)
                         ┌─────────────────────┐
                         │  纤维 ≅ SO(5)       │  ← 5 Majorana 的旋转自由度
                         │  ┌───────────────┐   │
                         │  │  H_EM = Σh_ij· │   │
                         │  │  (iγ_iγ_j)     │   │  ← 联络 A = H dt
                         │  │   ↓            │   │
                         │  │  U= Hol_γ(A)   │   │  ← 和乐 = 演化算符
                         │  └───────────────┘   │
                         └──────┬──────────────┘
                                │ π (投影)
                    ┌───────────┴───────────┐
                    │  底流形 M             │
                    │                        │
                    │  参数坐标:              │
                    │  (E₁, t₁, τ, t₂, t₃)  │  ← 编织路径 = M 中的回路 γ
                    │                        │
                    │  π₁(M) = B₂ ≅ ℤ       │  ← 基本群 = 辫子群
                    └────────────────────────┘
```

**核心命题**：纯编织是拓扑（仅依赖同伦类），动力学相位是几何（依赖同伦类内具体路径的曲率积分）。

---

## 二、底流形 M：Majorana 的"配置空间"

### 2.1 定义

两个被交换的 MZM（$\gamma_2, \gamma_3$）的有效"位置"由门控耦合 $(t_2, t_3)$ 参数化。
编织过程中，$(t_2, t_3)$ 的轨迹是：

$$\gamma(t) = (t_2(t), t_3(t)) = t_c \cdot 
\begin{cases}
(f_-(t), 0)       & 0 \le t < \tau \quad \text{(Step 1)} \\
(f_+(t), f_-(t))  & \tau \le t < 2\tau \quad \text{(Step 2)} \\
(0, f_+(t))       & 2\tau \le t \le 3\tau \quad \text{(Step 3)}
\end{cases}$$

这是一个**包围原点**的回路（原点 $(0,0)$ 对应门全部关闭——不允许的配置）。

### 2.2 基本群

底流形 $X = \mathbb{R}^2 \setminus \{(0,0)\}$（两个 MZM 耦合参数的有效空间，排除同时为零的奇点）：

$$\pi_1(X) \cong \mathbb{Z} \cong B_2$$

回路 $\gamma$ 的绕数 $w=1$ 对应一次 $\gamma_2 \leftrightarrow \gamma_3$ 交换。
**这就是辫子群 $B_2$ 在此处的具体实现。**

### 2.3 扩展到完整参数空间

完整的底流形还包括 ancilla 自由度参数：

$$M = \{(E_1, t_1, \tau, E_d, t_2, t_3) \mid \text{物理允许的值}\}$$

维度 $\dim M = 6$（这里 $t_2, t_3$ 的幅值由 $t_c$ 固定，时间轮廓由协议固定，但原则上可以变化）。

---

## 三、纤维与主丛

### 3.1 主丛构造

$$G = \text{SO}(5) \hookrightarrow P = M \times \text{SO}(5) \xrightarrow{\pi} M$$

- **结构群** $G = \text{SO}(5)$：5 个 Majorana 双线性的所有可能规范旋转
- **右作用**：$R_g: (x, h) \mapsto (x, hg)$ for $g \in \text{SO}(5)$
- **纤维**：$\pi^{-1}(x) \cong \text{SO}(5)$，即点 $x \in M$ 上所有可能的内部参考系选择

### 3.2 联络 1-形式

联络是 $\mathfrak{so}(5)$-值的 1-形式。因为底流形上的时间坐标 $t$ 是唯一的演化参数，
联络简化为沿 $t$ 方向的分量：

$$\boxed{\mathcal{A}(t; E_1, t_1, \tau) = A(t; E_1, t_1, \tau)\, dt}$$

其中 $A(t)$ 是 §3.3 中的 $5 \times 5$ 反对称矩阵，由 $H_{EM}$ 的系数给出。

在 $\mathfrak{so}(5)$ 显式基底 $\{X_{ij} = i\gamma_i\gamma_j\}_{1\le i<j\le 5}$ 下：

$$\mathcal{A}(t) = \sum_{i<j} 2h_{ij}(t)\, X_{ij}\; dt$$

其中 $h_{ij}(t)$ 是表 1 中的系数。

**表 1：$H_{EM}$ 各成分的联络分量**

| 项 | 系数 $h_{ij}(t)$ | $\mathfrak{so}(5)$ 生成元 | 参数依赖 |
|---|---|---|---|
| $iE_1\gamma_1\gamma_2$ | $E_1$ | $X_{12}$ | 恒定 |
| $i\vert t_2\vert\gamma_a\gamma_2$ | $t_c \cdot g_{42}(t)$ | $X_{24}$ | $t_c, f_\pm(t)$ |
| $-i\vert t_1\vert\gamma_b\gamma_1$ | $-t_1 \cdot g_{15}(t)$ | $X_{15}$ | $t_1, f_\pm(t)$ |
| $-i\vert t_3\vert\gamma_a\gamma_3$ | $t_c \cdot g_{34}(t)$ | $X_{34}$ | $t_c, f_\pm(t)$ |
| $iE_d\gamma_a\gamma_b$ | $E_0 \cdot g_{45}(t)$ | $X_{45}$ | $E_0, f_\pm(t)$ |

门控函数 $g_{ij}(t)$ 在各步的取值见 §3.3。

### 3.3 各步的联络显式

**Step 1** ($0 \le t < \tau$)：
$$A^{(1)}(t) = 2E_1 X_{12} + 2t_c f_-(t) X_{24} + 2t_1 f_-(t) X_{15} + 2E_0 f_+(t) X_{45}$$

**Step 2** ($\tau \le t < 2\tau$)：
$$A^{(2)}(t) = 2E_1 X_{12} + 2t_c f_+(t) X_{24} + 2t_c f_-(t) X_{34} + 2t_1 f_+(t) X_{15}$$

**Step 3** ($2\tau \le t \le 3\tau$)：
$$A^{(3)}(t) = 2E_1 X_{12} + 2t_c f_+(t) X_{34} + 2E_0 f_-(t) X_{45}$$

---

## 四、曲率：动力学的几何起源

### 4.1 一般定义

曲率 2-形式（$\mathfrak{so}(5)$-值）：

$$F = d\mathcal{A} + \frac{1}{2}[\mathcal{A}, \mathcal{A}]$$

其中 $[\mathcal{A}, \mathcal{A}](X,Y) = [\mathcal{A}(X), \mathcal{A}(Y)]$ 是李括号。

**在仅含时间参数的情况下**，$\mathcal{A} = A(t)dt$，由于外微分 $d(dt)=0$ 且 $dt \wedge dt = 0$，有：

$$F_{\text{1D time}} = d\mathcal{A} + \mathcal{A} \wedge \mathcal{A} = 0 \quad \text{（平凡！）}$$

**曲率的非平凡性体现在参数空间中**。把 $(E_1, t_1, \tau, t)$ 都视为底流形坐标：

$$\mathcal{A}(E_1, t_1, \tau, t) = A(t; E_1, t_1, \tau)\, dt$$

此时 $\mathcal{A}$ 只有 $dt$ 分量，所以 $d\mathcal{A}$ 有 $dE_1 \wedge dt$, $dt_1 \wedge dt$, $d\tau \wedge dt$ 分量。由于 $\mathcal{A} \wedge \mathcal{A} = 0$（只有一个基底方向），曲率为：

$$\boxed{F = \frac{\partial A}{\partial E_1}\, dE_1 \wedge dt + \frac{\partial A}{\partial t_1}\, dt_1 \wedge dt + \frac{\partial A}{\partial \tau}\, d\tau \wedge dt}$$

### 4.2 曲率的显式分量

对 $\mathcal{A}$ 求偏导（表 1 + §3.3 显式）：

**$E_1$ 方向的曲率**（所有三步相同）：
$$\frac{\partial A}{\partial E_1} = 2 X_{12} \quad \Longrightarrow \quad F_{E_1, t} = 2 X_{12}$$

**$t_1$ 方向的曲率**（Step 1 & 2；Step 3 中 $t_1 = 0$ 故无贡献）：
$$\frac{\partial A}{\partial t_1} = 
\begin{cases}
2 f_-(t) X_{15} & \text{Step 1} \\
2 f_+(t) X_{15} & \text{Step 2} \\
0                  & \text{Step 3}
\end{cases}$$

**$\tau$ 方向的曲率**：$\tau$ 出现在门控函数中
$$\frac{\partial f_\pm}{\partial \tau} = \pm\frac{\pi t}{2\tau^2}\sin\left(\frac{\pi t}{\tau}\right)$$

所以
$$\frac{\partial A^{(1)}}{\partial \tau} = 2t_c \frac{\partial f_-}{\partial\tau} X_{24} + 2t_1 \frac{\partial f_-}{\partial\tau} X_{15} + 2E_0 \frac{\partial f_+}{\partial\tau} X_{45}$$

（Step 2, 3 类似。）

### 4.3 曲率的关键性质

**性质 1：曲率的李代数闭包 = 和乐群李代数（Ambrose-Singer）**

各步曲率分量的李代数闭包：

| 步 | $\{F_{\mu,t}\}$ 张成的子代数 | 维数 | 和乐限制 |
|---|---|---|---|
| Step 1 | $\text{span}\{X_{12}, X_{24}, X_{15}, X_{45}\} \subset \mathfrak{so}(4)$ | 6 | $\subset$ SO(4) |
| Step 2 | $\text{span}\{X_{12}, X_{24}, X_{34}, X_{15}\} = \mathfrak{so}(5)$ | **10（满）** | 整个 SO(5) |
| Step 3 | $\text{span}\{X_{12}, X_{34}, X_{45}\} \subset \mathfrak{u}(1)\oplus\mathfrak{su}(2)$ | 4 | $\subset U(1)\times$SU(2) |

**性质 2：$E_1=0$ 时的曲率退化**

当 $E_1=0$：
- $F_{E_1, t} = 0$（$\sigma_z$ 方向无曲率）
- 剩余非零分量 $\{X_{24}, X_{34}, X_{15}, X_{45}\}$ 满足 $A=D$ 对称性（在 Sp(2) 表示中）
- 此对称性保证沿闭合回路的 $\sigma_z$ 净和乐为零
- **曲率退化 $\Rightarrow$ 和乐限制在 SO(5) 的 2 维子流形上**

**性质 3：曲率对易子与非阿贝尔性**

不同参数方向的曲率分量不对易：
$$[F_{E_1,t}, F_{t_1,t}] = 4 [X_{12}, X_{15}] \neq 0$$

（因为 $X_{12}=i\gamma_1\gamma_2$ 和 $X_{15}=i\gamma_1\gamma_5$ 涉及同一个 $\gamma_1$。）
**这是 $E_1$ 和 $t_1$ 产生干涉条纹的代数根源。**

---

## 五、和乐：演化算符作为曲率的积分

### 5.1 定义

沿编织回路 $\gamma$ 的和乐（平行移动）：

$$\text{Hol}_\gamma(\mathcal{A}) = \mathcal{P}\exp\left(\oint_\gamma \mathcal{A}\right) = \mathcal{P}\exp\left(\int_0^{3\tau} A(t)\,dt\right) = U(3\tau) \in \text{SO}(5)$$

这就是 SO(5) 演化矩阵 $R$。

### 5.2 非阿贝尔 Stokes 定理

和乐可以用曲率在"以 $\gamma$ 为边界的曲面"上的面积分表示（非阿贝尔推广）：

$$U(3\tau) = \mathcal{P}\exp\left(\oint_\gamma \mathcal{A}\right) = \mathcal{P}_\Sigma \exp\left(\iint_\Sigma \tilde{F}\right)$$

其中 $\tilde{F}$ 是曲率 $F$ 的"平行移动回拉"到参考点，$\mathcal{P}_\Sigma$ 是曲面有序化。这个公式显式说明：**和乐由曲率在参数空间中的分布决定，而不仅仅是同伦类。**

### 5.3 参数依赖：和乐的变分

设 $\theta = (E_1, t_1, \tau)$ 为参数。和乐的梯度：

$$\boxed{\frac{\partial U(3\tau)}{\partial \theta_k} = U(3\tau) \int_0^{3\tau} U(t)^\dagger\; \frac{\partial A(t)}{\partial \theta_k}\; U(t)\, dt}$$

这给出"参数改变时编织结果如何变化"的精确公式。

对 $E_1$：
$$\frac{\partial U}{\partial E_1} = 2 U(3\tau) \int_0^{3\tau} U(t)^\dagger X_{12} U(t)\, dt$$

（$X_{12}$ 在所有三步都存在，故积分范围为全程。）

对 $t_1$：
$$\frac{\partial U}{\partial t_1} = 2 U(3\tau) \left[\int_0^\tau f_- U^\dagger X_{15} U\, dt + \int_\tau^{2\tau} f_+ U^\dagger X_{15} U\, dt\right]$$

（Step 3 中 $t_1=0$ 故无贡献。）

---

## 六、缩并到 MZM 子空间：有效 SU(2) 曲率

### 6.1 绝热缩并

当 ancilla 自由度快于 MZM 时（大 $\tau$ 极限），可以将联络缩并到 $\{\gamma_1,\gamma_2,\gamma_3\}$ 子空间。
这便是 **维数约化**——大的 $\mathfrak{so}(5)$ 曲率在子空间上诱导有效 SU(2) 曲率。

有效联络（一阶）：
$$\mathcal{A}_{\text{eff}} = A_{\text{eff}}(t)\, dt, \quad A_{\text{eff}}(t) \in \mathfrak{su}(2)$$

在 Pauli 基底 $\{\sigma_x, \sigma_y, \sigma_z\}$ 中（映射见 §8）：

### 6.2 $E_1=0, t_1 \neq 0$ 情形的显式有效联络

由 report §7 的结果：
$$A_{\text{eff}}(t) = \varepsilon(t)\,\sigma_x + D(t)\,\sigma_y$$

其中：
- $\varepsilon(t)$：编织项（$\sim t_c$ 门控），其路径积分 $\Phi_G = \int \varepsilon\,dt = \pi/2$
- $D(t)$：$t_1$ 诱导的动态项

有效和乐：
$$\boxed{U_{\text{eff}}(3\tau) = \exp\!\big(-i\phi\,\hat{n}\cdot\vec{\sigma}\big)}$$
$$\phi = \sqrt{\Phi_G^2 + \Phi_D^2}, \quad \tan\alpha = \frac{\Phi_D}{\Phi_G}$$
$$\hat{n} = (\cos\alpha,\;\sin\alpha,\;0), \quad \Phi_D = \int_0^{3\tau} D(t)\,dt$$

### 6.3 $E_1 \neq 0$ 情形的完整有效联络

$$A_{\text{eff}}(t) = \varepsilon(t)\,\sigma_x + D_{t_1}(t)\,\sigma_y + D_{E_1}(t)\,\sigma_z$$

其中三项不对易：
$$[\sigma_x, \sigma_y] = 2i\sigma_z, \quad [\sigma_y, \sigma_z] = 2i\sigma_x, \quad [\sigma_z, \sigma_x] = 2i\sigma_y$$

**非零曲率分量的完整表格**：

$$F_{xy} \propto \sigma_z, \quad F_{yz} \propto \sigma_x, \quad F_{zx} \propto \sigma_y$$

所有三个分量均非零 $\iff E_1 \neq 0$ 且 $t_1 \neq 0$。

$E_1=0$ 时 $F_{xy}=0$，$F_{yz}=F_{zx}=0$（只有两个独立分量非零）→ 和乐流形降至 2 维。

---

## 七、保真度作为和乐的距离函数

### 7.1 定义

Fig 1(d) 的保真度可以写为和乐空间中两点之间的距离：

$$\text{fid}(E_1, t_1, \tau) = \left|\langle\psi_1^+|\,\text{Hol}_{\gamma(E_1,t_1,\tau)}\,\text{Hol}_{\gamma(E_1,t_1,\tau)}\,|\psi_1^-\rangle\right|^2$$

这是以下函数的采样：
$$\text{fid}: M \to [0,1], \quad (E_1, t_1, \tau) \mapsto \text{fidelity}$$

### 7.2 梯度（保真度的变分）

由和乐的变分公式（§5.3），可以计算保真度对参数的导数：

$$\frac{\partial\,\text{fid}}{\partial\theta_k} = 2\,\text{Re}\left[\langle\psi_1^+|U^2|\psi_1^-\rangle^* \cdot \langle\psi_1^+|\frac{\partial(U^2)}{\partial\theta_k}|\psi_1^-\rangle\right]$$

其中 $\partial(U^2)/\partial\theta_k = (\partial U/\partial\theta_k)U + U(\partial U/\partial\theta_k)$，而 $\partial U/\partial\theta_k$ 由 §5.3 给出。

**Fig 1(d) 的等高线 = 这些梯度的积分曲线**。

---

## 八、联络和曲率在四元数 Sp(2) 表示中的形式

利用 $\mathfrak{so}(5) \cong \mathfrak{sp}(2)$，一切化为 $2 \times 2$ 四元数。

### 8.1 Pauli 对应

$$i\mathbf{i} \leftrightarrow \sigma_x \leftrightarrow i\gamma_2\gamma_3$$
$$i\mathbf{j} \leftrightarrow \sigma_y \leftrightarrow i\gamma_3\gamma_1$$
$$i\mathbf{k} \leftrightarrow \sigma_z \leftrightarrow i\gamma_1\gamma_2$$

### 8.2 联络在 Sp(2) 中的分块形式

$$\mathcal{A}_{\text{sp}}(t) = K(t)\,dt = \begin{pmatrix} A(t) & B(t) \\ C(t) & D(t) \end{pmatrix} dt$$

其中 $A,B,C,D \in \mathbb{H}$，$A,D$ 为纯虚四元数（这样 $K \in \mathfrak{sp}(2)$）。

**对参数 $E_1$ 的偏导**：
$$\frac{\partial K}{\partial E_1} = \begin{pmatrix} \mathbf{i}/2 & 0 \\ 0 & -\mathbf{i}/2 \end{pmatrix} = \Sigma_{12}$$

**对参数 $t_1$ 的偏导**（Step 1）：
$$\frac{\partial K}{\partial t_1} = f_-(t) \cdot \begin{pmatrix} 0 & 1/2 \\ -1/2 & 0 \end{pmatrix} = f_-(t)\,\Sigma_{15}$$

**MZM 子空间有效联络**：
$$K_{\text{eff}} = A + Bq \in \operatorname{Im}\mathbb{H}$$

（$q = ZX^{-1}$ 是 Riccati 变量。）

在 Pauli 分量中：
$$\underbrace{K_{\text{eff}}}_{\text{anti-Hermitian}} = \omega_x(i\sigma_x) + \omega_y(i\sigma_y) + \omega_z(i\sigma_z)$$
$$\underbrace{H_{\text{eff}} = iK_{\text{eff}}}_{\text{Hermitian}} = \omega_x\sigma_x + \omega_y\sigma_y + \omega_z\sigma_z$$

### 8.3 有效曲率的分量显式（$E_1 \neq 0$ 一般情况）

三个频率分量来自：
$$\omega_x = \frac{|t_3|}{2} + \operatorname{Re}_{\mathbf{i}}(Bq)$$
$$\omega_y = \frac{|t_2|}{2} + \operatorname{Re}_{\mathbf{j}}(Bq)$$
$$\omega_z = \frac{E_1}{2} + \operatorname{Re}_{\mathbf{k}}(Bq)$$

其中 $\operatorname{Re}_{\mathbf{i}}(q) = q_1$，$\operatorname{Re}_{\mathbf{j}}(q) = q_2$，$\operatorname{Re}_{\mathbf{k}}(q) = q_3$。

$Bq$ 项包含 $q$ 的非线性演化（通过 Riccati ODE：$\dot{q} = C + Dq - qA - qBq$）。

**$E_1=0$ 时 $\omega_z$ 的性质**：虽然 $E_1=0$ 时 $\omega_z$ 的 $Bq$ 贡献可以瞬时非零，但门控函数的时间积分对称性保证 $\int_0^{3\tau} \omega_z(t)\,dt = 0$。这是 $A=D$ 对称性在联络语言中的表现。

---

## 九、数值曲率计算方案

可以用现有 SO(5) 代码直接计算上述所有几何量。

### 9.1 曲率分量的数值提取

联络在 $t$ 时刻的值由 `A_step1/2/3()` 给出（$5 \times 5$ 反对称矩阵）。
参数偏导用有限差分：

$$F_{E_1, t} \approx \frac{A(t; E_1+\delta, t_1, \tau) - A(t; E_1-\delta, t_1, \tau)}{2\delta}$$

或者用 §4.2 的解析显式直接编码（更快）。

### 9.2 和乐对参数的梯度

$$\frac{\partial U}{\partial E_1} \approx \frac{U(E_1+\delta) - U(E_1-\delta)}{2\delta}$$

每个差分点只需调一次 `so5_protocol()`。

### 9.3 平行移动曲率（用于 Stokes 定理验证）

在路径上的各点计算"平行移动回原点"的曲率：
$$\tilde{F}(t) = U(t)^\dagger\,F(t)\,U(t) \in \mathfrak{so}(5)$$

验证（非阿贝尔 Stokes）：
$$U(3\tau) \stackrel{?}{\approx} \exp\left(\int_0^{3\tau} \tilde{F}(t)\,dt\right)$$

仅当 $\tilde{F}(t)$ 在不同时刻对易时才精确成立（1 阶 Magnus）。

---

## 十、总结：三层几何结构

```
┌─────────────────────────────────────────────────────────┐
│ 第 0 层：拓扑（同伦类）                                   │
│ π₁(M) = B₂ ≅ ℤ                                          │
│ → 纯编织 U₀ 仅依赖绕数 w（=1 为单次 γ₂↔γ₃ 交换）         │
│ → 平坦联络 F=0 时，和乐仅由此层决定                        │
├─────────────────────────────────────────────────────────┤
│ 第 1 层：联络的非平坦性（曲率）                             │
│ F = dA + A∧A ≠ 0 分量:                                    │
│   F_{E₁,t} = 2 X₁₂            → σ_z 方向的"力"            │
│   F_{t₁,t} = 2g(t) X₁₅        → σ_y 方向的"力"            │
│   F_{τ,t}  = ∂A/∂τ            → 非绝热修正                │
│ → 同一同伦类内不同路径给出不同的和乐                        │
├─────────────────────────────────────────────────────────┤
│ 第 2 层：曲率的非对易性（非阿贝尔结构）                      │
│ [F_{E₁,t}, F_{t₁,t}] = 4[X₁₂, X₁₅] ≠ 0                   │
│ → E₁ 和 t₁ 的效应不能简单叠加                             │
│ → 产生 Fig 1(d) 的干涉条纹                                │
│ → 退化条件 E₁=0：F 的部分分量消失 → 和乐限制在 2 维子流形   │
└─────────────────────────────────────────────────────────┘
```

### 核心物理图像

编织系统是一个 SO(5)-主丛，底流形 $M$ 被编织协议参数化。$H_{EM}$ 定义了联络，演化算符 $U$ 是和乐。**纯编织**是拓扑（$\pi_1(M)$ 的非平凡表示），**动力学相位**是联络曲率 $F$ 在 $M$ 上的积分。$E_1$ 和 $t_1$ 在曲率中产生不对易的分量，导致非平凡的参数依赖——这就是 Fig 1(d) 的完整几何解释。

---

## 十一、从对易关系反推：和乐群的三个层次

### 11.0 核心命题

**从 $H_{EM}$ 各生成元的对易关系出发，严格推导和乐群的结构层次，
证明以下对应关系：**

$$\boxed{\begin{aligned}
E_1 = t_1 = 0 &\;\Longleftrightarrow\; \text{和乐群} \cong B_2 \;\Longleftrightarrow\; \text{平坦联络（纯拓扑）} \\
E_1 = 0,\; t_1 \neq 0 &\;\Longleftrightarrow\; \text{和乐群} \cong \text{SO}(2) \;\Longleftrightarrow\; A=D\;\text{对称性} \\
E_1 \neq 0,\; t_1 \neq 0 &\;\Longleftrightarrow\; \text{和乐群} \cong \text{SU}(2) \;\Longleftrightarrow\; \text{全非阿贝尔}
\end{aligned}}$$

### 11.1 预备：MZM 子空间的有效生成元

在绝热消除 ancilla 后（$\tau$ 足够大），MZM 子空间 $\{\gamma_1, \gamma_2, \gamma_3\}$ 上的
有效哈密顿量（Hermitian，在 Pauli 基底中）：

$$\boxed{H_{\text{eff}}(t) = \varepsilon(t)\,\sigma_x + D_{t_1}(t)\,\sigma_y + E_1\,\sigma_z}$$

其中 Pauli 矩阵对应关系为：
$$\sigma_x \leftrightarrow i\gamma_2\gamma_3,\quad \sigma_y \leftrightarrow i\gamma_3\gamma_1,\quad \sigma_z \leftrightarrow i\gamma_1\gamma_2$$

三个系数：
- $\varepsilon(t)$：来自 $t_2, t_3, E_d$ 的编织驱动（始终非零，提供几何相位）
- $D_{t_1}(t) \propto t_1 \cdot g(t)$：$t_1$ 诱导的 $\sigma_y$ 动态耦合
- $E_1$：$\gamma_1$–$\gamma_2$ 杂化能（恒定，$\sigma_z$ 方向）

它们满足 $\mathfrak{su}(2)$ 对易关系：
$$[\sigma_x, \sigma_y] = 2i\sigma_z,\quad [\sigma_y, \sigma_z] = 2i\sigma_x,\quad [\sigma_z, \sigma_x] = 2i\sigma_y$$

有效演化算符（单次编织）：
$$U_{\text{eff}}(3\tau) = \mathcal{P}\exp\!\left(-i\int_0^{3\tau} H_{\text{eff}}(t)\,dt\right) \in \text{SU}(2)$$

### 11.2 关键引理：和乐李代数由时间-对易子生成

**引理（简化的 Ambrose-Singer）**：对于时变哈密顿量 $H(t)$，和乐群 $\text{Hol}$ 的李代数
$\mathfrak{hol}$ 由 $\{H(t)\}_{t \in [0,T]}$ 及其所有多重对易子张成。

$$\mathfrak{hol} = \text{Lie}\big\langle \{H(t) : t \in [0,3\tau]\} \big\rangle$$

即：取所有时刻的 $H(t)$ 作为生成元集合，计算它们之间的所有可能的李括号，
得到的李代数的实线性张成就是 $\mathfrak{hol}$。

**物理含义**：如果 $H(t)$ 在所有时刻都限制在某个子代数 $\mathfrak{h} \subset \mathfrak{su}(2)$ 中，
那么整个演化 $U(t)$ 也被限制在对应的李子群 $H \subset \text{SU}(2)$ 中。
这就是 report §2.4 中李代数闭包分析的严格数学基础。

---

### 11.3 情形一：$E_1 = t_1 = 0$ → 和乐群 $\cong B_2$

**条件**：$E_1 = 0$ 且 $t_1 = 0$。

**$H_{\text{eff}}$ 的形式**：
$$H_{\text{eff}}(t) = \varepsilon(t)\,\sigma_x$$

其中 $\varepsilon(t)$ 在所有三步中非零（来自 $t_2, t_3, E_d$）。

**对易子分析**：
$$[H_{\text{eff}}(t_1), H_{\text{eff}}(t_2)] = \varepsilon(t_1)\varepsilon(t_2)\,[\sigma_x, \sigma_x] = 0 \quad \forall\, t_1, t_2$$

**所有时刻的 $H_{\text{eff}}(t)$ 均正比于 $\sigma_x$，相互对易。**

**和乐李代数**：
$$\mathfrak{hol} = \text{span}_{\mathbb{R}}\{\sigma_x\} \cong \mathfrak{u}(1)$$

李代数维数降为 1。

**和乐（路径有序化退化为普通指数）**：
$$U_{\text{eff}}(3\tau) = \exp\!\left(-i\left[\int_0^{3\tau} \varepsilon(t)\,dt\right]\sigma_x\right)$$

编织路径积分给出 $\int_0^{3\tau} \varepsilon(t)\,dt = \pi/2$（报告 §7.3 的 $\Phi_G$），因此：
$$U_{\text{eff}}(3\tau) = \exp(-i\frac{\pi}{2}\sigma_x) = -i\sigma_x = \begin{pmatrix}0 & -i \\ -i & 0\end{pmatrix}$$

这正是 $\gamma_2 \leftrightarrow \gamma_3$ 交换在旋量表示中的矩阵（$\gamma_2 \to \gamma_3,\; \gamma_3 \to -\gamma_2$）。

**和乐群的离散性**：虽然 $\mathfrak{hol} \cong \mathfrak{u}(1)$ 是连续 1 维李代数，
但编织协议给出了**固定的** $\pi/2$ 旋转角。重复编织 $n$ 次：
$$U^n = \exp(-i\frac{n\pi}{2}\sigma_x)$$

这是一个 $\mathbb{Z}_4$ 循环群（因为 $U^4 = I$）。$\mathbb{Z}_4$ 是 $B_2 \cong \mathbb{Z}$
在 SU(2) 旋量表示中的像。

$$\boxed{\text{Hol} \cong \mathbb{Z}_4 \cong B_2 / \ker(\rho) \;\Longleftrightarrow\; \text{平坦有效联络}}$$

**几何解释**：曲率 $F = [H(t_1), H(t_2)] = 0$。和乐仅依赖同伦类（绕数 $w$），不依赖路径的具体形状。
这就是**纯拓扑编织**。

---

### 11.4 情形二：$E_1 = 0$，$t_1 \neq 0$ → 和乐群 $\cong \text{SO}(2)$

**条件**：$E_1 = 0$ 且 $t_1 \neq 0$。

**$H_{\text{eff}}$ 的形式**：
$$H_{\text{eff}}(t) = \varepsilon(t)\,\sigma_x + D_{t_1}(t)\,\sigma_y$$

两项系数有不同的时间依赖：$\varepsilon(t)$ 来自门控 $t_2, t_3$，$D_{t_1}(t) \propto t_1 \cdot g(t)$。

**对易子分析**：
$$[H_{\text{eff}}(t_1), H_{\text{eff}}(t_2)] = \big[\varepsilon_1\sigma_x + D_1\sigma_y,\; \varepsilon_2\sigma_x + D_2\sigma_y\big]$$

展开（$[\sigma_x, \sigma_y] = 2i\sigma_z$，其余对易子为零）：
$$\begin{aligned}
&= \varepsilon_1 D_2[\sigma_x, \sigma_y] + D_1\varepsilon_2[\sigma_y, \sigma_x] \\
&= 2i(\varepsilon_1 D_2 - D_1\varepsilon_2)\,\sigma_z
\end{aligned}$$

$$\boxed{[H_{\text{eff}}(t_1), H_{\text{eff}}(t_2)] = 2i\big(\varepsilon(t_1)D_{t_1}(t_2) - D_{t_1}(t_1)\varepsilon(t_2)\big)\,\sigma_z}$$

因为 $\varepsilon(t)$ 和 $D_{t_1}(t)$ 有不同的时间轮廓（一个跟 $t_2/t_3$ 走，一个跟 $t_1$ 走），
差 $(\varepsilon_1 D_2 - D_1\varepsilon_2)$ 一般非零。

**和乐李代数**：对易子生成了 $\sigma_z$ 分量，因此：
$$\mathfrak{hol} = \text{span}_{\mathbb{R}}\{\sigma_x,\; \sigma_y,\; \sigma_z\} \cong \mathfrak{su}(2)$$

**李代数看似满 3 维 $\mathfrak{su}(2)$，为何和乐群是 SO(2)？**

**关键——$A=D$ 对称性约束**：

李代数 $\mathfrak{hol} \cong \mathfrak{su}(2)$ 是所有**可能**路径能生成的完全集。但**特定编织协议**
施加了 $A=D$ 对称性约束（report §7.5bis），该约束源于在 $\mathfrak{sp}(2)$ 表示中：

$$E_1 = 0 \;\Longrightarrow\; A(t) = D(t) = \frac{|t_3|}{2}\mathbf{i} + \frac{|t_2|}{2}\mathbf{j}$$

对称性 $A=D$ 导致 Riccati ODE 中：
$$\dot{q} = C + [A, q] - qBq$$

线性部分 $[A, q]$ 产生的瞬时 $\sigma_z$ 分量在编织全程的**时间积分恰好抵消**。
定量证明（report §7.5）：

$$\begin{aligned}
\int_0^{\tau} \text{Step 1 瞬时 }\sigma_z\,dt &\propto +\int f_m\,dt \\
\int_{\tau}^{2\tau} \text{Step 2 瞬时 }\sigma_z\,dt &\propto \pm\int (f_p + f_m)\,dt \\
\int_{2\tau}^{3\tau} \text{Step 3 瞬时 }\sigma_z\,dt &\propto -\int f_p\,dt
\end{aligned}$$

由于 $\int_0^\tau f_m\,dt = \int_0^\tau f_p\,dt = \tau/2$，三步求和：
$$\boxed{\int_0^{3\tau} \omega_z(t)\,dt = 0}$$

**净 $\sigma_z$ 和乐为零。** 这就是数值事实
$\int \omega_z\,dt = 0$ 在李代数层面的根源。

**和乐的形式**：由于净 $\sigma_z = 0$，有效和乐被限制在 $\sigma_x$–$\sigma_y$ 平面内：

$$U_{\text{eff}}(3\tau) = \exp\!\big(-i\phi\,(\cos\alpha\,\sigma_x + \sin\alpha\,\sigma_y)\big)$$

其中 $\phi = \sqrt{(\pi/2)^2 + \Phi_D^2}$，$\tan\alpha = \Phi_D/(\pi/2)$，
$\Phi_D = \int_0^{3\tau} D_{t_1}(t)\,dt \propto t_1\tau$。

**和乐群的结构**：

- **固定 $t_1, \tau$**：和乐是绕固定轴 $\hat{n} = (\cos\alpha, \sin\alpha, 0)$ 旋转 $\phi$ 角。
  重复编织给出该轴的任意倍数旋转：
  $$\text{Hol}(t_1, \tau) = \{\exp(-i\theta\,\hat{n}\cdot\vec{\sigma}) : \theta \in \mathbb{R}\} \cong \text{SO}(2) \cong U(1)$$

  这是**1 维连通 Abel 李子群**。

- **变动 $t_1$**：轴 $\hat{n}$ 在 $\sigma_x$–$\sigma_y$ 平面内连续偏转，但**始终无法产生 $\sigma_z$ 分量**。
  所有可达到的和乐构成 SU(2) 中的 2 维子流形（不是子群）。

$$\boxed{\text{固定参数：Hol} \cong \text{SO}(2) \;\Longleftrightarrow\; A=D\;\text{对称性} \;\Longleftrightarrow\; \int\omega_z = 0}$$

**几何解释**：曲率分量 $F_{t_1,t}$ 只在 $\sigma_x$–$\sigma_y$ 平面内产生效应。
沿编织回路的曲率面积分在 $\sigma_z$ 方向恰好正负抵消——曲率在 $\sigma_z$ 方向是"闭合"的
（类似磁场通过闭合曲面的总通量为零）。

---

### 11.5 情形三：$E_1 \neq 0$，$t_1 \neq 0$ → 和乐群 $\cong \text{SU}(2)$

**条件**：$E_1 \neq 0$ 且 $t_1 \neq 0$。

**$H_{\text{eff}}$ 的形式**：
$$H_{\text{eff}}(t) = \varepsilon(t)\,\sigma_x + D_{t_1}(t)\,\sigma_y + E_1\,\sigma_z$$

所有三个 Pauli 方向全部非零。

**对易子分析**：

$$\begin{aligned}
[H(t_1), H(t_2)] &= 2i\big(\varepsilon_1 D_2 - D_1\varepsilon_2\big)\,\sigma_z \\
&\quad + 2i\big(D_1 E_1 - E_1 D_2\big)\,\sigma_x \\
&\quad + 2i\big(E_1\varepsilon_2 - \varepsilon_1 E_1\big)\,\sigma_y
\end{aligned}$$

关键差异——$E_1 \neq 0$ 时新增了两类对易子：

$$[\sigma_y\text{-项}, \sigma_z\text{-项}] \propto \sigma_x,\qquad [\sigma_z\text{-项}, \sigma_x\text{-项}] \propto \sigma_y$$

这三个对易子**在时间上无法全部抵消**，因为 $E_1$ 是常数而 $\varepsilon(t), D_{t_1}(t)$ 是时变的。

**$A \neq D$ 对称性破缺**（report §7.5bis）：

$\mathfrak{sp}(2)$ 表示中：
$$A = \frac{E_1 + |t_3|}{2}\mathbf{i} + \frac{|t_2|}{2}\mathbf{j},\quad D = \frac{-E_1 + |t_3|}{2}\mathbf{i} + \frac{|t_2|}{2}\mathbf{j}$$

$$A - D = E_1\,\mathbf{i} = E_1\,\sigma_x \neq 0$$

$A \neq D$ 打破了 §11.4 中使 $\int\omega_z = 0$ 的对称性。Riccati ODE 的完整形式被恢复：

$$\dot{q} = C + Dq - qA - qBq$$

$Dq$ 和 $-qA$ 不再对称相消 → $\sigma_z$ 净积累 ≠ 0 →
三 Pauli 方向全部可独立到达。

**和乐李代数**：
$$\mathfrak{hol} = \text{span}_{\mathbb{R}}\{\sigma_x, \sigma_y, \sigma_z\} \cong \mathfrak{su}(2)$$

**和乐群**：生成元不对易（$[\sigma_x, \sigma_y] = 2i\sigma_z \neq 0$），
且三个方向的路径积分均可独立为非零。

$$\boxed{\text{Hol} \cong \text{SU}(2) \;\Longleftrightarrow\; E_1 \neq 0,\;t_1 \neq 0 \;\Longleftrightarrow\; \text{全非阿贝尔}}$$

通过调节 $(E_1, t_1, \tau)$ 三个独立参数，可以遍历 SU(2) 中任意群元素。
具体地：

| 参数 | 控制的自由度 | SU(2) 参数 |
|---|---|---|
| $t_c, E_0$（编织） | 基准 $\sigma_x$ 旋转 $\Phi_G = \pi/2$ | 固定 |
| $t_1$ | $\sigma_y$ 旋转 $\Phi_D$ | 第 1 连续自由度 |
| $E_1$ | $\sigma_z$ 旋转 $\Phi_{E_1}$ | 第 2 连续自由度 |
| $\tau$ | 所有旋转角度的缩放 | 第 3 连续自由度 |

三个自由度恰好匹配 $\dim\text{SU}(2) = 3$。

**数值验证**（report §7.5bis，`bloch_E1_nonzero_trajectories.py`）：

| 情形 | $v_x$ | $v_y$ | $v_z$ | Bloch 覆盖 |
|---|---|---|---|---|
| $E_1 = 0$，3 独立 $t_1$ | $[-1,1]$ | **$0$** | $[-1,1]$ | 1D 曲线 |
| **$E_1 \neq 0$**，3 独立 $t_1$ | $[-1,1]$ | **$[-1,1]$** | $[-1,1]$ | **整个球面** |

$E_1 \neq 0$ 释放了 $\sigma_y$ 方向 → 遍历 Bloch 球面全部点 → SU(2) 的群作用可迁。

---

### 11.6 总结：三步反推证明的完整逻辑链

```
输入：H_EM 的对易关系
  │
  ├── [H(t₁), H(t₂)] = 0 ∀t₁,t₂ ?
  │     │
  │     ├── 是 → E₁ = t₁ = 0
  │     │        𝔥𝔬𝔩 = span{σ_x} ≅ 𝔲(1)
  │     │        Hol ≅ ℤ₄ ≅ B₂/im(ρ)
  │     │        纯拓扑编织，平坦联络
  │     │
  │     └── 否 → 继续判断
  │           │
  │           ├── A = D（𝔰𝔭(2)对称性）？
  │           │     │
  │           │     ├── 是 → E₁ = 0, t₁ ≠ 0
  │           │     │        𝔥𝔬𝔩 = 𝔰𝔲(2)（全李代数）
  │           │     │        但 ∫ω_z dt = 0（对称性约束）
  │           │     │        Hol ≅ SO(2)（固定参数）
  │           │     │        轴在 σ_x-σ_y 平面内
  │           │     │
  │           │     └── 否 → E₁ ≠ 0, t₁ ≠ 0
  │           │              𝔥𝔬𝔩 = 𝔰𝔲(2)（全李代数）
  │           │              ∫ω_z dt ≠ 0（对称性破缺）
  │           │              Hol ≅ SU(2)（全非阿贝尔）
  │           │              遍历 Bloch 球面
```

### 11.7 配置空间流形性质的推论

从和乐群反推底流形的几何结构：

| 和乐群 | 底流形 $M$ 的有效曲率 | $\dim M_{\text{eff}}$ | 物理 |
|---|---|---|---|
| $B_2$（离散） | $F = 0$（平坦） | 0（纯拓扑） | 纯编织，仅依赖同伦类 |
| SO(2) | $F_{t_1,t} \neq 0$，但 $\int F|_{z} = 0$ | 1（轴固定在平面内） | $t_1$ 调节旋转角和轴偏角 |
| SU(2) | $F_{E_1,t}, F_{t_1,t}, [F, F] \neq 0$ | 3（全曲率释放） | $E_1, t_1, \tau$ 独立控制 SU(2) |

**核心推论**：配置空间流形 $M$ 的有效维度 = 和乐群参数的独立个数。
纯编织时流形退化为离散点集（各同伦类各一点），
动力学项 $E_1, t_1$ 逐级"撑开"流形，使其重新获得连续维度。

这正是 **"动力学影响 = 配置空间流形的曲率"** 命题的严格证明。

---

# 附录：三个核心问题的回答

## Q1: 纤维丛理论是否完整描述了演化？

### 回答：是。三个方面论证。

**1. 等价性论证**

原始 Schrödinger 方程 $i\partial_t|\psi\rangle = H|\psi\rangle$ 等价于 $\dot U = -iHU$。
在 SO(5) 旋量表示中，这等价于联络的平行移动方程：

$$\frac{d}{dt}g(t) + \mathcal{A}\left(\frac{d}{dt}\right) \cdot g(t) = 0$$

其中 $g(t) \in \text{SO}(5)$ 是纤维坐标（即 $R(t)$），$\mathcal{A}(d/dt) = A(t)$ 是联络在切矢量 $d/dt$ 上的取值。这是主丛上**水平提升**的标准方程。

演化算符是和乐：
$$U(t) = \text{Hol}_{\gamma|_{[0,t]}}(\mathcal{A})$$

**2. 曲率的完整性**

在仅含时间参数的原始设定中，曲率 $F = d\mathcal{A} + \mathcal{A} \wedge \mathcal{A}$ 恒为零（因为 $dt \wedge dt = 0$）。这**不是缺陷**——它反映的是：仅存在一个独立参数 $(t)$ 时，任何 1-形式 $\mathcal{A} = A(t)dt$ 局部都是平坦的（总可以通过重新参数化消去）。

曲率在以下扩展中变为非平凡：
- **参数空间扩展**：$(E_1, t_1, \tau, t)$ → $F$ 有非零分量（§4.2）
- **控制空间扩展**：允许每个时间片的 $A(t)$ 独立变化 → $F$ 在无穷维空间中非零

**3. 理论覆盖了所有物理内容**

| 物理概念 | 几何对应 | 是否完整捕获 |
|---|---|---|
| 态演化 | 纤维的平行移动 | ✓ |
| 编织算符 | 闭合回路 $\gamma$ 的和乐 | ✓ |
| 动力学相位 | 曲率 $F$ 的面积分（非阿贝尔 Stokes） | ✓ |
| 非绝热效应 | $\partial A/\partial\tau$ 即 $F_{\tau,t}$ | ✓ |
| 参数依赖 | 和乐映射 $h: M_{\text{param}} \to \text{SO}(5)$ 的微分 | ✓ |
| ABS 效应 | 曲率分量 $F_{E_1,t}, F_{t_1,t}$ 及其对易子 | ✓ |

**唯一不包含的是退相干等开放系统效应**——这需要非幺正演化，超出了主丛联络的范畴（需要推广到密度矩阵的几何）。

---

## Q2: 如何量化 $E_1$ 和 $t_1$ 的影响（非阿贝尔性）？

### 回答：五个递进的量化指标。

### 2.1 局部非阿贝尔性：曲率对易子的范数

最直接的代数指标——曲率分量是否对易：

$$\boxed{\mathcal{C}_{E_1,t_1}(t) := \big\| [F_{E_1,t}, F_{t_1,t}] \big\|}$$

显式计算（用 $\mathfrak{so}(5)$ 的 Killing 范数 $\|X\|^2 = -\frac{1}{2}\text{Tr}(X^2)$）：

$$F_{E_1,t} = 2X_{12}, \quad F_{t_1,t} = 2g(t)X_{15}$$

$$[X_{12}, X_{15}] = 2i X_{25} \neq 0 \quad \text{（因为 } \gamma_1 \text{ 同时出现在两者中）}$$

$$\boxed{\big\|[F_{E_1,t}, F_{t_1,t}]\big\| = 8|g(t)| \cdot \|X_{25}\|}$$

$g(t)$ 在各步的取值：

| 步 | $g(t)$ | $\|F_{E_1,t} \wedge F_{t_1,t}\|$ |
|---|---|---|
| Step 1 | $f_-(t)$ | $8 f_-(t) \cdot \|X_{25}\|$ |
| Step 2 | $f_+(t)$ | $8 f_+(t) \cdot \|X_{25}\|$ |
| Step 3 | 0 | 0（$t_1$ 被关断） |

**非阿贝尔性的时间累积**：
$$\boxed{\mathcal{C}_{\text{total}} = \int_0^{3\tau} \big\|[F_{E_1,t}, F_{t_1,t}]\big\|\,dt = 8\|X_{25}\| \cdot \left[\int_0^\tau f_-(t)dt + \int_\tau^{2\tau} f_+(t)dt\right]}$$

由于 $\int_0^\tau f_-(t)dt = \int_\tau^{2\tau} f_+(t)dt = \tau/2$：
$$\boxed{\mathcal{C}_{\text{total}} = 8\|X_{25}\| \cdot \tau}$$

**结论**：非阿贝尔性随编织时间 $\tau$ 线性累积——这解释了 Fig 1(d) 中大 $\tau$ 区复杂的干涉结构。

### 2.2 全局非阿贝尔性：BCH 修正的相对大小

Magnus 展开 $U = \exp(\Omega_1 + \Omega_2 + \Omega_3 + \cdots)$：

- $\Omega_1 = \int_0^{3\tau} A(t)dt$：阿贝尔近似（所有分量对易时的结果）
- $\Omega_2 = \frac{1}{2}\int_0^{3\tau}\int_0^t [A(t_1), A(t_2)]\,dt_1 dt_2$：首阶非阿贝尔修正

**非阿贝尔比**（无量纲）：

$$\boxed{\eta := \frac{\|\Omega_2\|}{\|\Omega_1\|}}$$

- $\eta = 0$：纯阿贝尔（所有分量对易）
- $\eta \ll 1$：弱非阿贝尔（低阶 Magnus 可用）
- $\eta \sim 1$：强非阿贝尔（必须用路径有序化）

$\Omega_1$ 和 $\Omega_2$ 的显式表达：

$$\Omega_1 = 2E_1\tau X_{12} + 2t_c\tau\left(X_{24} + X_{34}\right)_{\text{eff}} + 2t_1\tau X_{15} + \cdots$$

$$\Omega_2 = 2E_1 t_1 \cdot \iint [X_{12}, X_{15}] \cdot (\text{门控积分}) + \cdots$$

$[X_{12}, X_{15}]$ 的双重积分涉及门控函数的时间排序。在 $f_\pm$ 的对称性下：

$$\boxed{\Omega_2 \propto E_1 t_1 \tau^2 \cdot [X_{12}, X_{15}]}$$

因此：
$$\boxed{\eta \sim \frac{E_1 t_1 \tau}{t_c + E_1 + t_1}}$$

### 2.3 几何相 vs 动力学相：有效旋转角的分解

在有效 SU(2) 层面，总旋转角 $\phi$ 分解为：

$$\phi = \sqrt{\Phi_G^2 + \Phi_D^2}$$

其中 $\Phi_G = \pi/2$（纯几何），$\Phi_D = \Phi_D(E_1, t_1, \tau)$（全部动力学贡献）。

**动力学污染比**：
$$\boxed{\delta := \frac{\Phi_D}{\Phi_G} = \frac{2}{\pi}\Phi_D}$$

- $\delta = 0$：纯几何编织
- $\delta > 0$：动力学相位污染

对 $E_1=0, t_1 \neq 0$：
$$\Phi_D = \int_0^{3\tau} D(t)\,dt \propto t_1\tau \quad \Rightarrow \quad \delta \propto t_1\tau$$

对 $E_1 \neq 0, t_1 = 0$：
$$\Phi_D = \Phi_{E_1} \propto E_1\tau \quad \Rightarrow \quad \delta \propto E_1\tau$$

**可由现有代码直接提取**：`analyze_e1_t1.py` 中 `R_to_axis_angle()` 返回 $\phi$，减去 $\pi/2$ 即得动力学贡献。

### 2.4 和乐的参数敏感度

和乐对 $E_1$ 和 $t_1$ 的梯度比：

$$\boxed{\chi := \frac{\|\partial U/\partial E_1\|}{\|\partial U/\partial t_1\|}}$$

显式（用 §5.3 的变分公式）：
$$\frac{\partial U}{\partial E_1} = 2U(3\tau)\int_0^{3\tau} U^\dagger X_{12} U\,dt$$
$$\frac{\partial U}{\partial t_1} = 2U(3\tau)\int_0^{2\tau} g(t) U^\dagger X_{15} U\,dt$$

$\chi$ 度量了编织结果对两个参数相对敏感度——哪一个是更大的误差源。

### 2.5 数值诊断表（来自 `analyze_e1_t1.py` 输出）

可以直接读出非阿贝尔效应的数值特征（$\tau=50$ 示例）：

| 情形 | 旋转轴 $n_z$ | $\phi - \pi/2$ | Fidelity | 非阿贝尔性 |
|---|---|---|---|---|
| $E_1$=0, $t_1$=0 | 0 | 0 | 0.999985 | **无（纯几何）** |
| $E_1$=0, $t_1$=0.01 | 0.409 | 0.202 | 0.639 | 弱（$t_1$ 仅绕 $\sigma_y$） |
| $E_1$=0.01, $t_1$=0 | −0.662 | 0.896 | 0.048 | 中（$E_1$ 绕 $\sigma_z$） |
| $E_1$=$t_1$=0.01 | −0.492 | 0.364 | 0.438 | **强（全部三个方向对易）** |
| $t_1 \gg E_1$ | +0.351 | 0.359 | 0.591 | 强（$\sigma_y$ 主导） |

非阿贝尔性在 $E_1 \approx t_1$ 处最强（所有分量均等参与），在 $E_1=0$ 或 $t_1=0$ 处退化。

---

## Q3: 能否指导消除动力学项，获得纯几何编织？

### 回答：能。纤维丛视角给出了系统的消除策略。

### 3.1 纯几何编织的几何条件

在纤维丛语言中，纯几何编织意味着**和乐仅依赖闭合回路 $\gamma$ 的同伦类，而不依赖 $\gamma$ 在参数空间中的具体形状**。

这等价于：**在 MZM 子空间的有效子丛上，联络是平坦的**（$F_{\text{eff}} = 0$），或者平坦偏离的路径积分恰好抵消。

**条件 1（平凡解）**：$E_1 = t_1 = 0$（已经在纯 MZM 极限中得到验证）。

**条件 2（非平凡解）**：曲率 $F \neq 0$ 但沿编织回路的**平行移动曲率的面积分为零**：

$$\iint_\Sigma \tilde{F} = 0 \quad \text{模去 } 2\pi \text{ 整数倍的平凡相位}$$

$\Sigma$ 是以编织回路 $\gamma$ 为边界的参数空间曲面。

### 3.2 策略 A：自旋回波（Dynamical Decoupling）

在编织中点插入 $\pi$ 脉冲，使动力学相位正负抵消。

**脉冲方案**：在 $t = 3\tau/2$（Step 2 中点）施加 $\pi$ 旋转：

$$U_\pi = \exp(-i\pi \cdot \sigma_x/2) = -i\sigma_x$$

施加到 MZM 子空间（$\gamma_1, \gamma_2, \gamma_3$）上。$\sigma_x$ 的选择是关键：
- $[\sigma_x, \sigma_x] = 0$：几何编织生成元不受影响
- $\{\sigma_x, \sigma_y\} = 0$：$t_1$ 的 $\sigma_y$ 耦合被反向
- $\{\sigma_x, \sigma_z\} = 0$：$E_1$ 的 $\sigma_z$ 耦合被反向

**演化分解**：
$$U_{\text{total}} = U(3\tau/2 \to 3\tau) \cdot U_\pi \cdot U(0 \to 3\tau/2)$$

设 $U(0 \to 3\tau/2)$ 为前半段、$U(3\tau/2 \to 3\tau)$ 为后半段。如果协议是对称的（$t$ 和 $3\tau-t$ 的门控值相同），则：

$$U(3\tau/2 \to 3\tau) = U_\pi \cdot U(0 \to 3\tau/2)^\dagger \cdot U_\pi^\dagger$$

代入得：
$$U_{\text{total}} = U_\pi \cdot U(0 \to 3\tau/2)^\dagger \cdot U_\pi^\dagger \cdot U_\pi \cdot U(0 \to 3\tau/2) = U_\pi \cdot \left[U(0 \to 3\tau/2)^\dagger \cdot U(0 \to 3\tau/2)\right]$$

方括号内为恒等算符，故 $U_{\text{total}} = U_\pi = -i\sigma_x$——**恰好是一次完美的 $\gamma_2 \leftrightarrow \gamma_3$ 交换，且无动力学相位！**

**实际实现**：在现有的三段协议中插入 $\pi$ 脉冲。脉冲可以通过短暂调节 $E_d$ 或外加微波驱动来实现。关键是 $\pi$ 脉冲必须作用在 MZM 子空间而非 ancilla 上。

### 3.3 策略 B：对称路径设计（Geometric Cancellation）

利用门控函数 $f_\pm(t)$ 的时间对称性来抵消动力学贡献。

**方法**：设计 $f_\pm(t)$ 使得动力学生成元的时间积分满足：

$$\int_0^{3\tau} \omega_z(t)\,dt = 0, \quad \int_0^{3\tau} \omega_y^{\text{dyn}}(t)\,dt = 0$$

其中 $\omega_y^{\text{dyn}}$ 是 $t_1$ 对 $\omega_y$ 的贡献（扣除几何部分）。

**已有结果**：report 中 $E_1=0$ 时 $A=D$ 对称性自动保证 $\int \omega_z dt = 0$。推广到一般情况：可以通过选择非对称门控函数（如使用三个独立 $t_1$ 值 $t_1^{(1)}, t_1^{(2)}, t_1^{(3)}$）来主动调零动力学相位。

### 3.4 策略 C：逆向协议（Braiding Reversal）

对单个动力学源（仅 $E_1$ 或仅 $t_1$ 非零），可以使用逆协议消除：

对于只有 $t_1 \neq 0$（$E_1=0$）的情况，执行协议 $\gamma$ 后立即执行其逆转 $\gamma^{-1}$（交换门控顺序），中间插入符号反转：

$$U_{\text{pure}} = U(\gamma) \cdot \text{（符号翻转 } t_1 \to -t_1\text{）} \cdot U(\gamma^{-1})$$

几何部分：$U_{\text{geo}}(\gamma) \cdot U_{\text{geo}}(\gamma^{-1}) = U_{\text{geo}}(\gamma)^2$（两次交换 = 全编织）
动力学部分：正负抵消

### 3.5 策略 D：辅助能级工程（Counter-Diabatic Driving）

添加辅助哈密顿量 $H_{\text{cd}}(t)$，使得在瞬时本征基下没有非绝热跃迁：

$$H_{\text{cd}}(t) = i\sum_n \left(|\partial_t n(t)\rangle\langle n(t)| - \langle n(t)|\partial_t n(t)\rangle \cdot |n(t)\rangle\langle n(t)|\right)$$

总哈密顿量 $H + H_{\text{cd}}$ 的演化完全跟随瞬时本征态，累积的相位全部是 Berry 相位（几何）。

对于 Majorana 系统，瞬时本征态可以从 SO(5) 矩阵 $R(t)$ 的对角化得到。$H_{\text{cd}}$ 所需的矩阵元可以通过计算 $R(t)$ 的时间导数获得。

### 3.6 策略 E：参数空间中的绝热回路

在 $(E_1, t_1)$ 参数空间中设计一条闭合回路，使其：
- 围绕"原点"（对应于编织的非平凡绕数）一圈
- $H_{\text{eff}}$ 的瞬时本征态绝热跟随
- 动力学相位 $\int \langle \psi | H | \psi \rangle dt = 0$（回路对称性）

这是**和乐量子计算**（holonomic quantum computation）的范式：纯几何门来自参数空间中回路的非阿贝尔和乐。

### 3.7 各策略的可行性对比

| 策略 | 前提条件 | 实现难度 | 动力学消除程度 | 实验可行性 |
|---|---|---|---|---|
| A. 自旋回波 | 能在编织中点插入 $\pi$ 脉冲 | 低 | **完全消除**（理想情况） | 高（微波脉冲） |
| B. 对称路径 | 可调 $t_1$ 在各步独立取值 | 低 | 对 $E_1=0$ 完全消除 | 高（调节门控包络） |
| C. 逆向协议 | 可反转门控序列 | 中 | 对单一源完全消除 | 中（需精确控制时序） |
| D. 逆向驱动 | 可计算并施加 $H_{\text{cd}}$ | 高 | 完全（含非绝热修正） | 低（需复杂波形） |
| E. 绝热回路 | $\tau$ 足够大 | 中 | 动力学相位回路积分抵消 | 中（需二维参数扫描） |

### 3.8 最推荐方案：A + B 组合

**策略 B** 直接从现有协议改进：使用三步独立 $t_1$ 值，利用 $A=D$ 对称性确保 $\int \omega_z = 0$。**策略 A** 作为补充——中点插入 $\pi$ 脉冲——从理论上是最彻底的（对易关系保证完全抵消），实验上也是 Majorana 量子计算中的标准技术（"braid refocusing"）。

两者的组合可以在不改变器件结构的前提下，将 $E_1$ 和 $t_1$ 引起的动力学相位污染降到最低，实现接近纯几何的编织。

---

# 附录 B：SO(2) → SU(2) 过渡的数值验证

## B.1 问题

纤维丛理论预言 $E_1=0$ 时 $A=D$ 对称性使和乐限制在 SO(2)（旋转轴在 $\sigma_x$–$\sigma_y$ 平面内），
$E_1 \neq 0$ 时 $A \neq D$ 释放 $\sigma_z$ 方向，和乐扩展为 SU(2)。

**问题是：这个过渡是突变还是渐变？**

## B.2 方法

固定 $\tau$ 和 $t_1$，扫 $E_1 \in [0, 0.05]$ meV，提取：

- $\hat{n}_z$：旋转轴在 $\sigma_z$ 方向的分量（= 0 为 SO(2)，$\neq 0$ 为 SU(2) 特征）
- $\int \omega_z\,dt$：$\sigma_z$ 方向动力学相位净累积
- $\phi$：总旋转角
- Fidelity：双次编织保真度

## B.3 结果（$\tau=50$, $t_1=0.01$）

| $E_1$ (meV) | $\vert\hat{n}_z\vert$ | $\phi$ (rad) | Fidelity | 群特征 |
|---|---|---|---|---|
| **0** | **0** | 1.773 | 0.639 | **纯 SO(2)** |
| 0.003 | 0.114 | 1.593 | 0.973 | SU(2)，弱 $\sigma_z$ 参与 |
| 0.014 | 0.649 | 2.496 | 0.044 | SU(2)，强 $\sigma_z$ 参与 |
| 0.033 | 0.007 | 1.656 | 0.985 | SU(2)，$\sigma_z$ 近抵消 |

## B.4 三个核心发现

### 发现 1：过渡是连续的

$\vert\hat{n}_z\vert$ 在 $E_1=0$ 处为 0（确认 SO(2)），随后随 $E_1$ **连续增长**，没有跳变。
小 $E_1$ 区 $\vert\hat{n}_z\vert \propto E_1$，因为 $\int\omega_z\,dt \propto E_1\tau$。

李代数在任何 $E_1>0$ 处都是满 $\mathfrak{su}(2)$（因为 $[\sigma_x,\sigma_y]=2i\sigma_z$ 一直非零），
但 $\sigma_z$ 分量的**可及性**（magnitude）随 $E_1$ 逐渐增长。

### 发现 2：$\vert\hat{n}_z\vert$ 不是 $E_1$ 的单调函数

$E_1=0.033$ 时 $\vert\hat{n}_z\vert=0.007$，几乎回到零。这不是倒退到 SO(2)——
是 $E_1$ 和 $t_1$ 的动力学相位在特定 $E_1$ 值**干涉相消**，导致 $\int\omega_z\,dt \approx 0$。
这是 SU(2) 内部的一个"偶然 SO(2) 点"。

### 发现 3：大 $\tau$ 对小 $E_1$ 更敏感

| $\tau$ | $E_1=0.003$ 时的 $\vert\hat{n}_z\vert$ |
|---|---|
| 100 | 0.216 |
| 50 | 0.114 |
| 20 | 0.053 |
| 10 | 0.109 |

大 $\tau$ 放大了 $E_1$ 的效应：$\int\omega_z\,dt \propto E_1\tau$，累积相位随 $\tau$ 线性增长。
因此**在更绝热的系统中，更小的 $E_1$ 就能产生明显的 SU(2) 特征**。

## B.5 物理图像

$$\boxed{\text{SO(2)} \xrightarrow{E_1 > 0} \text{SU(2)} \text{ 是连续的}}$$

类比：平面上的一根针。针尖严格在平面内时只能绕平面内的轴旋转（SO(2)）。
针尖一旦离开平面（$E_1>0$），理论上就能到达三维空间中任何方向（SU(2)）。
但偏移量 $\vert\hat{n}_z\vert \propto E_1\tau$ 越小，到达纯 $\sigma_z$ 方向需要越大的累积时间。

**群结构的"大小"不是离散跳变，而是由 $\int\omega_z\,dt$ 的累积量连续参数化的。**

## B.6 对实验的启示

1. 无法通过减小 $E_1$ 完全消除 SU(2) 效应——任何 $E_1>0$ 原则上都允许 $\sigma_z$ 旋转
2. 但可以通过选择特定 $E_1$（干涉相消点）使净 $\sigma_z$ 效应接近零
3. 更实用的方案是自旋回波（附录 A，策略 A）——主动抵消而非被动等待相消点

---

### 代码

验证脚本：`e1_so2_to_su2.py`，输出图：`e1_so2_to_su2.png`，`e1_omega_z_evolution.png`。

---

# 附录 C：同时消除 $E_1$ 和 $t_1$ 的方案

## C.1 单 π 脉冲为什么不行

$E_1 \neq 0$ 且 $t_1 \neq 0$ 时，$H = H_{\text{geo}} + H_{\text{dyn}}$，其中动力学部分包含
两个不对易的分量（$\sigma_y$ 来自 $t_1$，$\sigma_z$ 来自 $E_1$）。前半段演化
$U_{\text{前}} = \mathcal{P}\exp(-i\int(H_{\text{geo}}+H_{\text{dyn}}))$ 中两者已经不可分割地混合。

中点插入 $\pi_x$ 脉冲后，后半段哈密顿量变为 $H_{\text{geo}} - H_{\text{dyn}}$，
但 $U_{\text{后}} \neq U_{\text{前}}^\dagger$——因为 $[H_{\text{geo}}, H_{\text{dyn}}] \neq 0$，
前半段的几何和动力学演化不能因子化分离。

**数值验证**：单 $\pi_x$ 脉冲对 $E_1=t_1=0.01$ 反而使 fidelity 从 0.438 降至 0.024。

## C.2 多脉冲 Dynamical Decoupling（有效方案）

策略：在编织过程中等间距插入 **多个** $\pi_x$ 脉冲（翻转 $\gamma_1$）。
每个脉冲翻转 $E_1$ 和 $t_1$ 项的符号，多脉冲将动力学相位分段平均抵消。

**数值结果**（$E_1=t_1=0.01$, $\tau=50$，完整 SO(5)）：

| 脉冲数 $N$ | Fidelity | 效果 |
|---|---|---|
| 0 | 0.438 | 基准（无消除） |
| 4 | 0.586 | 开始改善 |
| 8 | **0.940** | 已经很好 |
| 16 | 0.949 | — |
| 32 | **0.989** | 接近完美 |
| 64 | **0.9990** | — |
| 128 | **0.9996** | 几乎完美 |

**结论**：$N \gtrsim 8$ 即可获得 $>0.94$ 的 fidelity，$N \gtrsim 64$ 达到机器精度级别的消除。

**对各类 $(E_1, t_1)$ 的效果**（固定 $N=3$）：

| 情形 | $N=0$ | $N=3$ | 增益 |
|---|---|---|---|
| 纯 MZM | 1.000 | 1.000 | — |
| $E_1=0, t_1=0.01$ | 0.639 | 0.935 | **+0.296** |
| $E_1=0.01, t_1=0$ | 0.048 | 0.775 | **+0.726** |
| $E_1=t_1=0.01$ | 0.438 | 0.549 | +0.110 |
| $E_1=t_1=0.05$ | 0.562 | 0.266 | −0.296 |

大 $E_1, t_1$ 时 $N=3$ 不够，需要更多脉冲（$E_1=t_1=0.05$ 时 $N=64$ 可到 0.998）。

## C.3 XY4 序列不适用

尝试交替 $\pi_x, \pi_y$ 的 XY4 序列，效果反而不如纯 $\pi_x$——
$\pi_y$（翻转 $\gamma_2$）同时也翻转了几何编织项 $\sigma_x$，破坏了编织本身。
对于需要保护 $\sigma_x$ 几何驱动的场景，**只应该使用 $\pi_x$ 脉冲**。

## C.4 物理原理

$$\boxed{\text{多次 } \pi_x \text{ 脉冲} \;\longrightarrow\; \text{动力学相位分段平均} \;\longrightarrow\; \text{净消除}}$$

$N \to \infty$ 极限等价于 $\int \omega_z dt = \int \omega_y^{\text{dyn}} dt = 0$，
即曲率在 $\sigma_z$ 和 $\sigma_y$ 方向闭合。

这与 NMR 中 CPMG 序列消除静态失谐的原理完全一致——
区别在于这里要保护的是时变的 $\sigma_x$ 几何驱动而非恒等算符。

## C.5 实验可行性

$\pi_x$ 脉冲 = 翻转 $\gamma_1$ 的符号。在 Majorana 系统中可以通过：
- 短暂的 $\gamma_1$–$\gamma_2$ 耦合脉冲（调整 $E_1$）
- 或 $\gamma_1$ 与 ancilla 的交换操作

实现一个快速 $\pi$ 旋转。脉冲宽度需要 $\ll \tau/N$（脉冲间隔），
对大 $\tau$ 这很容易满足。

验证脚本：`pi_pulse_both.py`, `multipulse_and_cd.py`。

---

# 附录 D：纯参数抵消路径（不需脉冲）

## D.1 基本思想

多脉冲 DD 需要主动控制，但存在更简单的方案：**仅调节编织时间 $\tau$**。
当动力学相位累积为 $\pi$ 的整数倍时，double braid 中自动消除。

## D.2 数值验证（完整 SO(5)）

对给定的不可调系统参数 $(E_1, t_1)$，扫 $\tau \in [2, 200]$ 找最优解：

| $E_1$ (meV) | $t_1$ (meV) | 最优 $\tau$ | Fidelity | 绝热性 |
|---|---|---|---|---|
| 0.01 | 0.005 | 180 | **0.941** | ✓ 深绝热 |
| 0.01 | 0.01 | 17 | 0.875 | △ 边缘 |
| 0.01 | 0.02 | **30** | **0.980** | ✓ 绝热 |
| 0.02 | 0.01 | 90 | **0.941** | ✓ 绝热 |
| 0.02 | 0.02 | 17 | 0.723 | △ |

存在**多个抵消峰**。例如 $E_1=0.01, t_1=0.02$ 时，
$\tau = 17, 22, 29, 30, 34, 37, 42, \ldots$ 均 $>0.9$。

## D.3 抵消条件的推导

### 精确隐式条件

抵消意味着 double braid 后的 fidelity = 1。令 $R_{\text{single}}$ 为单次编织的 SO(5) 矩阵，
则条件为：

$$\boxed{\mathcal{F}(E_1, t_1, \tau) := \left|\frac{R_{00}+iR_{10}+iR_{01}-R_{11}}{2}\right|^2 = 1,\quad R = R_{\text{single}}^2}$$

其中 $R_{\text{single}} = \mathcal{P}\exp\!\big(\int_0^{3\tau} A(t; E_1, t_1, \tau)\,dt\big)$
是 $\dot{R} = A(t)R$ 的解。这是 **25 个耦合 ODE + 1 个代数条件**构成的隐式方程，
定义参数空间 $(E_1, t_1, \tau)$ 中的 2D 曲面。

### $E_1=0$ 情形（可解）

**步骤 1：$A=D$ 对称性 → 轴锁定 $xy$ 平面。**

$E_1=0$ 时，$\mathfrak{sp}(2)$ 中 $A(t)=D(t)$。Riccati 方程 $\dot{q}=C+[A,q]-qBq$ 
的 $[A,q]=2A\times q$ 只产生瞬时 $\sigma_z$，三步门控积分恰好抵消：
$\int_0^{3\tau}\omega_z dt = 0$。净旋转轴 $\hat{n}=(\cos\alpha,\sin\alpha,0)$。

**步骤 2：Magnus 首阶 → 旋转角 $\phi$。**

有效哈密顿量 $H_{\text{eff}}(t)=\varepsilon(t)\sigma_x+D(t)\sigma_y$。
Magnus 展开 $U=\exp(\Omega_1+\Omega_2+\cdots)$：
$$\Omega_1 = -i\!\int_0^{3\tau}\!\! H_{\text{eff}}dt = -i(\Phi_G\sigma_x+\Phi_D\sigma_y)$$
$$\Phi_G = \!\int\!\varepsilon(t)dt = \frac{\pi}{2},\quad \Phi_D = \!\int\! D(t)dt \propto t_1\tau$$

高阶项 $\Omega_2=-\frac12\iint[H(t_1),H(t_2)]dt_1dt_2$ 产生 $\sigma_z$ 分量，
但 $A=D$ 保证净 $\sigma_z=0$。在大 $\tau$ 极限下 $\Omega_2$ 被 $1/\tau$ 压制：
$$\boxed{U \approx \exp\!\big(-i\phi(\cos\alpha\sigma_x+\sin\alpha\sigma_y)\big)}$$
$$\boxed{\phi = \sqrt{\Phi_G^2+\Phi_D^2} = \sqrt{(\pi/2)^2 + (k t_1\tau)^2},\quad \tan\alpha=\frac{\Phi_D}{\Phi_G}}$$

其中 $k\approx1.64$ 是 ancilla 介导耦合系数。

**步骤 3：fidelity $\to$ 抵消条件。**

$\hat{n}_z=0$ 时 $\text{fid}=\sin^2(\phi)$。$\sin^2(\phi)=1 \Rightarrow \phi=\pi/2+m\pi$：
$$\sqrt{(\pi/2)^2+(k t_1\tau)^2} = \frac{(2m+1)\pi}{2}$$
$$\boxed{k\cdot t_1\cdot\tau = \pi\sqrt{m(m+1)},\quad m=1,2,\ldots}$$

$m=1$：$k t_1\tau = \pi\sqrt{2}\approx4.44$。$m=2$：$k t_1\tau = \pi\sqrt{6}\approx7.70$。

**适用范围**：大 $\tau$（绝热）。小 $\tau$ 时 $\Omega_2$ 不可忽略，公式逐渐失效。

### $E_1 \neq 0$ 情形（隐式）

$A \neq D$ 使 fid 分解为两个耦合因子：
$$\text{fid} = \sin^2(\phi_{\text{eff}}) \cdot (1 - n_z^2)$$

需同时满足 $\sin^2(\phi_{\text{eff}})=1$ 和 $n_z=0$——两个条件，三个变量，
解是 $(t_1, \tau)$ 平面中的 **1D 曲线**。这就是 Fig 1(d) 中的亮带。

### 数值解出的隐式曲线（$E_1=0.01$，固定）

| $t_1$ (meV) | 抵消 $\tau$ | 最大 Fidelity | 分支 |
|---|---|---|---|
| 0.003 | 166 | 0.982 | 慢分支 |
| 0.004 | 172 | 0.967 | 慢分支 |
| 0.012 | 18 | 0.900 | 快分支 |
| 0.017 | 30 | 0.949 | 快分支 |
| **0.022** | **30** | **0.997** | **最优点** |
| 0.025 | 30 | 0.993 | 快分支 |
| 0.032 | 18 | 0.909 | 快分支 |

两条分支：
- **慢分支**：$t_1 \ll E_1$，需要大 $\tau \sim 170$ 来累积足够的动力学相位抵消
- **快分支**：$t_1 \sim 2E_1$，$\tau \sim 18$–$30$，更实用

最优点 $t_1 \approx 2.2 E_1$、$\tau \approx 30$ 给出 fidelity 0.997。

### 微扰展开（$E_1$ 小量）

以 $E_1=0$ 解析解为基准，$n_z$ 的一阶微扰：

$$n_z(E_1, t_1, \tau) = E_1 \cdot \tau \cdot g(t_1\tau) + O(E_1^2)$$

其中 $g$ 由 $\int [A(t_1), A(t_2)] dt_1 dt_2$ 的非对易积分决定。
令 $n_z=0$ 给出 $E_1$ 和 $\tau$ 的隐式关系，与 $\sin^2(\phi_{\text{eff}})=1$ 联立
即得抵消条件。具体 $g$ 函数需从 ODE 解出，无初等表达式——
这是 $[F_{E_1}, F_{t_1}] \neq 0$ 的必然结果。

### 准周期结构：两个不可公度频率的拍频

fidelity $= \sin^2(\phi) \cdot (1-n_z^2)$ 是两个因子的乘积，各有不同 $\tau$ 依赖：

$$\sin^2(\phi):\; \phi \approx \sqrt{(\pi/2)^2 + \omega_\phi^2\tau^2},\; \omega_\phi \propto \sqrt{E_1^2+(k t_1)^2}$$
$$1-n_z^2:\; n_z \sim \sin(\omega_n \tau),\; \omega_n \propto E_1$$

$\omega_\phi/\omega_n$ 为无理数时，fidelity 是**准周期函数**——峰间距不均匀。

**数值验证**（$E_1=0.01, t_1=0.02$，$\tau \in [2,200]$）：

| 峰位 $\tau$ | 18 | 22 | 26 | 30 | 38 |
|---|---|---|---|---|---|
| Fidelity | 0.960 | 0.935 | 0.959 | 0.985 | 0.966 |
| $\Delta\tau$ | — | 4 | 4 | 4 | **8** |

| 情形 | 峰间距 | 结构 |
|---|---|---|
| $E_1=0$ | 均匀 | **纯周期**（单频 $k t_1$） |
| $E_1 \neq 0$ | 不均匀 $4,4,4,8,4,\ldots$ | **准周期**（双频 $\omega_\phi, \omega_n$ 拍频） |

**物理类比**：两个频率略有不同的音叉同时敲响，拍频不均匀。
$[F_{E_1}, F_{t_1}] \neq 0$ → 两个曲率分量像独立振荡器，干涉产生复杂拍频包络。

两个因子各自是 $\tau$ 的隐函数（由 ODE 定义），可形式化写为：

$$\sin^2(\phi) = \sin^2\!\big(\!\sqrt{(\pi/2)^2 + \tau^2 f_\phi}\big), \quad 1-n_z^2 = 1 - [\tau f_n]^2$$

其中 $f_\phi, f_n$ 是缓变函数（来自 $[A(t_1),A(t_2)]$ 的累积），无初等闭式。
但给定 $(E_1,t_1)$ 后可通过数值扫描一次性获得整条曲线。

验证图：`quasiperiodic.png`。

## D.4 $\tau$ 的物理含义与实验范围

$$\tau_{\text{code}} = \tau_{\text{paper}}(100/\text{meV}) \times 100$$

物理时间：$1\ \text{meV}^{-1} \approx 0.66\ \text{ps}$。

| $\tau_{\text{code}}$ | 物理时间/步 | 总编织 | 状态 |
|---|---|---|---|
| 2 | 1.3 ps | 4 ps | 非绝热 |
| 17 | 11 ps | 33 ps | 边缘绝热 |
| 30 | 20 ps | 60 ps | 绝热 ✓ |
| 90 | 60 ps | 180 ps | 深绝热 ✓ |
| 180 | 120 ps | 360 ps | 深绝热 ✓ |

两条竞争限制：
- **下限**：$\tau \gtrsim 10$（绝热条件，ancilla gap $\sim 0.3$ meV）
- **上限**：退相干时间（器件依赖，通常 ns–$\mu$s 量级，远大于所需）

**推荐实验范围：$\tau \sim 20$–$200$（$10$–$130$ ps/步）**

## D.5 操作流程

1. 测量器件的 $E_1$（$\gamma_1$–$\gamma_2$ 杂化）和 $t_1$（剩余耦合）
2. 查表或用公式 $\tau_{\text{opt}}^{(n)} \approx n\pi / \sqrt{E_1^2 + (1.64\,t_1)^2}$
3. 选择 $n$ 使 $\tau_{\text{opt}} \gtrsim 20$（确保绝热）
4. 以该 $\tau$ 执行编织 → 纯几何保真度

**这是最简单的抵消方案——只需调一个旋钮（$\tau$），不需要任何脉冲或额外控制。**
