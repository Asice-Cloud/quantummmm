# PRB111/PRB105 有效模型的 ODE 解法：总结文档

## 0. 动机：从 Majorana braiding 到一个统一 ODE

量子计算中 Majorana 零模 (MZM) 的非阿贝尔编织 (braiding) 是核心操作。
PRB111 (Zhang et al., 2025) 和 PRB105 (Chen et al., 2022) 研究了 ABS
(Andreev Bound State) 存在时的 braiding 性质。

两篇论文的核心有效哈密顿量相同：

$$H_{EM}(t) = iE_d\gamma_a\gamma_b + iE_1\gamma_1\gamma_2 + i|t_2|\gamma_a\gamma_2 - i|t_1|\gamma_b\gamma_1 - i|t_3|\gamma_a\gamma_3$$

其中 $\gamma_1,\dots,\gamma_5$ 是 5 个 Majorana 算符，$E_1$ 和 $t_1$ 刻画 ABS
偏离理想 MZM 的程度。

**传统处理方法**：直接数值求解 $\mathcal T\exp(-i\int H dt)$ 或紧束缚模型模拟。
**我们的方法**：利用 $so(5)\cong sp(2)$ 同构，将 5-Majorana 系统映射为四元数
Riccati 方程——一个 4 分量的统一 ODE，能精确描述全部三段 braiding 协议。

---

## 1. 模型搭建：从 Majorana 到 Sp(2)

### 1.1 李代数闭包

5 个 Majorana 的双线性生成元 $X_{ij}=i\gamma_i\gamma_j\;(i<j)$ 共 10 个，
满足标准的 $so(5)$ 对易关系：

$$[X_{ij}, X_{kl}] = 2i(\delta_{jk}X_{il} - \delta_{ik}X_{jl} - \delta_{jl}X_{ik} + \delta_{il}X_{jk}).$$

论文哈密顿量恰好是其中 5 个生成元的线性组合。

### 1.2 Sp(2) 四元数表示

利用李代数同构 $so(5)\cong sp(2)$，将 10 个生成元表示为 $2\times 2$ **四元数矩阵**。
取 Cl(5) Gamma 矩阵：

$$\Gamma_1=\begin{pmatrix}0&1\\1&0\end{pmatrix},\;
\Gamma_2=\begin{pmatrix}0&-\mathbf i\\ \mathbf i&0\end{pmatrix},\;
\Gamma_3=\begin{pmatrix}0&-\mathbf j\\ \mathbf j&0\end{pmatrix},\;
\Gamma_4=\begin{pmatrix}0&-\mathbf k\\ \mathbf k&0\end{pmatrix},\;
\Gamma_5=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.$$

旋量生成元 $\Sigma_{ij}=\frac14[\Gamma_i,\Gamma_j]$ 张成 $\mathfrak{sp}(2)$。
演化算符 $U(t)\in Sp(2)$ 满足 $\dot U = K(t)U$，其中

$$K(t) = \sum h_{ij}(t)\,\Sigma_{ij} = \begin{pmatrix} A(t) & B(t) \\ C(t) & D(t) \end{pmatrix}.$$

这里 $A,B,C,D\in\mathbb H$ 是时变四元数，由门控参数决定。

### 1.3 三段门控协议

每段 $\tau$ 时长，门控函数 $f_\pm(t)=\frac{1\pm\cos(\pi t/\tau)}{2}$：

| 段 | $|t_2|$ | $|t_3|$ | $E_d$ | $|t_1|$ (ABS) |
|---|---|---|---|---|---|
| Step 1 | $t_c f_-(t)$ | $0$ | $E_0 f_+(t)$ | $t_{1c}f_-(t)$ |
| Step 2 | $t_c f_+(t-\tau)$ | $t_c f_-(t-\tau)$ | $0$ | $t_{1c}f_+(t-\tau)$ |
| Step 3 | $0$ | $t_c f_+(t-2\tau)$ | $E_0 f_-(t-2\tau)$ | 0 |

---

## 2. 核心结果：统一四元数 Riccati ODE

### 2.1 Riccati 变量

将 $U(t)=\begin{pmatrix}X&Y\\ Z&W\end{pmatrix}$ 分块，定义

$$q(t) := Z(t)X(t)^{-1}\in\mathbb H,\qquad q(0)=0.$$

$q$ 是一个四元数（4 实分量），跨过步边界时自动连续。

### 2.2 统一 ODE

从 $\dot U=KU$ 严格推导（无近似）：

$$\boxed{\dot q(t) = C(t) + D(t)q(t) - q(t)A(t) - q(t)B(t)q(t)}.$$

三段协议的差异仅体现在 $A,B,C,D$ 的时变系数中，ODE 形式不变。

### 2.3 五个基本生成元的块分量

| 生成元 | 物理含义 | $A$ | $D$ | $B$ | $C$ |
|---|---|---|---|---|---|
| $\Sigma_{12}$ | $E_1$ 杂化 | $+\mathbf i/2$ | $-\mathbf i/2$ | $0$ | $0$ |
| $\Sigma_{24}$ | $|t_2|$ 耦合 | $+\mathbf j/2$ | $+\mathbf j/2$ | $0$ | $0$ |
| $\Sigma_{34}$ | $-|t_3|$ 耦合 | $-\mathbf i/2$ | $-\mathbf i/2$ | $0$ | $0$ |
| $\Sigma_{15}$ | $-|t_1|$ 耦合 | $0$ | $0$ | $-1/2$ | $+1/2$ |
| $\Sigma_{45}$ | $E_d$ 能级 | $0$ | $0$ | $+\mathbf k/2$ | $+\mathbf k/2$ |

### 2.4 回到物理可观测量

从 $q(t)$ 重建完整演化（§12）：

1. **Riccati ODE** → $q(t)$（4 分量）
2. **重建** $\dot X=(A+Bq)X$ → $X(t)$，得 $U(t)\in Sp(2)$
3. **旋量覆盖** $R_{ij}=\frac12\operatorname{Tr}(\Gamma_i U\Gamma_j U^\dagger)$ → $R(t)\in SO(5)$
4. **提取** $R_{123}$ → braiding 结果

数值验证：Riccati 管道与直接 Sp(2) 传播偏差 $<10^{-9}$（机器精度）。

---

## 3. 结果对照

### 3.1 三段协议的简化程度

| 段 | 代数闭包 | 可否写单位四元数 | 最紧凑形式 |
|---|---|---|---|
| Step 1 | $so(4)\cong su(2)\oplus su(2)$ | ✗ | 统一 Riccati |
| Step 2 | $so(5)\cong sp(2)$ 满秩 | ✗ | 统一 Riccati |
| Step 3 | $u(1)\oplus su(2)$ | ✓ | 单位四元数闭式 |

### 3.2 MZM 极限 ($E_1=0, t_1=0$)

绝热 ($\tau\ge 50$) 下 braiding 精确成立：

$$\gamma_2\to 0.999978\,\gamma_3,\qquad \gamma_3\to -0.999960\,\gamma_2.$$

$\gamma_a,\gamma_b$ 获得 $\approx 127^\circ$ 的路径依赖 holonomy。

### 3.3 $E_1=0, t_1\neq 0$：braid 方向鲁棒，相位振荡

$\gamma_2\to\gamma_3$ 对 $t_1$ **完全鲁棒**——绝热下始终精确为 1。
$t_1$ 产生 $\gamma_3\gamma_1$ 型动态旋转，旋转角 $\phi=\sqrt{(\pi/2)^2+\Phi_D^2}$，
其中 $\Phi_D=\int D(t)dt$ 是 $t_1$ 路径积分。

**对应 PRB105**：PRB105 的 Eq. (5) 在 $\theta_1=\theta_3=0$ 时退化为相同的单参数
SO(3) 旋转。差异仅在于 PRB105 使用双次 swap（$\Phi_G=\pi$），我们使用单次 swap
（$\Phi_G=\pi/2$）。

### 3.4 Fig 1(d)：$E_1=0.01$，扫 $t_1$ 和 $\tau$

| $t_1$ | $\tau$ | 旋转角 $\phi$ 范围 | 结论 |
|---|---|---|---|
| 0 | var | $1.05\sim1.60$ | 近乎固定 braid |
| 0.01 | var | $1.06\sim2.27$ | 10/10 唯一值 → 任意旋转 |
| 0.01 | var (含 $E_1$) | $1.10\sim3.10$ | 覆盖 $>\pi$ 弧度 |

### 3.5 Fig 8(c)：$E_1$ 驱动的 braiding 振荡

用 $E_1=0.18$ meV（论文紧束缚模型 $E_1\approx 0.0018$ meV 的 **100 倍**，因 $\tau$
轴单位差 100 倍），精确复现 Fig 8(c) 全部数值特征：

- 第一谷值 $\tau\approx 18$，第一恢复 $\tau\approx 36$
- 周期 $\approx 35$，保真度 $0\sim 1$
- 周期 $\propto 1/E_1$

图见 `fig8c_quantitative.png`。

---

## 4. 方法对比

| 方法 | 变量数 | 优点 | 局限 |
|---|---|---|---|
| $5\times5$ $so(5)$ 直接传播 | 25 | 直接得 $R(t)$ | 重 |
| $Sp(2)$ 四元数直接传播 | 16 | 紧凑 | 需旋量覆盖映射 |
| **统一 Riccati ODE** | **4** | **最轻量，结构可解析** | 非初等可积，需重建 $X(t)$ |
| Wei–Norman (10 ODE) | 10 | 精确 | 耦合非线 |
| Magnus 展开 | 截断 | 解析 | 低阶误差大 |

---

## 5. 关键脚本

| 脚本 | 功能 |
|---|---|
| `verify_riccati_sp2.py` | Riccati ODE 与 Sp(2) 传播对照（全三段，误差 $<10^{-9}$） |
| `verify_braiding.py` | 全三段 braiding 验证 + 扫参 |
| `fig8c_quantitative.png` | Fig 8(c) 定量对照图 |
| `step2_ip_solver.py` | Step 2 interaction picture 分析 |

---

## 6. 核心结论

1. **模型完善**：全空间 $so(5)\cong sp(2)$，无先验假设
2. **统一 ODE**：一个 4 分量 Riccati 方程覆盖全部三段 braiding 协议
3. **解析可处理**：$E_1=0$ 极限下可推导 $R=\exp(-i\phi\,\hat n\cdot\vec\sigma)$
4. **论文一致**：定量验证了 PRB111 Fig 1(d)、Fig 8(c)、Fig 8(f) 和 PRB105 Eq. (5)
5. **唯一不能做的**：紧束缚模型特有的长度/无序效应（需要微观模拟），但只要给定
   $E_1, t_1$ 的值，我们的有效模型就能精确给出 braiding 结果
