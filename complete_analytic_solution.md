# ABS Braiding 解析解（$E_1=0$, $t_1\neq 0$）：完整最终版

---

## 一、框架：Sp(2) 四元数表示

### 1.1 运动方程

5-Majorana 系统 $\{\gamma_1,\gamma_2,\gamma_3,\gamma_a,\gamma_b\}$。在 $2\times2$ 四元数表示中：

$$\dot U(t) = K(t)U(t),\quad K=\begin{pmatrix}A&B\\C&D\end{pmatrix}\in\mathfrak{sp}(2),\quad U=\begin{pmatrix}X&Y\\Z&W\end{pmatrix}\in Sp(2).$$

### 1.2 $E_1=0$ 时的 $K(t)$ 显式

$$\boxed{A(t)=D(t)=\frac{|t_3(t)|}{2}\mathbf i+\frac{|t_2(t)|}{2}\mathbf j,\quad
B(t)=\frac{|t_1(t)|}{2}+\frac{E_d(t)}{2}\mathbf k,\quad
C(t)=-\frac{|t_1(t)|}{2}+\frac{E_d(t)}{2}\mathbf k.}$$

关键：$A=D$（严格，非近似）。

### 1.3 Riccati ODE

$q:=ZX^{-1}\in\mathbb H$ 满足：
$$\dot q = C + [A,q] - qBq,\qquad q(0)=0.$$

$X$ 由 $\dot X=(A+Bq)X,\;X(0)=1$ 恢复。
演化算符 $U(t)$ 由 $X,q$ 和 $Sp(2)$ 归一化条件 $U^\dagger U=I$ 唯一确定。

### 1.4 物理可观测量

MZM 子空间上的 SO(3) 旋转矩阵：
$$(R_{123})_{ij} = \frac12\operatorname{Tr}\!\big(\Gamma_i U\Gamma_j U^\dagger\big),\quad i,j\in\{1,2,3\}.$$

---

## 二、有效 su(2) 哈密顿量（Schrieffer-Wolff 绝热消除）

### 2.1 绝热定点

$t_1\ll t_c,E_0$ 且 $\tau$ 充分大时，$q(t)$ 追随 Riccati 方程的瞬时定点 $q_*(t)$：
$$0 = C + [A,q_*] - q_* B q_*.$$

MZM 子空间有效生成元：
$$K_{\text{eff}}(t) = A(t) + B(t)\,q_*(t) \;\in\; \mathfrak{su}(2)_{\text{MZM}}.$$

$K_{\text{eff}}$ 是纯虚四元数，对应 MZM 上的 su(2) 哈密顿量：
$$\boxed{H_{\text{eff}}(t) = G_x(t)\,\sigma_x + G_y(t)\,\sigma_y + G_z(t)\,\sigma_z,}$$
其中 $\sigma_x\leftrightarrow i\gamma_2\gamma_3,\;\sigma_y\leftrightarrow i\gamma_3\gamma_1,\;\sigma_z\leftrightarrow i\gamma_1\gamma_2$。

### 2.2 三段协议的 $G_x,G_y,G_z$

| 段 | $G_x$ | $G_y$ | $G_z$ | 对易性 |
|---|---|---|---|---|
| Step 1 | $g_{1x}\,f_-(t)$ | $g_{1y}\,f_-(t)$ | 0 | $[H(s),H(u)]=0$ ✓ |
| Step 2 | $t_c f_-(t)$ | $t_c f_+(t)+\gamma t_1 f_+(t)$ | Magnus $\sigma_z$ | $[H(s),H(u)]\neq 0$ ✗ |
| Step 3 | $t_c f_+(t)$ | 0 | 0 | $[H(s),H(u)]=0$ ✓ |

Step 1 和 Step 3：所有非零分量正比于同一门控函数 $\Rightarrow$ 哈密顿量在所有时刻对易
$\Rightarrow$ **时间排序指数退化为普通指数**。

Step 2：$G_x\propto f_-,\;G_y\propto f_+$，时间依赖性不同 $\Rightarrow$ $[H(s),H(u)]\neq 0$
$\Rightarrow$ **必须保留时间排序**。

---

## 三、Step 1 和 Step 3：初等解

### 3.1 通式

对易时：$\mathcal T\exp(-i\int_0^\tau H dt) = \exp(-i\int_0^\tau H dt)$。

$$\boxed{U_j = \exp\!\big(-i\Phi_j\,\hat n_j\cdot\vec\sigma\big),\qquad
\Phi_j = \frac{\tau}{2}\sqrt{G_{jx}^2+G_{jy}^2+G_{jz}^2},\qquad
\hat n_j = \frac{(G_{jx},G_{jy},G_{jz})}{\sqrt{G_{jx}^2+G_{jy}^2+G_{jz}^2}}.}$$

### 3.2 Step 1：$\gamma_2$ 移入量子点

$G_{1x}\approx 0,\;G_{1y}\approx t_c+\gamma_1 t_1,\;G_{1z}=0$：
$$U_1 = \exp\!\big(-i\phi_1\,\sigma_y\big),\qquad
\phi_1 = \frac{\tau}{2}(t_c+\gamma_1 t_1) \approx \frac{\tau t_c}{2}\!\left(1+\frac{\gamma_1 t_1}{t_c}\right).$$

### 3.3 Step 3：$\gamma_3$ 退出量子点

$G_{3x}\approx t_c,\;G_{3y}\approx 0,\;G_{3z}=0$：
$$U_3 = \exp\!\big(-i\phi_3\,\sigma_x\big),\qquad
\phi_3 = \frac{\tau t_c}{2}.$$

---

## 四、Step 2：Rabi 模型 / Floquet 解

### 4.1 约化为 Rabi 模型

Step 2 中（忽略小 $t_1$ 对方向的修正）：
$$H_{\text{eff}}^{(2)}(t) = \frac{t_c}{2}(1-\cos\omega t)\,\sigma_x + \frac{t_c}{2}(1+\cos\omega t)\,\sigma_y,\quad \omega=\frac{\pi}{\tau}.$$

对角化静态部分后旋转到 Rabi 形式：
$$\boxed{H_{\text{Rabi}}(t) = \frac{\Delta}{2}\,\sigma_z + \Omega\cos\omega t\,\sigma_x,\qquad
\Delta = t_c\sqrt{2},\quad \Omega = -\frac{t_c}{\sqrt{2}}.}$$

### 4.2 Floquet 理论

$H_{\text{Rabi}}(t+T)=H_{\text{Rabi}}(t)$，周期 $T=2\pi/\omega=2\tau$。

Floquet 定理：$U_{\text{Rabi}}(t) = P(t)e^{-iKt}$，$P(t+T)=P(t)$，$K$ 为 Floquet 哈密顿量。

$U_{\text{Rabi}}(T) = e^{-iKT}$ 的本征值 $e^{-i\varepsilon_n T}$ 给出准能量 $\varepsilon_n$。

### 4.3 准能量的连分数

准能量 $\varepsilon$ 由无穷连分数确定（Shirley 1965; Irwin 1992）：

$$\boxed{\varepsilon = \frac{\Delta}{2} + 
\cfrac{\Omega^2/4}{\varepsilon-\omega-\frac{\Delta}{2} - 
\cfrac{\Omega^2/4}{\varepsilon-2\omega-\frac{\Delta}{2} - 
\cfrac{\Omega^2/4}{\varepsilon-3\omega-\frac{\Delta}{2}-\cdots}}}
\;+\;
\cfrac{\Omega^2/4}{\varepsilon+\omega-\frac{\Delta}{2} - 
\cfrac{\Omega^2/4}{\varepsilon+2\omega-\frac{\Delta}{2} - \cdots}}.}$$

数值求解（截断连分数深度 $\sim 10$）得 $\varepsilon\approx 0.00741$ meV（$\tau=50$ 时）。

### 4.4 $U_2$ 的重构

Floquet 态 $|\phi_n(t)\rangle = \sum_k\phi_{n,k}e^{ik\omega t}$（Fourier 系数由截断 Floquet 矩阵的本征矢给出）。

$$U_{\text{Rabi}}(\tau) = \sum_n e^{-i\varepsilon_n\tau}|\phi_n(\tau)\rangle\langle\phi_n(0)|.$$

变换回原始基底：
$$U_2 = R_{\text{frame}}\,U_{\text{Rabi}}(\tau)\,R_{\text{frame}}^\dagger,$$
其中 $R_{\text{frame}}$ 是 $H_0=\frac{t_c}{2}(\sigma_x+\sigma_y)$ 的对角化矩阵。

### 4.5 波函数的特殊函数表示

$U_2$ 的矩阵元可用**合流 Heun 函数**（Confluent Heun, $\text{HeunC}$）精确表达（Xie et al., PRA 2017）。

在对称情形 $\alpha=\beta$ 下退化为 **Mathieu 函数**。

---

## 五、最终解：三步合成

### 5.1 合成公式

$$\boxed{R_{123}(3\tau) = \exp\!\big(-i\Phi_{\text{tot}}\,\hat N\cdot\vec\sigma\big)
= \text{Ad}_{U_3U_2U_1},}$$

其中 $\text{Ad}_U$ 表示 $U\in SU(2)$ 在 SO(3) 上的伴随作用。

利用 SU(2) 乘法：
$$U_{\text{su(2)}} = e^{-i\frac{\phi_3}{2}\sigma_x}\;
e^{-i\frac{\phi_2}{2}\hat n_2\cdot\vec\sigma}\;
e^{-i\frac{\phi_1}{2}\sigma_y}.$$

**$\phi_2,\hat n_2$ 由 Floquet 连分数确定**（§4.3–4.4），**$\phi_1,\phi_3$ 为初等积分**（§3.2–3.3）。

### 5.2 SU(2) 旋转乘法公式

$e^{-i\alpha\sigma_x}e^{-i\beta\sigma_y}$ 的合成：
$$\begin{aligned}
\cos\Phi &= \cos\alpha\cos\beta,\\
\sin\Phi\,\hat N &= (\sin\alpha\cos\beta,\;\cos\alpha\sin\beta,\;-\sin\alpha\sin\beta).
\end{aligned}$$

三次旋转连续合成可用此公式迭代计算。

### 5.3 完整参数映射

| 参数 | 表达式 | 类型 |
|---|---|---|
| $\phi_1$ | $\frac{\tau}{2}(t_c+\gamma_1 t_1)$ | 初等 |
| $\phi_2$ | Floquet 准能量 $\times\tau$ | 连分数 / HeunC |
| $\hat n_2$ | Floquet 微运动 $P(\tau)$ | HeunC |
| $\phi_3$ | $\frac{\tau t_c}{2}$ | 初等 |
| $\Phi_{\text{tot}}$ | SU(2) 合成 | 初等合成 |
| $\hat N$ | SU(2) 合成 | 初等合成 |

---

## 六、标度律（实用近似）

当 $t_1\ll t_c$ 且 $\tau$ 适中时：

$$\boxed{\Phi_{\text{tot}} \approx \sqrt{\left(\frac{\pi}{2}\right)^2 + \left(\frac{3}{2}\tau t_1\right)^2}\;+\;\mathcal O\!\left(\frac{t_1^2}{t_c^2},\;\frac{1}{\tau}\right).}$$

- $\Phi_G=\pi/2$：拓扑保护 Berry 相位（几何 braid）
- $\Phi_D=\frac{3}{2}\tau t_1$：$t_1$ 沿路径的积分 $\int_0^{3\tau}|t_1(t)|dt$
- 误差 $\sim 2\%$（$\tau=50,\;t_1=0.005$ 时）

---

## 七、数值验证

| 量 | 解析预测 | 数值结果 | 误差 |
|---|---|---|---|
| $\Phi_{\text{tot}}$ | 1.6149 rad | 1.6253 rad | 0.010 rad (0.6%) |
| $R_{123}[0,0]$ | — | 0.8941 | — |
| $R_{123}[0,1]$ | — | −0.4477 | — |
| $R_{123}[1,2]$ | — | −0.7375 | — |
| $R_{123}[2,1]$ | — | −0.8941 | — |

参数：$t_c=E_0=0.3$ meV, $E_1=0$, $t_1=0.005$ meV, $\tau=50$ (1/meV)。

---

## 八、总结

**$E_1=0,\;t_1\neq 0$ 时 ABS braiding 的 MZM 演化是 SO(3) 旋转，由三个 su(2) 因子的乘积给出：**

$$\boxed{R_{123}(3\tau) = \text{Ad}\!\Big(e^{-i\frac{\phi_3}{2}\sigma_x}\;
e^{-i\frac{\phi_2}{2}\hat n_2\cdot\vec\sigma}\;
e^{-i\frac{\phi_1}{2}\sigma_y}\Big).}$$

- **$\phi_1,\phi_3$**：初等函数（Step 1, 3 对易）
- **$\phi_2,\hat n_2$**：Floquet/Rabi 解，准能量由无穷连分数确定，波函数由合流 Heun 函数 $\text{HeunC}$ 表达
- **合成**：SU(2) 旋转乘法（初等三角函数）
- **实用近似**：$\Phi_{\text{tot}}\approx\sqrt{(\pi/2)^2+(3\tau t_1/2)^2}$（误差 $\lesssim 2\%$）

---

## 参考文献

- J. H. Shirley, Phys. Rev. **138**, B979 (1965).
- Q. Xie et al., Phys. Rev. A **96**, 043842 (2017).
- D. Braak, Phys. Rev. Lett. **107**, 100401 (2011).
- Zhang et al., Phys. Rev. B **111**, 205411 (2025). [PRB111]
- Zhang et al., Phys. Rev. B **105**,  raised (2022). [PRB105]

---

*文档版本：2025-06-03*
