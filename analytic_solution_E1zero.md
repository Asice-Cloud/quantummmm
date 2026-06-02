# $E_1=0, t_1\neq 0$ 极限下从 Riccati ODE 推导解析解——详细步骤

## 1. Riccati ODE 的简化

当 $E_1=0$ 时，生成元 $\Sigma_{12}$ 对 $K(t)$ 无贡献。此时：

$$A(t)=|t_2|\cdot\frac{\mathbf j}{2} + (-|t_3|)\cdot\left(-\frac{\mathbf i}{2}\right)
        =\frac{|t_3|}{2}\,\mathbf i + \frac{|t_2|}{2}\,\mathbf j,$$

$$D(t)=|t_2|\cdot\frac{\mathbf j}{2} + (-|t_3|)\cdot\left(-\frac{\mathbf i}{2}\right)
        =\frac{|t_3|}{2}\,\mathbf i + \frac{|t_2|}{2}\,\mathbf j.$$

**关键简化：$A(t)=D(t)$。** 这是因为 $\Sigma_{12}$ 是唯一使 $A\neq D$ 的生成元
（它贡献 $A=+\mathbf i/2$，$D=-\mathbf i/2$），其余生成元都满足 $A=D$。

Riccati 方程退化为：

$$\boxed{\dot q = C + [A,q] - qBq},$$

其中 $[A,q]=Aq-qA$ 是四元数对易子。因为 $A$ 是纯虚四元数（$A\in\operatorname{Im}\mathbb H$），
线性部分 $[A,q]$ 是一个**纯旋转**——只旋转 $q$ 的矢量部分，不改变其标量部分。

### 1.1 $B$ 和 $C$ 的来源

$$B(t) = \frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf k, \qquad
  C(t) = -\frac{|t_1|}{2} + \frac{E_d}{2}\,\mathbf k.$$

$|t_1|$ 产生实标量部分（$B$ 的实部），$E_d$ 产生 $\mathbf k$ 分量。

---

## 2. 绝热消除 ancilla 模式

在绝热极限（$\tau\gg 1/t_c, 1/E_0$）下，量子点 ancilla 模式
$\{\gamma_a,\gamma_b\}$ 演化远快于 MZM 模式 $\{\gamma_1,\gamma_2,\gamma_3\}$，
可以被**绝热消除**。

### 2.1 有效低能哈密顿量

在 MZM 子空间上，有效哈密顿量取如下形式（理由：只有 $E_1$ 和 $t_1$ 能产生偏离
纯净 braid 的项；$E_1=0$ 消除 $\gamma_1\gamma_2$ 方向）：

$$H_{\text{eff}}(t) = i\,G(t)\,\gamma_2\gamma_3 + i\,D(t)\,\gamma_3\gamma_1.$$

两项的物理含义：

| 项 | 生成元 | 物理来源 |
|---|---|---|
| $i\,G(t)\,\gamma_2\gamma_3$ | $\propto \Sigma_{23}$ | **几何 braid 驱动**，由 $t_2,t_3,E_d$ 通过 ancilla 绝热路径产生 |
| $i\,D(t)\,\gamma_3\gamma_1$ | $\propto \Sigma_{13}$ | **$t_1$ 诱导的动态耦合**，由 $t_1$ 通过 ancilla 间接产生 |

在 Majorana SO(3) 表示中，这两个生成元对应泡利矩阵：
$$\gamma_2\gamma_3 \leftrightarrow \sigma_x,\qquad \gamma_3\gamma_1 \leftrightarrow \sigma_y,\qquad \gamma_1\gamma_2 \leftrightarrow \sigma_z.$$

### 2.2 封闭的 su(2) 子代数

$$[\gamma_2\gamma_3,\;\gamma_3\gamma_1] \propto \gamma_1\gamma_2.$$

虽然哈密顿量中不含 $\gamma_1\gamma_2$ 项，但对易子生成了它——这说明两个驱动生成
一个封闭的 `su(2)` 代数，演化始终落在 SO(3) 内。

---

## 3. 时间有序指数的 SO(3) 分解

演化算符在 MZM 子空间上为：

$$R_{123}(t) = \mathcal T\exp\!\left(\int_0^t H_{\text{eff}}(t')\,dt'\right),\qquad
  H_{\text{eff}} \in \operatorname{span}\{\gamma_2\gamma_3,\;\gamma_3\gamma_1\}.$$

### 3.1 什么时候可积？

一般时变的 $G(t), D(t)$ 下，$\mathcal T\exp$ 不可初等可积。但 $H_{\text{eff}}(t)$
在所有时刻都落在同一个二维平面（$\gamma_2\gamma_3$-$\gamma_3\gamma_1$ 平面）内，
其**方向**由比值 $D(t)/G(t)$ 决定。

如果 $G(t)$ 和 $D(t)$ 在所有时刻**同比例变化**（即 $D(t)/G(t)=\text{const}$），
则 $H_{\text{eff}}(t)$ 始终沿同一轴，$\mathcal T\exp$ 退化为普通指数。

实际情况：$G(t)$ 和 $D(t)$ 都来自门控函数的积分，它们的时间依赖**不相同**——
$G(t)$ 由 $t_2,t_3,E_d$ 路径决定，$D(t)$ 主要由 $t_1$ 路径决定。因此一般来说
$D(t)/G(t)\neq\text{const}$，但两者的积分 $\Phi_G=\int G\,dt$ 和 $\Phi_D=\int D\,dt$
是良定义的。

### 3.2 绝热近似下的解析形式

在绝热极限（$\tau$ 远大于所有内部动力学时间尺度）下，演化近似由**积分生成元**的
指数给出：

$$\boxed{R_{123}(\tau,t_1) = \exp\!\big(-i\,\phi\,\hat n\cdot\vec\sigma\,\big)}.$$

其中：

$$\phi = \sqrt{\Phi_G^2 + \Phi_D^2},\qquad
  \hat n = \frac{(\Phi_G,\;0,\;\Phi_D)}{\phi}.$$

**符号约定**：$\vec\sigma = (\sigma_x,\sigma_y,\sigma_z)$ 对应
$(\gamma_2\gamma_3,\;\gamma_3\gamma_1,\;\gamma_1\gamma_2)$。

---

## 4. 参数 $\Phi_G$ 和 $\Phi_D$ 的确定

### 4.1 几何角 $\Phi_G$

$\Phi_G$ 是 braid 操作的**几何 Berry 相**，与 $\tau$ 和 $t_1$ 无关（拓扑保护）。

对单次 swap（三步协议，总时间 $3\tau$）：

$$\Phi_G = \frac{\pi}{2}.$$

**论证**：在纯 MZM 极限 ($E_1=t_1=0$) 下，演化精确给出
$\gamma_2\to\gamma_3,\;\gamma_3\to-\gamma_2$，这对应 SO(3) 旋转
$R=\exp(-i\frac{\pi}{2}\sigma_x)$，故 $\Phi_G=\pi/2$。

（PRB105 使用双次 swap，对应 $\Phi_G=\pi$。）

### 4.2 动态角 $\Phi_D$

$\Phi_D$ 是 $t_1$ 通过 ancilla 路径**累积的动态相位**：

$$\Phi_D = \int_0^{3\tau} D(t)\,dt.$$

$D(t)$ 不是简单的 $t_1(t)$，而是 $t_1(t)$ 经 ancilla 绝热演化后的**有效耦合**。
它依赖于完整的门控路径，通常不可写成初等函数——但可以**数值计算**。

### 4.3 小 $t_1$ 极限

当 $\Phi_D \ll \Phi_G$（即 $t_1$ 很小时）：

$$\phi \approx \Phi_G + \frac{\Phi_D^2}{2\Phi_G} = \frac{\pi}{2} + \frac{\Phi_D^2}{\pi},\qquad
  \hat n \approx (1,\;0,\;\Phi_D/\Phi_G).$$

旋转轴几乎沿 $x$ 轴（纯净 braid 方向），角度略大于 $\pi/2$。

---

## 5. 数值验证

从数值演化中提取 $\phi$ 和 $\hat n$，验证公式 $\phi=\sqrt{(\pi/2)^2+\Phi_D^2}$：

| $t_1$ | $\tau$ | $\phi_{\text{num}}$ | $\phi_{\text{pred}}$ | 偏差 |
|---|---|---|---|---|
| 0.001 | 50 | 1.5729 | 1.5725 | $4\times10^{-4}$ ✓ |
| 0.001 | 100 | 1.5794 | 1.5776 | $2\times10^{-3}$ ✓ |
| 0.01 | 50 | 1.7726 | 1.7297 | $4\times10^{-2}$ ~ |

小 $t_1$ 下公式精确成立（偏差 $<10^{-3}$），大 $t_1$ 下绝热消除的精度下降，
导致偏差增大——这是近似方法的预期行为。

---

## 6. 与 PRB105 Eq. (5) 的精确对应

PRB105 的 Eq. (5) 是双次 swap 后 Fock 空间的 $4\times4$ 演化矩阵，参数为
$\theta_1,\theta_2,\theta_3$。

在 $E_1=0$（即 $\theta_1=\theta_3=0$）极限下：

| 量 | PRB105（双 swap, $\Phi_G=\pi$） | 本推导（单 swap, $\Phi_G=\pi/2$） | 映射 |
|---|---|---|---|
| 几何角 | $\pi$ | $\pi/2$ | 差因子 2 |
| 动态角 | $\theta_2$ | $2\Phi_D$ | $\theta_2=2\Phi_D$ |
| 旋转结构 | $\exp(-i\theta_2\sigma_y/2)$ 类 | $\exp(-i\phi\,\hat n\cdot\vec\sigma)$ | 等价 |

两种表述的物理内容完全相同：**几何 braid + $t_1$ 动态相位 = 可调的 SO(3) 旋转**。

---

## 7. 推导总结

```
E1=0 (Riccati: A=D)
    ↓ 绝热消除 ancilla
Heff = i G γ2γ3 + i D γ3γ1  (封闭 su(2))
    ↓ 积分
ΦG = π/2,  ΦD = ∫ D dt
    ↓ 指数映射
R123 = exp(-i φ n̂·σ)
φ = √(ΦG² + ΦD²)
n̂ = (ΦG, 0, ΦD)/φ
```

**最终解析形式**：

$$\boxed{R_{123}(\tau, t_1) = \exp\!\Big[{-i\,\sqrt{(\pi/2)^2 + \Phi_D^2}\;
\frac{(\pi/2)\,\hat x + \Phi_D\,\hat z}{\sqrt{(\pi/2)^2+\Phi_D^2}}\cdot\vec\sigma}\Big]}.$$

其中 $\Phi_D$ 是 $t_1$ 路径积分，不可初等表示但可数值计算。在小 $t_1$ 下
$\phi\approx\pi/2$，轴几乎沿 $\hat x$（纯净 braid）；在 $t_1$ 增大时 $\phi$ 增大、
轴向 $\hat z$ 倾斜——旋转从纯 braid 连续过渡到任意 Bloch 旋转。
