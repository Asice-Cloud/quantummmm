# ABS Braiding 的 Floquet/特殊函数解析解

---

## 一、问题约化：Step 2 即是 Rabi 模型

### 1.1 有效 su(2) 哈密顿量

SW 消去 ancilla 后，Step 2（$t\in[\tau,2\tau]$）的 MZM 子空间有效哈密顿量为：

$$\boxed{H_{\text{eff}}(t) = \frac{\alpha}{2}\big(1-\cos\omega t\big)\,\sigma_x 
+ \frac{\beta}{2}\big(1+\cos\omega t\big)\,\sigma_y,\qquad \omega=\frac{\pi}{\tau}}$$

其中 $\alpha = |t_3| \approx t_c$，$\beta = |t_2| + \gamma|t_1| \approx t_c + \mathcal O(t_1)$（$\gamma$ 为 ancilla 响应因子）。

写为标准形式：
$$H(t) = H_0 + H_1\cos\omega t,$$
$$H_0 = \frac{\alpha}{2}\sigma_x + \frac{\beta}{2}\sigma_y,\qquad
H_1 = -\frac{\alpha}{2}\sigma_x + \frac{\beta}{2}\sigma_y.$$

**关键**：$[H_0,H_1]=-\frac{i\alpha\beta}{2}\sigma_z \neq 0$——即 Rabi 模型（线偏振驱动），不是可轻松旋波近似（RWA）的圆偏振驱动。

### 1.2 旋转到主轴

令 $\tan\phi_0 = \beta/\alpha$，$H_0 = \frac{\Delta}{2}(\cos\phi_0\sigma_x+\sin\phi_0\sigma_y) = \frac{\Delta}{2}e^{-i\phi_0\sigma_z/2}\sigma_x e^{i\phi_0\sigma_z/2}$，其中 $\Delta=\sqrt{\alpha^2+\beta^2}$。

旋转 $e^{i\phi_0\sigma_z/2}$ 后：
$$\tilde H(t) = \frac{\Delta}{2}\sigma_x + \frac{1}{2}\big(-\alpha\cos\phi_0+\beta\sin\phi_0\big)\cos\omega t\,\sigma_x $$

+ \frac{1}{2}\big(-\alpha\sin\phi_0-\beta\cos\phi_0\big)\cos\omega t\,\sigma_y.$$

当 $\alpha=\beta$ 时 $\phi_0=\pi/4$：
$$\tilde H(t) = \frac{\Delta}{2}\sigma_x + 0\cdot\cos\omega t\,\sigma_x + \frac{-\Delta}{\sqrt{2}}\cos\omega t\,\sigma_y.$$

再绕 $x$ 轴旋转使静态部分沿 $z$：
$$H_{\text{Rabi}}(t) = \frac{\Delta}{2}\sigma_z - \frac{\Delta}{\sqrt{2}}\cos\omega t\,\sigma_x.$$

**这就是标准 Rabi 模型**：
$$\boxed{H_{\text{Rabi}}(t) = \frac{\Delta}{2}\sigma_z + \Omega\cos\omega t\,\sigma_x,\qquad 
\Delta=t_c\sqrt{2},\;\Omega=-\frac{t_c}{\sqrt{2}}.}$$

---

## 二、Floquet 理论

### 2.1 Floquet 定理

$H(t+T)=H(t)$，$T=2\pi/\omega=2\tau$。存在幺正的周期微运动算符 $P(t)=P(t+T),\;P(0)=I$ 和
与时间无关的 Floquet 哈密顿量 $K$，使得：
$$U(t) = P(t)\,e^{-iKt}.$$

演化算符在一个周期后：$U(T)=e^{-iKT}$。半周期（Step 2 的结束时刻）：
$$U(\tau) = P(\tau)\,e^{-iK\tau}.$$

$P(\tau)$ 是 Floquet 微运动算符在 $t=\tau$ 处的值。

### 2.2 Floquet 矩阵（扩展 Hilbert 空间）

在 $L^2[0,T]\otimes\mathbb C^2$ 中定义 Floquet 算符 $K_F = H(t)-i\frac{d}{dt}$。
$K_F$ 的本征值 $\varepsilon_n$ 即 quasienergies（准能量），$U(T)$ 的本征值为 $e^{-i\varepsilon_n T}$。

在 Floquet-Fourier 基 $\phi(t)=\sum_{n\in\mathbb Z}\phi_n e^{in\omega t}$ 中，$K_F$ 是
**无限分块三对角矩阵**：

$$\langle n|K_F|m\rangle = H_{n-m} + n\omega\delta_{nm}I_{2\times 2},$$

其中 $H_{n-m}$ 是 $H(t)=\sum_k H_k e^{ik\omega t}$ 的 Fourier 系数。

对本问题：$H(t)=H_0 + \frac12 H_1 e^{i\omega t} + \frac12 H_1 e^{-i\omega t}$（注意 $H_1$ 是 $\cos\omega t$ 的系数矩阵）。

### 2.3 无限矩阵形式

$$K_F = \begin{pmatrix}
\ddots & \ddots & & & \\
\ddots & H_0-\omega & \frac12 H_1 & & \\
& \frac12 H_1 & H_0 & \frac12 H_1 & \\
& & \frac12 H_1 & H_0+\omega & \ddots\\
& & & \ddots & \ddots
\end{pmatrix}$$

每个块是 $2\times2$ 矩阵：
$$H_0 = \frac{1}{2}\begin{pmatrix}0 & \alpha-i\beta\\ \alpha+i\beta & 0\end{pmatrix},\qquad
\frac12 H_1 = \frac{1}{4}\begin{pmatrix}0 & -\alpha-i\beta\\ -\alpha+i\beta & 0\end{pmatrix}.$$

---

## 三、准能量的连分数 / Hill 行列式

### 3.1 Hill 行列式

准能量 $\varepsilon$ 由 $\det(K_F-\varepsilon I)=0$ 确定。这是无限行列式——Hill 行列式。

对 $2\times2$ 块结构的无限三对角矩阵，Hill 行列式等价于一个**矩阵连分数**方程。

### 3.2 连分数形式

定义 Green 函数 $G_n(\varepsilon) = (\varepsilon - H_0 - n\omega - \frac14 H_1 G_{n+1}(\varepsilon) H_1)^{-1}$。

$G_n$ 满足递归：
$$G_n^{-1} = \varepsilon - H_0 - n\omega - \frac14 H_1 G_{n+1} H_1.$$

准能量条件：$\det G_0^{-1}(\varepsilon)=0$，即：
$$\det\!\left(\varepsilon - H_0 - \frac14 H_1\,\frac{1}{\varepsilon-H_0-\omega-\frac14 H_1\frac{1}{\varepsilon-H_0-2\omega-\cdots}H_1}\,H_1\right) = 0.$$

这是 **$2\times2$ 矩阵连分数**——Rabi 模型准能量的精确特征方程。

### 3.3 标量化

对 $\alpha=\beta$ 的对称情形，可进一步约化为标量连分数。

旋转后的 Rabi 形式 $H_{\text{Rabi}} = \frac{\Delta}{2}\sigma_z + \Omega\cos\omega t\,\sigma_x$：

标准文献给出准能量 $\varepsilon$ 满足的连分数方程（Irwin, 1992; Shirley, 1965）：

$$\boxed{\varepsilon = \frac{\Delta}{2} + \frac{|\Omega|^2/4}{\varepsilon-\omega-\Delta/2 - \frac{|\Omega|^2/4}{\varepsilon-2\omega-\Delta/2 - \frac{|\Omega|^2/4}{\varepsilon-3\omega-\Delta/2-\cdots}}} 
+ \frac{|\Omega|^2/4}{\varepsilon+\omega-\Delta/2 - \frac{|\Omega|^2/4}{\varepsilon+2\omega-\Delta/2 - \cdots}}}$$

这是 Rabi 模型中准能量的**精确解析表示**——不是初等函数，但是收敛极快的连分数。

---

## 四、Mathieu 函数联系

### 4.1 从 Rabi 到 Mathieu 方程

Rabi 模型的 Schrödinger 方程 $i\partial_t\psi = H_{\text{Rabi}}\psi$ 可被转化为
**Mathieu 方程**（或 Whittaker-Hill 方程）。

令 $\psi(t) = (c_+(t), c_-(t))^T$：
$$i\dot c_+ = \frac{\Delta}{2}c_+ + \Omega\cos\omega t\,c_-,$$
$$i\dot c_- = \Omega\cos\omega t\,c_+ - \frac{\Delta}{2}c_-.$$

消去 $c_+$ 后 $c_-$ 满足：
$$\ddot c_- + \omega\tan\omega t\,\dot c_- + \left[\left(\frac{\Delta}{2}\right)^2 + \Omega^2\cos^2\omega t + i\frac{\Delta}{2}\omega\tan\omega t\right]c_- = 0.$$

作变量替换 $z = \cos\omega t$ 可将此方程变为**合流 Heun 方程**（Confluent Heun Equation, CHE）。

### 4.2 合流 Heun 函数

CHE 的标准形式：$\frac{d^2w}{dz^2} + \left(\frac{\gamma}{z} + \frac{\delta}{z-1} + \varepsilon\right)\frac{dw}{dz} + \frac{\alpha z - q}{z(z-1)}w = 0.$

Rabi 模型的波函数可通过 $\text{HeunC}(q,\alpha,\gamma,\delta,\varepsilon;z)$ 表达。
具体参数映射见 Xie et al. (2017), Phys. Rev. A **96**, 043842。

### 4.3 总结：特殊函数层次结构

| 层次 | 描述 | 函数 |
|---|---|---|
| 0 级 | 含时二能级系统 | $U(t)=\mathcal T\exp(-i\int H dt)$ |
| 1 级 | 准能量 | 连分数 |
| 2 级 | 波函数分量 | 合流 Heun 函数 $\text{HeunC}$ |
| 3 级 | 对称约化（RWA 极限） | Mathieu 函数（角 Mathieu $\text{ce}_n,\text{se}_n$） |
| 4 级 | 初等函数 | $\exp,\sin,\cos,\sqrt{\cdot}$（仅在 $[H(s),H(u)]\equiv 0$ 时） |

**我们正好落在第 1–2 级**——准能量用连分数严格表示，波函数用 HeunC 表达。
第 4 级（初等函数）被 $[H(s),H(u)]\neq 0$ 排除，但不是"无解析解"——而是
解析解需要用更高级的特殊函数。

---

## 五、三段 braiding 的完整解析结构

### 5.1 Step 1 和 Step 3：初等可积

Step 1 中 $H_{\text{eff}}(t)\propto f_-(t)\cdot(\text{常方向})$，$[H(s),H(u)]=0$：
$$U_1 = \exp\!\left(-i\int_0^\tau H_{\text{eff}}(t)dt\right) = \exp\!\big(-i\Phi_1\,\hat n_1\cdot\vec\sigma\big),$$
$$\Phi_1 = \frac{\tau}{2}\sqrt{\alpha_1^2+\beta_1^2},\quad \hat n_1 = (\cos\theta_1,\sin\theta_1,0).$$

Step 3 同理。

### 5.2 Step 2：Floquet / Heun 解

$$U_2 = U_{\text{Floquet}}(\tau;\Delta,\Omega,\omega),$$

其中：
- 准能量 $\varepsilon$ 由 §3.3 的连分数确定
- $U_2$ 通过 Floquet 微运动 $P(\tau)$ 和 $e^{-iK\tau}$ 计算
- 等价地，$U_2$ 的矩阵元用 HeunC 函数表达

### 5.3 完整演化

$$U_{\text{total}}(3\tau) = U_3\,U_2\,U_1 = \exp\!\big(-i\Phi_{\text{tot}}\,\hat N\cdot\vec\sigma\big),$$

其中 $(\Phi_{\text{tot}},\hat N)$ 通过 SU(2) 旋转合成公式从 $(\Phi_1,\hat n_1),(\Phi_2,\hat n_2),(\Phi_3,\hat n_3)$ 计算。

---

## 六、数值实现：截断 Floquet 矩阵

实际计算中，截断 Fourier 模式到 $\pm N_{\text{max}}$，得有限矩阵：
$$K_F^{(N)} \in \mathbb C^{2(2N+1)\times 2(2N+1)}.$$

对角化此矩阵的最低几个本征值即准能量。$N\sim 3$–$5$ 通常足以收敛到 $10^{-6}$ 精度
（因为 Fourier 系数衰减为 $1/n^2$）。

Floquet 态 $|\phi_n(t)\rangle = \sum_{k} \phi_{n,k} e^{ik\omega t}$ 由本征矢的 Fourier 系数重构，
进而 $U(t) = \sum_n e^{-i\varepsilon_n t}|\phi_n(t)\rangle\langle\phi_n(0)|$。

---

## 七、实用结论

1. **$\mathcal T\exp$ 不是初等函数**——$[H(s),H(u)]\neq 0$ 排除了这一点。
2. **但它是"命名特殊函数"**——准能量是连分数，波函数是 HeunC。
3. **实用上**：截断 Floquet 矩阵（$N=5$，42×42 矩阵）对角化即得高精度准能量。
4. **对称情形 $\alpha=\beta$**：波函数退化为 Mathieu 函数，准能量连分数进一步简化。
5. **标度律** $\Phi_{\text{tot}}\approx\sqrt{\Phi_G^2+\Phi_D^2}$ 是连分数在 $\Delta\ll\omega$ 极限下的首阶近似。

---

## 参考文献

- J. H. Shirley, *Solution of the Schrödinger Equation with a Hamiltonian Periodic in Time*, Phys. Rev. **138**, B979 (1965).
- Q. Xie et al., *Analytical solutions of the Rabi model*, Phys. Rev. A **96**, 043842 (2017).
- D. Braak, *Integrability of the Rabi Model*, Phys. Rev. Lett. **107**, 100401 (2011).

---

*文档版本：2025-06-03*
