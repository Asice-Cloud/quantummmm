# PRB 113 Supplementary Information 详细解读

## 文献信息

- **标题**: Supplementary information for "Adiabatic nonabelian braiding of imperfect Majoranas"
- **作者**: Maximilian Nitsch, Viktor Svensson, William Samuelson, Konstantin Nestmann, Jeroen Danon, Karsten Flensberg, Rubén Seoane Souto, Martin Leijnse
- **日期**: March 26, 2026
- **出处**: Physical Review B 113

## 核心问题

在有 MBS（Majorana Bound States）空间重叠（imperfect Majoranas）的情况下，如何实现绝热非阿贝尔编织（braiding），并获得解析解。

---

## 第 I 节：编织协议的实现——基于最小 Kitaev 链

### S1. 单个最小 Kitaev 链

一个最小 Kitaev 链由两个**自旋极化量子点 (QD)** 组成 ($l = 1, 2$)，各持有单个费米子态，通过窄超导体耦合：

$$H = \sum_{l=1,2} \xi_l n_l + t_{\text{cot}} c_1^\dagger c_2 + \Delta_{\text{car}} c_1^\dagger c_2^\dagger + H.c. \tag{S1}$$

- $\xi_l$：量子点的单粒子能量
- $t_{\text{cot}}$：弹性共隧穿振幅 (elastic co-tunneling)
- $\Delta_{\text{car}}$：交叉 Andreev 反射振幅 (crossed Andreev reflection)

当 $t_{\text{cot}} = \Delta_{\text{car}}$ 且 $\xi_1, \xi_2 = 0$ 时，系统有一个**简并基态**，描述为两个 MBSs $\gamma, \tilde{\gamma}$，各自局域在一个 QD 上。

偏离 $t_{\text{cot}} = \Delta_{\text{car}}$ 会导致 MBS 耦合：

$$H_c = \epsilon i\gamma\tilde{\gamma}, \quad \epsilon = \Delta_{\text{car}} - t_{\text{cot}} \tag{S2}$$

### S2. 编织装置

连接三个最小 Kitaev 链构成编织装置。内层 QD 通过隧道耦合连接：

$$H = \sum_{k=1}^{3} \epsilon_k i\gamma_k \tilde{\gamma}_k + t_{12} i\gamma_1 \gamma_2 + t_{13} i\gamma_1 \gamma_3 \tag{S3}$$

通过超导回路中的磁通 $\Phi_{12}, \Phi_{13}$ 使 $t_{12}, t_{13}$ 为实数。

编织协议通过时变参数 $\epsilon_1(t), t_{12}(t), t_{13}(t)$ 实现，保持 $\epsilon_2 = \epsilon_3 = 0$。

### S3. 空间重叠的引入（关键！）

**"非完美"** 体现在外 MBS $\tilde{\gamma}_k$ 部分位于内层 QD 上，与 $\gamma_k$ 有空间重叠。这通过变换描述：

$$\gamma_k \to \frac{1}{\sqrt{1+\zeta^2}} \gamma_k + i\frac{\zeta}{\sqrt{1+\zeta^2}} \tilde{\gamma}_k \tag{S4}$$

- $\zeta = 0$：完美局域化 MBS
- $\zeta = 1$：平庸态（trivial state）
- MBS 极化: $M = \frac{1-\zeta^2}{1+\zeta^2}$

> ⚠️ 文中假设 $\zeta$ 是常数系统参数，与 $\epsilon_k$ 独立。

### S4. 变换后的哈密顿量

变换使系统哈密顿量变为：

$$H \to H + \tilde{H} \tag{S5}$$

$$H = \frac{\epsilon_1}{1+\zeta^2} i\gamma_1 \tilde{\gamma}_1 + \frac{1}{\sqrt{1+\zeta^2}} t_{12} i\gamma_1 \gamma_2 + \frac{1}{\sqrt{1+\zeta^2}} t_{13} i\gamma_1 \gamma_3 \tag{S6}$$

$$\tilde{H} = \frac{\zeta^2}{1+\zeta^2} t_{12} i\tilde{\gamma}_1 \tilde{\gamma}_2 + \frac{\zeta^2}{1+\zeta^2} t_{13} i\tilde{\gamma}_1 \tilde{\gamma}_3$$

$\tilde{H}$ 的效应是**解除基态简并** ($\delta E \neq 0$)。如果协议相对于 $\delta E$ 是绝热的，结果是平凡的——系统保持在原始态。

### S5. 修正项恢复简并

通过调节左右 MBS 对的 $t_{\text{cot}}$ 和 $\Delta_{\text{car}}$ 来引入修正：

$$H_{\text{corr}} = \Lambda i\gamma_2 \tilde{\gamma}_2 + \Lambda i\gamma_3 \tilde{\gamma}_3 \tag{S7}$$

$\Lambda(t)$ 是时变参数，需要与协议同步调谐。为简化解析解，选择两个系统相同的 overlap。

全哈密顿量：

$$H_{\text{tot}} = H + \tilde{H} + H_{\text{corr}} \tag{S8}$$

### S6. 无量纲化

基态扇区与激发态之间的能隙：

$$\Delta = \sqrt{\epsilon_1^2 + \frac{1}{(1+\zeta^2)^2}(t_{12}^2 + t_{13}^2)} \tag{S9}$$

定义无量纲参数：

$$\rho_1 = \epsilon_1/\Delta, \quad \rho_{k=2,3} = \frac{1}{\sqrt{1+\zeta^2}} t_{1k}/\Delta, \quad \lambda = \Lambda/\Delta, \quad \eta = \zeta^2 \tag{S10}$$

得到**主论文 Eq.(1)** 的无量纲哈密顿量：

$$h = H_{\text{tot}}/\Delta = i\gamma_1 \gamma_\Delta + \eta(\rho_2 i\tilde{\gamma}_1 \tilde{\gamma}_2 + \rho_3 i\tilde{\gamma}_1 \tilde{\gamma}_3) + \lambda(i\gamma_2 \tilde{\gamma}_2 + i\gamma_3 \tilde{\gamma}_3) \tag{S11}$$

其中 $\gamma_\Delta = \rho_1 \tilde{\gamma}_1 + \rho_2 \gamma_2 + \rho_3 \gamma_3$，满足 $\sum \rho_i^2 = 1$。

---

## 第 II 节：哈密顿量对角化

### 步骤 1：旋转 $\tilde{\gamma}_1, \gamma_2, \gamma_3$ 基

引入 $\gamma_\Delta$：

$$\gamma_\Delta = \cos\theta \,\tilde{\gamma}_1 + \sin\theta\cos\phi \,\gamma_2 + \sin\theta\sin\phi \,\gamma_3 \tag{S12}$$

角度参数化：

$$\rho_1 = \cos\theta,\quad \rho_2 = \sin\theta\cos\phi,\quad \rho_3 = \sin\theta\sin\phi \tag{S13}$$

引入两个与 $\gamma_\Delta$ 正交的 MBSs：

$$\gamma_{\theta'} = -\sin\theta \,\tilde{\gamma}_1 + \cos\theta\cos\phi \,\gamma_2 + \cos\theta\sin\phi \,\gamma_3 \tag{S14}$$
$$\gamma_{\phi'} = -\sin\phi \,\gamma_2 + \cos\phi \,\gamma_3$$

### 步骤 2：旋转 $\tilde{\gamma}_2, \tilde{\gamma}_3$ 基

$$\gamma_\eta = \cos\phi \,\tilde{\gamma}_2 + \sin\phi \,\tilde{\gamma}_3 \tag{S15}$$
$$\gamma_{\eta'} = -\sin\phi \,\tilde{\gamma}_2 + \cos\phi \,\tilde{\gamma}_3$$

$\gamma_1$ 不旋转。此时哈密顿量变为：

$$h = i\gamma_1 \gamma_\Delta + i\gamma_\Delta \gamma_\eta (\lambda\sin\theta + \eta\cos\theta\sin\theta) + i\gamma_\eta \gamma_{\theta'} (\eta\sin^2\theta - \lambda\cos\theta) + i\gamma_{\phi'} \gamma_{\eta'} \lambda \tag{S16}$$

**已实现**: $\gamma_{\phi'}, \gamma_{\eta'}$ 对与其他 MBS 解耦。

### 步骤 3：通过 ansatz 对角化剩余部分

$$\gamma_1^D = \alpha\gamma_1 + \beta\gamma_\eta,\quad \gamma_\eta^D = \alpha\gamma_\eta - \beta\gamma_1 \tag{S17}$$
$$\gamma_\Delta^D = \mu\gamma_\Delta + \nu\gamma_{\theta'},\quad \gamma_{\theta'}^D = \mu\gamma_{\theta'} - \nu\gamma_\Delta$$

定义缩写：

$$\tilde{\lambda} = \lambda\sin\theta + \eta\cos\theta\sin\theta,\quad \tilde{\eta} = \eta\sin^2\theta - \lambda\cos\theta \tag{S19}$$

代入后哈密顿量变为：

$$h = i\gamma_1^D \gamma_\Delta^D (\alpha\mu - \tilde{\lambda}\beta\mu + \tilde{\eta}\beta\nu) + i\gamma_\eta^D \gamma_{\theta'}^D (\beta\nu + \tilde{\lambda}\alpha\nu + \tilde{\eta}\alpha\mu)$$
$$+ i\gamma_\Delta^D \gamma_\eta^D (\mu\beta + \tilde{\lambda}\alpha\mu - \tilde{\nu}\alpha\nu) + i\gamma_1^D \gamma_{\theta'}^D (-\alpha\nu + \tilde{\lambda}\beta\nu + \tilde{\nu}\beta\mu) + i\gamma_{\phi'} \gamma_{\eta'} \lambda \tag{S18}$$

### 解耦条件

要求 $\gamma_1^D, \gamma_\Delta^D$ 与 $\gamma_\eta^D, \gamma_{\theta'}^D$ 解耦：

$$\mu\beta + \tilde{\lambda}\alpha\mu - \tilde{\nu}\alpha\nu = 0 \tag{S20}$$
$$-\alpha\nu + \tilde{\lambda}\beta\nu + \tilde{\nu}\beta\mu = 0$$

引入角度 $\theta_\alpha, \theta_\mu$ 满足该条件：

$$\theta_\mu = -\frac{1}{2}\arctan\left(\frac{2\tilde{\lambda}\tilde{\eta}}{1+\tilde{\lambda}^2-\tilde{\eta}^2}\right) \tag{S21}$$
$$\theta_\alpha = -\arctan(-\tilde{\eta}\tan\theta_\mu + \tilde{\lambda})$$

$$\alpha = \cos\theta_\alpha,\quad \beta = \sin\theta_\alpha,\quad \mu = \cos\theta_\mu,\quad \nu = \sin\theta_\mu \tag{S22}$$

### 最终对角化形式

$$h = \tilde{\Delta} \, i\gamma_1^D \gamma_\Delta^D + \varepsilon \, i\gamma_\eta^D \gamma_{\theta'}^D + \lambda \, i\gamma_{\phi'} \gamma_{\eta'} \tag{S23}$$

其中：

$$\tilde{\Delta} = \alpha\mu - \tilde{\lambda}\beta\mu + \tilde{\eta}\beta\nu > 0 \tag{S24}$$
$$\varepsilon = \beta\nu + \tilde{\lambda}\alpha\nu + \tilde{\eta}\alpha\mu = \tilde{\eta}\alpha/\mu \tag{S26}$$

**低能扇区**: $i\gamma_1^D \gamma_\Delta^D = -1$（在整个协议中保持不变）。

---

## 第 III 节：协议对宇称扇区的依赖

### 总宇称守恒

MBS 基变换关系：

$$\gamma_1 \to \gamma_1^D,\quad \tilde{\gamma}_1 \to \gamma_\Delta^D,\quad \gamma_2 \to \gamma_{\theta'}^D,\quad \tilde{\gamma}_2 \to \gamma_\eta^D,\quad \gamma_3 \to \gamma_{\phi'},\quad \tilde{\gamma}_3 \to \gamma_{\eta'} \tag{S27}$$

总宇称 $\sigma = \pm 1$ 在任意协议点：

$$\sigma = i\gamma_1\tilde{\gamma}_1 \cdot i\gamma_2\tilde{\gamma}_2 \cdot i\gamma_3\tilde{\gamma}_3 = i\gamma_1^D \gamma_\Delta^D \cdot i\gamma_{\theta'}^D \gamma_\eta^D \cdot i\gamma_{\phi'} \gamma_{\eta'} \tag{S28}$$

利用 $i\gamma_1^D \gamma_\Delta^D = -1$：

$$i\gamma_{\phi'} \gamma_{\eta'} = \sigma \cdot i\gamma_\eta^D \gamma_{\theta'}^D \tag{S29}$$

代入哈密顿量：

$$h = \tilde{\Delta} \, i\gamma_1^D \gamma_\Delta^D + (\varepsilon + \sigma \cdot \lambda) \, i\gamma_\eta^D \gamma_{\theta'}^D \tag{S30}$$

### 简并条件（主论文 Eq.11）

$$\varepsilon(\eta, \lambda) + \sigma\lambda = 0 \tag{S31}$$

### 费米子极限分析 ($\eta = 1, \theta \ll 1$)

在协议的起点和终点处展开：

$$\varepsilon(\theta, \lambda) = -\lambda - \theta^2\frac{1}{(-1+\lambda)} + O(\theta^3) \tag{S32}$$

$$\Rightarrow 0 = (-1+\sigma)\lambda - \theta^2\frac{1}{(-1+\lambda)} + O(\theta^3)$$

- **奇宇称 ($\sigma = -1$)**: 要求 $\lambda = 0 + O(\theta^2)$，只需小幅修正
- **偶宇称 ($\sigma = 1$)**: 要求 $\lambda = 1 + O(\theta)$，需要大幅修正

偶宇称下：

$$E_{\text{even, low}}/\Delta = -1 + O(\theta) \tag{S33}$$
$$E_{\text{even, high}}/\Delta = 1 \pm (\varepsilon - \lambda) = 1 \pm 2 + O(\theta) \tag{S34}$$

> 🔑 **结论**: 修正项将高能态拉到 $E_{\text{even, high}} = -1$，对 $\eta = 1$ 完全闭合低能扇区能隙。对 $\eta < 1$，能隙保持有限但减小。因此存在两个不同协议，取决于总宇称。

### 符号选择的影响

- $\rho_2 \to -\rho_2$ 等效于 $(\gamma_2, \tilde{\gamma}_2) \to (-\gamma_2, -\tilde{\gamma}_2)$，**不改变物理**
- $\rho_3 \to -\rho_3$、$\eta$ 的符号也不影响
- $\lambda$ 的符号不重要，因其值由优化简并决定
- **$\rho_1 \to -\rho_1$** 会交换奇偶宇称，协议表现不同
- $\tilde{\gamma}_3 \to -\tilde{\gamma}_3$ 也会交换奇偶宇称

---

## 第 IV 节：两组 MBS 耦合的对角化（理论铺垫）

哈密顿量可写为：

$$H = i\vec{\gamma}_L^T M \vec{\gamma}_R \tag{S38}$$

其中 $M$ 是实矩阵。对 $M$ 做 SVD（奇异值分解）：$M = U^T S V$，$U, V$ 为实正交矩阵，$S$ 为对角矩阵。实正交矩阵保持正则反对易关系，因此：

$$\vec{\chi}_L = U\vec{\gamma}_L,\quad \vec{\chi}_R = V\vec{\gamma}_R \tag{S39-40}$$

是合法 MBS 基。在此基下：

$$H = i\vec{\chi}_L^T S \vec{\chi}_R \tag{S41}$$

**两组不混合**，绝热演化时各自独立演化。

---

## 第 V 节：Berry 相位的计算 ⭐

### 投影算符

$$P = \frac{1}{2}(1 - i\gamma_1^D \gamma_\Delta^D) \tag{S43}$$

Berry 相位（主论文 Eq.12）：

$$U(\eta) = \mathcal{T} \exp\left[-\oint [P, \partial_\omega P] d\omega\right] \tag{S42}$$

### 分三步计算

协议参数化分为三段：

$$\oint d\omega = \underbrace{\int_{\pi/2}^{0} d\theta}_{\text{Step 1: }\tau_{12}} + \underbrace{\int_{0}^{\pi/2} d\phi}_{\text{Step 2: }\tau_{23}} + \underbrace{\int_{0}^{\pi/2} d\theta}_{\text{Step 3: }\tau_{31}} \tag{S44}$$

---

### 第一步 $\tau_{12}$ ($\phi = 0$, $\theta: \pi/2 \to 0$)

$$U_{12} = \mathcal{T} \exp\left[-\int_0^{\pi/2} [P, \partial_\theta P] d\theta\right] \tag{S45}$$

计算 $\partial_\theta P$：

$$\partial_\theta P = -\frac{i}{2} \partial_\theta(\gamma_1^D \gamma_\Delta^D) = -\frac{i}{2}(\gamma_\eta^D \gamma_\Delta^D \cdot \partial_\theta\theta_\alpha + \gamma_1^D \gamma_{\theta'}^D \cdot (1 + \partial_\theta\theta_\mu)) \tag{S46-47}$$

在 $\phi = 0$ 时简化：

$$\gamma_1^D \gamma_\eta^D = \gamma_1 \tilde{\gamma}_2,\quad \gamma_\Delta^D \gamma_{\theta'}^D = \tilde{\gamma}_1 \gamma_2 \tag{S49}$$

$$[P, \partial_\theta P] = \frac{1}{2}(\gamma_1 \tilde{\gamma}_2 \cdot \partial_\theta\theta_\alpha + \tilde{\gamma}_1 \gamma_2 \cdot (1 + \partial_\theta\theta_\mu)) \tag{S50}$$

> 🔑 **关键简化**: 被积函数由两个作用在不同 MBS 对上的常量矩阵组成，被积函数与自身对易 — 可以**忽略时间排序** $\mathcal{T}$！

积分结果：

$$-\int_0^{\pi/2} [P, \partial_\theta P] d\theta = \tilde{\gamma}_2 \gamma_1 \,\theta_{\alpha,23}/2 + \gamma_2 \tilde{\gamma}_1 \,(\pi/4 + \theta_{\mu,23}/2) \tag{S51}$$

其中定义了**在 $\theta = \pi/2$ 处**的常数角：

$$\theta_{\alpha,23} \equiv \theta_\alpha(\theta = \pi/2, \phi),\quad \theta_{\mu,23} \equiv \theta_\mu(\theta = \pi/2, \phi) \tag{S52}$$

> 注意 $\theta_\alpha, \theta_\mu$ 在第二协议步中确实为常数！

$$U_{12} = \exp\left(\tilde{\gamma}_2 \gamma_1 \,\theta_{\alpha,23}/2\right) \exp\left(\gamma_2 \tilde{\gamma}_1 \,(\pi/4 + \theta_{\mu,23}/2)\right) \tag{S53}$$

### 第三步 $\tau_{31}$

利用对称性，交换 $2 \leftrightarrow 3$ 并取厄米共轭以反转时间排序：

$$U_{31} = \exp\left(\gamma_1 \tilde{\gamma}_3 \,\theta_{\alpha,23}/2\right) \exp\left(\tilde{\gamma}_1 \gamma_3 \,(\pi/4 + \theta_{\mu,23}/2)\right) \tag{S54}$$

---

### 第二步 $\tau_{23}$ ($\theta = \pi/2$, $\phi: 0 \to \pi/2$)

这是**最困难的部分**，因为 $\gamma_1 \gamma_{\eta'}$ 和 $\tilde{\gamma}_1 \gamma_{\phi'}$ 在积分期间会演化，不是常数算符。

由对角化基定义可得：

$$\partial_\phi \gamma_1^D = \beta_{23} \gamma_{\eta'},\quad \partial_\phi \gamma_\Delta^D = \mu_{23} \gamma_{\phi'} \tag{S58}$$

其中 $\mu_{23} \equiv \mu(\theta=\pi/2,\phi)$ 等。

$$[P, \partial_\phi P] = \frac{1}{2}(\mu^2 \gamma_2 \gamma_3 + \beta^2 \tilde{\gamma}_2 \tilde{\gamma}_3 + \alpha\beta \gamma_1 \gamma_{\eta'} - \mu\nu \tilde{\gamma}_1 \gamma_{\phi'}) \tag{S59}$$

**关键技巧**: 拆分为两组对易的 MBS — $(\tilde{\gamma}_1, \gamma_2, \gamma_3)$ 和 $(\gamma_1, \tilde{\gamma}_2, \tilde{\gamma}_3)$，分别求解后再合并。

用泡利矩阵表示第一组：

$$\gamma_3 \tilde{\gamma}_1 = i\sigma_x,\quad \gamma_2 \tilde{\gamma}_1 = -i\sigma_y,\quad \gamma_2 \gamma_3 = i\sigma_z \tag{S63}$$

化为微分方程：

$$\partial_\omega U_{23}^{\mu\nu}(\omega) = \frac{i}{2}(-\mu^2 \sigma_z - \mu\nu(\sin\omega\,\sigma_y + \cos\omega\,\sigma_x)) U_{23}^{\mu\nu}(\omega) \tag{S65}$$

初始条件: $U_{23}^{\mu\nu}(\omega=0) = \mathbb{1}_2$

解出：

$$U_{23}^{\mu\nu} = e^{-i\pi/4\,\sigma_z} \, e^{-i\pi/4\,\nu(\mu\sigma_y - \nu\sigma_z)} \tag{S69}$$

转回 MBS 语言：

$$U_{23}^{\mu\nu} = e^{\pi/4 \gamma_3 \gamma_2} \, e^{-\pi/4 \nu(\mu\gamma_3 \tilde{\gamma}_1 + \nu\gamma_3 \gamma_2)} \tag{S70}$$

第二组 ($\alpha\beta$ 部分) 类似求解。

### 完整编织结果

三部分相乘：

$$U(\eta) = U_{31} U_{23} U_{12} = e^{\frac{\pi}{4}(1+\nu_{23})\gamma_3 \gamma_2} \, e^{\frac{\pi}{4}(1-\alpha_{23})\tilde{\gamma}_3 \tilde{\gamma}_2} \tag{S72}$$

其中用到了 $\mu_{23} = \cos\theta_{\mu,23},\;\nu_{23} = \sin\theta_{\mu,23},\;\alpha_{23} = \cos\theta_{\alpha,23},\;\beta_{23} = \sin\theta_{\alpha,23}$。

---

## 第 VI 节：第二协议区间简并条件的解

在 $\theta = \pi/2$（即 $\tau_{23}$ 区间）时，简并条件 (S31) 简化为：

$$\frac{\alpha_{23}}{\mu_{23}} = -\sigma\frac{\lambda_{23}}{\eta} \tag{S74}$$

结合 (S21, S22) 解出具体的参数值：

$$\lambda_{23} = -\sigma\frac{\eta}{\sqrt{1+\eta^2}}$$

$$\nu_{23} = \frac{\sigma\eta^2}{\sqrt{1+\eta^2+\eta^4}},\quad \alpha_{23} = \frac{1}{\sqrt{1+\eta^2+\eta^4}} \tag{S76}$$

$$\mu_{23} = \frac{\sqrt{1+\eta^2}}{\sqrt{1+\eta^2+\eta^4}},\quad \beta_{23} = \frac{\sigma\eta\sqrt{1+\eta^2}}{\sqrt{1+\eta^2+\eta^4}}$$

最终得到主论文 Eq.(13) 的角度：

$$\phi = \frac{\pi}{4}(1+\nu_{23}) = \frac{\pi}{4}\left(1 + \frac{\sigma\eta^2}{\sqrt{1+\eta^2+\eta^4}}\right) \tag{S77}$$

$$\tilde{\phi} = \frac{\pi}{4}(1-\alpha_{23}) = \frac{\pi}{4}\left(1 - \frac{1}{\sqrt{1+\eta^2+\eta^4}}\right)$$

---

## 第 VII 节：宇称测量与 MBS 相似度的函数关系

### MBS 相似度（主论文 Eq.15）

$$S(U) = \text{Tr}[U^2(\eta)^\dagger U^2(0)]/4$$

推导过程：

$$S(U) = \left|\text{Tr}\left[\exp(-2(\phi - \tilde{\phi})\gamma_3\gamma_2) \exp(2\frac{\pi}{4}\gamma_3\gamma_2)\right]\right|^2/4$$
$$= \left|\text{Tr}\left[\exp\left(2\left(\frac{\pi}{4} - \phi + \tilde{\phi}\right)\gamma_3\gamma_2\right)\right]\right|^2/4$$
$$= \cos^2\left(2\left(\frac{\pi}{4} - \phi + \tilde{\phi}\right)\right)$$
$$= \sin^2(2(\phi - \tilde{\phi})) \tag{S78}$$

### 宇称测量

系统初始化为奇宇称 ($\sigma = -1$) 的基态扇区内的态 $|\psi\rangle$：

$$i\gamma_2\tilde{\gamma}_2|\psi\rangle = i\gamma_3\tilde{\gamma}_3|\psi\rangle = \pm|\psi\rangle \tag{S79}$$

计算 $U^2(\eta)$ 对宇称算符的变换：

$$U^2(\eta) i\gamma_2\tilde{\gamma}_2 U^2(\eta)^\dagger = i(\cos(4\phi)\gamma_2 + \sin(4\phi)\gamma_3)(\cos(4\tilde{\phi})\tilde{\gamma}_2 + \sin(4\tilde{\phi})\tilde{\gamma}_3) \tag{S80}$$

测量期望值（利用 $\langle\psi|i\gamma_2\tilde{\gamma}_3|\psi\rangle = \langle\psi|i\gamma_3\tilde{\gamma}_2|\psi\rangle = 0$）：

$$\langle\psi|U^2(\eta) i\gamma_2\tilde{\gamma}_2 U^2(\eta)^\dagger|\psi\rangle = \pm(1 - 2\sin^2(2(\phi - \tilde{\phi}))) = \pm(1 - 2S(U)) \tag{S81}$$

> 🎯 **核心结论**: 宇称测量的期望值直接与 MBS 相似度 $S(U)$ 线性相关！可以通过测量宇称来推断编织保真度。

---

## 第 VIII 节：数值模拟

求解含时薛定谔方程：

$$\partial_t U(\eta, t) = iH(t)U(\eta, t),\quad U(\eta, 0) = \mathbb{1} \tag{S82}$$

耦合参数使用光滑阶跃函数：

$$f(t) = 1/2 + \tanh(k\cos(2\pi t/T))/2 \tag{S83}$$

$$\rho_1 = f(t/T)/N,\quad \rho_2 = f(t/T - 1/3)/N,\quad \rho_3 = f(t/T - 2/3)/N \tag{S84}$$

$N$ 是归一化因子确保 $\sum_i \rho_i^2 = 1$。$k=10$ 控制阶跃的锐度。$\lambda$ 在每个时间步通过数值求根算法解 Eq.(S31) 来确定。

### 数值结果（图 S1）

- **修正后的协议**消除了动力学相位误差
- 短时间协议 ($T\Delta = 60 \pm 12$)：有小的 diabatic 误差，但修正后仍明显优于未修正
- 长时间协议 ($T\Delta = 200 \pm 40$)：修正后的数值结果与解析解几乎完美匹配
- 未修正协议有系统性的 $S(U)$ 偏离

---

## 第 IX 节：非对称 MBS 重叠

推广到三个链有不同的重叠参数 $\eta_1 \neq \eta_2 \neq \eta_3$：

$$H/\Delta = \rho_1 i\gamma_1\tilde{\gamma}_1 + \rho_2 i\gamma_1\gamma_2 + \rho_3 i\gamma_1\gamma_3$$
$$+ \rho_2\sqrt{\eta_1\eta_2}\,i\tilde{\gamma}_1\tilde{\gamma}_2 + \rho_3\sqrt{\eta_1\eta_3}\,i\tilde{\gamma}_1\tilde{\gamma}_3$$
$$+ \lambda_2 i\gamma_2\tilde{\gamma}_2 + \lambda_3 i\gamma_3\tilde{\gamma}_3 \tag{S85}$$

修正项现在包含两个独立参数 $\lambda_2, \lambda_3$，在每个时间步数值优化以确保基态简并。

### 经验有效重叠公式

通过测试不同的平均方式，发现以下经验公式效果最好：

$$\eta_{\text{eff}} = \sqrt{\eta_1\eta_2\eta_3} \tag{S86}$$

将 $\eta_{\text{eff}}$ 代入解析公式中替代 $\eta$。图 S2 显示数值解与使用 $\eta_{\text{eff}}$ 的解析曲线吻合得很好。对于过大的 $\eta_{\text{eff}}$，优化会失败。

---

## 总结：推导逻辑链

```
三个最小 Kitaev 链耦合 (Eq.S3)
        │
        ▼
引入 MBS 空间重叠 ζ (Eq.S4)
        │
        ▼
变换后出现 H̃ 项打破基态简并
        │
        ▼
引入修正项 H_corr 恢复简并 (Eq.S7)
        │
        ▼
无量纲化 → 主论文 Eq.(1) 的 h (Eq.S11)
        │
        ▼
对角化哈密顿量 — 三步基旋转 (Sec.II)
        │
        ▼
Berry 相位计算 — 三步分治 (Sec.V)
        │
        ▼
第二区间简并条件求解 (Sec.VI)
        │
        ▼
最终编织算符:
U(η) = exp(π/4·(1+ν₂₃)γ₃γ₂) · exp(π/4·(1-α₂₃)γ̃₃γ̃₂)
        │
        ▼
MBS 相似度: S(U) = sin²(2(φ-φ̃))
        │
        ▼
宇称测量期望 ∝ 1 - 2S(U)
```

## 核心贡献

这个补充材料的核心贡献是：**在 MBS 存在空间重叠的"非完美"情况下**，

1. 通过引入修正项 $H_{\text{corr}}$ 恢复基态简并
2. **解析地**对角化了哈密顿量
3. **解析地**计算了 Berry 相位
4. 得到了编织操作的**精确解析表达式**
5. 建立了宇称测量与编织保真度之间的**直接函数关系**
6. 推广到非对称 MBS 重叠情况，给出经验有效参数公式
