## 结合支撑材料的完整推导：从路径到哈密顿，再到 MZM / ABS 与我们的模型

这份支撑材料的主线非常清楚：先设计门控路径，再得到时变微观哈密顿量；随后把它投影到低能零模子空间，得到 Majorana 的有效耦合；最后区分出理想 MZM 的非阿贝尔 braid，以及 ABS 参与时的动力学相积累。下面把这条链完整写出来，并和我们自己的 Bloch / 有效单比特模型对应起来。

### 1. 从控制路径开始：门电压如何进入模型

补充材料里，四个门的控制变量记作 $g_1(t),g_2(t),g_3(t),g_4(t)$。它们决定了纳米线和中心 trivial superconductor 的隧穿耦合

$$
t_{ic}(t)=g_i(t)t_0.
$$

这里 $t_0$ 是统一的基准跃迁幅度。路径设计的本质，就是按时间调节 $g_i(t)$，让某些 Majorana 先杂化、再交换、再回到杂化。

braiding 前的初始配置是

$$
g_1=g_2=0,\qquad g_3=g_4=1,
$$

这意味着只有某些臂与中心区相连，三个低能 Majorana 对被局域在不同的 domain wall 上。

### 2. 低能有效哈密顿量：从路径到 Majorana 耦合

在这个初始构型下，补充材料给出的低能有效哈密顿量为

$$
H_eff = iE_1\gamma_1\gamma_2 + iE_2\gamma_3\gamma_4 + iE_3\gamma_5\gamma_{50} + iE_4\gamma_6\gamma_{60}
+ i t_c(\cdots) + h.c.
$$

其中后面的 $i t_c(\cdots)$ 展开成门控变量的两两耦合项，例如 $g_1g_2\gamma_2\gamma_3$、$g_1g_3\gamma_{50}\gamma_2$ 等。

当系统满足 $N_x \gg N_c$ 时，有限尺寸项 $E_j$ 相对更小，可以忽略，于是得到更简洁的低能模型

$$
H_eff = i t_c(各对 Majorana 之间的耦合) + h.c.
$$

这一步就是“路径 $g_i(t)$ 直接变成 Majorana 双线性耦合”的关键。

### 3. 三步 braid 的核心逻辑

braiding 通过三步完成 $\gamma_2$ 和 $\gamma_3$ 的交换。每一步都由先开某个门、再关某个门组成，物理含义是把 MZM 通过 trivial superconductor 传送到另一条纳米线。

#### Step 1

在第一步中，固定 $g_2=0$，于是哈密顿量简化为

$$
H_eff^{(1)} = i t_c( g_1 g_3\gamma_{50}\gamma_2 + g_1 g_4\gamma_2\gamma_{60} + g_3 g_4\gamma_{50}\gamma_{60}) + h.c.
$$

对角化后出现两个瞬时零模：

$$
\gamma_\theta = \frac{1}{C}\left(\frac{\gamma_2}{g_1}+\frac{\gamma_{50}}{g_3}+\frac{\gamma_{60}}{g_4}\right),
$$

$$
\gamma_\phi = \gamma_3.
$$

在这个步骤里，$G_1$ 关、$G_3$ 开，于是 $\gamma_2$ 被传送到 $W_3$，即与中心区附近的辅助零模发生连接。

#### Step 2

第二步固定 $g_3=0$，有效哈密顿量变成

$$
H_eff^{(2)} = i t_c( g_1 g_4\gamma_{50}\gamma_2 + g_2 g_4\gamma_3\gamma_2 + g_2 g_4\gamma_3\gamma_{60}) + h.c.
$$

对应的瞬时零模可写为

$$
\gamma_\theta = \gamma_{50},
$$

$$
\gamma_\phi = \frac{1}{C}\left(-\frac{\gamma_2}{g_1}+\frac{\gamma_3}{g_2}+\frac{\gamma_{60}}{g_4}\right).
$$

这一步把 $\gamma_3$ 传送到 $W_1$，并天然带上一个 $\pi$ 相位（补充材料明确指出这个负号来自对角化）。

#### Step 3

第三步固定 $g_1=0$，得到

$$
H_eff^{(3)} = i t_c( g_2 g_3\gamma_{50}\gamma_3 + g_2 g_4\gamma_3\gamma_{60} + g_3 g_4\gamma_{50}\gamma_{60}) + h.c.
$$

瞬时零模为

$$
\gamma_\theta = \frac{1}{C}\left(\frac{\gamma_3}{g_2}+\frac{\gamma_{50}}{g_3}+\frac{\gamma_{60}}{g_4}\right),
$$

$$
\gamma_\phi = -\gamma_2.
$$

这一步把 $\gamma_2$ 传送到 $W_2$。至此，$\gamma_2$ 和 $\gamma_3$ 的空间位置被交换。

### 4. braid 的本质：不是简单的直积，而是累积演化

补充材料特别强调一个点：每一步的展开都必须基于**前一步结束后的状态**，不能总是基于初始基底来写。也就是说，braid 是带累积效应的。

对第 1 步，杂化演化可以写成

$$
U^{(1)}=\exp\!\left(\frac{\theta^{(1)}}{2}\gamma_1\gamma_2\right),
$$

其中

$$
\frac{\theta^{(1)}}{2}=\int_0^T \epsilon^{(1)}(t)\,dt.
$$

它对 Majorana 的作用是一个连续旋转，例如

$$
\gamma_1\to \cos\theta^{(1)}\,\gamma_1+\sin\theta^{(1)}\,\gamma_2,
$$

$$
\gamma_2\to -\sin\theta^{(1)}\,\gamma_1+\cos\theta^{(1)}\,\gamma_2.
$$

在步骤 1 到 3 的组合里，braid 算符是

$$
B=\exp\!\left(\frac{\pi}{4}\gamma_2\gamma_3\right),
$$

而总演化不是单独的 braid，也不是单独的杂化，而是“杂化段 + braid 段 + 杂化段”串起来的结果。

把前 3 步写成一个整体，补充材料给出的结构是

$$
U^{(3)}\,U^{(2)}\,B\,U^{(1)}.
$$

后 3 步重复前 3 步，因此整个 6 步的 braid 过程是

$$
U_F = \big(U^{(3)}\,U^{(2)}\,B\,U^{(1)}\big)^2.
$$

这正是“理想 braid + 累积动力学相”最清楚的形式。

### 5. MZM 怎么在这里被描述出来

MZM 的物理特征是：它们在 braid 过程中保持零能附近的局域化，而且交换过程不破坏能隙。

在补充材料里，这一点表现为：

1. 初始时，三对 Majorana 局域在不同 domain wall；
2. 交换过程中，瞬时零模始终存在；
3. 最终交换完成后，$\gamma_2$ 和 $\gamma_3$ 的空间位置互换，但态的非阿贝尔结构保留。

从几何上看，补充材料把三条有效低能 Majorana 写成一个三维向量

$$
\vec\gamma=(\gamma_2,\gamma_3,\gamma_{50}),
$$

并把耦合参数写成一个“方向向量”

$$
\vec\delta = (t_1,t_2,t_3)=t_c(g_1g_4,g_2g_4,g_3g_4).
$$

这样低能哈密顿量可以概括成

$$
H = i\gamma_{60}(\vec\delta\cdot\vec\gamma).
$$

这说明瞬时零模就是与 $\vec\delta$ 垂直的方向；也就是说，MZM 可以看成“在一个缓慢变化的方向场 $\vec\delta$ 上保持垂直的零模向量”。

这和我们自己的 Bloch 球模型是一一对应的：

- $\vec\delta$ 对应我们的 $\vec d$；
- 垂直于 $\vec\delta$ 的零模对应 Bloch 球上的几何轨迹；
- 围绕原点的闭合路径对应非阿贝尔 braid。

### 6. ABS 怎么被描述出来

ABS 的关键区别是：它不是严格的孤立零模，而是一对弱耦合的 MZM。补充材料给出了两个非常重要的证据。

#### 6.1 相位 knob 和 domain wall

如果 $\phi_0=0$，在 $t=2T$ 的时候会形成 domain wall，能隙关闭，并在 bulk 中额外出现一对零模。这会导致信息泄漏，braiding 失败。

如果 $\phi_0\neq0$，domain wall 被平滑化，整个过程能隙保持打开，braid 的非阿贝尔特性才得以保留。

这说明 ABS 和 MZM 的区别不是“有没有零模”这么简单，而是：

- MZM 是拓扑保护的、能隙打开下的零模；
- ABS 更容易受 phase 和局域耦合影响，带来动力学相和态混合。

#### 6.2 LDOS 的 twin-peak 证据

补充材料的 LDOS 图显示：在 $t=T$ 时，ABS 的局域密度分裂成两个部分，权重相等；这说明一个 ABS 可以理解成两个耦合的 MZM。

随着 braid 时间 $T$ 变化，LDOS 在 $t=5T$、$t=6T$ 处发生振荡，这种对 $T$ 的正弦式依赖正是耦合能量不为零的信号。

因此，ABS 在这个论文里并不是“坏掉的 MZM”这么笼统，而是一个明确的物理对象：

> ABS = 两个弱耦合的 MZM，带着显著的动力学相积累。

### 7. 这些结果如何映射到我们的模型

我们的模型写成

$$
H_eff(t)=d_0(t)I+\vec d(t)\cdot\vec\sigma.
$$

和补充材料的对应关系是：

- $d_0(t)$ 对应累计动力学相；
- $\vec d(t)$ 对应路径在逻辑子空间中的几何方向；
- 轨迹绕原点闭合，对应理想 braid；
#### 9.2.3.3 4x4 数值矩阵流程

如果要把上面的结果直接喂给数值程序，最方便的写法是先把张量哈密顿量写成标准基
$\{|00\rangle,|01\rangle,|10\rangle,|11\rangle\}$ 下的 4x4 矩阵。令

$$
H_{\otimes}(t)=t_c
\begin{pmatrix}
-(g_1g_2+g_3g_4) & 0 & 0 & -(g_1g_4+g_2g_3)+i(g_1g_3+g_2g_4) \\
0 & -g_1g_2+g_3g_4 & (g_1g_4-g_2g_3)+i(g_1g_3-g_2g_4) & 0 \\
0 & (g_1g_4-g_2g_3)-i(g_1g_3-g_2g_4) & g_1g_2-g_3g_4 & 0 \\
-(g_1g_4+g_2g_3)-i(g_1g_3+g_2g_4) & 0 & 0 & g_1g_2+g_3g_4
\end{pmatrix},
$$

其中 $g_i=g_i(t)$，或者在第 $k$ 段路径里替换成 $g_i^{(k)}(u)$。若显式保留 $+\mathrm{h.c.}$ 且系数取实数，则整个矩阵整体再乘一个 2。

于是第 $k$ 段的数值积分可按网格 $u_m=m/M$ 离散为

$$
R_{\otimes}^{(k)}
\approx
\prod_{m=M-1}^{0}
\exp\!\left[-iT\,\Delta u\,H_{\otimes}^{(k)}(u_m)\right],
\qquad \Delta u=\frac{1}{M}.
$$

这里乘积从右向左累积，和路径有序指数 $\mathcal T_u$ 的方向一致。六步总门就是

$$
R_{\otimes,F}
\approx
\big(R_{\otimes}^{(3)}R_{\otimes}^{(2)}R_{\otimes}^{(1)}\big)^2.
$$

如果只关心逻辑子空间，就用投影算符

$$
P_L=|01\rangle\langle01|+|10\rangle\langle10|
$$

取出

$$
R_{L,F}=P_L R_{\otimes,F} P_L.
$$

然后与论文的理想 braid 门比较：

$$
U_\mathrm{braid}=\exp\!\left(\frac{\pi}{4}\gamma_a\gamma_b\right)
=\frac{1}{\sqrt2}(I+\gamma_a\gamma_b)
\sim\exp\!\left(-i\frac{\pi}{4}\sigma_c\right).
$$

数值检验时，最直接的残差是

$$
\epsilon_\mathrm{braid}
=
\min_{\phi\in\mathbb R}
\big\|e^{i\phi}R_{L,F}-U_\mathrm{braid}\big\|_F,
$$

而如果要检查 braid/YBE 形式本身，就把 $\check R_\mathrm{braid}$ 作为两体门嵌入三体空间，验证

$$
\check R_{12}\check R_{23}\check R_{12}-\check R_{23}\check R_{12}\check R_{23}=0.
$$

这就是把你的路径门和论文里的理想 braid 直接放到同一个数值框架里的具体做法。
$$
\mathbf g^{(3)}(t)=(0,1-s(t),s(t),1),\qquad t\in[2T,3T],
$$

并在 $t\in[3T,6T]$ 重复一次前 3 步。这样就把补充材料中的六步路径严格地写成一个分段控制轨迹。

门控电压与隧穿强度满足

$$
t_{ic}(t)=g_i(t)t_0.
$$

在初始点 $\mathbf g(0)=(0,0,1,1)$，在三步之后依次到达

$$
(1,0,0,1)\to(0,1,0,1)\to(0,0,1,1),
$$

因此 $\gamma_2$ 与 $\gamma_3$ 的空间位置被交换。

#### 9.2 从路径到微观 Majorana 哈密顿量

在零模与辅助 Majorana 的基底中，低能哈密顿量为

$$
H_M(t)= iE_1\gamma_1\gamma_2+iE_2\gamma_3\gamma_4+iE_3\gamma_5\gamma_{50}+iE_4\gamma_6\gamma_{60}
+it_c\sum_{m<n} C_{mn}(t)\gamma_m\gamma_n+h.c.,
$$

其中

$$
C_{12}=g_1g_2,\quad C_{25}=g_1g_3,\quad C_{2,60}=g_1g_4,\quad C_{35}=g_2g_3,\quad C_{3,60}=g_2g_4,\quad C_{5,60}=g_3g_4.
$$

当 $N_x\gg N_c$ 时，$E_j\to0$，于是可忽略有限尺寸项，得到

$$
H_M(t)=it_c\big(g_1g_2\gamma_2\gamma_3+g_1g_3\gamma_{50}\gamma_2+g_1g_4\gamma_2\gamma_{60}+g_2g_3\gamma_3\gamma_{50}+g_2g_4\gamma_3\gamma_{60}+g_3g_4\gamma_{50}\gamma_{60}\big)+h.c.
$$

这就是“路径数据”进入“有效哈密顿量”的第一步。

#### 9.2.1 从路径参数到张量系数

路径 $p(t)$ 本身不是 $h_{uv}(t)$，它先决定一组随时间变化的控制幅度 $\lambda_k(t)=\lambda_k(p(t))$。对我们的装置来说，这些控制幅度至少包括：

$$
t_{ic}(t)=g_i(t)t_0,
$$

以及由 $V_x(t)$、$V_D(t)$、$\phi_0(t)$ 生成的局域能级偏移和相位因子。于是微观哈密顿量可统一写成

$$
H_{\otimes}(t)=H_{\otimes}^{(0)}+\sum_k \lambda_k(p(t))\,O_k,
$$

其中 $O_k$ 是固定的张量积算符基。若改写到 Pauli 张量基

$$
H_{\otimes}(t)=\sum_{\mu,\nu\in\{0,x,y,z\}} h_{\mu\nu}(t)\,\sigma_\mu\otimes\sigma_\nu,
$$

那么系数就是对这组基底的投影：

$$
h_{\mu\nu}(t)=\frac14\operatorname{tr}\!\left[(\sigma_\mu\otimes\sigma_\nu)\,H_{\otimes}(t)\right]
=\sum_k C_{\mu\nu}^{(k)}\,\lambda_k(p(t)).
$$

这里 $C_{\mu\nu}^{(k)}$ 是由模型结构和基底选择决定的固定矩阵元。于是可以把路径依赖理解成：

$$
p(t)\Longrightarrow \lambda_k(p(t))\Longrightarrow h_{\mu\nu}(t).
$$

在我们的最小模型里，通常有三类最重要的对应关系：

- $g_i(t)$ 主要控制含有 $\gamma_m\gamma_n$ 的交换/隧穿项，因此对应 $h_{xx},h_{yy},h_{xy},h_{yx}$ 这些非对角张量系数；
- $V_x(t)$ 和 $V_D(t)$ 主要进入对角项，因此对应 $h_{z0},h_{0z},h_{zz}$；
- $\phi_0(t)$ 主要控制复相位，因此会把同一组非对角系数在 $x/y$ 方向上重新分解，等价于旋转 $h_{xx},h_{yy},h_{xy},h_{yx}$ 的相对权重。

所以，路径到张量系数的关系不是一一硬编码，而是“控制参数先决定耦合幅度，再由固定 Pauli 张量基读出系数”。

严格复现链：

$$
p(t)\rightarrow \lambda_k(t)\rightarrow h_{\mu\nu}(t)\rightarrow H_{\otimes}(t)\rightarrow R(u)\rightarrow \mathrm{YBE}.
$$

$$
h_{\mu\nu}(t)=\frac14\operatorname{tr}\!\left[(\sigma_\mu\otimes\sigma_\nu)H_{\otimes}(t)\right],
\qquad
H_{\otimes}(t)=\sum_{\mu,\nu}h_{\mu\nu}(t)\,\sigma_\mu\otimes\sigma_\nu.
$$

#### 9.2.2 从 Majorana 到 Pauli 串，再到逻辑 sigma

如果你想把上面的 $H_M(t)$ 严格变成 qubit 哈密顿量，最标准的做法是先做 Jordan--Wigner 映射。取

$$
\gamma_{2j-1}=\Big(\prod_{m<j}\sigma_m^z\Big)\sigma_j^x,
\qquad
\gamma_{2j}=\Big(\prod_{m<j}\sigma_m^z\Big)\sigma_j^y.
$$

于是同一费米模上的杂化项直接变成

$$
i\gamma_{2j-1}\gamma_{2j}=-\sigma_j^z.
$$

对最近邻的 Majorana 二次项，有

$$
i\gamma_{2j-1}\gamma_{2j+1}=\sigma_j^y\sigma_{j+1}^x,
$$

$$
i\gamma_{2j-1}\gamma_{2j+2}=\sigma_j^y\sigma_{j+1}^y,
$$

$$
i\gamma_{2j}\gamma_{2j+1}=-\sigma_j^x\sigma_{j+1}^x,
$$

$$
i\gamma_{2j}\gamma_{2j+2}=-\sigma_j^x\sigma_{j+1}^y.
$$

因此一般地，

$$
H_M(t)=i\sum_{a<b}\epsilon_{ab}(t)\gamma_a\gamma_b
\xrightarrow{\mathrm{JW}}
H_P(t)=\sum_\alpha c_\alpha(t)S_\alpha,
$$

其中 $S_\alpha$ 是 Pauli 串或两体张量积，$c_\alpha(t)$ 则由 $\epsilon_{ab}(t)$ 和 JW 产生的 $Z$ 串共同决定。若只关心逻辑比特，继续投影到编码子空间

$$
H_L(t)=P_L H_P(t)P_L
=d_0(t)I+d_x(t)\sigma_x+d_y(t)\sigma_y+d_z(t)\sigma_z.
$$

这时的系数可以直接从矩阵元读出：

$$
d_0(t)=\tfrac12\operatorname{tr}H_L(t),\qquad
d_a(t)=\tfrac12\operatorname{tr}\big(H_L(t)\sigma_a\big),\quad a=x,y,z.
$$

如果采用 Bravyi--Kitaev，则把上面的费米子哈密顿量先编码成 BK qubit 哈密顿量：

$$
H_{\mathrm{BK}}(t)=U_{\mathrm{BK}}\,H_M(t)\,U_{\mathrm{BK}}^\dagger
=\sum_\alpha \tilde c_\alpha(t)\,\tilde S_\alpha,
$$

其中 $\tilde S_\alpha$ 是 BK 编码下的 Pauli 串，典型权重满足 $\mathrm{wt}(\tilde S_\alpha)=O(\log N)$，比 JW 的链式 $Z$ 串更短。对应链条写成

$$
H_M(t)\xrightarrow{\mathrm{BK}}H_{\mathrm{BK}}(t)\xrightarrow{P_L}H_L(t)=P_LH_{\mathrm{BK}}(t)P_L.
$$

若只保留 tetron 逻辑子空间，则

$$
H_L(t)=d_0(t)I+d_x(t)\sigma_x+d_y(t)\sigma_y+d_z(t)\sigma_z,
$$

而 $d_0,d_x,d_y,d_z$ 仍由

$$
d_a(t)=\tfrac12\operatorname{tr}\big(H_L(t)\sigma_a\big),\quad a=x,y,z,\qquad d_0(t)=\tfrac12\operatorname{tr}H_L(t)
$$

读出。BK 只改变中间 qubit 表示，不改变最终的逻辑 Bloch 旋转。

#### 9.2.2.1 具体两 qubit Pauli 张量积表达式

取四个相关 Majorana 的固定顺序

$$
(\Gamma_1,\Gamma_2,\Gamma_3,\Gamma_4)=(\gamma_2,\gamma_3,\gamma_{50},\gamma_{60}),
$$

并用两 qubit 的 JW 表示

$$
\Gamma_1=X\otimes I,\qquad
\Gamma_2=Y\otimes I,\qquad
\Gamma_3=Z\otimes X,\qquad
\Gamma_4=Z\otimes Y.
$$

则六个双线性项变成

$$
i\Gamma_1\Gamma_2=-Z\otimes I,
$$

$$
i\Gamma_3\Gamma_1=-Y\otimes X,
$$

$$
i\Gamma_1\Gamma_4=Y\otimes Y,
$$

$$
i\Gamma_2\Gamma_3=-X\otimes X,
$$

$$
i\Gamma_2\Gamma_4=-X\otimes Y,
$$

$$
i\Gamma_3\Gamma_4=-I\otimes Z.
$$

因此支撑材料中的低能哈密顿量可写成

$$
H_{\otimes}(t)=t_c\big[
-g_1g_2\,Z\otimes I
-g_1g_3\,Y\otimes X
+g_1g_4\,Y\otimes Y
-g_2g_3\,X\otimes X
-g_2g_4\,X\otimes Y
-g_3g_4\,I\otimes Z
\big].
$$

若保留显式 $+\mathrm{h.c.}$ 且系数实数，则整体只差一个常数因子 2。这个式子就是从支撑材料直接落到 Pauli 张量积基的具体表达式。

#### 9.2.2.2 张量系数到逻辑系数

把上式写成

$$
H_{\otimes}(t)=\sum_{\mu,\nu\in\{0,x,y,z\}} h_{\mu\nu}(t)\,\sigma_\mu\otimes\sigma_\nu,
$$

则非零系数为

$$
h_{ZI}(t)=-t_c g_1g_2,
\qquad
h_{YX}(t)=-t_c g_1g_3,
\qquad
h_{YY}(t)=+t_c g_1g_4,
$$

$$
h_{XX}(t)=-t_c g_2g_3,
\qquad
h_{XY}(t)=-t_c g_2g_4,
\qquad
h_{IZ}(t)=-t_c g_3g_4,
$$

其余 $h_{\mu\nu}(t)=0$（若显式保留 $+\mathrm{h.c.}$ 且系数实数，则整体统一乘以 2）。

在逻辑基底 $\{|01\rangle,|10\rangle\}$ 上，使用投影字典

$$
X\otimes X\mapsto +\sigma_x,\qquad
Y\otimes Y\mapsto +\sigma_x,\qquad
X\otimes Y\mapsto -\sigma_y,\qquad
Y\otimes X\mapsto +\sigma_y,
$$

$$
Z\otimes I\mapsto +\sigma_z,\qquad
I\otimes Z\mapsto -\sigma_z,
$$

得到

$$
d_x(t)=h_{XX}(t)+h_{YY}(t)=t_c\big(g_1g_4-g_2g_3\big),
$$

$$
d_y(t)=-h_{XY}(t)+h_{YX}(t)=t_c\big(g_2g_4-g_1g_3\big),
$$

$$
d_z(t)=h_{ZI}(t)-h_{IZ}(t)=t_c\big(g_3g_4-g_1g_2\big),
$$

$$
d_0(t)=h_{II}(t)-h_{ZZ}(t)=0.
$$

于是该两 qubit Pauli 张量积哈密顿量在逻辑子空间中的严格结果为

$$
H_L(t)=d_x(t)\sigma_x+d_y(t)\sigma_y+d_z(t)\sigma_z
=t_c\Big[
(g_1g_4-g_2g_3)\sigma_x
+(g_2g_4-g_1g_3)\sigma_y
+(g_3g_4-g_1g_2)\sigma_z
\Big],
$$

或在保留全局标量偏置时写成

$$
H_L(t)=d_0(t)I+\vec d(t)\cdot\vec\sigma.
$$

#### 9.2.3 由 Pauli 张量哈密顿构造 YBE 算符

如果目标是 Yang--Baxter 方程，那么真正的起点不是逻辑子空间里的 $H_L$，而是张量积空间里的 $H_{\otimes}(t)$ 或其等价的 Pauli 串表达。先选一个 Pauli 张量基参数化

$$
H_{\otimes}(u)=\sum_{\mu,\nu\in\{0,x,y,z\}} h_{\mu\nu}(u)\,\sigma_\mu\otimes\sigma_\nu.
$$

然后在这个张量空间上定义对应的两体算符

$$
R(u)=\mathcal P\exp\!\left(-i\int_0^u H_{\otimes}(s)\,ds\right),
$$

若 $H_{\otimes}$ 在参数区间内近似常数，也可写成

$$
R(u)=e^{-iu H_{\otimes}}.
$$

这时 $R(u)$ 本身也可以展开成 Pauli 张量基：

$$
R(u)=\sum_{\mu,\nu} r_{\mu\nu}(u)\,\sigma_\mu\otimes\sigma_\nu,
\qquad
r_{\mu\nu}(u)=\frac14\operatorname{tr}\!\left[(\sigma_\mu\otimes\sigma_\nu)R(u)\right].
$$

Yang--Baxter 方程检验的就是这个两体算符是否满足

$$
R_{12}(u)R_{13}(u+v)R_{23}(v)=R_{23}(v)R_{13}(u+v)R_{12}(u).
$$

所以顺序必须是：

$$
p(t)\Longrightarrow H_{\otimes}(t)\Longrightarrow R(u)\Longrightarrow \text{YBE 检验}.
$$

$$
H_{\otimes}(u)\Longrightarrow R(u)=\mathcal P\exp\!\left(-i\int_0^u H_{\otimes}(s)\,ds\right)
\Longrightarrow R_{12}R_{13}R_{23}=R_{23}R_{13}R_{12}.
$$

如果把前面的六段控制路径直接代回去，那么这里的 $R(u)$ 还可以写成分段有序指数的乘积。设第 $k$ 段路径对应的张量哈密顿量为 $H_{\otimes}^{(k)}(t)$，则

$$
R^{(k)}=\mathcal T\exp\!\left(-i\int_{(k-1)T}^{kT}H_{\otimes}^{(k)}(t)\,dt\right),
\qquad
R_F=R^{(6)}R^{(5)}\cdots R^{(1)}.
$$

在每一段内若 $H_{\otimes}^{(k)}$ 近似常数，还可近似为

$$
R_F\approx e^{-iT H_{\otimes}^{(6)}}e^{-iT H_{\otimes}^{(5)}}\cdots e^{-iT H_{\otimes}^{(1)}}.
$$

这就是把“路径”真正积分进 $R(u)$ 之后得到的整体两体门；然后再投影到逻辑子空间，才得到下面的 $R_L(u)$ 和 Bloch 旋转图像。

$$
R_L(u)=P_LR(u)P_L\sim e^{-i\Phi(u)}\exp\!\left(-i\frac{\Theta(u)}{2}\hat n\cdot\vec\sigma\right).
$$

然后才是

$$
H_L(t)=P_L H_{\otimes}(t)P_L\quad\text{或}\quad R_L(u)=P_L R(u)P_L,
$$

用于把同一个模型投影到逻辑比特上，解释成 $SU(2)$ 的 Bloch 旋转。

注意：不是任意一个 Pauli 张量积哈密顿量都自动给出 YBE 解；它必须落在一个可积 ansatz 里，例如六顶点 / XXZ 型的 Pauli 张量结构。我们的写法只是把“可积候选空间”先搭出来，再去验 YBE，而不是反过来先投影到单比特再谈方程。

#### 9.2.3.1 按设计路径代入的精确 Pauli 哈密顿量

把 9.1 的三段路径写成统一的段内参数 $u\in[0,1]$：

$$
\mathbf g^{(1)}(u)=(u,0,1-u,1),\qquad
\mathbf g^{(2)}(u)=(1-u,u,0,1),\qquad
\mathbf g^{(3)}(u)=(0,1-u,u,1),
$$

并在 $t\in[3T,6T]$ 重复一次前 3 段。把它直接代入 9.2.2.1 的两比特 Pauli 基哈密顿量，得到每一段的精确微观形式：

$$
H_{\otimes}^{(1)}(u)=t_c\big[-u(1-u)\,Y\otimes X+u\,Y\otimes Y-(1-u)\,I\otimes Z\big],
$$

$$
H_{\otimes}^{(2)}(u)=t_c\big[-u(1-u)\,Z\otimes I+(1-u)\,Y\otimes Y-u\,X\otimes Y\big],
$$

$$
H_{\otimes}^{(3)}(u)=t_c\big[-u(1-u)\,X\otimes X-(1-u)\,X\otimes Y-u\,I\otimes Z\big].
$$

因此，整个设计路径对应的精确张量基哈密顿量就是这三个段 Hamiltonian 的周期拼接，而不是任何 XXZ 模板。

如果继续投影到逻辑子空间 $\mathcal H_L=\mathrm{span}\{|01\rangle,|10\rangle\}$，则三段的精确逻辑哈密顿量分别为

$$
H_L^{(1)}(u)=t_c\big[u\,\sigma_x-u(1-u)\,\sigma_y+(1-u)\,\sigma_z\big],
$$

$$
H_L^{(2)}(u)=t_c\big[(1-u)\,\sigma_x+u\,\sigma_y-u(1-u)\,\sigma_z\big],
$$

$$
H_L^{(3)}(u)=t_c\big[-u(1-u)\,\sigma_x+(1-u)\,\sigma_y+u\,\sigma_z\big].
$$

这就是按照支撑材料路径本身得到的精确 Pauli 基哈密顿量。

#### 9.2.3.2 分段积分后的 $R$

若每段的物理时间长度为 $T$，则段内演化的精确结果是路径有序指数

$$
R_L^{(k)}=\mathcal T_u\exp\!\left(-iT\int_0^1 H_L^{(k)}(u)\,du\right),\qquad k=1,2,3.
$$

六步总演化就是

$$
R_{L,F}=\big(R_L^{(3)}R_L^{(2)}R_L^{(1)}\big)^2.
$$

如果不投影、直接保留两比特张量空间，则同理有

$$
R_{\otimes}^{(k)}=\mathcal T_u\exp\!\left(-iT\int_0^1 H_{\otimes}^{(k)}(u)\,du\right),
\qquad
R_{\otimes,F}=\big(R_{\otimes}^{(3)}R_{\otimes}^{(2)}R_{\otimes}^{(1)}\big)^2.
$$

这里没有把它伪装成 XXZ 的闭式矩阵，因为这条路径下的生成元是显式的、但非交换的 Pauli 组合；能写出的严格结果就是上面的分段有序指数乘积。

如果你要一个真正“可直接写下”的闭式版本，可以把每段冻结到一阶 Magnus 平均：

$$
\bar H_{\otimes}^{(1)}=t_c\left[-\frac16\,Y\otimes X+\frac12\,Y\otimes Y-\frac12\,I\otimes Z\right],
$$

$$
\bar H_{\otimes}^{(2)}=t_c\left[-\frac16\,Z\otimes I+\frac12\,Y\otimes Y-\frac12\,X\otimes Y\right],
$$

$$
\bar H_{\otimes}^{(3)}=t_c\left[-\frac16\,X\otimes X-\frac12\,X\otimes Y-\frac12\,I\otimes Z\right].
$$

于是取 \(\omega=t_c\sqrt{19}/6\) 后，三段的闭式近似为

$$
R_{\otimes,\mathrm{avg}}^{(1)}=e^{-iT\bar H_{\otimes}^{(1)}}
=\cos(\omega T)I-i\sin(\omega T)\frac{-Y\otimes X+3Y\otimes Y-3I\otimes Z}{\sqrt{19}},
$$

$$
R_{\otimes,\mathrm{avg}}^{(2)}=e^{-iT\bar H_{\otimes}^{(2)}}
=\cos(\omega T)I-i\sin(\omega T)\frac{-Z\otimes I+3Y\otimes Y-3X\otimes Y}{\sqrt{19}},
$$

$$
R_{\otimes,\mathrm{avg}}^{(3)}=e^{-iT\bar H_{\otimes}^{(3)}}
=\cos(\omega T)I-i\sin(\omega T)\frac{-X\otimes X-3X\otimes Y-3I\otimes Z}{\sqrt{19}}.
$$

对应的逻辑子空间平均哈密顿量则是

$$
\bar H_L^{(1)}=t_c\left[\frac12\sigma_x-\frac16\sigma_y+\frac12\sigma_z\right],
\qquad
\bar H_L^{(2)}=t_c\left[\frac12\sigma_x+\frac12\sigma_y-\frac16\sigma_z\right],
\qquad
\bar H_L^{(3)}=t_c\left[-\frac16\sigma_x+\frac12\sigma_y+\frac12\sigma_z\right].
$$

如果把论文里的理想 braid 代入，正确检查的不是一般谱参数 YBE，而是 braid relation。理想 braid 生成元是

$$
U_\mathrm{braid}}=\exp\!\left(\frac{\pi}{4}\gamma_a\gamma_b\right)
=\frac{1}{\sqrt2}\big(I+\gamma_a\gamma_b\big)
\sim\exp\!\left(-i\frac{\pi}{4}\sigma_c\right).
$$

把它记成 \(\check R_\mathrm{braid}}\)，则三体上应满足

$$
\check R_{12}\check R_{23}\check R_{12}=\check R_{23}\check R_{12}\check R_{23}.
$$

若再用标准约定 \(R=P\check R\)，它就对应常数型 YBE

$$
R_{12}R_{13}R_{23}=R_{23}R_{13}R_{12}.
$$

因此，结论是：论文里的理想 braid 门本身满足 braid relation，也就是在 \(\check R\) 约定下的 YBE；但我们这条路径得到的 \(R_{\otimes,F}\) 只是一个一般的分段演化，除非数值上恰好收敛到 \(U_\mathrm{braid}}\)，否则它不会自动满足 YBE。

#### 9.2.3.3 4x4 数值矩阵流程

如果要把上面的结果直接喂给数值程序，最方便的写法是先把张量哈密顿量写成标准基
$\{|00\rangle,|01\rangle,|10\rangle,|11\rangle\}$ 下的 4x4 矩阵。令

$$
H_4(t)=t_c
\begin{pmatrix}
-(g_1g_2+g_3g_4) & 0 & 0 & -(g_1g_4+g_2g_3)+i(g_1g_3+g_2g_4) \\
0 & -g_1g_2+g_3g_4 & (g_1g_4-g_2g_3)+i(g_1g_3-g_2g_4) & 0 \\
0 & (g_1g_4-g_2g_3)-i(g_1g_3-g_2g_4) & g_1g_2-g_3g_4 & 0 \\
-(g_1g_4+g_2g_3)-i(g_1g_3+g_2g_4) & 0 & 0 & g_1g_2+g_3g_4
\end{pmatrix},
$$

其中 $g_i=g_i(t)$，或者在第 $k$ 段路径里替换成 $g_i^{(k)}(u)$。若显式保留 $+\mathrm{h.c.}$ 且系数取实数，则整个矩阵整体再乘一个 2。

于是第 $k$ 段的数值积分可按网格 $u_m=m/M$ 离散为

$$
R_4^{(k)}
\approx
\prod_{m=M-1}^{0}
\exp\!\left[-iT\,\Delta u\,H_4^{(k)}(u_m)\right],
\qquad \Delta u=\frac{1}{M}.
$$

这里乘积从右向左累积，和路径有序指数 $\mathcal T_u$ 的方向一致。六步总门就是

$$
R_{4,F}
\approx
\big(R_4^{(3)}R_4^{(2)}R_4^{(1)}\big)^2.
$$

如果只关心逻辑子空间，就用投影算符

$$
P_L=|01\rangle\langle01|+|10\rangle\langle10|
$$

取出

$$
R_{L,F}=P_L R_{4,F} P_L.
$$

然后与论文的理想 braid 门比较：

$$
U_\mathrm{braid}=\exp\!\left(\frac{\pi}{4}\gamma_a\gamma_b\right)
=\frac{1}{\sqrt2}(I+\gamma_a\gamma_b)
\sim\exp\!\left(-i\frac{\pi}{4}\sigma_c\right).
$$

数值检验时，最直接的残差是

$$
\epsilon_\mathrm{braid}
=
\min_{\phi\in\mathbb R}
\big\|e^{i\phi}R_{L,F}-U_\mathrm{braid}\big\|_F,
$$

而如果要检查 braid/YBE 形式本身，就把 \(\check R_\mathrm{braid}\) 作为两体门嵌入三体空间，验证

$$
\check R_{12}\check R_{23}\check R_{12}-\check R_{23}\check R_{12}\check R_{23}=0.
$$

这就是把你的路径门和论文里的理想 braid 直接放到同一个数值框架里的具体做法。

#### 9.3 三步交换的瞬时零模求解

对任意时刻，将零模写成线性组合

$$
\Gamma(t)=a_2(t)\gamma_2+a_3(t)\gamma_3+a_{50}(t)\gamma_{50}+a_{60}(t)\gamma_{60}.
$$

零模条件是

$$
[H_M(t),\Gamma(t)]=0.
$$

补充材料中由于 $\gamma_{60}$ 与 $\gamma_{50}$ 形成高能配对，真正参与传输的三维低能子空间可取

$$
\vec\gamma=(\gamma_2,\gamma_3,\gamma_{50}),\qquad \vec\delta=t_c(g_1g_4,g_2g_4,g_3g_4).
$$

于是传输部分的有效哈密顿量可写为

$$
H_{\mathrm{tr}}(t)=i\gamma_{60}(\vec\delta(t)\cdot\vec\gamma).
$$

设

$$
\Gamma_a=\vec a\cdot\vec\gamma.
$$

代入可得零模条件等价于

$$
\vec a\cdot\vec\delta=0.
$$

也就是说，零模就是三维向量空间中与 $\vec\delta$ 垂直的方向。取与 $\vec\delta$ 正交的两条单位矢量 $\hat e_\theta,\hat e_\phi$，便得到瞬时零模

$$
\gamma_\theta=\vec\gamma\cdot\hat e_\theta,\qquad \gamma_\phi=\vec\gamma\cdot\hat e_\phi.
$$

这正是补充材料中式 (S4)–(S9) 的统一几何表述。三步路径分别对应 $\vec\delta$ 在三个坐标平面中的转动，因此零模方向连续转移：

$$
\gamma_2\to\gamma_{50}\to\gamma_3\to-\gamma_2,
$$

从而完成交换。

具体到每一步，写成与补充材料一致的瞬时零模即为

Step 1: $g_2=0$, 
$H^{(1)}=it_c(g_1g_3\gamma_{50}\gamma_2+g_1g_4\gamma_2\gamma_{60}+g_3g_4\gamma_{50}\gamma_{60})+h.c.$

$$
\gamma_\theta=\frac{1}{C}\left(\frac{\gamma_2}{g_1}+\frac{\gamma_{50}}{g_3}+\frac{\gamma_{60}}{g_4}\right),\qquad \gamma_\phi=\gamma_3.
$$

Step 2: $g_3=0$, 
$H^{(2)}=it_c(g_1g_4\gamma_{50}\gamma_2+g_2g_4\gamma_3\gamma_2+g_2g_4\gamma_3\gamma_{60})+h.c.$

$$
\gamma_\theta=\gamma_{50},\qquad \gamma_\phi=\frac{1}{C}\left(-\frac{\gamma_2}{g_1}+\frac{\gamma_3}{g_2}+\frac{\gamma_{60}}{g_4}\right).
$$

Step 3: $g_1=0$, 
$H^{(3)}=it_c(g_2g_3\gamma_{50}\gamma_3+g_2g_4\gamma_3\gamma_{60}+g_3g_4\gamma_{50}\gamma_{60})+h.c.$

$$
\gamma_\theta=\frac{1}{C}\left(\frac{\gamma_3}{g_2}+\frac{\gamma_{50}}{g_3}+\frac{\gamma_{60}}{g_4}\right),\qquad \gamma_\phi=-\gamma_2.
$$

这里负号就是补充材料强调的 $\pi$ 相位来源。

#### 9.4 braid 算符与六步总演化

Majorana 交换由

$$
U_B=\exp\!\left(\frac{\pi}{4}\gamma_2\gamma_3\right)
$$

实现，因为它满足

$$
U_B^\dagger\gamma_2U_B=\gamma_3,\qquad U_B^\dagger\gamma_3U_B=-\gamma_2.
$$

更一般地，若

$$
U(\theta)=\exp\!\left(\frac{\theta}{2}\gamma_2\gamma_3\right),
$$

则

$$
U(\theta)^\dagger\gamma_2U(\theta)=\gamma_2\cos\theta+\gamma_3\sin\theta,
$$

$$
U(\theta)^\dagger\gamma_3U(\theta)=-\gamma_2\sin\theta+\gamma_3\cos\theta.
$$

取 $\theta=\pi/2$，便得到 $\gamma_2\to\gamma_3$、$\gamma_3\to-\gamma_2$ 的 braid 交换。

六步路径的总演化写成

$$
U_F=\big(U^{(3)}U^{(2)}BU^{(1)}\big)^2,
$$

其中每个 $U^{(k)}$ 都是对应一步的绝热演化算符

$$
U^{(k)}=\mathcal T\exp\!\left(-i\int_{(k-1)T}^{kT}H^{(k)}(t)\,dt\right).
$$

这说明 braid 不是一个单点操作，而是带有前后步累积的分段演化。

#### 9.5 ABS 的两 Majorana 解释与 LDOS 判据

若只看 ABS，本质上可写成一对弱耦合 Majorana

$$
H_{\mathrm{ABS}}=i\epsilon_{\mathrm{ABS}}\eta_1\eta_2,
$$

因此能级是

$$
E_\pm=\pm\epsilon_{\mathrm{ABS}}.
$$

这意味着 ABS 不是严格零模，而是具有有限劈裂的双峰结构。其 LDOS 可概括为

$$
\rho(x,E;t)=\sum_n |\psi_n(x,t)|^2\,\delta(E-E_n(t)).
$$

因此在 $t=T$ 处出现的双峰分裂，就是两个弱耦合 Majorana 的空间分离；而后续随 $T$ 的振荡，则对应

$$
\epsilon_{\mathrm{ABS}}\neq0
$$

时的动力学混合。

相位 knob 的作用则是控制 domain wall 是否闭合：

$$
\phi_0=0\ \Rightarrow\ \Delta_{\min}(2T)=0,
$$

$$
\phi_0\neq0\ \Rightarrow\ \Delta_{\min}(t)>0\ \text{for all }t.
$$

前者允许额外零模从 bulk 漏出，后者保持 gap 打开，因此前者破坏 braid，后者保护 braid。

#### 9.6 投影到我们的 Bloch 模型

严格复现链：

$$
H_M(t)\xrightarrow{\mathrm{JW}}H_P(t)=\sum_\alpha c_\alpha(t)S_\alpha\xrightarrow{P_L}H_L(t)=d_0(t)I+\vec d(t)\cdot\vec\sigma.
$$

$$
H_{\otimes}(t)=\sum_{\mu,\nu}h_{\mu\nu}(t)\,\sigma_\mu\otimes\sigma_\nu,
\qquad
H_L(t)=P_LH_{\otimes}(t)P_L.
$$

这里要按我们的模型链条来写：先由路径得到张量积形式的微观哈密顿量，再投影到编码逻辑子空间，最后才得到单个 $\sigma$。

先写微观层的 Pauli 张量积展开：

$$
H_{\otimes}(t)=\sum_{\mu,\nu\in\{0,x,y,z\}} h_{\mu\nu}(t)\,\sigma_\mu\otimes\sigma_\nu.
$$

这里的 $h_{\mu\nu}(t)$ 由路径参数 $p(t)$ 决定；若从 Majorana 双线性出发，则它们先经 JW 变成 Pauli 串或两体张量积，再整理成上式。选定逻辑子空间 $\mathcal H_L=\mathrm{span}\{|0_L\rangle,|1_L\rangle\}$，令投影算符为 $P_L$，则有效哈密顿量是

$$
H_L(t)=P_L H_{\otimes}(t)P_L.
$$

若在逻辑基底 $\{|0_L\rangle,|1_L\rangle\}$ 上展开，便有

$$
H_L(t)=
\begin{pmatrix}
\langle 0_L|H_{\otimes}(t)|0_L\rangle & \langle 0_L|H_{\otimes}(t)|1_L\rangle \\
\langle 1_L|H_{\otimes}(t)|0_L\rangle & \langle 1_L|H_{\otimes}(t)|1_L\rangle
\end{pmatrix}
=\begin{pmatrix}
\varepsilon_0(t) & \Delta(t) \\
\Delta^*(t) & \varepsilon_1(t)
\end{pmatrix}.
$$

把这个 $2\times2$ 厄米矩阵按 Pauli 基展开，就得到

$$
H_L(t)=\frac{\varepsilon_0(t)+\varepsilon_1(t)}{2}I
+\operatorname{Re}\Delta(t)\,\sigma_x
-\operatorname{Im}\Delta(t)\,\sigma_y
+\frac{\varepsilon_0(t)-\varepsilon_1(t)}{2}\sigma_z.
$$

因此

$$
d_0(t)=\frac{\varepsilon_0(t)+\varepsilon_1(t)}{2},\qquad
d_x(t)=\operatorname{Re}\Delta(t),\qquad
d_y(t)=-\operatorname{Im}\Delta(t),\qquad
d_z(t)=\frac{\varepsilon_0(t)-\varepsilon_1(t)}{2}.
$$

任意 $2\times2$ 厄米矩阵都可唯一写成

$$
H_L(t)=d_0(t)I+d_x(t)\sigma_x+d_y(t)\sigma_y+d_z(t)\sigma_z
=d_0(t)I+\vec d(t)\cdot\vec\sigma.
$$

其中

$$
d_0(t)=\tfrac12\operatorname{tr}H_L(t),\qquad d_a(t)=\tfrac12\operatorname{tr}\big(H_L(t)\sigma_a\big).
$$

于是总演化为

$$
U(T)=\mathcal T\exp\!\left(-i\int_0^T H_L(t)\,dt\right)
=e^{-i\Phi(T)}\,\mathcal T\exp\!\left(-i\int_0^T \vec d(t)\cdot\vec\sigma\,dt\right),
$$

其中

$$
\Phi(T)=\int_0^T d_0(t)\,dt.
$$

若每一段路径内 $\vec d(t)$ 近似固定，则对第 $k$ 段有

$$
U_k=e^{-i\tau_k d_{0,k}}\left[\cos(|\vec d_k|\tau_k)I-i\sin(|\vec d_k|\tau_k)\hat d_k\cdot\vec\sigma\right],
$$

从而

$$
U(T)=\prod_k U_k= e^{-i\Phi(T)}\,U_{\mathrm{SU(2)}}.
$$

当 $\vec d(t)$ 绕原点闭合时，$U_{\mathrm{SU(2)}}$ 就是 Bloch 球上的旋转门；当 $\Phi(T)$ 只贡献全局相位时，就得到论文所说的“理想 braid + 相位动力学累积 = Bloch 旋转”。

#### 9.7 复现标准

因此，按支撑材料复现我们的模型，实际要验证的等式链就是

$$
p(t)\Longrightarrow H_M(t)\Longrightarrow H_L(t)=d_0I+\vec d\cdot\vec\sigma
\Longrightarrow U(T)=e^{-i\Phi(T)}U_{\mathrm{SU(2)}}.
$$

只要数值结果满足：

$$
\Delta_{\min}(\phi_0\neq0)>0,\qquad \rho(x,t=T)\text{ split into two equal-weight peaks},
$$

$$
U_F=\big(U^{(3)}U^{(2)}BU^{(1)}\big)^2,
$$

以及

$$
U(T)=e^{-i\Phi(T)}\exp\!\left(-i\frac{\Theta(T)}{2}\hat n\cdot\vec\sigma\right),
$$

那么我们就不是在“解释”论文，而是在同一控制链上复现它。

### 10. 最后总结

这份支撑材料把论文主文中的结论完整地补成了一个严格链条：

路径 $g_i(t)$  ->  微观哈密顿量 $H_M(t)$  ->  低能零模 $H_eff(t)$  ->  瞬时 MZM / ABS 描述  ->  braid 门 + 累积相位。

对我们自己的模型来说，这意味着：只要我们能把路径设计成让 $\vec d(t)$ 在逻辑子空间里实现闭合绕行，并把 $d_0(t)$ 处理成全局相位或回波消去，那么就能在同一个框架下同时描述：

- MZM 的拓扑 braid；
- ABS 的动力学偏移；
- 最终的 Bloch 球单比特旋转。

这就是这份支撑材料对我们模型最直接的验证意义。
