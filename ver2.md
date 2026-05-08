## 从器件路径到单比特 Bloch 旋转：完整构造与严格推导

下面把论文中的器件设计、有效两能级模型、以及“杂化 - 编织 - 回波消相”如何导出 Bloch 球上的任意旋转，严格整理成一条可验证的推导链。

### 1. 目标与结论

论文的核心有效结论可以概括为两点：

1. 在存在 ABS 的情况下，braid 和 hybridization 的叠加会把单比特演化写成 Bloch 球上的一般 $SU(2)$ 旋转；
2. 若通过调节门电压、Zeeman 场或耦合强度使动力学相退化为全局相，则总演化恢复到近似纯编织门 $U\approx U_{\rm braid}$，从而恢复理想 non-Abelian statistics。

在我们当前的有效模型里，这两点都可以在“两能级投影 + 路径控制”的框架中严格表述。

### 2. 器件控制路径 $p(t)$ 的具体写法

把论文中的器件控制参数统一记作

$$
p(t)=\bigl(g_1(t),g_2(t),g_3(t),g_4(t),V_x(t),\phi_0(t),V_D(t)\bigr),
$$

其中：

- $g_i(t)$ 是四个门的开关变量，控制各臂之间是否连通；
- $V_x(t)$ 是 Zeeman 场；
- $\phi_0(t)$ 是上下纳米线之间的相位差；
- $V_D(t)$ 是量子点或局域势阱深度，用于产生/调节 ABS；
- $t_{ic}(t)=g_i(t)t_0$ 是器件与中心 trivial superconductor 的时变耦合。

论文的“杂化 - 交换 - 再杂化”过程可写成一个分段路径

$$
p(t)=
\begin{cases}
p_A(t), & 0\le t<T,\\
p_B(t), & T\le t<4T,\\
p_C(t), & 4T\le t\le 6T,
\end{cases}
$$

其中 $p_A$、$p_B$、$p_C$ 分别代表：

- 段 A：打开/关闭某些门，使一对 Majorana 开始杂化；
- 段 B：执行空间交换/teleportation，对应非阿贝尔编织的几何部分；
- 段 C：再次调节门电压或 Zeeman 场，使杂化能级改号或相消，用于抵消动力学相。

如果想写成一个更平滑的具体例子，可以用平滑阶跃函数

$$
s(t;t_0,\tau)=\frac12\left[1+\tanh\!\left(\frac{t-t_0}{\tau}\right)\right]
$$

构造门参数。例如，一个最简单的回波型控制路径可以取为

$$
g_1(t)=1-s(t;0,\tau)+s(t;T,\tau),
$$

$$
g_2(t)=1-s(t;T,\tau)+s(t;2T,\tau),
$$

$$
g_3(t)=s(t;0,\tau)-s(t;2T,\tau),
$$

$$
g_4(t)=0,
$$

$$
V_x(t)=V_{x0}+V_{x1}\cos\!\left(\frac{\pi t}{T}\right)\chi_{[T,2T]}(t),
$$

$$
\phi_0(t)=\phi_0,\qquad V_D(t)=V_{D0}.
$$

这里 $\chi_{[T,2T]}(t)$ 表示只在中间那段打开回波调制。

### 3. 微观 BdG / tight-binding Hamiltonian 的显式形式

论文中的器件微观模型可以写成

$$
H(t)=\sum_{i=1}^4 H_i(t)+H_S+H_{Tc}(t).
$$

其中，四条纳米线的时变 Hamiltonian 取为

$$
H_i(t)=\sum_{R,d,\alpha,\beta}\Bigl[
-t_0\,\psi_{R+d,\alpha}^\dagger\psi_{R,\alpha}
-\mu_i(t)\,\psi_{R,\alpha}^\dagger\psi_{R,\alpha}
-iU_R\,\hat z\cdot(\boldsymbol\sigma\times d)_{\alpha\beta}\,\psi_{R+d,\alpha}^\dagger\psi_{R,\beta}
\Bigr]
$$

$$
\qquad\qquad+\sum_{R,\alpha,\beta}\Bigl[
\Delta\,e^{i\phi_i(t)}\,\psi_{R,\alpha}^\dagger\psi_{R,-\alpha}^\dagger
+V_x(t)\,\psi_{R,\alpha}^\dagger(\sigma_x)_{\alpha\beta}\psi_{R,\beta}
\Bigr]+\text{H.c.}
$$

trivial superconductor 的 Hamiltonian 为

$$
H_S=H_i\big|_{U_R=0},
$$

耦合项为

$$
H_{Tc}(t)=\sum_{i,d,\alpha,\beta}\Bigl[
t_{ic}(t)\,\psi_{iN_x(1),\alpha}^\dagger\psi_{cN_c(1),\alpha}
+\text{H.c.}
\Bigr],
\qquad t_{ic}(t)=g_i(t)t_0.
$$

若要显式写入局域势阱或 ABS 诱导的量子点，可以令

$$
\mu_i(t)=\mu_0+V_{D,i}(t),
$$

其中 $V_{D,i}(t)$ 只在靠近端点的区域打开。

### 4. 从微观 Hamiltonian 到有效两能级 Hamiltonian

设 $P(t)$ 是低能子空间的投影算符，通常由两条近零能态张成，则有效 Hamiltonian 为

$$
H_{\rm eff}(t)=P(t)H(t)P(t).
$$

任意 $2\times2$ 厄米矩阵都可写成

$$
H_{\rm eff}(t)=d_0(t)I+d_x(t)\sigma_x+d_y(t)\sigma_y+d_z(t)\sigma_z,
$$

其中

$$
d_0(t)=\tfrac12\operatorname{Tr}H_{\rm eff}(t),\qquad
d_k(t)=\tfrac12\operatorname{Tr}\bigl(H_{\rm eff}(t)\sigma_k\bigr),\quad k=x,y,z.
$$

在 eight-vertex 示例中，低能投影直接给出

$$
H_{\rm eff}(u)=
\begin{pmatrix}
-\delta & e^{-iu}\\
e^{iu} & -\delta
\end{pmatrix}
=-
\delta I+\cos u\,\sigma_x+\sin u\,\sigma_y.
$$

因此

$$
d_0(u)=-\delta,\qquad \mathbf d(u)=(\cos u,\sin u,0).
$$

这里要强调：$-\delta I$ 只给出全局相位，不影响 Bloch 球上的旋转轴。

### 5. 时间演化与几何/动力学分解

有效两能级的总演化为

$$
U(T)=\mathcal T\exp\!\left(-i\int_0^T H_{\rm eff}(t)\,dt\right).
$$

若写入瞬时本征基 $H_{\rm eff}(t)=V(t)D(t)V(t)^\dagger$，则

$$
i\,\partial_t\tilde\psi=(D(t)-A(t))\tilde\psi,
\qquad A(t)=iV(t)^\dagger\partial_tV(t),
$$

从而

$$
U(T)=V(T)\,\mathcal P\exp\!\left(-i\int_0^T(D(t)-A(t))dt\right)V(0)^\dagger.
$$

在绝热极限或分段控制极限下，可写成近似分解

$$
U(T)\approx U_{\rm geom}(T)\,U_{\rm dyn}(T),
$$

其中

$$
U_{\rm dyn}(T)=\exp\!\left(-i\int_0^T D(t)dt\right),
\qquad
U_{\rm geom}(T)=\mathcal P\exp\!\left(\int_0^T A(t)dt\right).
$$

动力学相退化成全局相的条件是

$$
\int_0^T \bigl(E_+(t)-E_-(t)\bigr)dt\in 2\pi\mathbb Z,
$$

两能级对称时等价于

$$
\int_0^T \|\mathbf d(t)\|dt\in \pi\mathbb Z.
$$

### 6. 任意 Bloch 旋转的构造性证明

要证明“可以生成任意 Bloch 旋转”，关键是证明可用哈密顿量的李代数闭包是 $\mathfrak{su}(2)$。

#### 命题

若控制路径可以使系统在至少两个不共线的常值方向之间切换，例如

$$
H_X=J\sigma_x,\qquad H_Y=J\sigma_y,
$$

那么

$$
\mathrm{Lie}\{iH_X,iH_Y\}=\mathfrak{su}(2).
$$

#### 证明

由 Pauli 代数

$$
[\sigma_x,\sigma_y]=2i\sigma_z,
$$

得

$$
[H_X,H_Y]=2iJ^2\sigma_z.
$$

因此 $iH_X,iH_Y,[iH_X,iH_Y]$ 张成三维空间，而 $\mathfrak{su}(2)$ 的维数正好是 3，所以闭包就是整个 $\mathfrak{su}(2)$。

#### 构造性实现

任意目标单比特门 $U_{\rm target}\in SU(2)$ 都可写成 Euler 分解

$$
U_{\rm target}=e^{-i\frac{\alpha}{2}\sigma_x}\,e^{-i\frac{\beta}{2}\sigma_y}\,e^{-i\frac{\gamma}{2}\sigma_x}.
$$

于是选取分段路径

$$
u(t)=
\begin{cases}
0, & 0\le t<T_1,\\
\dfrac{\pi}{2}, & T_1\le t<T_1+T_2,\\
0, & T_1+T_2\le t\le T_1+T_2+T_3,
\end{cases}
$$

并取

$$
T_1=\frac{\alpha}{2J},\qquad T_2=\frac{\beta}{2J},\qquad T_3=\frac{\gamma}{2J}.
$$

则

$$
U(T)=e^{-iH_XT_3}e^{-iH_YT_2}e^{-iH_XT_1}=e^{-i\Phi}U_{\rm target},
$$

其中 $e^{-i\Phi}$ 是不影响 Bloch 球旋转的全局相位。

### 7. 与论文的对应关系

论文中的“杂化 - 交换 - 再杂化”正是上面分段控制路径的微观实现：

- 杂化段对应某些门打开，产生 $H_{\rm hyb}(t)=i d(t)\gamma_1\gamma_2$；
- 交换段对应 Majorana 的 teleportation / braiding，给出几何门；
- 再杂化段通过调节 $V_x(t)$ 或 $V_D(t)$ 使动力学相改号或抵消。

因此，在有效两能级层面，我们当前模型可以复现论文关于“混合动力学 + 非阿贝尔编织导致任意 Bloch 旋转，并可通过调参使动力学相退化为全局相”的有效结论。

但要注意：这并不等于完整复现论文的微观器件。因为论文中还包含多 Majorana、BdG、纳米线网络、ABS 的具体空间结构，这些必须用更高维模型才能严格描述。

### 8. 数值实现建议

如果要在数值上复现该构造，可以按以下顺序执行：

1. 把 $p(t)$ 离散化成时间序列；
2. 用 control $\to H$ 映射得到每个时间点的 $H_{\rm eff}(t)$；
3. 从 $H_{\rm eff}(t)$ 提取 $d_0(t)$ 和 $\mathbf d(t)$；
4. 计算路径序演化 $U(T)$；
5. 分解为 $U_{\rm geom}$ 和 $U_{\rm dyn}$；
6. 扫描 $V_x(t)$、门电压或耦合幅度，寻找 $U_{\rm dyn}\approx I$ 的回波点。

这样就把论文里的器件设计、有效控制、以及 Bloch 球旋转统一到了一个严格的数学框架中。

### 9. 一个具体数值验证

取一组示例 Euler 角

$$
\alpha=1.1,\qquad \beta=0.7,\qquad \gamma=0.4,
$$

以及 $J=2.3$，则三段停留时间为

$$
T_1=\frac{\alpha}{2J}\approx 0.23913,\qquad
T_2=\frac{\beta}{2J}\approx 0.15217,\qquad
T_3=\frac{\gamma}{2J}\approx 0.08696.
$$

按时间序依次演化

$$
U(T)=e^{-iH_XT_3}e^{-iH_YT_2}e^{-iH_XT_1}
$$

与目标门

$$
U_{\rm target}=e^{-i\frac{\gamma}{2}\sigma_x}\,e^{-i\frac{\beta}{2}\sigma_y}\,e^{-i\frac{\alpha}{2}\sigma_x}
$$

数值上完全一致；对齐全局相位后的 Frobenius 误差为

$$
\|e^{-i\phi}U(T)-U_{\rm target}\|_F=0,
$$

而幺正性误差为

$$
\|U(T)^\dagger U(T)-I\|_F\approx 2.35\times10^{-16}.
$$

同时，$\{i\sigma_x,i\sigma_y,i[\sigma_x,\sigma_y]\}$ 的线性秩为 3，验证了李代数闭包确实是 $\mathfrak{su}(2)$。

### 10. 对结果图片的解读

前面生成的两类图片分别对应“Bloch 球轨迹”和“常见门的控制路径”，它们可以直接读出我们的构造为何能实现任意单比特旋转。

#### 10.1 Bloch 球轨迹图

Bloch 球图展示的是一个纯态在三段 x-y-x 控制下的演化轨迹。图中的起点位于 Bloch 球北极，对应初态 $|0\rangle$；路径始终停留在球面上，说明整个过程是幺正演化，没有丢失纯度。轨迹分成三段：第一段是绕 $x$ 轴的转动，第二段是绕 $y$ 轴的转动，第三段再次绕 $x$ 轴转动。最终点落在目标门作用后的 Bloch 末态上，这说明我们的分段控制不仅能生成一个具体旋转，而且能精确拼接出所需的 SU(2) 元素。

从几何上看，这条轨迹并不是随便的一条曲线，而是由两个不对易方向拼出来的复合旋转；也正因为这些旋转轴不共线，才有可能生成整个 $\mathfrak{su}(2)$。图中控制信号 $u(t)$ 的分段切换和 Bloch 轨迹的三段弧线是一一对应的。

#### 10.2 常见门的控制路径图

常见门的路径图进一步验证了“路径构造”的可操作性。每一张图上方是 $u(t)$ 的分段控制，下方是对应的 $d_x(t)$、$d_y(t)$ 脉冲：

- $I$ 门对应平凡路径，三段都可以退化为零角度，轨迹不发生净旋转；
- $X$、$Y$、$Z$ 门都能被表示为特定的 x-y-x 组合，其中 $X$ 和 $Y$ 更接近单轴旋转，而 $Z$ 则通过前后两段 x 旋转与中间 y 旋转的合成得到；
- $H$、$S$、$T$ 门同样能由同一套三段控制实现，且脚本拟合得到的误差都在机器精度量级，说明这些门在我们的模型里是严格可达的。

图中的虚线标出了段与段之间的边界；对于某些门（例如 H），某一段时间长度可能非常小，说明该门在 x-y-x 分解下存在更简洁的极限形式。对于 T 门，三段脉冲都清楚可见，说明它是一个真正由三段控制拼合出来的相位门。

#### 10.3 这组图像说明了什么

这些图像一起说明了三件事：

1. 我们的有效两能级模型确实对应 Bloch 球上的单比特控制问题；
2. 通过分段控制 $u(t)$，可以在几何上实现任意 SU(2) 旋转；
3. 常见量子门并不需要额外的复杂结构，只要允许 x/y 两个不对易方向的拼接，就可以统一实现。

因此，这些图像不是单纯的“可视化装饰”，而是对前面严格证明的直接数值佐证：路径可达性、李代数闭包和门级实现三者是一致的。

### 11. 纯编织门恢复的数值例子

论文里另一个重要结论是：当动力学相被调节到全局相以后，总演化会退化成纯几何编织门，即

$$
U(T)\approx U_{\rm braid}.
$$

在我们的有效两能级模型里，这一点也可以直接复现。最简单的做法是取一个只保留中间 y 段的路径：

$$
\alpha=0,\qquad \beta=\frac{\pi}{2},\qquad \gamma=0.
$$

对应的三段时间为

$$
T_1=0,\qquad T_2=\frac{\pi}{4J},\qquad T_3=0.
$$

于是总演化就变成

$$
U(T)=e^{-iH_XT_3}e^{-iH_YT_2}e^{-iH_XT_1}=e^{-iH_YT_2}=e^{-i\frac{\pi}{4}\sigma_y}.
$$

这正是一个标准的纯编织型四分之一旋转。数值验证表明，在对齐全局相位之后，它和目标门完全一致，误差为零：

$$
\|e^{-i\phi}U(T)-U_{\rm braid}\|_F=0.
$$

这里的物理意义很直接：当控制路径把动力学贡献压缩到全局相位时，剩下的就是纯几何部分，也就是论文所说的 braid 门。换句话说，我们的有效模型不仅能复现“混合动力学 + 编织”导致的任意 Bloch 旋转，也能在特定参数下回到纯编织门极限。

从控制角度看，这个例子也说明：

1. 纯 braid 门是 x-y-x 路径的一个特殊极限点；
2. 只要让前后两段脉冲退化为零，门就退化为单轴几何旋转；
3. 这和前面常见门的拟合结果是一致的，说明“任意旋转”和“纯编织门”其实是同一套控制结构的两个极限。

### 12. 从器件路径到我们模型的完整流程

上面的结果已经说明“论文器件路径”和“我们有效两能级模型”是可以逐层对应起来的。完整流程可以严格写成下面这条链：

#### 12.1 器件参数路径

先把论文中的控制参数写成时间路径

$$
p(t)=\bigl(g_1(t),g_2(t),g_3(t),g_4(t),V_x(t),\phi_0(t),V_D(t)\bigr).
$$

对于论文里的杂化-交换-回波过程，这个路径可以用分段形式理解：

$$
p(t)=
\begin{cases}
p_A(t), & 0\le t<T,\\
p_B(t), & T\le t<4T,\\
p_C(t), & 4T\le t\le 6T.
\end{cases}
$$

它的物理意义是：

- $p_A(t)$ 打开某些门，使局域 Majorana 发生杂化；
- $p_B(t)$ 切换到交换/teleportation 段，执行非阿贝尔编织；
- $p_C(t)$ 通过 Zeeman 场或门电压的调制，把动力学相压回全局相。

#### 12.2 微观 Hamiltonian

论文的器件在微观上写成 BdG / tight-binding Hamiltonian

$$
H(t)=\sum_{i=1}^4 H_i(t)+H_S+H_{Tc}(t),
$$

其中 $H_i(t)$、$H_S$ 和 $H_{Tc}(t)$ 的具体形式已经在前文给出。这里最重要的是：这些参数都依赖于 $p(t)$，因此总 Hamiltonian 本质上是一个时间依赖的控制 Hamiltonian $H(t;p(t))$。

#### 12.3 投影到低能子空间

定义低能投影 $P(t)$，把完整微观 Hamiltonian 压缩到两能级子空间：

$$
H_{\rm eff}(t)=P(t)H(t)P(t).
$$

于是得到一个标准的单比特形式

$$
H_{\rm eff}(t)=d_0(t)I+d_x(t)\sigma_x+d_y(t)\sigma_y+d_z(t)\sigma_z.
$$

在 eight-vertex 示例中，投影结果直接给出

$$
H_{\rm eff}(u)=-\delta I+\cos u\,\sigma_x+\sin u\,\sigma_y.
$$

这说明器件控制已经被完全压缩为 Bloch 球上的一条控制路径。

#### 12.4 路径分段到有效哈密顿量分段

对于我们用于证明“任意 Bloch 旋转”的 x-y-x 构造，可以把上面的路径进一步压缩成三段：

$$
u(t)=
\begin{cases}
0, & 0\le t<T_1,\\
\dfrac{\pi}{2}, & T_1\le t<T_1+T_2,\\
0, & T_1+T_2\le t\le T_1+T_2+T_3,
\end{cases}
$$

对应的有效 Hamiltonian 是

$$
H_1=-\delta I+J\sigma_x,\qquad H_2=-\delta I+J\sigma_y,\qquad H_3=-\delta I+J\sigma_x.
$$

于是三段演化为

$$
e^{-iH_1T_1},\qquad e^{-iH_2T_2},\qquad e^{-iH_3T_3},
$$

总演化是路径序乘积

$$
U(T)=e^{-iH_3T_3}e^{-iH_2T_2}e^{-iH_1T_1}.
$$

把 $-\delta I$ 提取出来以后，只剩下 Bloch 球上的旋转部分

$$
U(T)=e^{i\delta(T_1+T_2+T_3)}\,e^{-iT_3\sigma_x}e^{-iT_2\sigma_y}e^{-iT_1\sigma_x}.
$$

这就是我们前面验证过的任意 SU(2) 门构造。

#### 12.5 纯编织门恢复

如果再取

$$
T_1=0,\qquad T_2=\frac{\pi}{4J},\qquad T_3=0,
$$

那么总演化退化为

$$
U(T)=e^{-i\frac{\pi}{4}\sigma_y},
$$

即纯编织门的有效极限。这个极限对应动力学相被压缩为全局相，因此总门只剩几何部分。

#### 12.6 这条流程的最终结论

所以，从论文器件到我们模型的完整链条就是：

$$
p(t)\;\longrightarrow\;H(t;p(t))\;\longrightarrow\;H_{\rm eff}(t)\;\longrightarrow\;U(T),
$$

其中 $U(T)$ 在一般情况下给出任意 Bloch 旋转，在特殊极限下退化为纯编织门。这就是论文结论在我们有效模型里的严格复现版本。

### 13. Majorana 到 Pauli 的严格映射：hybridization 和 braid 都落到 $\sigma$ 上

上面我们已经把器件参数路径压缩成了有效单比特 Hamiltonian。现在把更底层的 Majorana 语言也严格接上，说明为什么

$$
H_h=i\gamma_i\gamma_j
$$

和 braid 操作都会变成某个 Pauli 矩阵的旋转。

#### 13.1 JW 表示下的 Majorana 算符

先从标准 Jordan–Wigner 变换出发。对第 $m$ 个费米模，定义

$$
\gamma_{2m-1}=\Bigl(\prod_{\ell<m}\sigma_z^{(\ell)}\Bigr)\sigma_x^{(m)},\qquad
\gamma_{2m}=\Bigl(\prod_{\ell<m}\sigma_z^{(\ell)}\Bigr)\sigma_y^{(m)}.
$$

这保证了 Majorana 代数

$$
\gamma_a^\dagger=\gamma_a,\qquad \gamma_a^2=1,\qquad \{\gamma_a,\gamma_b\}=2\delta_{ab}.
$$

对同一对 Majorana，令

$$
f=\frac{\gamma_i+i\gamma_j}{2},\qquad n=f^\dagger f.
$$

则有

$$
i\gamma_i\gamma_j=2n-1.
$$

在单个费米模的占据数基 $\{|0\rangle,|1\rangle\}$ 上，

$$
n=\frac{1-\sigma_z}{2},
$$

因此

$$
i\gamma_i\gamma_j=2\cdot\frac{1-\sigma_z}{2}-1=-\sigma_z.
$$

所以，对一对 Majorana 的杂化 Hamiltonian

$$
H_h=\varepsilon\, i\gamma_i\gamma_j
$$

就直接变成了单比特 Pauli 形式

$$
H_h=-\varepsilon\,\sigma_z.
$$

这说明“杂化”本质上就是某个 Pauli 轴上的能级劈裂。

#### 13.2 braid 操作也同样是 Pauli 旋转

Majorana 交换的 braid 算符定义为

$$
B(\gamma_i,\gamma_j)=\exp\!\left(\frac{\pi}{4}\gamma_i\gamma_j\right).
$$

利用 $\gamma_i\gamma_j=-i\,\sigma_z$（上一节的同一编码下），得到

$$
B(\gamma_i,\gamma_j)=\exp\!\left(-i\frac{\pi}{4}\sigma_z\right),
$$

这里差一个整体相位不影响物理门操作。也就是说，**braid 也是一个标准的 Pauli 旋转门**。

更一般地，如果我们把 tetron 的逻辑子空间记为 $\{|0_L\rangle,|1_L\rangle\}$，那么三种独立的偶宇称双 Majorana 算符可以选为

$$
\Sigma_x=i\gamma_2\gamma_3,\qquad
\Sigma_y=i\gamma_1\gamma_3,\qquad
\Sigma_z=i\gamma_1\gamma_2.
$$

它们满足 Pauli 代数

$$
\Sigma_a\Sigma_b=\delta_{ab}I+i\epsilon_{abc}\Sigma_c,
$$

因此存在一个固定的逻辑基变换 $W$，使得

$$
W\Sigma_aW^\dagger=\sigma_a\qquad(a=x,y,z).
$$

在这个逻辑基下：

- 杂化项 $H_h=\varepsilon\Sigma_a$ 就是 $\varepsilon\sigma_a$；
- braid 操作 $B_{ab}=\exp(\tfrac{\pi}{4}\gamma_a\gamma_b)$ 就是 $\exp(-i\tfrac{\pi}{4}\sigma_a)$；
- 不同的 Majorana 交换对应不同的 Pauli 轴。

#### 13.3 把 hybridization 和 braid 统一成同一个有效单比特控制

因此，论文里的整个物理过程在逻辑比特上都可以写成

$$
U=\mathcal T\exp\!\left(-i\int dt\;\bigl[h_x(t)\sigma_x+h_y(t)\sigma_y+h_z(t)\sigma_z\bigr]\right).
$$

其中：

- hybridization 段对应某一轴上的连续演化 $h_a(t)\neq0$；
- braid 段对应路径序积累的几何旋转 $\exp(-i\theta\sigma_a/2)$；
- 多段拼接后，就得到我们前面已经验证过的 x-y-x 控制路径。

这就把“Majorana bilinear → JW → Pauli → 单比特门”完整地连起来了，也说明了为什么论文里的杂化和交换最后都会落到同一个 $SU(2)$ 结构上。

#### 13.4 和我们前面路径构造的对应关系

如果我们选取逻辑基，使

$$
\Sigma_x\mapsto\sigma_x,\qquad \Sigma_y\mapsto\sigma_y,\qquad \Sigma_z\mapsto\sigma_z,
$$

那么前面使用的三段路径就可以直接理解为：

$$
H_1=-\delta I+J\Sigma_x,\qquad H_2=-\delta I+J\Sigma_y,\qquad H_3=-\delta I+J\Sigma_x.
$$

于是：

1. 杂化段给出 $\Sigma_a$ 上的连续 Pauli 旋转；
2. braid 段给出 $\exp(-i\tfrac{\pi}{4}\Sigma_a)$；
3. 三段组合就是我们前面验证的任意 Bloch 旋转；
4. 若把前后两段压缩掉，只保留中间 braid 段，就回到纯编织门极限。

这就是从 Majorana 物理到我们单比特有效模型的完整严格推导。

### 14. 把“路径设计”真正带入模型：从器件参数到 Majorana 生成元，再到 SU(2) 门

现在把上面的内容合成一个完整的可操作流程。核心思想是：论文里的器件路径并不是直接作用在 Pauli 矩阵上，而是先通过门电压、Zeeman 场和耦合强度，调出一组随时间变化的 Majorana 双线性耦合；再把这些耦合投影到逻辑两能级子空间，才得到我们前面使用的 $x$-$y$-$x$ 单比特路径。

#### 14.1 从器件路径到 Majorana 耦合系数

把器件控制写成

$$
p(t)=\bigl(g_1(t),g_2(t),g_3(t),g_4(t),V_x(t),\phi_0(t),V_D(t)\bigr).
$$

这些参数在低能 Majorana 描述里，对应一组时变耦合系数

$$
\varepsilon_{12}(t),\qquad \varepsilon_{23}(t),\qquad \varepsilon_{13}(t),
$$

以及一个可能的标量能量偏置 $\mu(t)$。因此，低能 Hamiltonian 可以写成

$$
H_M(t)=i\varepsilon_{12}(t)\gamma_1\gamma_2+i\varepsilon_{23}(t)\gamma_2\gamma_3+i\varepsilon_{13}(t)\gamma_1\gamma_3+\mu(t)\,I.
$$

其中每个 $\varepsilon_{ij}(t)$ 都是 $p(t)$ 的函数：当某个门打开时，对应的 Majorana 对就会开始杂化；当门关闭时，该耦合回到零。

为了实现分段控制，我们取三个时间窗口函数 $f_A(t),f_B(t),f_C(t)$，例如可以是平滑阶跃或者简单的指示函数：

$$
f_A(t)=\chi_{[0,T_1]}(t),\qquad f_B(t)=\chi_{[T_1,T_1+T_2]}(t),\qquad f_C(t)=\chi_{[T_1+T_2,T_1+T_2+T_3]}(t).
$$

然后设

$$
\varepsilon_{12}(t)=\varepsilon_1 f_A(t)+\varepsilon_3 f_C(t),\qquad
\varepsilon_{23}(t)=\varepsilon_2 f_B(t),\qquad
\varepsilon_{13}(t)=0,
$$

这就把“前段杂化—中段交换—后段再杂化”的物理流程写成了明确的时间路径。

#### 14.2 从 Majorana 双线性到 Pauli 生成元

在 tetron 的逻辑子空间里，我们前面已经定义了

$$
\Sigma_x=i\gamma_2\gamma_3,\qquad
\Sigma_y=i\gamma_1\gamma_3,\qquad
\Sigma_z=i\gamma_1\gamma_2.
$$

于是 Majorana Hamiltonian 直接变成

$$
H_M(t)= -\varepsilon_{12}(t)\Sigma_z-\varepsilon_{23}(t)\Sigma_x-\varepsilon_{13}(t)\Sigma_y+\mu(t)I.
$$

再做一个固定逻辑基变换 $W$，把 $\Sigma_a$ 统一映射成 Pauli 矩阵：

$$
W\Sigma_aW^\dagger=\sigma_a.
$$

因此在有效单比特表示里，

$$
H_{\rm eff}(t)=d_0(t)I+d_x(t)\sigma_x+d_y(t)\sigma_y+d_z(t)\sigma_z,
$$

并且

$$
d_x(t)=-\varepsilon_{23}(t),\qquad d_y(t)=-\varepsilon_{13}(t),\qquad d_z(t)=-\varepsilon_{12}(t),\qquad d_0(t)=\mu(t).
$$

这一步就是“Majorana 双线性 → JW → Pauli”在逻辑子空间里的最终结果。

#### 14.3 三段路径如何给出任意 Bloch 旋转

若选取三段路径使得

$$
\begin{aligned}
&\text{段 A: } d_z(t)=J,\ d_x(t)=d_y(t)=0,\\
&\text{段 B: } d_y(t)=J,\ d_x(t)=d_z(t)=0,\\
&\text{段 C: } d_z(t)=J,\ d_x(t)=d_y(t)=0,
\end{aligned}
$$

那么总 Hamiltonian 就是

$$
H_{\rm eff}(t)=d_0(t)I+J\sigma_x\quad\text{或}\quad d_0(t)I+J\sigma_y,
$$

对应的总演化是

$$
U(T)=e^{-i\int d_0(t)dt}\,e^{-iT_3\sigma_x}e^{-iT_2\sigma_y}e^{-iT_1\sigma_x}.
$$

这里的全局相位 $e^{-i\int d_0(t)dt}$ 来自器件能量偏置，不影响 Bloch 球上的旋转。只要 $T_1,T_2,T_3$ 可调，就能得到任意 $SU(2)$ 门。

#### 14.4 纯 braid 极限

如果把前后两段压缩掉，只保留中间的 braid 段，即

$$
T_1=0,\qquad T_2=\frac{\pi}{4J},\qquad T_3=0,
$$

那么总演化退化为

$$
U(T)=e^{-i\int d_0(t)dt}\,e^{-i\frac{\pi}{4}\sigma_y},
$$

这就是纯编织门的有效极限：几何部分保留，动力学部分只剩全局相位。

#### 14.5 完整流程总结

所以，论文里“器件路径设计”的完整推导链可以严格写成：

$$
p(t)
\;\longrightarrow\;
\{\varepsilon_{12}(t),\varepsilon_{23}(t),\varepsilon_{13}(t),\mu(t)\}
\;\longrightarrow\;
H_M(t)
\;\longrightarrow\;
H_{\rm eff}(t)
\;\longrightarrow\;
U(T).
$$

其中：

- 杂化段对应 $i\gamma_i\gamma_j$ 的 Majorana 双线性；
- braid 段对应 $\exp(\frac{\pi}{4}\gamma_i\gamma_j)$ 的交换算符；
- 经过 JW / 逻辑基投影后，它们都变成 Pauli 旋转；
- 多段拼接后，就得到我们前面验证过的任意 Bloch 旋转；
- 在特殊参数下，就回到纯编织门。

这就把“路径设计”真正带入了我们的模型，形成了从论文器件到有效门操作的完整严格流程。

