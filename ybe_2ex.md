# Eight-Vertex 模型 → Pauli 表达 → MZM/ABS 与非阿贝尔控制（完整示例）

 

# 一、目标

本例展示如何从 **eight-vertex (R(u))** 出发：

$$
R(u)\ \longrightarrow\ H(u)\ \longrightarrow\ \text{Pauli 表达}\ \longrightarrow\ \vec d(u)
$$

并得到：

* MZM vs ABS 的统一判据
* 非阿贝尔演化的来源
* 可实现的控制路径

 

# 二、Eight-Vertex (R(u)) 结构


## 15.0.1 论文器件对应的时间依赖 Hamiltonian

如果把论文里的纳米线器件直接写成时间依赖 BdG Hamiltonian，可以记成

$$
H(t)=\sum_{i=1}^4 H_i(t)+H_S+H_{Tc}(t),
$$

其中四条纳米线的单条 Hamiltonian 为

$$
H_i(t)=\sum_{R,d,\alpha,\beta}\Bigl[
-t_0\,\psi_{R+d,\alpha}^\dagger\psi_{R,\alpha}
-\mu_i(t)\,\psi_{R,\alpha}^\dagger\psi_{R,\alpha}
-iU_R\,\hat z\cdot(\boldsymbol\sigma\times d)_{\alpha\beta}\,\psi_{R+d,\alpha}^\dagger\psi_{R,\beta}
\Bigr]
$$

$$
\qquad\qquad+
\sum_{R,\alpha,\beta}\Bigl[
\Delta\,e^{i\phi_i(t)}\,\psi_{R,\alpha}^\dagger\psi_{R,-\alpha}^\dagger
+V_x(t)\,\psi_{R,\alpha}^\dagger (\sigma_x)_{\alpha\beta}\psi_{R,\beta}
\Bigr]+\text{H.c.}
$$

Trivial superconductor 一般写成

$$
H_S=H_i\big|_{U_R=0},
$$

而耦合项写成

$$
H_{Tc}(t)=\sum_{i,d,\alpha,\beta}\Bigl[t_{ic}(t)\,\psi_{iN_x(1),\alpha}^\dagger\psi_{cN_c(1),\alpha}+\text{H.c.}\Bigr],
$$

并且

$$
t_{ic}(t)=g_i(t)\,t_0.
$$

如果要把门电压/局域势阱显式写进去，也可以令

$$
\mu_i(t)=\mu_0+V_{D,i}(t),
$$

其中 $V_{D,i}(t)$ 只在靠近端点的区域打开，用来模拟论文里的 QD confinement 或局域势阱。
标准形式：

$$
R(u)=
\begin{pmatrix}
a(u)&0&0&d(u)\\
0&b(u)&c(u)&0\\
0&c(u)&b(u)&0\\
d(u)&0&0&a(u)
\end{pmatrix}
$$

 

## 简化参数化（用于计算）

选取一个可解析子族：

$$
\begin{aligned}
a(u) &= \cos u \\
b(u) &= \sin u \\
c(u) &= \cos u \\
d(u) &= i\delta \sin u
\end{aligned}
$$

 

👉 得到：

$$
R(u)=
\begin{pmatrix}
\cos u & 0 & 0 & i\delta \sin u \\
0 & \sin u & \cos u & 0 \\
0 & \cos u & \sin u & 0 \\
i\delta \sin u & 0 & 0 & \cos u
\end{pmatrix}
$$

 

# 三、从 (R(u)) 构造 Hamiltonian

定义：

$$
H(u)= i\partial_u R(u)R^{-1}(u)
$$

 

## 计算结果（主导结构）

小 δ 展开（关键简化）

为了抓物理本质，我们取：

$$
R^{-1} \approx R^\dagger(u)+ \mathcal O(\delta^2)
$$
（因为接近 unitary）,因此得到
$$
\boxed{
H(u)\approx
\underbrace{\cos u(\sigma_x\otimes\sigma_x)
+
\sin u(\sigma_y\otimes\sigma_y)}_{geometric rotation}
+
\underbrace{\delta(\sigma_z\otimes\sigma_z)}_{coupling, ->abs}
}
$$

 

## 物理含义

| 项       | 作用              |
|   --     |     -  |
| (XX, YY) | 交换/配对（几何旋转）     |
| (ZZ)     | 能量劈裂（MZM → ABS） |

 

# 四、低能子空间投影(具体的量子比特)

选择：

$$
\mathcal S = {|01\rangle,\ |10\rangle}
$$

 

## 有效 Hamiltonian

$$
H_{\text{eff}}(u)=
\begin{pmatrix}
-\delta & e^{-iu}\\
e^{iu} & -\delta
\end{pmatrix}
$$

 

# 五、Bloch 球表示
写成Bloch向量形式
$$
H_{\text{eff}}(u)= d_0(u) I + \vec d(u)\cdot \vec \sigma
$$

$$
\boxed{
d_0(u)=-\delta,\qquad \vec d(u)=(\cos u,\ \sin u,\ 0)
}
$$

 

# 六、几何解释（核心结果）

 

## 1️⃣ MZM 极限

$$
\delta = 0
$$

* 路径：赤道圆
* 穿过原点
* 能量 $E \approx 0$

 

👉 特点：

* 强 Berry 相
* 几何主导
* 接近理想 braid

 

 

## 2️⃣ ABS 情况

$$
\delta \neq 0
$$

* 路径：偏移圆（纬度圈）
* 不经过原点
* 有有限能隙,$E=\pm \sqrt{1+\delta^2}$

 

👉 特点：

* 动力学相显著
* 几何效应减弱
* 偏离 braid / YBE

 

 

## 🔥 核心图像

| 状态  | Bloch 路径 |
|  - |   -- |
| MZM | 围绕原点     |
| ABS | 偏移轨道     |

 

# 七、Berry 相与动力学相

 

## Berry 相
两能级系统:
$$
A = \frac{1}{2}(1-\cos \theta)d\phi\\
\cos \theta = \frac{\delta}{\sqrt{1+\delta^2}}
$$

一圈的berry phase:
$$
\gamma = \pi\left(1 - \frac{\delta}{\sqrt{1+\delta^2}}\right)
$$

$\delta=0:\gamma=\pi$,最大几何;

$\delta \rightarrow \infty:\gamma \rightarrow 0$

 

## 动力学相

$$
\phi_{\text{dyn}}=\int E(u)du = 2\pi\sqrt{1+\delta^2}
$$

 

定义一个比值判据

$$
\mathcal{Q} = \frac{\gamma}{\phi_{\text{dyn}}} = \frac{\pi(1-\delta/\sqrt{1+\delta^2})}{2\pi \sqrt{1+\delta^2}}
$$



| 区域         | 物理        |
|    - |  -   |
| $\delta=0$ | MZM       |
|  $\delta \ll 1$ | quasi-MZM |
|  $\delta \gg 1$ | ABS       |

总结来说：eight-vertex R(u)，在 Pauli 投影后，本质上就是一个 Bloch 球上的路径问题；MZM 和 ABS 的区别就是这条路径是否穿过原点。
 

# 八、非阿贝尔性的来源


## 判据

$$
[H(u_1),H(u_2)] \neq 0
$$

 

在本模型中：

$$
[XX,YY] \neq 0
$$

 

👉 结论：

> **非阿贝尔性来自 Pauli 项之间的不对易**

 

# 九、构造合成 braid 操作

 

## Hamiltonian

$$
H(u)= J[\cos u,XX + \sin u,YY] + \epsilon XY
$$

 

## 演化

$$
U=\mathcal{T}\exp\left(-i\int_0^{2\pi}H(u)du\right)
$$

 

## 得到门

$$
U_{\text{braid}} \approx
\exp\left(-i\frac{\pi}{4}\sigma_y\right)
$$

 

# 十、控制路径设计

 

## 三段路径

* 段 A：(XX \to YY)
* 段 B：旋转方向改变
* 段 C：回到初始

 

## 关键条件

* 路径闭合
* Hamiltonian 不对易
* 子空间隔离

 

# 十一、实验对应

 

| Pauli 项 | 实现方式 |
|   - |  - |
| XX, YY  | 耦合控制 |
| ZZ      | 能量偏置 |
| XY      | 相位驱动 |

 

## 控制步骤

1. 设置耦合
2. 扫描相位 (u)
3. 执行闭路径
4. 得到非阿贝尔门

 

# 十二、统一理解

 

## 从 eight-vertex 到物理：

$$
R(u)
\longrightarrow
H(u)
\longrightarrow
\vec d(u)
\longrightarrow
\text{几何路径}
$$

 

## MZM vs ABS：

$$
\boxed{
\text{是否穿过 }\vec d=0
}
$$

 

 

# 十三、最重要结论

> **Eight-vertex 模型在 Pauli 投影后，本质上是一个 Bloch 球路径问题：**
>
> * MZM：路径围绕原点（几何主导）
> * ABS：路径远离原点（动力学主导）
>
> 非阿贝尔性由路径的不对易结构决定。


# 十四、任意 Bloch 旋转的构造性实现

如果我们把可控参数 $u(t)$ 设计成分段常值路径，并允许在若干个取值上停留，那么该模型可以实现论文所说的“任意 Bloch 旋转”。其严格依据是：只要可用的两个有效哈密顿量方向不共线，它们的李代数闭包就生成 $\mathfrak{su}(2)$。

## 14.1 我们模型里的 $H_{\rm eff}(u)$ 到底是什么

对 eight-vertex 这个例子，投影到 $\{|01\rangle,|10\rangle\}$ 后得到的有效两能级 Hamiltonian 是

$$
H_{\rm eff}(u)=
\begin{pmatrix}
-\delta & e^{-iu}\\
e^{iu} & -\delta
\end{pmatrix}.
$$

把 $e^{-iu}=\cos u-i\sin u$ 展开后，它等价于

$$
H_{\rm eff}(u)= -\delta I + \cos u\,\sigma_x + \sin u\,\sigma_y.
$$

因此 Bloch 向量就是

$$
\mathbf d(u)=(\cos u,\sin u,0),\qquad d_0(u)=-\delta.
$$

这里要注意：$-\delta I$ 只贡献全局相位，不改变 Bloch 球上的旋转；真正决定旋转轴的是 $\cos u\,\sigma_x+\sin u\,\sigma_y$。

## 14.2 把路径分成三段以后，每一段的 Hamiltonian

我们选三个控制区间：

$$
u(t)=
\begin{cases}
0, & 0\le t<T_1,\\
\dfrac{\pi}{2}, & T_1\le t<T_1+T_2,\\
0, & T_1+T_2\le t\le T_1+T_2+T_3.
\end{cases}
$$

于是三段上分别有

$$
H_1=H_{\rm eff}(0)= -\delta I+\sigma_x,
$$

$$
H_2=H_{\rm eff}(\pi/2)= -\delta I+\sigma_y,
$$

$$
H_3=H_{\rm eff}(0)= -\delta I+\sigma_x.
$$

因为 $I$ 和 $\sigma_{x,y}$ 对易，所以每一段的演化可以拆成“全局相位”乘“Bloch 旋转”：

$$
e^{-iH_1T_1}=e^{+i\delta T_1}\,e^{-iT_1\sigma_x},
$$

$$
e^{-iH_2T_2}=e^{+i\delta T_2}\,e^{-iT_2\sigma_y},
$$

$$
e^{-iH_3T_3}=e^{+i\delta T_3}\,e^{-iT_3\sigma_x}.
$$

## 14.3 总演化结果怎么得到

因此总演化是三段乘积

$$
U(T)=e^{-iH_3T_3}e^{-iH_2T_2}e^{-iH_1T_1}.
$$

把上面的分解代入，得到

$$
U(T)=e^{i\delta(T_1+T_2+T_3)}\,
e^{-iT_3\sigma_x}e^{-iT_2\sigma_y}e^{-iT_1\sigma_x}.
$$

也就是说，真正作用在 Bloch 球上的部分就是

$$
U_{\rm Bloch}=e^{-iT_3\sigma_x}e^{-iT_2\sigma_y}e^{-iT_1\sigma_x}.
$$

如果我们把目标门写成 Euler 形式

$$
U_{\rm target}=e^{-i\frac{\alpha}{2}\sigma_x}e^{-i\frac{\beta}{2}\sigma_y}e^{-i\frac{\gamma}{2}\sigma_x},
$$

那么只要选

$$
T_1=\frac{\alpha}{2},\qquad T_2=\frac{\beta}{2},\qquad T_3=\frac{\gamma}{2},
$$

就有

$$
U_{\rm Bloch}=U_{\rm target},
$$

而前面的 $e^{i\delta(T_1+T_2+T_3)}$ 只是全局相位，不影响 Bloch 旋转。

对于 eight-vertex 对应的有效两能级模型，可取两个控制点

$$
u_X=0,\qquad u_Y=\frac{\pi}{2},
$$

对应的 traceless 有效哈密顿量为

$$
H_X=J\sigma_x,\qquad H_Y=J\sigma_y,
$$

因此

$$
[H_X,H_Y]=2iJ^2\sigma_z\neq 0.
$$

所以 $\{H_X,H_Y\}$ 的李代数闭包就是 $\mathfrak{su}(2)$。任意目标门 $U_{\rm target}\in SU(2)$ 都可写成 Euler 分解

$$
U_{\rm target}=e^{-i\frac{\alpha}{2}\sigma_x}\,e^{-i\frac{\beta}{2}\sigma_y}\,e^{-i\frac{\gamma}{2}\sigma_x}.
$$

于是定义分段路径

$$
u(t)=
\begin{cases}
0, & 0\le t<T_1,\\
\dfrac{\pi}{2}, & T_1\le t<T_1+T_2,\\
0, & T_1+T_2\le t\le T_1+T_2+T_3,
\end{cases}
$$

并令

$$
T_1=\frac{\alpha}{2J},\qquad T_2=\frac{\beta}{2J},\qquad T_3=\frac{\gamma}{2J}.
$$

则总演化满足

$$
U(T)=e^{-iH_XT_3}e^{-iH_YT_2}e^{-iH_XT_1}=e^{-i\Phi}U_{\rm target},
$$

其中 $e^{-i\Phi}$ 是全局相位，不影响 Bloch 球上的旋转。若希望路径更平滑，可以把阶跃切换替换为短时间插值，但李代数生成性质和可达性不变。

 

# 十五、把论文器件设计写成路径，并判断我们的模型能否复现

论文里的器件控制可以抽象为时间路径

$$
p(t)=\bigl(g_1(t),g_2(t),g_3(t),g_4(t),V_x(t),\phi_0(t),V_D(t)\bigr),
$$

其中 $g_i(t)$ 表示四个门的开关，$V_x(t)$ 表示 Zeeman 场，$\phi_0(t)$ 表示上下纳米线相位差，$V_D(t)$ 表示量子点/局域势阱深度。论文中的三类核心过程可以写成下面三段路径：

## 15.0 一个可直接使用的具体路径

为了把“杂化—交换—再杂化”写得更具体，可以取一个平滑阶跃函数

$$
s(t;t_0,\tau)=\frac{1}{2}\left[1+\tanh\!\left(\frac{t-t_0}{\tau}\right)\right],
$$

用它来定义门电压路径。一个最简单的 braiding 回波例子是：

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
\phi_0(t)=\phi_0\;\text{(常数)},\qquad V_D(t)=V_{D0}\;\text{(常数)}.
$$

这里 $\chi_{[T,2T]}(t)$ 表示只在中间那段打开 Zeeman 回波调制。这个路径的物理意义是：

- $0\le t<T$：打开/关闭 $G_1,G_3$，让一个 Majorana 在两臂之间开始杂化；
- $T\le t<2T$：切换到 $G_2,G_1$ 的组合，同时调制 $V_x(t)$，实现交换并部分抵消动力学相；
- $2T\le t\le 3T$：切回 $G_3,G_2$ 的组合，完成再杂化/闭合。

如果只想写成最简的“分段常值”路径，也可以直接写为

$$
p(t)=
\begin{cases}
(1,1,0,0,V_{x0},\phi_0,V_{D0}), & 0\le t<T,\\
(0,1,1,0,V_{x0}+V_{x1}\cos(\pi t/T),\phi_0,V_{D0}), & T\le t<2T,\\
(1,0,1,0,V_{x0},\phi_0,V_{D0}), & 2T\le t\le 3T.
\end{cases}
$$

这个例子不是唯一的，但足够说明论文里的器件设计可以被写成明确的时间路径。

1. **杂化段**：打开某些门，让一对 Majorana 混合，等效成
$$
	H_{\rm hyb}(t)=i\,d(t)\,\gamma_1\gamma_2.
$$
2. **交换段**：切换门的组合，使 Majorana 在不同臂之间交换，等效成编织操作
$$
	B(\gamma_i,\gamma_j)=\exp\!\left(\frac{\pi}{4}\gamma_i\gamma_j\right).
$$
3. **回波/再杂化段**：通过调节 $V_x(t)$ 或 $V_D(t)$，让杂化能级改号或相消，从而抵消动力学相。

因此，论文的器件过程可以写成一个分段控制路径

$$
p(t)=
\begin{cases}
p_A(t), & 0\le t<T,\\
p_B(t), & T\le t<4T,\\
p_C(t), & 4T\le t\le 6T,
\end{cases}
$$

其中 $p_A,p_B,p_C$ 分别代表“杂化—交换—杂化”三类过程。投影到有效两能级后，对应

$$
H_{\rm eff}(t)=d_0(t)I+\mathbf d(t)\cdot\boldsymbol\sigma,
$$

更具体地，若把论文里的“杂化/交换”压缩到两能级里，可以写成

$$
H_{\rm eff}(t)=d_0(t)I+d_x(t)\sigma_x+d_y(t)\sigma_y+d_z(t)\sigma_z,
$$

其中一个最常见的路径化选择是

$$
d_x(t)=J\cos\theta(t),\qquad d_y(t)=J\sin\theta(t),\qquad d_z(t)=\delta(t),
$$

这里 $\theta(t)$ 对应门电压/耦合切换引起的旋转角，$\delta(t)$ 对应 ABS 杂化或 Zeeman 偏置引起的能级劈裂。

总演化为路径序

$$
U(T)=\mathcal T\exp\!\left(-i\int_0^T H_{\rm eff}(t)\,dt\right).
$$

## 15.1 带入我们的模型能复现什么

把论文的器件路径压缩到我们当前的有效模型后，能够复现的是**几何旋转和动力学相竞争**，也就是 Bloch 球上的任意单比特旋转。原因是我们已经具备：

- $H_{\rm eff}=d_0 I+\mathbf d\cdot\sigma$ 这样的单比特生成形式；
- 通过分段选择 $u(t)$ 可以得到至少两个不共线方向的哈密顿量；
- 通过调节各段停留时间，可以让不同段的动力学相相互抵消或累积。

因此，**在有效两能级层面，我们可以复现论文“混合动力学 + 非阿贝尔编织导致任意 Bloch 旋转”的结论**。


## 15.3 如果要在数值上复现，怎么做

1. 先把 $p(t)$ 离散化成时间序列；
2. 用 control $\to H$ 映射得到每个时间点的 $H_{\rm eff}(t)$；
3. 从 $H_{\rm eff}(t)$ 提取 $\mathbf d(t)$；
4. 用路径序演化得到 $U(T)$，再分解为 $U_{\rm geom}$ 和 $U_{\rm dyn}$；
5. 扫描 $V_x(t)$ 或门电压的幅度，寻找 $U_{\rm dyn}\approx I$ 的回波点。

这条流程和论文的物理顺序是对应的，只是我们把复杂器件压缩成了一个可计算的有效控制模型。


# 十六、可扩展方向

 

1. 数值绘制 $\vec d(u)$ 轨迹
2. 计算 Berry 曲率分布
3. 分析 YBE 偏离
4. 优化路径以实现稳定量子门

 
