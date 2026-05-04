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

 

# 四、低能子空间投影

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
H_{\text{eff}}(u)= \vec d(u)\cdot \vec \sigma
$$

$$
\boxed{
\vec d(u) = (\cos u,\ \sin u,\ \delta)
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

 

# 十四、可扩展方向

 

1. 数值绘制 $\vec d(u)$ 轨迹
2. 计算 Berry 曲率分布
3. 分析 YBE 偏离
4. 优化路径以实现稳定量子门

 
