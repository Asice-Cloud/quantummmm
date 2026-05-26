# Majorana SO(2N) Closure Phase Theory

# 从 Clifford Algebra 到 Dynamical Lie Phase 的统一框架

---

# 一、理论目标

本文建立一个完全一般化的 Majorana 动力学框架。

核心思想：

不再：

- 预设 qubit
- 预设 Bloch sphere
- 预设 SU(2)
- 预设低能投影
- 预设 quaternion dynamics

而是：

从：

\[
2N
\]

个 Majorana 的一般 Hamiltonian 出发，

系统研究：

---

## 1. Clifford algebra 的一般结构

## 2. Quadratic Majorana Hamiltonian 的普适形式

## 3. SO(2N) dynamical algebra

## 4. Dynamical closure

## 5. Magnus hierarchy

## 6. Lie subalgebra classification

## 7. Holonomy enlargement

## 8. Fidelity geometry

## 9. PRB105 与 PRB111 的统一解释

---

# 二、Majorana 与 Clifford Algebra

定义：

\[
\gamma_i^\dagger=\gamma_i
\]

满足：

\[
\{\gamma_i,\gamma_j\}=2\delta_{ij}
\tag{1}
\]

对于：

\[
2N
\]

个 Majorana：

\[
\gamma_1,\gamma_2,\dots,\gamma_{2N}
\]

它们生成 Clifford algebra：

\[
Cl(2N)
\tag{2}
\]

---

# 三、Quadratic Majorana Hamiltonian 的普适形式

任意 quadratic Majorana Hamiltonian：

\[
\boxed{
H(t)
=
\frac i2
\sum_{i,j}
A_{ij}(t)\gamma_i\gamma_j
}
\tag{3}
\]

其中：

\[
A_{ij}=-A_{ji}
\]

是真反对称矩阵。

定义 generator：

\[
L_{ij}=i\gamma_i\gamma_j
\tag{4}
\]

则：

\[
H(t)=\sum_{i<j}A_{ij}(t)L_{ij}
\tag{5}
\]

---

# 四、SO(2N) 的生成

利用 Clifford algebra：

\[
\{\gamma_i,\gamma_j\}=2\delta_{ij}
\]

计算：

\[
[L_{ij},L_{kl}]
\]

得到：

\[
\boxed{
[L_{ij},L_{kl}]
=
2i(
\delta_{jk}L_{il}
-
\delta_{ik}L_{jl}
-
\delta_{jl}L_{ik}
+
\delta_{il}L_{jk}
)
}
\tag{6}
\]

这正是：

\[
\mathfrak{so}(2N)
\]

的结构常数。

因此：

\[
\boxed{
\{L_{ij}\}
\text{ 生成 }
\mathfrak{so}(2N)
}
\tag{7}
\]

维数：

\[
\dim\mathfrak{so}(2N)
=
N(2N-1)
\tag{8}
\]

---

# 五、真实物理系统的本质

真实系统并不会激活全部：

\[
L_{ij}
\]

而只会激活某个子集：

\[
\mathcal S
=
\{L_{i_1j_1},L_{i_2j_2},\dots\}
\tag{9}
\]

因此：

真正 physics 的核心不是：

\[
\mathfrak{so}(2N)
\]

本身。

而是：

\[
\boxed{
\mathfrak g
=
\mathrm{Lie}(\mathcal S)
}
\tag{10}
\]

即：

Hamiltonian generator 的 dynamical closure。

---

# 六、Dynamical Lie Algebra

定义：

\[
\boxed{
\mathfrak g
=
\mathrm{Lie}\{H(t)\}
}
\tag{11}
\]

即：

最小满足：

\[
[X,Y]\in\mathfrak g
\]

的 closure。

---

# 七、Closure Algorithm

给定：

\[
H(t)=\sum_ah_a(t)G_a
\]

其中：

\[
G_a\in\mathfrak{so}(2N)
\]

---

## Step 1

定义初始集合：

\[
S_0=\{G_a\}
\tag{12}
\]

---

## Step 2

计算所有 commutator：

\[
[G_i,G_j]
\]

---

## Step 3

若产生新的线性独立 generator：

加入集合。

得到：

\[
S_1
\]

---

## Step 4

继续计算：

\[
[S_1,S_1]
\]

---

## Step 5

重复直到：

\[
S_{n+1}=S_n
\]

于是：

\[
\boxed{
\mathfrak g
=
\mathrm{span}(S_n)
}
\tag{13}
\]

---

# 八、Closure 的真正物理意义

closure 不是数学游戏。

它真正描述的是：

\[
\boxed{
\text{Dynamics 实际能够访问的 transport manifold}
}
\]

即：

系统真正的 holonomy geometry。

---

# 九、最小 SU(2) Closure

考虑：

\[
H
=
\lambda_1L_{12}
+
\lambda_2L_{23}
\tag{14}
\]

初始：

\[
S_0=\{L_{12},L_{23}\}
\]

计算：

\[
[L_{12},L_{23}]
=
2iL_{13}
\tag{15}
\]

于是：

\[
S_1
=
\{L_{12},L_{23},L_{13}\}
\]

继续：

\[
[L_{12},L_{13}]
=-2iL_{23}
\]

\[
[L_{23},L_{13}]
=2iL_{12}
\]

closure 停止。

因此：

\[
\boxed{
\mathfrak g
\simeq
\mathfrak{su}(2)
}
\tag{16}
\]

---

# 十、SU(2) 的 Rodrigues / Quaternion Structure

定义：

\[
J_i=\frac12L_i
\]

满足：

\[
[J_i,J_j]=i\epsilon_{ijk}J_k
\tag{17}
\]

若：

\[
H(t)=\vec B(t)\cdot\vec J
\]

则：

\[
U(T)
=
\mathcal P
e^{-i\int dt\,\vec B(t)\cdot J}
\tag{18}
\]

若 closure 不扩张，

则：

\[
\boxed{
U=e^{-i\Theta\hat n\cdot J}
}
\tag{19}
\]

进一步：

\[
\boxed{
U
=
\cos\frac\Theta2
-
2i(\hat n\cdot J)
\sin\frac\Theta2
}
\tag{20}
\]

这对应：

- Bloch rotation
- Rodrigues formula
- quaternion transport

---

# 十一、重要更新：SU(2) 不是普适结构

SU(2)：

不是 Majorana 的本质。

而只是：

\[
\mathfrak{so}(2N)
\]

内部的一个特殊 closure。

因此：

Bloch sphere：

不是 Majorana 的普适几何。

Quaternion：

也不是先验结构。

它们仅对应：

\[
\mathfrak g\simeq\mathfrak{su}(2)
\]

这一 dynamical phase。

---

# 十二、SO(4) Enlargement Example

考虑：

\[
H
=
L_{12}+L_{23}+L_{34}
\tag{21}
\]

closure：

\[
L_{13},L_{24},L_{14}
\]

全部生成。

最终：

\[
\boxed{
\mathfrak g
=
\mathfrak{so}(4)
}
\tag{22}
\]

注意：

\[
\mathfrak{so}(4)
\simeq
\mathfrak{su}(2)\oplus\mathfrak{su}(2)
\tag{23}
\]

因此：

系统不再对应单 quaternion rotation。

而对应：

双 SU(2) transport。

---

# 十三、Magnus Expansion

时间演化：

\[
U(T)=\mathcal P e^{-i\int H(t)dt}
\tag{24}
\]

定义 Magnus expansion：

\[
U=e^{\Omega(T)}
\tag{25}
\]

其中：

---

## 一阶

\[
\Omega_1=-i\int Hdt
\tag{26}
\]

---

## 二阶

\[
\Omega_2
=-\frac12
\int dt_1dt_2
[H(t_1),H(t_2)]
\tag{27}
\]

---

## 三阶

\[
\Omega_3
\sim
[H,[H,H]]
\tag{28}
\]

---

# 十四、Magnus Hierarchy 的真正意义

由于：

\[
H(t)\subseteq\mathfrak g
\]

且：

\[
\mathfrak g
\]

对 commutator closure。

因此：

\[
\boxed{
\Omega_n\in\mathfrak g
}
\tag{29}
\]

于是：

\[
\boxed{
\Omega(T)\in\mathfrak g
}
\tag{30}
\]

因此：

Magnus hierarchy：

实际上就是：

\[
\boxed{
\text{Dynamical Lie algebra growth hierarchy}
}
\]

---

# 十五、Closure Enlargement

考虑：

理想 braid Hamiltonian：

\[
H_0(t)
\]

对应 closure：

\[
\mathfrak g_0
\]

加入 perturbation：

\[
H(t)=H_0(t)+\delta H(t)
\tag{31}
\]

重新定义：

\[
\boxed{
\mathfrak g
=
\mathrm{Lie}\{H_0,\delta H\}
}
\tag{32}
\]

---

# 十六、内部扰动

若：

\[
[X_b,\mathfrak g_0]
\subseteq
\mathfrak g_0
\tag{33}
\]

则：

closure 不扩张。

即：

\[
\mathfrak g=\mathfrak g_0
\tag{34}
\]

此时：

仅发生：

- rotation angle renormalization
- axis deformation
- phase shift

transport manifold 不变。

---

# 十七、真正的 Enlargement

若：

\[
[X_b,\mathfrak g_0]
\not\subseteq
\mathfrak g_0
\tag{35}
\]

则：

commutator 会生成新的 independent generator。

于是：

\[
\boxed{
\mathfrak g
>
\mathfrak g_0
}
\tag{36}
\]

这称为：

\[
\boxed{
\text{Holonomy Enlargement}
}
\]

---

# 十八、PRB105 的统一解释

PRB105：

对应：

\[
\mathfrak g\simeq\mathfrak{su}(2)
\]

因此：

- arbitrary Bloch rotation
- exact Rodrigues solution
- quaternion reduction

成立。

注意：

这并不是 Majorana 的普适性质。

而只是：

一个低维 closure fixed point。

---

# 十九、PRB111 的统一解释

PRB111：

本质上对应：

closure dimension growth。

即：

不同时间段 Hamiltonian：

\[
[H(t_1),H(t_2)]\neq0
\]

Magnus hierarchy 开始真正生成新的 transport direction。

因此：

- fidelity oscillation 更复杂
- path ordering 不再 collapse
- quaternion reduction 失效
- transport manifold enlargement

但：

重要的是：

Quadratic Majorana Hamiltonian 始终 closure 于：

\[
\mathfrak{so}(2N)
\]

因此：

PRB111 的复杂性：

不是离开 quadratic algebra。

而是：

\[
\boxed{
\mathfrak{so}(2N)
\text{ 内部的 closure enlargement}
}
\]

---

# 二十、Dynamical Lie Phase

定义：

\[
\boxed{
\mathcal P(H)
=
\mathrm{Lie}\{H(t)\}
}
\tag{37}
\]

称为：

系统的 dynamical Lie phase。

---

# 二十一、不同 Dynamical Phase

---

## Phase I

\[
\mathfrak g\simeq\mathfrak{su}(2)
\]

对应：

- Bloch sphere
- quaternion transport
- exact solvability

---

## Phase II

\[
\mathfrak g\simeq\mathfrak{so}(4)
\]

对应：

- coupled SU(2)
- multi-axis transport
- generalized rotation manifold

---

## Phase III

\[
\mathfrak g\simeq\mathfrak{so}(6)
\]

对应：

- full six-Majorana transport
- strongly non-Abelian dynamics
- complex fidelity geometry

---

# 二十二、Fidelity Geometry

定义理想演化：

\[
U_0
\]

实际演化：

\[
U
\]

定义：

\[
\Delta U=U_0^\dagger U
\tag{38}
\]

---

## SU(2) 情况

若：

\[
\Delta U=e^{-i\delta\theta\hat n\cdot J}
\]

则：

\[
F=\cos^2\frac{\delta\theta}2
\tag{39}
\]

fidelity 仅对应 rotation mismatch。

---

## Enlargement 情况

若：

\[
\Delta U\notin SU(2)
\]

则：

fidelity 不再由单 rotation angle 控制。

而对应：

更高维 holonomy geometry。

---

# 二十三、最终统一理论

整个理论最终统一为：

\[
Cl(2N)
\]

↓

\[
\mathfrak{so}(2N)
\]

↓

\[
\mathfrak g
=
\mathrm{Lie}\{H(t)\}
\]

↓

分类：

- SU(2)
- SO(4)
- SO(6)
- ...

↓

transport manifold

↓

Magnus hierarchy

↓

fidelity geometry

---

# 二十四、最终核心观点

Bloch sphere：

不是 Majorana 的本质。

Quaternion：

不是先验结构。

SU(2)：

不是普适 algebra。

真正普适的是：

\[
\boxed{
H(t)\subseteq\mathfrak{so}(2N)
}
\]

而真正 physics 的核心是：

\[
\boxed{
\mathfrak g
=
\mathrm{Lie}\{H(t)\}
}
\]

即：

Hamiltonian closure 所定义的 dynamical geometry。

---

# 二十五、真正的研究方向

最终真正需要研究的是：

给定：

\[
H(t)=\sum_{ij}A_{ij}(t)L_{ij}
\]

分类：

\[
\mathfrak g
=
\mathrm{Lie}\{H(t)\}
\]

研究：

- 哪些 graph closure 于 SU(2)
- 哪些 closure 于 SO(4)
- 哪些 closure 于完整 SO(2N)
- 哪些 perturbation 导致 enlargement
- fidelity 如何反映 enlargement
- Magnus expansion 如何诊断 transport geometry
- 是否存在 generalized quaternion structure
- 是否存在新的 exact-solvable Lie phase

---

# 二十六、总结

我们最终得到的，不再只是：

“Majorana braid 的 Bloch sphere 描述”。

而是：

\[
\boxed{
\text{SO(2N) Dynamical Closure Phase Theory}
}
\]

即：

一个从 Clifford algebra 出发，

通过 closure、Magnus hierarchy 与 Lie subalgebra classification，

统一描述：

- braid dynamics
- perturbation
- fidelity
- holonomy
- enlargement
- exact solvability

的完整理论框架。

