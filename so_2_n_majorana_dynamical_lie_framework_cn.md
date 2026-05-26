# 从 SO(2N) 出发的 Majorana Dynamical Lie Algebra 统一框架

# Clifford Algebra → Dynamical Closure → Holonomy Enlargement

---

# 一、研究目标

本文建立一个完全一般化的 Majorana 动力学框架：

不再：

- 预设 qubit
- 预设 Bloch sphere
- 预设 SU(2)
- 预设低能投影

而是：

从：

\[
2N
\]

个 Majorana 的一般 Hamiltonian 出发，

系统地研究：

---

## 1. Clifford algebra 的一般结构

## 2. SO(2N) 的生成

## 3. Dynamical Lie algebra

## 4. Magnus expansion 与 closure

## 5. Quaternion reduction 条件

## 6. Holonomy enlargement

## 7. Fidelity geometry

## 8. PRB105 / PRB111 的统一解释

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
\]

---

# 三、二次 Majorana Hamiltonian 的一般形式

任意 quadratic Majorana Hamiltonian：

\[
\boxed{
H(t)
=
\frac i2
\sum_{i,j}
A_{ij}(t)\gamma_i\gamma_j
}
\tag{2}
\]

其中：

\[
A_{ij}=-A_{ji}
\]

是真反对称矩阵。

定义 generator：

\[
L_{ij}=i\gamma_i\gamma_j
\tag{3}
\]

则：

\[
H(t)=\sum_{i<j}A_{ij}(t)L_{ij}
\tag{4}
\]

---

# 四、SO(2N) 的生成

计算 commutator：

\[
[L_{ij},L_{kl}]
\]

利用 Clifford algebra：

\[
\{\gamma_i,\gamma_j\}=2\delta_{ij}
\]

可得：

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
\tag{5}
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
\tag{6}
\]

维数：

\[
\dim\mathfrak{so}(2N)
=
N(2N-1)
\tag{7}
\]

---

# 五、真实 Hamiltonian 的本质

真实系统：

并不会激活全部：

\[
L_{ij}
\]

而只会激活某个子集：

\[
\mathcal S
=
\{L_{i_1j_1},L_{i_2j_2},\dots\}
\tag{8}
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
\tag{9}
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
\tag{10}
\]

即：

最小满足：

\[
[X,Y]\in\mathfrak g
\]

的 closure。

---

# 七、Closure Algorithm

给定初始集合：

\[
S_0=\{G_a\}
\tag{11}
\]

其中：

\[
H(t)=\sum_ah_a(t)G_a
\]

closure 算法：

---

## Step 1

计算：

\[
[G_i,G_j]
\]

---

## Step 2

若产生新的线性独立 generator：

加入集合。

得到：

\[
S_1
\]

---

## Step 3

继续计算：

\[
[S_1,S_1]
\]

---

## Step 4

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
\tag{12}
\]

---

# 八、最小 SU(2) Closure

考虑：

\[
H
=
\lambda_1L_{12}
+
\lambda_2L_{23}
\tag{13}
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
\tag{14}
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
\tag{15}
\]

---

# 九、SU(2) 的 Rodrigues Structure

定义：

\[
J_i=\frac12L_i
\]

满足：

\[
[J_i,J_j]=i\epsilon_{ijk}J_k
\tag{16}
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
\tag{17}
\]

若 closure 不扩张，

则：

\[
\boxed{
U=e^{-i\Theta\hat n\cdot J}
}
\tag{18}
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
\tag{19}
\]

这就是：

- Bloch rotation
- Rodrigues formula
- unit quaternion

的统一来源。

---

# 十、Quaternion Structure

对于 spin-1/2 representation：

\[
J_i=\frac12\sigma_i
\]

利用：

\[
(\hat n\cdot\sigma)^2=1
\tag{20}
\]

指数自动 collapse：

\[
U
=
\cos\frac\Theta2
-
i(\hat n\cdot\sigma)\sin\frac\Theta2
\tag{21}
\]

这正是单位四元数：

\[
q=a+b\mathbf i+c\mathbf j+d\mathbf k
\]

的标准形式。

因此：

\[
\boxed{
SU(2)
\simeq
\text{unit quaternion}
}
\tag{22}
\]

---

# 十一、Holonomy Enlargement

若 perturbation 满足：

\[
[X_b,\mathfrak g_0]
\subseteq
\mathfrak g_0
\tag{23}
\]

则：

closure 不扩张。

即：

\[
\mathfrak g=\mathfrak g_0
\]

系统仍 confined 于同一 manifold。

---

# 十二、真正的 Enlargement

若：

\[
[X_b,\mathfrak g_0]
\not\subseteq
\mathfrak g_0
\tag{24}
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
\tag{25}
\]

这称为：

\[
\boxed{
\text{Holonomy Enlargement}
}
\]

---

# 十三、SO(4) Enlargement Example

考虑：

\[
H
=
L_{12}+L_{23}+L_{34}
\tag{26}
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
\tag{27}
\]

注意：

\[
\mathfrak{so}(4)
\simeq
\mathfrak{su}(2)\oplus\mathfrak{su}(2)
\tag{28}
\]

因此：

系统不再对应单 quaternion rotation。

而是：

双 SU(2) transport。

---

# 十四、Magnus Expansion 与 Closure

时间演化：

\[
U(T)=\mathcal P e^{-i\int H(t)dt}
\tag{29}
\]

Magnus expansion：

\[
U=e^{\Omega(T)}
\tag{30}
\]

其中：

---

## 一阶

\[
\Omega_1=-i\int Hdt
\tag{31}
\]

---

## 二阶

\[
\Omega_2
=-\frac12
\int dt_1dt_2
[H(t_1),H(t_2)]
\tag{32}
\]

---

## 三阶

\[
\Omega_3
\sim
[H,[H,H]]
\tag{33}
\]

因此：

\[
\boxed{
\text{Magnus expansion 自动生成 Lie closure}
}
\tag{34}
\]

---

# 十五、Closure 的真正物理意义

closure 并不是数学游戏。

它真正描述的是：

\[
\boxed{
\text{Dynamics 真正能够访问的 transport direction}
}
\]

即：

系统 holonomy manifold。

---

# 十六、Fidelity Geometry

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
\tag{35}
\]

---

## SU(2) 情况

若：

\[
\Delta U
=
e^{-i\delta\theta\hat n\cdot J}
\]

则：

\[
F=\cos^2\frac{\delta\theta}2
\tag{36}
\]

fidelity 仅对应：

rotation mismatch。

---

## Enlargement 情况

若 closure 扩张：

则：

\[
\Delta U
\notin
SU(2)
\]

此时：

fidelity 不再由单 rotation angle 控制。

而对应：

更高维 holonomy geometry。

---

# 十七、PRB105 与 PRB111 的统一解释

---

## PRB105

对应：

\[
\mathfrak g\simeq\mathfrak{su}(2)
\]

因此：

- arbitrary Bloch rotation
- quaternion reduction
- exact Rodrigues solution

成立。

---

## PRB111

对应：

partial enlargement。

因此：

- fidelity oscillation 更复杂
- Magnus 高阶真正重要
- path ordering 不再 collapse
- transport manifold 扩张

---

# 十八、真实物理系统中的一般情况

真实 Majorana 系统通常包含：

---

## 1. 长程 overlap

\[
L_{14},L_{15},\dots
\]

---

## 2. 多 ABS coupling

---

## 3. disorder-induced hybridization

---

## 4. noisy time dependence

---

## 5. nonadiabatic correction

这些通常导致：

\[
\boxed{
\mathfrak g
\rightarrow
\mathfrak{so}(2N)
}
\]

即：

full holonomy enlargement。

---

# 十九、最终统一理论

因此：

Majorana dynamics 的真正结构为：

\[
\boxed{
Cl(2N)
\rightarrow
\mathfrak{so}(2N)
\rightarrow
\mathfrak g
\rightarrow
\text{Holonomy manifold}
}
\tag{37}
\]

其中：

\[
\mathfrak g
=
\mathrm{Lie}\{H(t)\}
\]

决定：

- exact solvability
- quaternion reduction
- Bloch sphere existence
- fidelity geometry
- transport topology
- holonomy enlargement

---

# 二十、真正的研究方向

最终真正需要研究的是：

---

## 给定：

\[
H(t)=\sum_{ij}A_{ij}(t)L_{ij}
\]

---

## 分类：

\[
\mathfrak g
=
\mathrm{Lie}\{H(t)\}
\]

---

## 研究：

- 哪些 graph closure 于 SU(2)
- 哪些 closure 于 SO(4)
- 哪些 closure 于完整 SO(2N)
- 哪些 perturbation 导致 enlargement
- fidelity 如何反映 enlargement
- Magnus expansion 如何诊断 transport geometry

---

# 二十一、核心结论

Bloch sphere：

不是 Majorana 的本质。

Quaternion：

不是先验假设。

SU(2)：

不是普适结构。

真正普适的是：

\[
\boxed{
H(t)\subset\mathfrak{so}(2N)
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

