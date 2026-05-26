# Full Majorana–Magnus–Holonomy Model V2

# 从 \(2N\) 个 Majorana 出发的完整统一框架（扩展版）

---

# 一、研究目标

我们希望建立一个完全统一的框架，用于研究：

- Majorana braiding
- 动力学扰动
- ABS/MZM 混合
- fidelity decay
- leakage
- non-Abelian holonomy
- dynamical Lie algebra enlargement

核心要求：

- 不预设 qubit
- 不预设 Bloch sphere
- 不预设 SU(2)
- 直接从 \(2N\) 个 Majorana 出发
- 用 Hamiltonian closure 自动决定系统的有效群结构

最终目标是建立：

\[
H(t)
\rightarrow
\mathfrak g
\rightarrow
G=e^{\mathfrak g}
\rightarrow
U(T)
\rightarrow
F(T)
\]

之间的完整对应关系。

---

# 二、\(2N\) 个 Majorana 与 Clifford algebra

定义：

\[
\gamma_i,
\quad i=1,\dots,2N
\]

满足：

\[
\{\gamma_i,\gamma_j\}=2\delta_{ij}
\tag{1}
\]

它们生成 Clifford algebra：

\[
Cl(2N)
\]

---

# 三、Majorana bilinear generator

定义：

\[
L_{ij}
=
\frac{i}{2}\gamma_i\gamma_j
\quad (i<j)
\tag{2}
\]

这些 generator 满足：

\[
[L_{ij},L_{kl}]
=
i(
\delta_{jk}L_{il}
-\delta_{ik}L_{jl}
-\delta_{jl}L_{ik}
+\delta_{il}L_{jk}
)
\tag{3}
\]

这正是：

\[
\mathfrak{so}(2N)
\]

的 Lie algebra。

因此：

\[
\boxed{
\mathrm{span}\{L_{ij}\}
\simeq
\mathfrak{so}(2N)
}
\tag{4}
\]

对应群：

\[
Spin(2N)
\]

---

# 四、完整 quadratic Majorana Hamiltonian

一般 Majorana Hamiltonian：

\[
H(t)
=
\sum_{i<j}
A_{ij}(t)L_{ij}
\tag{5}
\]

其中：

\[
A_{ij}(t)=-A_{ji}(t)
\]

因此：

\[
A(t)\in\mathfrak{so}(2N)
\]

于是：

\[
H(t)\in\mathfrak{so}(2N)
\]

---

# 五、时间演化与 Spin(2N)

系统时间演化：

\[
U(T)
=
\mathcal P
\exp\left(
-i\int_0^Tdt\,H(t)
\right)
\tag{6}
\]

由于：

\[
H(t)\in\mathfrak{so}(2N)
\]

因此：

\[
\boxed{
U(T)\in Spin(2N)
}
\tag{7}
\]

因此：

> Majorana dynamics 本质上是 Spin(2N) 上的 holonomy transport。

---

# 六、dynamical Lie algebra

定义：

\[
\boxed{
\mathfrak g
=
\mathrm{Lie}\{-iH(t)\}
}
\tag{8}
\]

即：

不断计算：

\[
[H(t_1),H(t_2)]
\]

以及更高 nested commutator，直到 closure。

这决定：

- 系统真正访问的 transport manifold
- 有效群结构
- 是否发生 holonomy enlargement
- 是否仍可 reduction 到 SU(2)

注意：

\[
\mathfrak g
\subseteq
\mathfrak{so}(2N)
\]

但一般不等于整个 \(\mathfrak{so}(2N)\)。

---

# 七、Magnus expansion

定义：

\[
U(T)=e^{\Omega(T)}
\tag{9}
\]

其中：

\[
\Omega
=
\Omega_1+\Omega_2+\Omega_3+\cdots
\]

这是 exact expansion。

---

# 八、一阶 Magnus

\[
\Omega_1
=
-i\int_0^Tdt\,H(t)
\tag{10}
\]

代入：

\[
H(t)
=
\sum_a h_a(t)G_a
\]

其中：

\[
G_a=L_{ij}
\]

得到：

\[
\Omega_1
=
-i\sum_a
\left(
\int_0^Tdt\,h_a(t)
\right)
G_a
\tag{11}
\]

关键：

> 一阶 Magnus 不生成新 generator。

它只保留 Hamiltonian 中已有方向。

---

# 九、二阶 Magnus：closure 的真正来源

\[
\Omega_2
=
-\frac12
\int_0^Tdt_1
\int_0^{t_1}dt_2
[H(t_1),H(t_2)]
\tag{12}
\]

代入：

\[
H(t)
=
\sum_a h_a(t)G_a
\]

得到：

\[
\Omega_2
=
-\frac12
\sum_{ab}
\int dt_1dt_2
h_a(t_1)h_b(t_2)
[G_a,G_b]
\tag{13}
\]

因此：

> Magnus 二阶开始自动生成 dynamical Lie closure。

若：

\[
[G_a,G_b]
\notin
\mathrm{span}\{G_a\}
\]

则：

> 系统开始访问新的 transport direction。

---

# 十、三阶 Magnus

\[
\Omega_3
=
\frac16
\int dt_1dt_2dt_3
\Big(
[H_1,[H_2,H_3]]
+
[[H_1,H_2],H_3]
\Big)
\tag{14}
\]

这里：

nested commutator 会继续生成更高 closure。

因此：

\[
\boxed{
\text{Magnus expansion}
\iff
\text{continuous-time Lie closure}
}
\]

---

# 十一、扰动的统一定义

设：

\[
H(t)
=
H_0(t)+\delta H(t)
\tag{15}
\]

其中：

\[
H_0(t)
=
\sum_a\lambda_a(t)G_a
\]

定义原 closure：

\[
\mathfrak g_0
=
\mathrm{Lie}\{G_a\}
\tag{16}
\]

扰动：

\[
\delta H(t)
=
\sum_b\epsilon_b(t)X_b
\tag{17}
\]

其中：

\[
X_b\in\mathfrak{so}(2N)
\]

---

# 十二、完整 algebraic 判据

## 情况 A：internal deformation

若：

\[
X_b\in\mathfrak g_0
\tag{18}
\]

则：

\[
\boxed{
\mathfrak g=\mathfrak g_0
}
\tag{19}
\]

因此：

- 不发生群扩张
- 只发生内部 transport deformation
- 只改变 effective rotation axis
- holonomy confined 于原 manifold

若：

\[
\mathfrak g_0=\mathfrak{su}(2)
\]

则系统仍可 reduction 到 Bloch sphere。

---

## 情况 B：commuting extension

若：

\[
[X_b,\mathfrak g_0]=0
\tag{20}
\]

则：

\[
\boxed{
\mathfrak g
=
\mathfrak g_0\oplus\mathfrak u(1)
}
\tag{21}
\]

因此：

- perturbation 不改变 braid sector
- 只增加额外 dynamical phase
- 不发生 leakage

---

## 情况 C：holonomy enlargement

若：

\[
[X_b,\mathfrak g_0]
\not\subseteq
\mathfrak g_0
\tag{22}
\]

则：

Magnus 二阶开始生成新的 generator。

于是：

\[
\boxed{
\mathfrak g>\mathfrak g_0
}
\tag{23}
\]

因此：

- 系统离开原 transport manifold
- trajectory explores larger group manifold
- fidelity decay 不再只是 phase mismatch
- 出现真正 holonomy leakage

---

# 十三、holonomy enlargement 的严格定义

定义：

原 transport manifold：

\[
G_0=e^{\mathfrak g_0}
\]

实际 transport manifold：

\[
G=e^{\mathfrak g}
\]

若：

\[
\mathfrak g>\mathfrak g_0
\]

则：

\[
G_0\subset G
\]

因此：

trajectory 从原 subgroup orbit 漂移到更大 manifold。

这就是：

\[
\boxed{
\text{holonomy enlargement}
}
\]

---

# 十四、fidelity 的统一定义

定义 ideal braid：

\[
U_B
\]

实际演化：

\[
U(T)=e^{\Omega(T)}
\]

定义 deviation：

\[
\Delta U
=
U_B^\dagger U(T)
\tag{24}
\]

fidelity：

\[
F(T)
=
|\langle\psi_0|\Delta U|\psi_0\rangle|^2
\tag{25}
\]

---

# 十五、SU(2) closure 情况

若：

\[
\Omega(T)\in\mathfrak{su}(2)
\tag{26}
\]

则存在：

\[
\delta\theta(T),\hat n(T)
\]

使：

\[
\Delta U
=
\exp\left(
-i\delta\theta\,\hat n\cdot\sigma
\right)
\tag{27}
\]

因此：

\[
\boxed{
F(T)
=
\cos^2\frac{\delta\theta(T)}2
}
\tag{28}
\]

---

# 十六、non-SU(2) enlargement 情况

若：

\[
\Omega(T)
\notin
\mathfrak{su}(2)
\tag{29}
\]

则不存在单 angle reduction。

因此：

- fidelity 不再只是 rotation mismatch
- 系统离开单 Bloch manifold
- fidelity decay 对应 transport manifold leakage

---

# 十七、PRB105 与本框架

PRB105 的核心结果：

- arbitrary Bloch rotation
- ABS–MZM coupling
- fidelity oscillation

在本框架中：

- effective closure confined 于 SU(2)
- perturbation 只 renormalize effective rotation axis
- Magnus expansion 不生成新的 independent transport direction

因此：

\[
\mathfrak g=\mathfrak{su}(2)
\]

observable 表现为：

\[
F(T)=\cos^2\frac{\delta\theta(T)}2
\]

---

# 十八、PRB111 与本框架

PRB111 的核心现象：

- fidelity oscillation
- decay
- revival
- ABS fluctuation induced degradation

在本框架中：

- near-zero ABS 对应 commuting 或 near-commuting perturbation
- strong ABS fluctuation 会生成新的 Clifford direction
- Magnus commutator 开始 enlarges closure
- fidelity decay 来源于 holonomy drift

因此：

\[
\boxed{
\text{perturbation-induced holonomy enlargement}
}
\]

---

# 十九、full algebra 与 dynamical algebra 的严格区别

这是整个框架最关键的结构之一。

必须严格区分：

- full algebra
- dynamical algebra

---

## 1. full algebra

对于：

\[
2N
\]

个 Majorana，

所有 bilinear：

\[
L_{ij}
=
\frac{i}{2}\gamma_i\gamma_j
\]

共同张成：

\[
\mathfrak{so}(2N)
\]

这对应：

> 系统允许的全部 quadratic Majorana direction。

即整个 Clifford bilinear 空间。

---

## 2. dynamical algebra

真实 Hamiltonian：

\[
H_0(t)
=
\sum_a\lambda_a(t)G_a
\]

通常只包含部分 generator。

因此 dynamics 实际访问的 manifold 不是整个：

\[
\mathfrak{so}(2N)
\]

而是：

\[
\boxed{
\mathfrak g_0
=
\mathrm{Lie}\{H_0(t)\}
}
\tag{30}
\]

即 Hamiltonian 的 commutator closure。

---

# 二十、最直观例子

考虑 4 个 Majorana：

\[
\gamma_1,\gamma_2,\gamma_3,\gamma_4
\]

full algebra：

\[
\mathfrak{so}(4)
\]

维数为：

\[
6
\]

对应 generator：

\[
L_{12},L_{13},L_{14},L_{23},L_{24},L_{34}
\]

现在只取：

\[
H_0
=
\lambda_1L_{12}
+
\lambda_2L_{23}
\]

则：

\[
[L_{12},L_{23}]
=
iL_{13}
\]

因此：

\[
\mathfrak g_0
=
\mathrm{span}
\{L_{12},L_{23},L_{13}\}
\simeq
\mathfrak{su}(2)
\tag{31}
\]

注意：

- full algebra 仍然是：

\[
\mathfrak{so}(4)
\]

- 但 dynamics 只访问：

\[
\mathfrak{su}(2)
\]

因此：

> effective geometry 由 dynamical closure 决定，而不是由 Majorana 数量直接决定。

---

# 二十一、为什么需要讨论 \(X_b\in\mathfrak g_0\)

虽然：

\[
X_b\in\mathfrak{so}(2N)
\]

说明它是合法 quadratic Majorana operator。

但：

并不意味着：

\[
X_b\in\mathfrak g_0
\]

若：

\[
X_b\in\mathfrak g_0
\]

则：

> perturbation 只是原 transport manifold 内部 deformation。

不会生成新的 transport direction。

---

# 二十二、为什么需要讨论 \([X_b,\mathfrak g_0]\)

真正时间演化：

\[
U(T)
=
\mathcal P e^{-i\int Hdt}
\]

不是简单线性叠加。

由于 time-ordering：

Magnus expansion 自动生成：

\[
\Omega_2
\sim
[H(t_1),H(t_2)]
\]

因此：

即使：

\[
X_b
\notin
\mathfrak g_0
\]

也不一定立即导致 enlargement。

真正关键是：

\[
[X_b,\mathfrak g_0]
\]

是否生成新的 independent direction。

---

# 二十三、真正要检查的对象

真正需要计算的是：

\[
\boxed{
\mathfrak g
=
\mathrm{Lie}
\{\mathfrak g_0,X_b\}
}
\tag{32}
\]

即：

将 perturbation 加入后，closure 是否扩大。

---

# 二十四、几何解释

\[
\mathfrak{so}(2N)
\]

对应整个 transport manifold 的全部 tangent direction。

而：

\[
\mathfrak g_0
\]

对应原 dynamics 的 tangent distribution。

perturbation：

\[
X_b
\]

对应一个新的 vector field。

若：

\[
X_b\in\mathfrak g_0
\]

则：

新的 vector field 仍位于原 distribution 内。

不会 enlargement。

若：

\[
[X_b,\mathfrak g_0]
\not\subseteq
\mathfrak g_0
\]

则：

Lie bracket 会生成新的 tangent direction。

于是 distribution 扩张。

因此：

> holonomy enlargement 本质上是 reachable manifold enlargement。

---

# 二十五、Magnus 与 Lie closure 的本质关系

Magnus expansion：

\[
\Omega_2
\sim
[H_1,H_2]
\]

\[
\Omega_3
\sim
[H,[H,H]]
\]

会自动生成 dynamical Lie closure。

因此：

\[
\boxed{
\text{Lie enlargement}
\iff
\text{time-ordering induced closure generation}
}
\tag{33}
\]

这不是人为定义。

而是：

> 时间演化自身自动生成的 transport geometry。

---

# 二十六、最终统一结构

完整理论链条：

\[
Cl(2N)
\rightarrow
\mathfrak{so}(2N)
\rightarrow
H(t)
\rightarrow
\mathfrak g_0
\rightarrow
X_b
\rightarrow
\mathrm{Lie}\{\mathfrak g_0,X_b\}
\rightarrow
\Omega(T)
\rightarrow
U(T)
\rightarrow
F(T)
\]

其中：

- \(\mathfrak{so}(2N)\)：全部允许 generator
- \(\mathfrak g_0\)：当前 protocol 实际访问 subgroup
- perturbation：新的 vector field
- commutator closure：决定是否出现新 reachable direction
- Magnus expansion：自动生成这些新的 transport direction
- fidelity decay：对应 actual holonomy 离开原 subgroup manifold

---

# 二十七、下一步真正需要做的事情

1. 对具体 braid protocol 写出 explicit \(H(t)\)
2. 做 Magnus expansion 到二阶/三阶
3. 显式求 closure enlargement
4. 建立 physical perturbation dictionary
5. 建立：

\[
\text{physical noise}
\rightarrow
\text{Lie enlargement}
\rightarrow
\text{fidelity decay}
\]

的严格映射
6. 对照 PRB105 / PRB111 数值结果
7. 研究 Spin(2N) manifold 上的几何控制理论

