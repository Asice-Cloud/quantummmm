
# PRB111 / Cl(5) → so(5) → Sp(2) → Quaternion Riccati（完全展开严格版）

> 本文档目标：  
> 在**不省略定义、不隐含投影、不跳步**的前提下，从 Cl(5) + so(5) 全代数动力学出发，  
> 明确区分：
> - 必然结构（algebraic identities）
> - 表示选择（representation choice）
> - 坐标选择（chart / Riccati projection）

---

# 0. 基本对象与目标

我们从最原始对象开始：

\[
H(t)\in \mathfrak{so}(5)\subset \mathrm{Cl}(5)
\]

目标是系统说明：

1. so(5) 全空间动力学
2. 是否以及如何出现 Sp(2) 表示
3. quaternion 结构从何而来
4. Riccati 变量何时合法出现

---

# 1. Clifford代数 Cl(5)

## 1.1 Gamma矩阵定义

取 Cl(5) 生成元：

\[
\{\Gamma_i\}_{i=1}^5,\quad
\{\Gamma_i,\Gamma_j\}=2\delta_{ij}
\]

---

## 1.2 so(5)生成元

定义：

\[
L_{ij}=\frac14[\Gamma_i,\Gamma_j],\quad 1\le i<j\le 5
\]

因此：

- 维数：10
- 基底：\(L_{12}, L_{13}, ..., L_{45}\)

---

## 1.3 Lie代数结构

满足：

\[
[L_{ij},L_{kl}]
= \delta_{jk}L_{il}-\delta_{ik}L_{jl}
-\delta_{jl}L_{ik}+\delta_{il}L_{jk}
\]

👉 这是 so(5) 的完整结构，不依赖任何表示。

---

# 2. PRB111 驱动项（抽象层）

考虑论文选取的五个驱动：

\[
\mathcal{D}=
\{L_{12}, L_{2a}, L_{1b}, L_{3a}, L_{ab}\}
\]

其中：

- \(a=4, b=5\)

即：

\[
L_{2a}=L_{24},\quad
L_{1b}=L_{15},\quad
L_{3a}=L_{34},\quad
L_{ab}=L_{45}
\]

---

# 3. 全空间动力学（无任何分解）

## 3.1 Hamiltonian

\[
H(t)=\sum_{i<j} h_{ij}(t)L_{ij}
\]

---

## 3.2 演化方程

\[
\dot U(t)=-iH(t)U(t),\quad U(t)\in \mathrm{Spin}(5)
\]

---

## 3.3 adjoint形式（纯代数）

对任意 \(X\in \mathfrak{so}(5)\)：

\[
\dot X=[-iH(t),X]
\]

---

## ✔ 第一层结论

此层：

- 无 block
- 无 quaternion
- 无 projection

👉 完全全空间

---

# 4. Sp(2)结构：不是必然，而是同构选择

## 4.1 关键事实（Lie代数同构）

\[
\mathfrak{so}(5)\cong \mathfrak{sp}(2)
\]

但注意：

> 这是抽象同构，不是坐标变换

---

## 4.2 Sp(2)定义

\[
\mathfrak{sp}(2)=
\{X\in \mathbb{C}^{4\times4}\mid X^\dagger \Omega + \Omega X=0\}
\]

其中：

\[
\Omega=
\begin{pmatrix}
0 & I\\
-I & 0
\end{pmatrix}
\]

---

## 4.3 block结构（关键但仍为表示）

任意元素可写：

\[
X=
\begin{pmatrix}
A & B\\
-B^\dagger & C
\end{pmatrix}
\]

其中：

- \(A,C\) ：2×2 anti-Hermitian
- \(B\) ：2×2 complex matrix

---

## ✔ 重要说明

这一结构：

- 不是 so(5) intrinsic decomposition
- 是 Sp(2) representation choice

---

# 5. quaternion结构来源（第二次选择）

## 5.1 quaternion基底嵌入

定义：

\[
\mathbf{i} = i\sigma_x,\quad
\mathbf{j} = i\sigma_y,\quad
\mathbf{k} = i\sigma_z
\]

---

## 5.2 2×2复矩阵 ↔ quaternion

任意：

\[
q = a_0 I + a_1 (i\sigma_x)+a_2(i\sigma_y)+a_3(i\sigma_z)
\]

对应 quaternion：

\[
q = a_0 + a_1\mathbf{i}+a_2\mathbf{j}+a_3\mathbf{k}
\]

---

## ✔ 结论

quaternion 不是 Cl(5) 必然结构，而是：

> Sp(2) 表示 + 2×2矩阵 ↔ H 同构

---

# 6. so(5) → Sp(2) 的一般嵌入（关键严格点）

## 6.1 任意生成元展开

\[
H(t)=\sum h_{ij}(t)L_{ij}
\]

在 Sp(2) 表示中：

\[
H(t)\mapsto
\begin{pmatrix}
A(h) & B(h)\\
-B(h)^\dagger & C(h)
\end{pmatrix}
\]

---

## 6.2 映射性质（非常重要）

这里：

- \(A(h)\)：由某些 \(h_{ij}\) 线性组合
- \(B(h)\)：由其余线性组合
- \(C(h)=-A(h)\)（sp(2)条件）

👉 但**具体对应关系依赖 representation choice**

---

## ✔ 关键结论

不能唯一写出：

\[
L_{ij} \rightarrow \text{fixed quaternion matrix}
\]

除非指定 Gamma representation

---

# 7. Riccati变量：唯一非线性来源

## 7.1 block演化

设：

\[
U=
\begin{pmatrix}
X & Y\\
Z & W
\end{pmatrix}
\]

---

## 7.2 定义 Riccati变量（局部）

\[
q = YX^{-1}
\]

条件：

- \(X(t)\) 可逆
- 局部 big cell

---

## 7.3 推导（标准结果）

由：

\[
\dot U=-iHU
\]

得到：

\[
\boxed{
\dot q = B + A q - q C - q B^\dagger q
}
\]

---

## ✔ 重要说明

这一公式：

- 是纯 block algebra 推导
- 不依赖 quaternion
- 不依赖 so(5) 特性

---

# 8. quaternion Riccati出现的最后一步

## 8.1 quaternion识别

若进一步选：

\[
B(h)\in \mathbb{R}^4 \cong \mathbb{H}
\]

则：

\[
q \in \mathbb{H}
\]

---

## 8.2 得到 quaternion形式

\[
\boxed{
\dot q
=
B
+ A q - q A
- q B^\dagger q
}
\]

（利用 sp(2) 中 \(C=-A\)）

---

## 8.3 展开为实四元数形式

\[
q = q_0 + q_1 \mathbf{i}+q_2 \mathbf{j}+q_3 \mathbf{k}
\]

非线性项：

\[
q q^\dagger q = |q|^2 q
\]

---

# 9. 三层结构总结（核心）

## 第一层（必然）

- so(5) Lie algebra
- Spin(5) flow

---

## 第二层（表示选择）

- so(5) ≅ sp(2)
- block matrix form

---

## 第三层（坐标选择）

- Riccati变量
- quaternion识别

---

# 10. 最终核心结论

\[
\boxed{
\text{Quaternion Riccati 不是 Cl(5) intrinsic结构}
}
\]

而是：

\[
\boxed{
\text{Cl(5) flow + Sp(2) representation + local chart}
}
\]

的组合结果。

