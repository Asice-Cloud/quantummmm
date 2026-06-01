# 从 2N Majorana 到 Lie closure、Wei–Norman 分解与 so(5)（PRB111框架对照）

> 本文档整理一个从 **2N Majorana 算符体系 → Lie 代数闭包 → Wei–Norman 分解 → so(5)结构 → 四元数关系 → ODE构建** 的完整推导链条，并说明如何对照 PRB111 类论文中的参数化方式进行映射。

---

# 1. 2N Majorana 体系的基本结构

## 1.1 Majorana算符定义
考虑 \(2N\) 个 Majorana 算符：

\[
\gamma_i = \gamma_i^\dagger, \quad \{\gamma_i,\gamma_j\}=2\delta_{ij}, \quad i=1,...,2N
\]

它们生成 Clifford 代数：Cl(2N)

---

## 1.2 二次型哈密顿量（一般形式）

所有 Majorana 二次型 Hamiltonian：
\[
H(t)=\frac{i}{4}\sum_{i,j=1}^{2N} A_{ij}(t)\,\gamma_i\gamma_j
\]

其中：

- \(A(t) = -A^T(t)\)
- 实反对称矩阵

---

## 1.3 生成李代数结构
定义生成元：

\[
X_{ij} = \frac{i}{2}\gamma_i\gamma_j, \quad i<j
\]

则：

\[
[X_{ij},X_{kl}] = i(\delta_{jk}X_{il} - \delta_{ik}X_{jl} - \delta_{jl}X_{ik} + \delta_{il}X_{jk})
\]

因此：

\[
\mathfrak{g} \cong so(2N)
\]

---

# 2. Lie closure（李闭包）

## 2.1 闭包定义
给定生成集合：

\[\mathcal{S} = \{X_{ij}\}\]

闭包定义：
\[
\mathrm{Lie}(\mathcal{S}) = \text{所有反复对易生成的线性空间}
\]

---

## 2.2 在 Majorana 系统中的意义

- 任意二次 Majorana Hamiltonian → so(2N)
- 若加限制耦合结构 → 子代数（如 so(5)）

---

# 3. 从 so(2N) 降维到 so(5)（关键步骤）

## 3.1 选择子空间
在 PRB111 类模型中，通常选择：

- 5 个有效生成元（或 5 个独立 Majorana bilinear组合）

记为：

\[
\{B_1, B_2, B_3, B_4, B_5\}
\]

---

## 3.2 闭合条件
要求：

\[
[B_i, B_j] = \sum_k c_{ij}^k B_k
\]

即闭合在 5 维空间

---

## 3.3 so(5) 判定
若结构常数满足 so(5) 李代数结构，则：

\[
\mathfrak{g}_{eff} \cong so(5)
\]

---

# 4. so(5) 与 Clifford / 四元数结构

## 4.1 Clifford 代数关系

\[
so(5) \cong sp(2) \cong usp(4)
\]

而：

- quaternion algebra \(\mathbb{H}\)
- Clifford algebra Cl(0,4)

之间存在同构关系

---

## 4.2 四元数表示
定义四元数基：

\[
1, \; i, j, k
\]

扩展到 5 维时：

- 一个标量 + 四元数向量
- 或 SU(2) × SU(2) 表示结构

---

## 4.3 物理意义

so(5) 常用于：

- SO(5) superspin
- Majorana零模耦合
- 拓扑态有效描述

---

# 5. Wei–Norman 分解

## 5.1 基本思想
时间演化算符：

\[
U(t) = \mathcal{T} e^{-i \int H(t) dt}
\]

写成李群分解：

\[
U(t)=\prod_{k=1}^{n} e^{x_k(t) B_k}
\]

---

## 5.2 so(5) 参数化
设：

\[
U(t)=e^{x_1(t)B_1} e^{x_2(t)B_2} \cdots e^{x_5(t)B_5}
\]

---

## 5.3 导数结构

\[
\dot{U}U^{-1} = \sum_k \dot{x}_k \; \mathrm{Ad}_{e^{x_1 B_1}...e^{x_{k-1} B_{k-1}}}(B_k)
\]

---

## 5.4 得到 ODE 系统

与哈密顿量匹配：

\[
\dot{U}U^{-1} = -i H(t)
\]

得到：

\[
\sum_k M_k(x) \dot{x}_k = h_k(t)
\]

其中：

- \(M_k(x)\)：由 Adjoint representation 给出
- \(h_k(t)\)：Hamiltonian 投影

---

# 6. Adjoint 表示与 ad(Bi)

## 6.1 定义

\[
\mathrm{ad}_{B_i}(X) = [B_i, X]
\]

---

## 6.2 矩阵表示
在基 \(\{B_i\}\) 上：

\[
\mathrm{ad}(B_i)_{jk} = c_{ij}^k
\]

---

## 6.3 Wei–Norman核心

\[
M(x) = \prod e^{x_i \mathrm{ad}(B_i)}
\]

---

# 7. PRB111 框架映射（通用结构说明）

> 注意：以下为**通用 Majorana-driven so(5) 模型映射方式**，用于对照 PRB111 类论文变量结构。

## 7.1 Hamiltonian结构映射
一般论文写作形式：

\[
H(t) = \sum_i E_i(t) B_i + \sum_j t_j(t) B_j
\]

对应：

- \(E_i(t)\)：局域场 / onsite terms
- \(t_j(t)\)：耦合项

---

## 7.2 生成元选取策略
PRB111类结构通常：

- 从 Majorana bilinears中选5个独立组合
- 通过对易验证闭合

---

## 7.3 so(5)识别步骤

1. 计算所有 \([B_i,B_j]\)
2. 写成线性组合
3. 提取结构常数 tensor
4. 与 so(5)标准结构比对

---

# 8. ODE最终形式（Wei–Norman核心结果）

最终系统：

\[
\dot{x}(t) = M^{-1}(x) h(t)
\]

其中：

- \(x(t) = (x_1,...,x_5)\)
- \(h(t)\) 来自 Hamiltonian 投影

---

# 9. 数值计算流程（实现路径）

## Step 1
构造 Majorana 表示

## Step 2
生成 so(2N) bilinears

## Step 3
选取闭合子代数（so(5)）

## Step 4
计算 ad(B_i)

## Step 5
构造 M(x)

## Step 6
解 ODE

---

# 10. 总结

完整链条：

\[
\text{Majorana} \rightarrow so(2N) \rightarrow \text{closure} \rightarrow so(5) \rightarrow \text{quaternion structure} \rightarrow \text{Wei–Norman ODE}
\]

---

# 如果下一步继续
可以进一步严格展开：

- 明确 5 个 \(B_i\) 的具体 Majorana 表达
- 完整计算 ad(B_i) 10×10 或 5×5 矩阵
- 按 PRB111 原式逐项匹配 \(E(t), t(t)\)
- 写出显式 M(x)
- 推导可计算 ODE 系统

