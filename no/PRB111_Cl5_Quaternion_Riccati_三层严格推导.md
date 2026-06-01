
# PRB111 / Cl(5) → so(5) → Quaternion Riccati 三层严格推导（无歧义版本）

## 0. 总体目标

从
\[
H(t)\in \mathfrak{so}(5)\subset \mathrm{Cl}(5)
\]
出发，在**不引入任何先验投影假设**的前提下，分三层构造：

- 第一层：全空间 Spin(5) 动力学
- 第二层：Lie 代数结构层（不选坐标）
- 第三层：在结构选择下出现 quaternion Riccati

---

# 第一层：全空间动力学（无分解）

## 1.1 Cl(5) Clifford 结构

取 Gamma 矩阵：
\[
\{\Gamma_i\}, i=1,\dots,5,\quad \{\Gamma_i,\Gamma_j\}=2\delta_{ij}
\]

so(5) 生成元：
\[
L_{ij}=\frac{1}{4}[\Gamma_i,\Gamma_j]
\]

---

## 1.2 PRB111 哈密顿全展开

\[
H(t)=\sum_{i<j} h_{ij}(t)L_{ij}
\]

关注驱动项：
\[
L_{12}, L_{2a}, L_{1b}, L_{3a}, L_{ab}
\]

---

## 1.3 Spin(5) 全局动力学

\[
\dot U(t)=-iH(t)U(t),\quad U(t)\in \mathrm{Spin}(5)
\]

等价 adjoint 表达：
\[
\dot X=[-iH(t),X],\quad X\in \mathfrak{so}(5)
\]

---

## ✔ 第一层结论

无 quaternion、无 Riccati、无投影。

---

# 第二层：Lie 代数结构层（不选表示）

## 2.1 Lie algebra 同构

\[
\mathfrak{so}(5)\cong \mathfrak{sp}(2)
\]

但此处仅为抽象同构，不引入坐标。

---

## 2.2 Cartan 分解（未指定）

\[
\mathfrak{so}(5)=\mathfrak{k}\oplus\mathfrak{p}
\]

但：
- 未指定 \(\mathfrak{k}\)
- 未指定 embedding
- 未指定 representation

---

## 2.3 生成元结构分类

任意 \(L_{ij}\) 仅满足：

- compact-like:
\[
[\mathfrak{k},\mathfrak{k}]\subset \mathfrak{k}
\]

- mixed:
\[
[\mathfrak{k},\mathfrak{p}] \subset \mathfrak{k}\oplus\mathfrak{p}
\]

---

## ✔ 第二层结论

仍无 quaternion / Riccati。

---

# 第三层：结构选择 → Quaternion Riccati

## 3.1 Sp(2) 表示选择

选择 symplectic representation：

\[
X=
\begin{pmatrix}
A & B\\
-B^\dagger & C
\end{pmatrix}
\]

---

## 3.2 Quaternion 识别

\[
\mathbb{R}^4 \cong \mathbb{H}
\]

令：
\[
B \leftrightarrow q\in \mathbb{H}
\]

---

## 3.3 Big cell Riccati 变量

局部定义：
\[
q = YX^{-1}
\]

（需 \(X\) 可逆）

---

## 3.4 动力学推导

由：
\[
\dot U=-iHU
\]

得到 block 方程：

\[
\boxed{
\dot q = B + A q - q C - q B^\dagger q
}
\]

---

## 3.5 quaternion 化形式

取：
\[
A=\alpha \mathbf{k}+\beta \mathbf{i},\quad B=q
\]

得到：

\[
\boxed{
\dot q
=
q
+(\alpha \mathbf{k}+\beta \mathbf{i})q
+q(\alpha \mathbf{k}+\beta \mathbf{i})
-|q|^2 q
}
\]

---

# 总结

## 核心结构链

\[
\text{Cl(5) flow}
\rightarrow
\text{Sp(2) polarization}
\rightarrow
\text{local chart}
\rightarrow
q=YX^{-1}
\rightarrow
\text{quaternion Riccati}
\]

---

## 核心结论

- Riccati **不是 Cl(5) 必然结构**
- Riccati **依赖 symplectic 分裂 + 局部坐标**
- quaternion 来自 representation choice

