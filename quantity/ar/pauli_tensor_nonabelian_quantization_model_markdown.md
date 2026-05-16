# Pauli张量积表示下的 non-Abelian 量化模型

---

# 1. 基础Hilbert空间

定义单体空间：

\[
V=\mathbb C^2
\]

定义全空间：

\[
\mathcal H_N=(\mathbb C^2)^{\otimes N}
\]

Pauli矩阵：

\[
\sigma^0=I,
\quad
\sigma^x,
\quad
\sigma^y,
\quad
\sigma^z
\]

二元算符空间：

\[
\mathrm{End}(V\otimes V)
\]

的标准基：

\[
\boxed{
\{\sigma^\alpha\otimes\sigma^\beta\}
}
\]

其中：

\[
\alpha,\beta\in\{0,x,y,z\}
\]

因此任意二元Hamiltonian都可展开为：

\[
\boxed{
h(u)
=
\sum_{\alpha,\beta}
h_{\alpha\beta}(u)
\sigma^\alpha\otimes\sigma^\beta
}
\]

这里：

- \(u\)：路径参数 / braid parameter
- \(h_{\alpha\beta}(u)\)：路径依赖耦合

---

# 2. 局域路径演化

定义局域演化算符：

\[
\boxed{
R(u)
=
\mathcal T
\exp
\left(
-i\int_0^u h(s)ds
\right)
}
\]

其中：

\[
R(u)\in \mathrm{End}(V\otimes V)
\]

即：

\[
R(u)\in \mathrm{Mat}(4,\mathbb C)
\]

---

# 3. 三体嵌入

在三体空间：

\[
V^{\otimes3}
\]

定义：

\[
\boxed{
h_{12}
=
\sum_{\alpha,\beta}
h_{\alpha\beta}
\sigma^\alpha\otimes\sigma^\beta\otimes I
}
\]

\[
\boxed{
h_{23}
=
\sum_{\mu,\nu}
h_{\mu\nu}
I\otimes\sigma^\mu\otimes\sigma^\nu
}
\]

于是：

\[
R_{12}=e^{-i\delta u h_{12}}
\]

\[
R_{23}=e^{-i\delta u h_{23}}
\]

---

# 4. Yang–Baxter结构

定义Yang–Baxter deviation：

\[
\boxed{
\Delta
=
R_{12}R_{23}R_{12}
-
R_{23}R_{12}R_{23}
}
\]

若：

\[
\Delta=0
\]

则满足Yang–Baxter关系。

---

# 5. 小步长展开

取：

\[
\delta u\ll1
\]

则：

\[
R=e^{-ih\delta u}
\]

展开：

\[
\boxed{
R
=
I
-i\delta u h
-\frac{\delta u^2}{2}h^2
+O(\delta u^3)
}
\]

因此：

\[
R_{12}=I-i\delta u h_{12}+O(\delta u^2)
\]

\[
R_{23}=I-i\delta u h_{23}+O(\delta u^2)
\]

---

# 6. YBE deviation展开

将展开代入：

\[
\Delta
=
R_{12}R_{23}R_{12}
-
R_{23}R_{12}R_{23}
\]

一阶项完全抵消。

二阶项也完全抵消。

最低非平凡阶出现在三阶：

\[
\boxed{
\Delta
=
-\delta u^3
[h_{12},h_{23}]
+
O(\delta u^4)
}
\]

因此：

\[
\boxed{
\text{non-Abelian性由局域Hamiltonian的不对易性决定}
}
\]

---

# 7. Pauli张量积下的具体计算

现在计算：

\[
[h_{12},h_{23}]
\]

代入：

\[
h_{12}
=
\sum_{\alpha,\beta}
h_{\alpha\beta}
\sigma^\alpha\otimes\sigma^\beta\otimes I
\]

\[
h_{23}
=
\sum_{\mu,\nu}
h_{\mu\nu}
I\otimes\sigma^\mu\otimes\sigma^\nu
\]

得到：

\[
[h_{12},h_{23}]
=
\sum
h_{\alpha\beta}
h_{\mu\nu}
[
\sigma^\alpha\otimes\sigma^\beta\otimes I,
I\otimes\sigma^\mu\otimes\sigma^\nu
]
\]

由于只有第二个site重叠：

\[
\boxed{
[h_{12},h_{23}]
=
\sum
h_{\alpha\beta}
h_{\mu\nu}
\sigma^\alpha
\otimes
[\sigma^\beta,\sigma^\mu]
\otimes
\sigma^\nu
}
\]

---

# 8. Pauli代数

使用：

\[
[\sigma^a,\sigma^b]
=
2i\epsilon_{abc}\sigma^c
\]

于是：

\[
\boxed{
[h_{12},h_{23}]
=
2i
\sum
h_{\alpha\beta}
h_{\mu\nu}
\epsilon_{\beta\mu\gamma}
\sigma^\alpha
\otimes
\sigma^\gamma
\otimes
\sigma^\nu
}
\]

因此：

\[
\boxed{
\Delta
\sim
2i\delta u^3
\sum
h_{\alpha\beta}
h_{\mu\nu}
\epsilon_{\beta\mu\gamma}
\sigma^\alpha
\otimes
\sigma^\gamma
\otimes
\sigma^\nu
}
\]

---

# 9. non-Abelian量化

定义Frobenius non-Abelian measure：

\[
\boxed{
\mathcal N
=
\sqrt{
\mathrm{Tr}(\Delta^\dagger\Delta)
}
}
\]

利用Pauli正交关系：

\[
\mathrm{Tr}(P_aP_b)=2^N\delta_{ab}
\]

得到：

\[
\boxed{
\mathcal N
\propto
\delta u^3
\sqrt{
\sum
|h_{\alpha\beta}|^2
|h_{\mu\nu}|^2
\epsilon_{\beta\mu\gamma}^2
}
}
\]

由于：

\[
\epsilon_{\beta\mu\gamma}^2=1
\]

仅当：

\[
\beta\neq\mu
\]

因此最终：

\[
\boxed{
\mathcal N
\sim
\delta u^3
\sqrt{
\sum_{\beta\neq\mu}
|h_{\alpha\beta}|^2
|h_{\mu\nu}|^2
}
}
\]

---

# 10. 物理意义

因此：

\[
\boxed{
\text{non-Abelian性来自不同Pauli通道之间的不对易耦合}
}
\]

---

# 11. 特殊情况

## 11.1 单Pauli通道

若Hamiltonian只有：

\[
\sigma^x\otimes\sigma^x
\]

则：

\[
[h_{12},h_{23}]=0
\]

因此：

\[
\boxed{
\mathcal N=0
}
\]

系统为Abelian。

---

## 11.2 多Pauli竞争

若同时存在：

\[
\sigma^x\otimes\sigma^x
\]

与：

\[
\sigma^y\otimes\sigma^y
\]

则：

\[
[\sigma^x,\sigma^y]\neq0
\]

因此：

\[
\boxed{
\mathcal N>0
}
\]

出现non-Abelian braid structure。

---

# 12. 与Schur complement联系

定义：

\[
H=
\begin{pmatrix}
H_{PP} & H_{PQ}\\
H_{QP} & H_{QQ}
\end{pmatrix}
\]

Green函数：

\[
G(z)=(z-H)^{-1}
\]

Schur complement：

\[
\boxed{
\Sigma(z)
=
H_{PQ}(z-H_{QQ})^{-1}H_{QP}
}
\]

于是有效Hamiltonian：

\[
\boxed{
H_{eff}
=
H_{PP}+\Sigma(z)
}
\]

因此：

\[
h_{\alpha\beta}
\rightarrow
h_{\alpha\beta}^{eff}
}
\]

从而：

\[
\boxed{
\mathcal N
\rightarrow
\mathcal N(\Sigma)
}
\]

即：

\[
\boxed{
ABS leakage
\Rightarrow
Pauli-channel mixing
\Rightarrow
non-Abelian deformation
}
\]

---

# 13. 最终统一结构

整个模型最终闭合为：

\[
\boxed{
h_{\alpha\beta}
\Rightarrow
[h_{12},h_{23}]
\Rightarrow
\Delta_{YB}
\Rightarrow
\mathcal N
}
\]

并且：

\[
\boxed{
\mathcal N
\sim
\sqrt{
\sum_{\beta\neq\mu}
|h_{\alpha\beta}|^2
|h_{\mu\nu}|^2
}
}
\]

---

# 14. 最终物理图像

整个理论的核心为：

\[
\boxed{
\text{局域Pauli张量积Hamiltonian}
\Rightarrow
\text{路径演化}
\Rightarrow
\text{braid algebra}
\Rightarrow
\text{non-Abelian统计}
}
\]

而：

\[
\boxed{
ABS
=
\Sigma(z)
\text{引入的额外Pauli mixing}
}
\]

最终导致：

\[
\boxed{
\text{braid deformation}
}
\]

即：

\[
\boxed{
\text{spectral topology}
\leftrightarrow
\text{operator algebra geometry}
}
\]

