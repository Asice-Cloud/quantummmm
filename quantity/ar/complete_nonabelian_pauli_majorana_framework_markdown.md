# 从 Pauli 张量积到 Majorana / Schur 补 / non-Abelian 量化的统一理论

---

# 1. 微观自旋模型

考虑开放边界 XY 链：

\[
H
=
\sum_{j=1}^{N-1}
\left(
J_x\sigma_j^x\sigma_{j+1}^x
+
J_y\sigma_j^y\sigma_{j+1}^y
\right)
+
h\sum_{j=1}^{N}\sigma_j^z
\]

Hilbert 空间：

\[
\mathcal H
=
\bigotimes_{j=1}^N \mathbb C^2
\]

---

# 2. Jordan–Wigner 变换

定义费米子：

\[
c_j
=
\left(
\prod_{k<j}\sigma_k^z
\right)
\frac{\sigma_j^x-i\sigma_j^y}{2}
\]

\[
c_j^\dagger
=
\left(
\prod_{k<j}\sigma_k^z
\right)
\frac{\sigma_j^x+i\sigma_j^y}{2}
\]

满足：

\[
\{c_i,c_j^\dagger\}=\delta_{ij}
\]

定义：

\[
t=J_x+J_y
\]

\[
\Delta=J_x-J_y
\]

\[
\mu=2h
\]

得到 Kitaev 链：

\[
H=
-t\sum_j(c_j^\dagger c_{j+1}+h.c.)
+
\Delta\sum_j(c_j c_{j+1}+h.c.)
-
\mu\sum_j c_j^\dagger c_j
\]

---

# 3. Majorana 表示

定义：

\[
\gamma_{A,j}=c_j+c_j^\dagger
\]

\[
\gamma_{B,j}=-i(c_j-c_j^\dagger)
\]

满足：

\[
\gamma^\dagger=\gamma
\]

\[
\{\gamma_i,\gamma_j\}=2\delta_{ij}
\]

哈密顿量：

\[
H
=
\frac{i}{2}
\sum_j
\Big[
-\mu\gamma_{A,j}\gamma_{B,j}
\]

\[
+(t+\Delta)\gamma_{B,j}\gamma_{A,j+1}
\]

\[
+(-t+\Delta)\gamma_{A,j}\gamma_{B,j+1}
\Big]
\]

---

# 4. 拓扑点与零模

取：

\[
t=\Delta,
\quad
\mu=0
\]

则：

\[
H_0
=
it\sum_{j=1}^{N-1}
\gamma_{B,j}\gamma_{A,j+1}
\]

边界 Majorana：

\[
\gamma_L=\gamma_{A,1}
\]

\[
\gamma_R=\gamma_{B,N}
\]

完全不进入哈密顿量。

定义复费米子：

\[
f=\frac{\gamma_L+i\gamma_R}{2}
\]

零模空间：

\[
\mathcal H_P
=
\mathrm{span}
\{
|0\rangle,
 f^\dagger|0\rangle
\}
\]

bulk 模空间：

\[
\mathcal H_Q
\]

整体：

\[
\mathcal H
=
\mathcal H_P
\oplus
\mathcal H_Q
\]

---

# 5. 偏离拓扑点与 ABS

定义：

\[
\delta=t-\Delta
\]

扰动：

\[
V
=
i\delta
\sum_j
\gamma_{A,j}\gamma_{B,j+1}
\]

总哈密顿量：

\[
H=H_0+V
\]

有限长度与 pairing imbalance 会导致：

\[
H_{PQ}\neq0
\]

出现 edge-bulk mixing。

---

# 6. P/Q 分块

写成：

\[
H
=
\begin{pmatrix}
H_{PP} & H_{PQ}\\
H_{QP} & H_{QQ}
\end{pmatrix}
\]

其中：

- \(H_{PP}\)：边界零模
- \(H_{QQ}\)：bulk Bogoliubov 模
- \(H_{PQ},H_{QP}\)：边界-bulk 耦合

---

# 7. Schur 补

定义 Green 函数：

\[
G(z)=(z-H)^{-1}
\]

块矩阵恒等式：

\[
G_{PP}(z)
=
(z-H_{PP}-\Sigma(z))^{-1}
\]

定义：

\[
\Sigma(z)
=

H_{PQ}(z-H_{QQ})^{-1}H_{QP}
\]

有效哈密顿量：

\[
\boxed{
H_{eff}(z)
=
H_{PP}+\Sigma(z)
}
\]

---

# 8. Schur 补展开

设：

\[
H_{QQ}|n\rangle=\epsilon_n|n\rangle
\]

则：

\[
\boxed{
\Sigma(z)
=
\sum_n
\frac{
H_{PQ}|n\rangle\langle n|H_{QP}
}{
z-\epsilon_n
}
}
\]

表示：

\[
P\to Q\to P
\]

的虚跃迁过程。

---

# 9. 最小链模型（N=3）

Majorana：

\[
(\gamma_{A1},\gamma_{B1},
\gamma_{A2},\gamma_{B2},
\gamma_{A3},\gamma_{B3})
\]

理想部分：

\[
H_0
=
it(
\gamma_{B1}\gamma_{A2}
+
\gamma_{B2}\gamma_{A3}
)
\]

扰动：

\[
V
=
i\delta(
\gamma_{A1}\gamma_{B2}
+
\gamma_{A2}\gamma_{B3}
)
\]

---

# 10. P 空间

边界零模：

\[
\gamma_L=\gamma_{A1}
\]

\[
\gamma_R=\gamma_{B3}
\]

定义：

\[
f=\frac{\gamma_{A1}+i\gamma_{B3}}2
\]

P 空间：

\[
\{
|0\rangle,
 f^\dagger|0\rangle
\}
\]

---

# 11. Q 空间

定义 bulk fermion：

\[
d_1=\frac{\gamma_{B1}+i\gamma_{A2}}2
\]

\[
d_2=\frac{\gamma_{B2}+i\gamma_{A3}}2
\]

Q 空间：

\[
\{
d_1^\dagger|0\rangle,
 d_2^\dagger|0\rangle
\}
\]

---

# 12. Hamiltonian 矩阵

基底：

\[
(
|0\rangle,
 f^\dagger|0\rangle,
 d_1^\dagger|0\rangle,
 d_2^\dagger|0\rangle
)
\]

理想部分：

\[
H_0
=
\begin{pmatrix}
0&0&0&0\\
0&0&0&0\\
0&0&2t&0\\
0&0&0&2t
\end{pmatrix}
\]

耦合：

\[
H_{PQ}
=
\delta
\begin{pmatrix}
1&1\\
1&1
\end{pmatrix}
\]

\[
H_{QP}=H_{PQ}^\dagger
\]

bulk：

\[
H_{QQ}
=
\begin{pmatrix}
2t&0\\
0&2t
\end{pmatrix}
\]

---

# 13. Schur 补具体结果

\[
\Sigma(z)
=
H_{PQ}(z-H_{QQ})^{-1}H_{QP}
\]

\[
(z-H_{QQ})^{-1}
=
\begin{pmatrix}
\frac1{z-2t}&0\\
0&\frac1{z-2t}
\end{pmatrix}
\]

得到：

\[
\boxed{
\Sigma(z)
=
\frac{2\delta^2}{z-2t}
\begin{pmatrix}
1&1\\
1&1
\end{pmatrix}
}
\]

有效哈密顿量：

\[
\boxed{
H_{eff}(z)
=
\frac{2\delta^2}{z-2t}
\begin{pmatrix}
1&1\\
1&1
\end{pmatrix}
}
\]

---

# 14. ABS 能级

特征值：

\[
\lambda_1=0
\]

\[
\lambda_2
=
\frac{4\delta^2}{z-2t}
\]

低能近似：

\[
z\approx0
\]

得到：

\[
\boxed{
E_{ABS}
\approx
-\frac{2\delta^2}{t}
}
\]

---

# 15. braid operator

定义：

\[
R(u)=e^{-iuH}
\]

投影：

\[
R_{eff}(u)
=
Pe^{-iuH}P
=
e^{-iuH_{eff}}
\]

展开：

\[
R_{eff}
=
I-iuH_{eff}+O(H_{eff}^2)
\]

---

# 16. braid deformation

定义：

\[
\Delta_{braid}
=
R_1R_2R_1
-
R_2R_1R_2
\]

设：

\[
R_i=R_0+\delta R_i
\]

其中：

\[
\delta R_i=-iu\Sigma_iR_0
\]

展开到一阶：

\[
\boxed{
\Delta_{braid}
\sim
[R_0,\Sigma]
}
\]

---

# 17. non-Abelian measure

定义：

\[
\boxed{
\mathcal N
=
\|\Delta_{braid}\|
}
\]

因此：

\[
\boxed{
\mathcal N
\sim
\|\Sigma\|
}
\]

代入：

\[
\boxed{
\mathcal N
\sim
\frac{(t-\Delta)^2}{t}
}
\]

---

# 18. MZM 与 ABS

## MZM

\[
t=\Delta
\]

则：

\[
\Sigma=0
\]

\[
\Delta_{braid}=0
\]

存在严格 non-Abelian braid representation。

---

## ABS

\[
t\neq\Delta
\]

则：

\[
\Sigma\neq0
\]

出现：

- edge-bulk mixing
- finite-energy ABS
- braid deformation
- YBE defect

---

# 19. 与 PRB 111.205411 对应

PRB 定义：

\[
\mathcal N_{PRB}
=
\|U_iU_jU_i-U_jU_iU_j\|
\]

本模型：

\[
\boxed{
\mathcal N_{PRB}
\sim
\|\Sigma\|
\sim
\frac{(t-\Delta)^2}{t}
}
\]

对应关系：

| PRB | 本模型 |
|---|---|
| ABS subspace | \(\mathcal H_P\) |
| hybridization | \(H_{PQ},H_{QP}\) |
| braid deviation | \(\Delta_{braid}\) |
| non-Abelian measure | \(\Sigma\) |

---

# 20. Pauli 张量积空间

定义 Pauli operator basis：

\[
\mathcal P_N
=
\{
\sigma^{\alpha_1}\otimes\cdots\otimes\sigma^{\alpha_N}
\}
\]

其中：

\[
\alpha_i\in\{0,x,y,z\}
\]

维度：

\[
\dim\mathcal P_N=4^N
\]

Pauli algebra：

\[
[\sigma^x,\sigma^y]=2i\sigma^z
\]

因此：

\[
\boxed{
\mathcal P_N
\text{ 本身是 non-Abelian operator algebra}
}
\]

---

# 21. Pauli braid generator

定义：

\[
\boxed{
B_j
=
\exp\left(
 i\theta\,
 \sigma_j^\alpha\sigma_{j+1}^\beta
\right)
}
\]

例如：

\[
B_j
=
 e^{i\frac\pi4\sigma_j^x\sigma_{j+1}^x}
\]

由于：

\[
(\sigma_i^\alpha\sigma_j^\beta)^2=I
\]

因此：

\[
B_j
=
\cos\theta
+
i\sin\theta\,
\sigma_j^\alpha\sigma_{j+1}^\beta
\]

特别：

\[
\theta=\frac\pi4
\]

得到：

\[
\boxed{
B_j
=
\frac1{\sqrt2}
(I+i\sigma_j^\alpha\sigma_{j+1}^\beta)
}
\]

---

# 22. braid relation

定义：

\[
A_j=\sigma_j^x\sigma_{j+1}^x
\]

满足：

\[
A_j^2=I
\]

以及：

\[
A_jA_{j+1}=-A_{j+1}A_j
\]

则：

\[
\boxed{
B_jB_{j+1}B_j
=
B_{j+1}B_jB_{j+1}
}
\]

即 braid group representation。

---

# 23. Majorana 与 Pauli 的对应

定义：

\[
\gamma_{2j-1}
=
\left(
\prod_{k<j}\sigma_k^z
\right)\sigma_j^x
\]

\[
\gamma_{2j}
=
\left(
\prod_{k<j}\sigma_k^z
\right)\sigma_j^y
\]

则：

\[
i\gamma_j\gamma_{j+1}
\Longleftrightarrow
\sigma^\alpha\sigma^\beta
\]

Majorana braid：

\[
U_j
=
 e^{\frac\pi4\gamma_j\gamma_{j+1}}
\]

对应：

\[
\boxed{
U_j
=
 e^{i\frac\pi4\sigma_j^\alpha\sigma_{j+1}^\beta}
}
\]

---

# 24. Pauli 空间中的 non-Abelian 量化

## (1) braid defect

\[
\boxed{
\mathcal N_1
=
\|B_iB_{i+1}B_i-B_{i+1}B_iB_{i+1}\|
}
\]

---

## (2) commutator curvature

\[
\boxed{
\mathcal N_2
=
\sum_{i<j}
\|[B_i,B_j]\|
}
\]

---

## (3) Berry curvature

定义 connection：

\[
A_\mu
=
B^{-1}\partial_\mu B
\]

曲率：

\[
F_{\mu\nu}
=
\partial_\mu A_\nu
-
\partial_\nu A_\mu
+
[A_\mu,A_\nu]
\]

量化：

\[
\boxed{
\mathcal N_3
=
\int
\mathrm{Tr}(F\wedge *F)
}
\]

---

## (4) operator entanglement

Schmidt 分解：

\[
B_j
=
\sum_\alpha
 s_\alpha
 A_\alpha\otimes C_\alpha
\]

定义：

\[
\boxed{
\mathcal N_4
=
-\sum_\alpha
 s_\alpha^2
 \log s_\alpha^2
}
\]

---

# 25. Schur 补与 Pauli algebra deformation

有效理论：

\[
R_{eff}
=
e^{-i(H_{PP}+\Sigma)}
\]

展开：

\[
H_{PP}
=
\sum c_{\alpha\beta}
\sigma^\alpha\otimes\sigma^\beta
\]

\[
\Sigma
=
\sum d_{\alpha\beta}
\sigma^\alpha\otimes\sigma^\beta
\]

因此：

\[
\boxed{
\Sigma
\text{ 本身是 Pauli algebra deformation}
}
\]

---

# 26. 最终统一结构

\[
\sigma_i^\alpha\sigma_j^\beta
\]

↓

Jordan–Wigner

↓

Kitaev chain

↓

Majorana zero modes

↓

P/Q decomposition

↓

\[
\Sigma
=
H_{PQ}(z-H_{QQ})^{-1}H_{QP}
\]

↓

\[
R_{eff}
=
e^{-i(H_{PP}+\Sigma)}
\]

↓

\[
\Delta_{braid}
\sim
[R_0,\Sigma]
\]

↓

\[
\mathcal N
\sim
\|\Sigma\|
\]

同时：

\[
B_j
=
 e^{i\theta\sigma^\alpha\sigma^\beta}
\]

直接在 Pauli operator space 上形成 braid representation。

---

# 27. 最终结论

\[
\boxed{
\text{non-Abelianity degradation}
=
\text{Majorana-bulk virtual mixing induced by Pauli tensor interactions}
}
\]

并且：

\[
\boxed{
\text{Pauli tensor-product operator space itself carries a non-Abelian braid representation}
}
\]

因此：

\[
\boxed{
\text{non-Abelian statistics can be quantified both:}
}
\]

1. 在 Majorana / Schur complement 有效理论中；
2. 在 Pauli operator algebra 中直接定义。

