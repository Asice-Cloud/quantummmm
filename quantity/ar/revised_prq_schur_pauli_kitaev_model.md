# 修正后的完整 PRQ–Schur–Pauli–Kitaev 非阿贝尔传播模型

---

# 0. 理论目标

本文给出修正后的完整理论结构。

核心修正：

原模型中：

\[
h_{\alpha\beta}(u)\sigma^\alpha\otimes\sigma^\beta
\]

被视为“原初 Hamiltonian”。

但在新的 PRQ 传播框架中：

\[
h_{\alpha\beta}^{\mathrm{eff}}
\]

并不是输入参数，而是由：

\[
P_iRQ_iRP_j
\]

诱导出的 effective propagation kernel。

因此：

- Pauli tensor 不再是理论本体；
- propagation structure 才是本体；
- non-Abelian 性来自传播路径的不交换；
- Kitaev chain 只是该传播理论的一个特殊实例。

---

# 1. 理论本体：投影传播结构

定义 Hilbert 空间：

\[
\mathcal H = \mathcal H_P \oplus \mathcal H_Q
\]

定义投影：

\[
P^2=P,\qquad Q^2=Q,\qquad P+Q=I
\]

传播算符：

\[
R
\]

可以是：

## (a) 时间演化

\[
R(t)=e^{-itH}
\]

## (b) Green 函数

\[
R(\omega)=(\omega-H)^{-1}
\]

## (c) 更一般 transfer/scattering operator

---

理论真正研究对象：

\[
P_iRQ_j
\]

以及：

\[
P_iRQ_jRP_k
\]

即：

“子空间之间的传播”。

---

# 2. Schur 补的严格出现

Hamiltonian 分块：

\[
H=
\begin{pmatrix}
H_{PP}&H_{PQ}\\
H_{QP}&H_{QQ}
\end{pmatrix}
\]

定义 resolvent：

\[
G(\omega)=(\omega-H)^{-1}
\]

则：

\[
G_{PP}
=
P(\omega-H)^{-1}P
\]

利用 block inversion：

\[
\boxed{
G_{PP}
=
(\omega-H_{PP}-\Sigma(\omega))^{-1}
}
\]

其中：

\[
\boxed{
\Sigma(\omega)
=
H_{PQ}
(\omega-H_{QQ})^{-1}
H_{QP}
}
\]

这是 exact identity。

---

# 3. PRQ 有效传播核

定义：

\[
\boxed{
K_{ij}(\omega)
=
P_iRQ_jRP_j
}
\]

在 resolvent 表示中：

\[
\boxed{
K_{ij}(\omega)
=
H_{P_iQ_j}
(\omega-H_Q)^{-1}
H_{Q_jP_i}
}
\]

因此：

effective propagation：

\[
P_i
\to
Q_j
\to
P_k
\]

由 Schur kernel 精确生成。

---

# 4. Kitaev 链：严格对应

标准 Kitaev Hamiltonian：

\[
H
=
-\mu\sum_i c_i^\dagger c_i
-
t\sum_i(c_i^\dagger c_{i+1}+h.c.)
+
\Delta\sum_i(c_ic_{i+1}+h.c.)
\]

---

## 4.1 Majorana 分解

定义：

\[
\gamma_{iA}=c_i+c_i^\dagger
\]

\[
\gamma_{iB}=-i(c_i-c_i^\dagger)
\]

满足：

\[
\{\gamma_{i\alpha},\gamma_{j\beta}\}
=
2\delta_{ij}\delta_{\alpha\beta}
\]

---

Hamiltonian：

\[
H
=
\frac{i}{2}\sum_i
\Big[
-\mu\gamma_{iA}\gamma_{iB}
+
(t+\Delta)\gamma_{iB}\gamma_{i+1,A}
+
(-t+\Delta)\gamma_{iA}\gamma_{i+1,B}
\Big]
\]

---

## 4.2 Kitaev 拓扑点

取：

\[
t=\Delta
\]

则：

\[
\boxed{
H
=
\frac{i}{2}
\sum_i
\Big[
-\mu\gamma_{iA}\gamma_{iB}
+
2t\gamma_{iB}\gamma_{i+1,A}
\Big]
}
\]

---

# 5. Pi 与 Qi 的严格定义

定义：

\[
\boxed{
P_i=\mathrm{span}\{\gamma_{iB}\}
}
\]

\[
\boxed{
Q_i=\mathrm{span}\{\gamma_{iA}\}
}
\]

因此：

每个 site：

\[
Q_i=\gamma_{iA},
\qquad
P_i=\gamma_{iB}
\]

---

# 6. Kitaev 链中的真实传播

Hamiltonian 中：

## onsite term

\[
-\frac{i\mu}{2}\gamma_{iA}\gamma_{iB}
\]

对应：

\[
P_i\leftrightarrow Q_i
\]

---

## nearest-neighbor term

\[
it\gamma_{iB}\gamma_{i+1,A}
\]

对应：

\[
P_i\leftrightarrow Q_{i+1}
\]

---

因此真实传播：

\[
\boxed{
P_i
\to
Q_{i+1}
\to
P_{i+1}
}
\]

或者等价：

\[
Q_i
\to
P_i
\to
Q_{i+1}
\]

---

# 7. Schur 传播核

消去 Q-sector：

\[
Q_i=\gamma_{iA}
\]

定义：

\[
G(\omega)=(\omega-H)^{-1}
\]

Schur kernel：

\[
\boxed{
T_{i\to i+1}(\omega)
=
H_{P_iQ_{i+1}}
(\omega-H_Q)^{-1}
H_{Q_{i+1}P_{i+1}}
}
\]

---

在 Kitaev 点：

\[
H_{QQ}=0
\]

因此：

\[
(\omega-H_Q)^{-1}
=
\frac1\omega
\]

---

于是：

\[
H_{P_iQ_{i+1}}=it
\]

\[
H_{Q_{i+1}P_{i+1}}
=
-\frac{i\mu}{2}
\]

得到：

\[
\boxed{
T_{i\to i+1}(\omega)
=
\frac{\mu t}{2\omega}
}
\]

---

# 8. 边界 Majorana 的严格推导

零模：

\[
H\psi=0
\]

写：

\[
\psi
=
\sum_i
(a_i\gamma_{iA}+b_i\gamma_{iB})
\]

代入 Hamiltonian。

得到：

\[
-\mu a_i +2t a_{i+1}=0
\]

因此：

\[
\boxed{
a_{i+1}
=
\frac{\mu}{2t}a_i
}
\]

递推：

\[
a_i
=
\left(\frac{\mu}{2t}\right)^{i-1}a_1
\]

若：

\[
\left|\frac{\mu}{2t}\right|<1
\]

则：

\[
a_i\sim e^{-i/\xi}
\]

指数衰减。

因此：

\[
\boxed{
\psi_L
=
\sum_i
\left(\frac{\mu}{2t}\right)^{i-1}
\gamma_{iA}
}
\]

为左边界 Majorana。

---

# 9. 传播与 Majorana 的严格关系

定义 transfer operator：

\[
\mathcal T:
\psi_i\mapsto\psi_{i+1}
\]

则：

\[
\boxed{
\psi_{i+1}
=
\frac{\mu}{2t}\psi_i
}
\]

即：

\[
\boxed{
\psi_i
=
\mathcal T^{i-1}\psi_1
}
\]

若：

\[
\rho(\mathcal T)<1
\]

则 bulk propagation 收缩。

因此：

边界 Majorana：

\[
\boxed{
\text{是 bulk propagation 的可正规化固定点}
}
\]

---

# 10. Jordan–Wigner 与 Pauli 表示

JW 变换：

\[
\gamma_{iA}
=
\left(\prod_{j<i}\sigma_j^z\right)\sigma_i^x
\]

\[
\gamma_{iB}
=
\left(\prod_{j<i}\sigma_j^z\right)\sigma_i^y
\]

---

因此：

\[
\gamma_{iB}\gamma_{i+1,A}
\sim
\sigma_i^y\sigma_{i+1}^x
\]

---

于是：

\[
\boxed{
h_{i,i+1}^{\mathrm{eff}}
=
K_{i,i+1}(\omega)
\sigma_i^y\otimes\sigma_{i+1}^x
}
\]

其中：

\[
K_{i,i+1}
=
T_{i\to i+1}
\]

---

# 11. 关键修正：Pauli Hamiltonian 不再 fundamental

旧模型：

\[
h_{\alpha\beta}(u)
\]

被视为输入。

---

新模型：

\[
\boxed{
h_{\alpha\beta}^{\mathrm{eff}}
=
\langle\alpha_i|
H_{PQ}
(\omega-H_Q)^{-1}
H_{QP}
|\beta_j\rangle
}
\]

---

即：

Pauli tensor coupling：

\[
\boxed{
\text{由 propagation 自动生成}
}
\]

---

# 12. 修正后的 Dyson 展开

定义 effective Hamiltonian：

\[
\boxed{
h^{\mathrm{eff}}(u,\omega)
=
\sum_{\alpha,\beta}
h_{\alpha\beta}^{\mathrm{eff}}(u,\omega)
\sigma^\alpha\otimes\sigma^\beta
}
\]

---

定义路径演化：

\[
\boxed{
R(u)
=
\mathcal T
\exp
\left(
-i\int_0^u h^{\mathrm{eff}}(s,\omega)ds
\right)
}
\]

---

Dyson 展开：

\[
R(u)=I+R^{(1)}+R^{(2)}+\cdots
\]

---

## 一阶

\[
R^{(1)}
=
-i\int_0^u ds\,h^{\mathrm{eff}}(s)
\]

---

## 二阶

\[
R^{(2)}
=
(-i)^2
\int_0^u ds_1
\int_0^{s_1}ds_2
h^{\mathrm{eff}}(s_1)
h^{\mathrm{eff}}(s_2)
\]

---

# 13. 三体嵌入

定义：

\[
V^{\otimes3}
\]

---

定义：

\[
h_{12}^{\mathrm{eff}}
=
\sum_{\alpha,\beta}
h_{\alpha\beta}^{\mathrm{eff}}
\sigma^\alpha\otimes\sigma^\beta\otimes I
\]

---

以及：

\[
h_{23}^{\mathrm{eff}}
=
\sum_{\mu,\nu}
h_{\mu\nu}^{\mathrm{eff}}
I\otimes\sigma^\mu\otimes\sigma^\nu
\]

---

# 14. Yang–Baxter deviation

定义：

\[
\boxed{
\Delta
=
R_{12}(u)
R_{23}(u+v)
R_{12}(v)
-
R_{23}(v)
R_{12}(u+v)
R_{23}(u)
}
\]

---

# 15. 最低非平凡阶

严格 Dyson 展开后：

\[
\boxed{
\Delta^{(3)}
\sim
[h_{12}^{\mathrm{eff}},h_{23}^{\mathrm{eff}}]
}
\]

---

代入 effective propagation kernel：

\[
\boxed{
\Delta^{(3)}
\sim
[
P_1RQ_1RP_2,\;
P_2RQ_2RP_3
]
}
\]

---

# 16. non-Abelian 来源

旧模型：

\[
[\sigma^\alpha,\sigma^\beta]\neq0
\]

---

新模型：

non-Abelian 性有两层来源：

---

## (A) propagation noncommutativity

\[
[K(s_1),K(s_2)]\neq0
\]

---

## (B) Pauli algebra noncommutativity

\[
[\sigma^\alpha,\sigma^\beta]\neq0
\]

---

因此：

\[
\boxed{
\text{non-Abelian性}
=
\text{传播路径不交换}
+
\text{Pauli代数不交换}
}
\]

---

# 17. Pauli 交换子

使用：

\[
[\sigma^a,\sigma^b]
=
2i\epsilon_{abc}\sigma^c
\]

因此：

\[
\boxed{
[h_{12}^{\mathrm{eff}},h_{23}^{\mathrm{eff}}]
=
2i
\sum
h_{\alpha\beta}^{\mathrm{eff}}
h_{\mu\nu}^{\mathrm{eff}}
\epsilon_{\beta\mu\gamma}
\sigma^\alpha\otimes\sigma^\gamma\otimes\sigma^\nu
}
\]

---

# 18. Frobenius non-Abelian measure

定义：

\[
\boxed{
\mathcal N
=
\sqrt{
\mathrm{Tr}(\Delta^\dagger\Delta)
}
}
\]

---

利用：

\[
\mathrm{Tr}(P_aP_b)=2^N\delta_{ab}
\]

得到：

\[
\boxed{
\mathcal N
\sim
\sqrt{
\sum
|h_{\alpha\beta}^{\mathrm{eff}}|^2
|h_{\mu\nu}^{\mathrm{eff}}|^2
\epsilon_{\beta\mu\gamma}^2
}
}
\]

---

# 19. 修正后的完整理论闭环

最终：

\[
\boxed{
P_iRQ_jRP_k
\Rightarrow
\Sigma(\omega)
\Rightarrow
h_{\alpha\beta}^{\mathrm{eff}}
\Rightarrow
R(u)
\Rightarrow
[h_{12}^{\mathrm{eff}},h_{23}^{\mathrm{eff}}]
\Rightarrow
\Delta_{YB}
\Rightarrow
\mathcal N
}
\]

---

# 20. 理论最终物理图像

整个理论本体：

不是 bare Hamiltonian。

而是：

\[
\boxed{
\mathcal N_{\mathrm{prop}}
=
\{P_iRQ_jRP_k\}
}
\]

形成的 propagation network。

其中：

- \(P_i\)：边界/低能 fiber；
- \(Q_j\)：bulk/junction continuum；
- \(R\)：传播；
- Schur 补：消去中间传播；
- Pauli tensor：传播的 operator representation；
- Majorana：传播固定点；
- topology：传播网络无法平凡闭合。

---

# 21. 最终核心结论

修正后的理论不再是：

\[
\text{Pauli path algebra}
\]

而是：

\[
\boxed{
\text{Schur-generated propagation geometry}
}
\]

Kitaev 链只是该传播几何中的一个特殊可解实例。

真正 fundamental 的对象是：

\[
\boxed{
P_iRQ_jRP_k
}
\]

而：

\[
h_{\alpha\beta}^{\mathrm{eff}}
\]

只是其 Pauli representation。
