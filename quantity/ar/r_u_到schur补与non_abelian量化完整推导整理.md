# R(u) 到 non-Abelian 量化：PRP、Schur补、YBE与实验观测统一推导

---

# 0. 基本对象与目标

考虑算符：

\[
R(u)=U(u)e^{-iK(u)}
\]

Hilbert空间分解：

\[
\mathcal H = \mathcal H_P \oplus \mathcal H_Q
\]

定义投影算符：

\[
P^2=P,\quad Q=I-P
\]

研究目标：

- PRP是否封闭
- non-Abelian结构如何产生
- ABS / MZM如何区分
- 与PRB 111.205411中non-Abelian measure的关系

---

# 1. PRP投影结构

定义有效演化：

\[
R_P(u)=PR(u)P
\]

封闭性条件：

\[
QR(u)P=0
\]

---

# 2. 指数展开（严格）

\[
e^{-iK} = \sum_{n=0}^{\infty} \frac{(-i)^n}{n!} K^n
\]

因此：

\[
R = U \sum_{n\ge0} \frac{(-i)^n}{n!} K^n
\]

投影：

\[
QRP = \sum_{n\ge0} \frac{(-i)^n}{n!} Q U K^n P
\]

封闭性等价于：

\[
QUK^nP=0,\quad \forall n
\]

---

# 3. 一阶展开（核心结构）

插入单位：

\[
I=P+Q
\]

对K展开：

\[
K = PKP + PKQ + QKP + QKQ
\]

---

计算：

\[
QUKP = QU(PKP + PKQ + QKP + QKQ)P
\]

展开：

\[
QUKP = QU(PKP) + QU(PKQ)P + QU(QKP)P + QU(QKQ)P
\]

---

# 4. 结论1：非封闭来源

PRP不封闭 ⇔ 存在：

\[
PKQ \neq 0 \quad \text{或} \quad QKP \neq 0
\]

即K在P/Q分解下存在非对角块。

---

# 5. Schur complement严格推导

考虑 resolvent：

\[
G(z) = (z-K)^{-1}
\]

K块结构：

\[
K =
\begin{pmatrix}
K_{PP} & K_{PQ} \\
K_{QP} & K_{QQ}
\end{pmatrix}
\]

---

# 5.1 Schur补公式（严格恒等式）

\[
G_{PP}(z)
=
(z - K_{PP} - K_{PQ}(z-K_{QQ})^{-1}K_{QP})^{-1}
\]

---

# 5.2 定义自能（Schur补）

\[
\Sigma(z)=K_{PQ}(z-K_{QQ})^{-1}K_{QP}
\]

有效哈密顿量：

\[
K_{\mathrm{eff}}(z)=K_{PP}+\Sigma(z)
\]

---

# 6. Schur补物理意义

\[
\Sigma = P \to Q \to P \text{虚跃迁回路}
\]

展开：

\[
\Sigma(z)=\sum_n \frac{K_{PQ}|n\rangle\langle n|K_{QP}}{z-\epsilon_n}
\]

---

# 7. PRP封闭条件

\[
PRP \text{封闭} \iff \Sigma(z)=0
\]

等价于：

\[
K_{PQ}=0,\quad K_{QP}=0
\]

---

# 8. MZM / ABS分类

## 8.1 MZM

\[
\Sigma=0
\]

- 无混合
- zero mode严格存在
- braid exact

---

## 8.2 ABS

\[
\Sigma\neq0
\]

- bulk混合
- 能级移动
- braid变形

---

# 9. 有效演化

\[
R_{\text{eff}}(u)=e^{-i(K_{PP}+\Sigma)}
\]

---

# 10. Connection定义

\[
A(u)=R^{-1}_{\text{eff}}\partial_u R_{\text{eff}}
\]

展开：

\[
A = A_0 + A_\Sigma
\]

其中：

\[
A_0=R_{PP}^{-1}\partial_u R_{PP}
\]

\[
A_\Sigma = \text{Schur修正项}
\]

---

# 11. 曲率（non-Abelian来源）

\[
F = \partial_u A + A^2
\]

展开：

\[
F = F_0 + F_\Sigma + [A_0,A_\Sigma]
\]

---

# 12. 几何判据

## MZM：

\[
\Sigma=0 \Rightarrow F=0
\]

flat connection

## ABS：

\[
\Sigma\neq0 \Rightarrow F\neq0
\]

non-flat connection

---

# 13. braid algebra deformation

定义：

\[
\Delta_{braid}=R_iR_j-R_jR_i
\]

展开：

\[
\Delta_{braid}=\Delta_0 + [R_0,\Sigma] + [\Sigma,\Sigma]
\]

---

# 14. YBE变形

\[
R_{12}R_{13}R_{23}-R_{23}R_{13}R_{12}
\sim \mathcal D(\Sigma)
\]

---

# 15. Drinfeld twist等价

存在F：

\[
R_{\text{eff}} = F_{21}R_0F^{-1}
\]

小变换：

\[
\Sigma \leftrightarrow \mathcal D_{R_0}(f)
\]

结论：

\[
\boxed{\Sigma \equiv \text{Drinfeld twist generator}}
\]

---

# 16. 量子群变形

\[
U_q(\mathfrak g) \to U_{q+\delta q}
\]

\[
\delta q \sim \Sigma
\]

---

# 17. non-Abelian统一量

\[
\boxed{
\mathcal N
=\|\Sigma\|
+\|\Delta_{braid}\|
+\|F\|
}
\]

---

# 18. 实验可观测量

## 18.1 能级偏移

\[
E_{ABS} \approx \mathrm{Re}\,\Sigma(0)
\]

---

## 18.2 寿命（泄漏）

\[
\Gamma = -\mathrm{Im}\,\Sigma
\]

---

## 18.3 braid phase error

\[
\Delta\theta \sim f(\Sigma)
\]

---

## 18.4 YBE违背

\[
\Delta_{YBE} \sim \Sigma + \Sigma^2
\]

---

# 19. PRB 111.205411对应关系

| PRB | 本模型 |
|-----|--------|
| ABS | P子空间 |
| braid deviation | Δ_braid |
| hybridization | K_PQ,K_QP |
| non-Abelian measure | Σ-induced deformation |

---

# 20. 最终统一结论

\[
\boxed{
\text{non-Abelianity} = f(\Sigma)
}
\]

且

\[
\boxed{
\Sigma = P\to Q\to P \text{虚跃迁回路}
}
\]

---

# 21. 总结一句话

PRP投影不是截断，而是：

> full system在P子空间上的Schur补有效理论，其非阿贝尔性完全由Σ决定。

