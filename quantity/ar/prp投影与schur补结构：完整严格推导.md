# PRP投影、动力学泄漏与Schur补结构（严格推导整理）

---

# 0. 基本对象

考虑算符

\[
R(u)=U(u)e^{-iK(u)}
\]

定义投影算符：

\[
P^2=P,\quad Q=I-P
\]

目标研究：

\[
R_P(u):=PR(u)P
\]

以及封闭性条件：

\[
QR(u)P=0
\]

---

# 1. PRP封闭性的严格展开

## 1.1 指数展开

\[
e^{-iK} = \sum_{n=0}^{\infty} \frac{(-i)^n}{n!} K^n
\]

代入 R：

\[
R = U \sum_{n=0}^{\infty} \frac{(-i)^n}{n!} K^n
\]

左乘 Q，右乘 P：

\[
QRP = \sum_{n=0}^{\infty} \frac{(-i)^n}{n!} Q U K^n P
\]

---

## 1.2 封闭性等价条件

\[
QRP=0 \iff QUK^nP=0,\ \forall n \ge 0
\]

---

# 2. 一阶结构完整展开

考虑 n=1：

\[
QUKP
\]

插入单位：

\[
I=P+Q
\]


## 2.1 对 K 的完全展开

\[
K = (P+Q)K(P+Q)
\]

展开得：

\[
K = PKP + PKQ + QKP + QKQ
\]

---

## 2.2 代入计算

\[
QUKP = QU(PKP + PKQ + QKP + QKQ)P
\]

逐项展开：

\[
QUKP = QU(PKP) + QU(PKQ)P + QU(QKP)P + QU(QKQ)P
\]

---

# 3. 封闭性条件的真实结构

即使满足纯 braid 封闭：

\[
QUP = 0
\]

仍有：

\[
QUKP = QU(PKQ)P + QU(QKP)P + QU(QKQ)P
\]

---

# 4. 结论1：动力学破坏封闭的来源

PRP不封闭 ⇔ 存在以下任一结构：

\[
PKQ \neq 0 \quad \text{或} \quad QKP \neq 0
\]

或更一般：

\[
K \text{在 } P/Q \text{分解下存在 off-diagonal block}
\]

---

# 5. 从动力学到Schur complement

## 5.1 从 resolvent 出发

考虑：

\[
G(z) = (z-K)^{-1}
\]

在 P/Q 分解下：

\[
K =
\begin{pmatrix}
K_{PP} & K_{PQ} \\
K_{QP} & K_{QQ}
\end{pmatrix}
\]

---

## 5.2 Block matrix inversion 定理

严格恒等式：

\[
G_{PP}(z)
=
(z - K_{PP} - K_{PQ}(z-K_{QQ})^{-1}K_{QP})^{-1}
\]

---

## 5.3 定义有效哈密顿量

\[
K_{\mathrm{eff}}(z)
=
K_{PP} + \Sigma(z)
\]

其中

\[
\Sigma(z)=K_{PQ}(z-K_{QQ})^{-1}K_{QP}
\]

---

# 6. Schur补物理意义

## 6.1 结构解释

\[
\Sigma(z)
=
P \to Q \to P \text{ 虚跃迁通道}
\]

---

## 6.2 展开形式

\[
(z-K_{QQ})^{-1} = \sum_n \frac{|n\rangle\langle n|}{z-\epsilon_n}
\]

因此：

\[
\Sigma(z)=\sum_n \frac{K_{PQ}|n\rangle\langle n|K_{QP}}{z-\epsilon_n}
\]

---

# 7. PRP封闭与Schur结构的等价性

## 7.1 封闭 ⇔ 无Schur补

\[
PRP \text{封闭} \iff \Sigma(z)=0
\]

等价于：

\[
K_{PQ}=0,\quad K_{QP}=0
\]

---

# 8. MZM / ABS的严格判据

---

## 8.1 MZM（拓扑保护）

\[
K_{PQ}=K_{QP}=0
\]

\[
\Rightarrow \Sigma(z)=0
\]

\[
\Rightarrow P \text{子空间严格不混合}
\]

\[
\Rightarrow E=0 \text{零模稳定}
\]

---

## 8.2 ABS（混合态）

\[
K_{PQ}\neq0 \quad \text{或} \quad K_{QP}\neq0
\]

\[
\Rightarrow \Sigma(z)\neq0
\]

\[
\Rightarrow E=0 \to E=\pm\epsilon
\]

---

# 9. 在R(u)中的对应结构

\[
R(u)=U(u)e^{-iK(u)}
\]

投影：

\[
PR(u)P
\]

等价于：

\[
R_{PP}(u) - R_{PQ}(u)R_{QQ}^{-1}(u)R_{QP}(u)
\]

---

# 10. YBE在投影后的含义

原始YBE：

\[
R_{12}R_{13}R_{23}=R_{23}R_{13}R_{12}
\]

投影后：

\[
R^{\mathrm{eff}} = R_{PP} - \Sigma
\]

---

## 10.1 结论

YBE成立 ⇔ Schur补满足一致性条件（cocycle-like constraint）

否则：

\[
\Sigma \text{导致YBE变形}
\]

---

# 11. 最终统一图像

---

## MZM

\[
\Sigma=0
\]

- PRP严格封闭
- YBE严格成立
- pure braid representation

---

## ABS

\[
\Sigma\neq0
\]

- PRP仅近似封闭
- YBE变形
- dynamical mixing存在

---

# 12. 最终一句话总结

\[
\boxed{
PRP封闭性 = P子空间在K的Schur补意义下是否不与Q发生虚跃迁耦合
}
\]

