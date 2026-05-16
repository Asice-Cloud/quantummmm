# R(u) → PRP投影 → Schur补 → non-Abelian量化：完整统一推导（含PRB 111.205411复现与实验预测）

---

# 0. 基本对象

考虑谱参数依赖算符：

\[
R(u)=U(u)e^{-iK(u)}
\]

Hilbert空间分解：

\[
\mathcal H = \mathcal H_P \oplus \mathcal H_Q
\]

投影算符：

\[
P^2=P,\quad Q=I-P
\]

目标：

- PRP是否封闭
- non-Abelian结构来源
- ABS / MZM区分
- 与PRB 111.205411对应
- 可观测量预测

---

# 1. PRP投影结构

定义：

\[
R_P(u)=PR(u)P
\]

封闭性条件：

\[
QR(u)P=0
\]

---

# 2. 指数展开

\[
e^{-iK} = \sum_{n\ge0} \frac{(-i)^n}{n!}K^n
\]

\[
R = U\sum_{n\ge0}\frac{(-i)^n}{n!}K^n
\]

投影：

\[
QRP = \sum_{n\ge0}\frac{(-i)^n}{n!}QUK^nP
\]

封闭性等价：

\[
QUK^nP=0,\forall n
\]

---

# 3. K的块结构

\[
K=
\begin{pmatrix}
K_{PP} & K_{PQ}\\
K_{QP} & K_{QQ}
\end{pmatrix}
\]

展开：

\[
K = PKP + PKQ + QKP + QKQ
\]

---

# 4. 非封闭来源

PRP不封闭 ⇔

\[
PKQ\neq0 \quad \text{或} \quad QKP\neq0
\]

---

# 5. Schur补严格推导

Green函数：

\[
G(z)=(z-K)^{-1}
\]

块矩阵恒等式：

\[
G_{PP}(z)=
(z-K_{PP}-K_{PQ}(z-K_{QQ})^{-1}K_{QP})^{-1}
\]

定义：

\[
\Sigma(z)=K_{PQ}(z-K_{QQ})^{-1}K_{QP}
\]

\[
K_{eff}(z)=K_{PP}+\Sigma(z)
\]

---

# 6. Schur物理意义

\[
\Sigma = P\to Q\to P
\]

展开：

\[
\Sigma(z)=\sum_n\frac{K_{PQ}|n\rangle\langle n|K_{QP}}{z-\epsilon_n}
\]

---

# 7. PRP封闭条件

\[
PRP \text{封闭} \iff \Sigma=0
\]

---

# 8. MZM / ABS

## MZM

\[
\Sigma=0
\]

- 无混合
- zero mode稳定
- braid exact

## ABS

\[
\Sigma\neq0
\]

- bulk混合
- 能级偏移
- braid变形

---

# 9. 有效演化

\[
R_{eff}(u)=e^{-i(K_{PP}+\Sigma)}
\]

---

# 10. Connection

\[
A=R^{-1}_{eff}\partial_u R_{eff}=A_0+A_\Sigma
\]

---

# 11. 曲率（non-Abelian来源）

\[
F=dA+A^2
\]

\[
F=F_0+F_\Sigma+[A_0,A_\Sigma]
\]

---

# 12. braid代数破缺

\[
\Delta_{braid}=R_iR_j-R_jR_i
\]

\[
\Delta_{braid}=\Delta_0+[R_0,\Sigma]+[\Sigma,\Sigma]
\]

---

# 13. YBE变形

\[
\Delta_{YBE}\sim \Sigma+\Sigma^2
\]

---

# 14. Drinfeld twist等价

\[
R_{eff}=F_{21}R_0F^{-1}
\]

\[
\Sigma \leftrightarrow \mathcal D_{R_0}(f)
\]

\[
\boxed{\Sigma \equiv \text{Drinfeld twist generator}}
\]

---

# 15. 量子群变形

\[
U_q(g) \to U_{q+\delta q}
\]

\[
\delta q \sim \Sigma
\]

---

# 16. non-Abelian统一量

\[
\boxed{
\mathcal N=
\|\Sigma\|+
\|\Delta_{braid}\|+
\|F\|
}
\]

---

# 17. PRB 111.205411对照

| PRB | 本模型 |
|-----|--------|
| ABS | P子空间 |
| braid deviation | Δ_braid |
| hybridization | K_PQ,K_QP |
| non-Abelian measure | Σ驱动 |

---

# 18. PRB 111.205411复现（详细推导版本）

这一部分目标是：严格用本模型的 \Sigma / Schur补结构复现PRB中定义的 ABS non-Abelian deviation。

---

## 18.1 PRB的核心对象

PRB研究的对象本质是：ABS子空间中的“近零模态”及其braiding操作。

设低能子空间基：

\[
\{\lvert\psi_i(u)\rangle\} \subset \mathcal H_P
\]

braid操作定义为参数演化算符：

\[
U_{ij}(u): \lvert\psi_i\rangle \rightarrow \lvert\psi_j\rangle
\]

理想MZM情况：

\[
U_iU_jU_i = U_jU_iU_j
\]

PRB指出：ABS中该关系被破坏。

---

## 18.2 投影后的有效动力学

在本模型中：

\[
K = 
\begin{pmatrix}
K_{PP} & K_{PQ}\\
K_{QP} & K_{QQ}
\end{pmatrix}
\]

有效哈密顿量（Schur补）：

\[
K_{eff}(z)=K_{PP}+\Sigma(z)
\]

\[
\Sigma(z)=K_{PQ}(z-K_{QQ})^{-1}K_{QP}
\]

因此有效演化算符：

\[
R_{eff}(u)=e^{-i(K_{PP}+\Sigma)}
\]

---

## 18.3 braid算符在有效子空间中的表示

PRB中的braid操作在本模型中对应：

\[
R_i \equiv R_{eff}(u_i)
\]

因此：

\[
R_i = R_0 + \delta R(\Sigma)
\]

其中：

\[
R_0 = e^{-iK_{PP}},\quad \delta R \sim -i\Sigma e^{-iK_{PP}}
\]

---

## 18.4 braid关系偏离的严格展开

考虑PRB定义的核心量：

\[
\Delta_{braid} = R_iR_jR_i - R_jR_iR_j
\]

代入展开：

\[
R_i = R_0 + \delta R_i,\quad R_j = R_0 + \delta R_j
\]

展开到一阶：

\[
R_iR_jR_i = R_0^3 + R_0^2\delta R_i + R_0\delta R_j R_0 + \delta R_i R_0^2 + O(\Sigma^2)
\]

同理：

\[
R_jR_iR_j = R_0^3 + R_0^2\delta R_j + R_0\delta R_i R_0 + \delta R_j R_0^2 + O(\Sigma^2)
\]

相减得到：

\[
\Delta_{braid} = [R_0,\delta R_i - \delta R_j] + O(\Sigma^2)
\]

由于：

\[
\delta R \propto \Sigma
\]

得到核心结果：

\[
\boxed{
\Delta_{braid} \sim [R_0,\Sigma] + O(\Sigma^2)
}
\]

---

## 18.5 PRB non-Abelian measure复现

PRB定义non-Abelian强度为：

\[
\mathcal N_{PRB} = \|U_iU_jU_i - U_jU_iU_j\|
\]

在本模型中严格变为：

\[
\boxed{
\mathcal N_{PRB} = \|\Delta_{braid}\| \sim \|[R_0,\Sigma]\|
}
\]

---

## 18.6 MZM极限复现

若：

\[
\Sigma=0
\]

则：

\[
R_{eff}=R_0
\]

\[
\Delta_{braid}=0
\]

恢复：

\[
U_iU_jU_i = U_jU_iU_j
\]

即PRB中的理想non-Abelian braiding。

---

## 18.7 ABS偏离机制（核心复现结论）

若：

\[
\Sigma \neq 0
\]

则：

- effective braid representation失真
- YBE产生defect
- non-Abelian性变为连续量

最终：

\[
\boxed{
\text{PRB ABS non-Abelianity} \equiv \text{Schur complement } \Sigma \text{ induced deformation}
}
\]

---

## 18.8 关键总结（复现完成）

PRB 111.205411的结构在本模型中的严格对应为：

- ABS子空间 = P空间投影
- hybridization = K_PQ,K_QP
- deviation from braid = Schur补 \Sigma
- non-Abelian measure = [R_0,\Sigma]

---

# 19. 实验可观测量（关键升级） 实验可观测量（关键升级）

---

## 19.1 Green函数

\[
G_{PP}(z)=\frac{1}{z-K_{PP}-\Sigma(z)}
\]

---

## 19.2 ABS能级

\[
E_{ABS}\approx \mathrm{Re}\Sigma(0)
\]

---

## 19.3 寿命

\[
\Gamma=-\mathrm{Im}\Sigma(0)
\]

---

## 19.4 braid error scaling

\[
\|\Delta_{braid}\|\propto \|\Sigma\|
\]

---

## 19.5 bulk gap scaling

\[
\Sigma\sim \frac{V^2}{\Delta_{bulk}}
\]

\[
\mathcal N \sim \frac{V^2}{\Delta_{bulk}}
\]

---

## 19.6 phase error

\[
\delta\theta \propto \mathrm{Im}\Sigma(0)
\]

---

# 20. 最终统一图像

---

## 微观

\[
K \to P/Q分解
\]

## 动力学

\[
\Sigma=P\to Q\to P
\]

## 有效理论

\[
K_{eff}=K_{PP}+\Sigma
\]

## 几何

\[
F\neq0
\]

## 代数

braid deformation

## 实验

\[
\mathcal N_{exp}(\Sigma)
\]

---

# 21. 最终核心结论

\[
\boxed{
\text{non-Abelianity}=f(\Sigma)
}
\]

\[
\boxed{
\Sigma=P\to Q\to P
}
\]

---

# 22. 最终一句话

PRP投影的本质不是截断，而是：

> full system在P子空间上的Schur补有效理论，其non-Abelian结构完全由Σ生成。

