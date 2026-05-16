# ABS谱函数与non-Abelian相图计算方法（Σ → PRB复现完整流程）

---

# 0. 总目标

建立从 Green函数出发的完整可计算链：

\[
\boxed{
H \rightarrow G(z) \rightarrow \Sigma(z) \rightarrow A(\omega) \rightarrow \mathcal{N} \rightarrow \text{phase diagram}
}
\]

用于复现 PRB 111.205411 中 non-Abelian 结构。

---

# 1. 起点：Green函数

给定哈密顿量：

\[
G(z)=(z-H)^{-1}
\]

P/Q分解：

\[
H=
\begin{pmatrix}
H_{PP} & H_{PQ} \\
H_{QP} & H_{QQ}
\end{pmatrix}
\]

投影Green函数：

\[
G_{PP}(z)=P(z-H)^{-1}P
\]

---

# 2. Schur补（核心量）

\[
G_{PP}(z)=(z-H_{PP}-\Sigma(z))^{-1}
\]

其中：

\[
\boxed{
\Sigma(z)=H_{PQ}(z-H_{QQ})^{-1}H_{QP}
}
\]

展开形式：

\[
\Sigma(z)=\sum_n \frac{V_n V_n^\dagger}{z-\epsilon_n}
\]

---

# 3. Σ(z)的物理意义

- Re Σ：能级重整化
- Im Σ：衰减 / broadening
- P→Q→P虚传播

定义耦合强度：

\[
\Gamma(\omega)=-\mathrm{Im}\,\Sigma(\omega)
\]

---

# 4. ABS谱函数

## 定义

\[
\boxed{
A(\omega)=-\frac{1}{\pi}\mathrm{Im}\,\mathrm{Tr}\,G_{PP}(\omega+i\eta)
}
\]

代入Schur补：

\[
A(\omega)= -\frac{1}{\pi}\mathrm{Im}\,\mathrm{Tr}\frac{1}{\omega-H_{PP}-\Sigma(\omega)}
\]

---

# 5. ABS判据

## (1) 极点条件

\[
\det(\omega-H_{PP}-\Sigma(\omega))=0
\]

## (2) 有限寿命

\[
\mathrm{Im}\,\Sigma(\omega)\neq 0
\]

## (3) ABS定义

谱函数中的宽峰：

- 非δ函数
- 有finite width

---

# 6. non-Abelian量化指标

PRB定义：

\[
\mathcal{N}_{PRB}=\|U_iU_jU_i-U_jU_iU_j\|
\]

在本模型中：

\[
U_i = R_i(u)
\]

---

## 6.1 effective映射

\[
R(u)=\mathcal{T}\exp\left(-i\int_0^u H_{eff}(s)ds\right)
\]

\[
H_{eff}(z)=H_{PP}+\Sigma(z)
\]

---

## 6.2 non-Abelian指标

\[
\boxed{
\mathcal{N}=\|R_iR_jR_i-R_jR_iR_j\|
}
\]

---

# 7. ABS ↔ non-Abelian关系

\[
\boxed{
\mathcal{N}(\omega) \sim \Gamma(\omega)=-\mathrm{Im}\Sigma(\omega)
}
\]

---

# 8. phase diagram构造

选择参数空间：

- 横轴：\(\mu\)（化学势）
- 纵轴：\(\lambda\)（bulk耦合强度）

---

# 9. phase indicator

定义：

\[
P(\mu,\lambda)=\int d\omega\;A(\omega)\Gamma(\omega)
\]

---

# 10. 相结构

## (1) MZM phase

\[
\Sigma(\omega)\approx 0
\Rightarrow P\approx 0
\]

---

## (2) ABS phase

\[
\Sigma(\omega)\neq 0
\Rightarrow P>0
\]

---

## (3) crossover region

\[
P \sim O(1)
\]

---

# 11. phase boundary

由以下条件确定：

\[
\boxed{
|\Sigma(\omega)| \sim \Delta_{topo}
}
\]

或：

\[
\mathrm{Im}\det(\omega-H_{PP}-\Sigma)=0
\]

---

# 12. PRB风格图（计算输出）

## 图1：谱函数

\[
A(\omega)
\]

- MZM：δ峰
- ABS：展宽峰

---

## 图2：non-Abelian强度

\[
\mathcal{N}(\mu,\lambda)
\]

---

## 图3：phase diagram

颜色标度：

- 蓝：MZM（Σ≈0）
- 红：ABS（Σ≠0）
- 过渡区：混合态

---

# 13. 完整计算流程（可实现）

## Step 1
构造 H

## Step 2
计算 G(z)

## Step 3
提取 Σ(z)

## Step 4
计算 A(ω)

## Step 5
构造 R(u)

## Step 6
计算 braid deviation

\[
R_iR_jR_i - R_jR_iR_j
\]

## Step 7
计算 non-Abelian指标

\[
\mathcal{N}
\]

## Step 8
扫描参数 (μ, λ)

## Step 9
绘制 phase diagram

---

# 14. 最终统一结论

\[
\boxed{
H \rightarrow \Sigma(z) \rightarrow A(\omega) \rightarrow \mathcal{N} \rightarrow \text{phase diagram}
}
\]

---

# 15. 一句话总结

non-Abelian性在本模型中等价于 Green函数 Schur补 Σ(z) 在谱函数中诱导的宽化结构，并可直接映射为 PRB 111.205411 中的 braid deviation 量。

