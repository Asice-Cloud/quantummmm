# Green函数与路径演化统一结构（R(u)修正版本）

---

# 1. 两类核心对象

在该体系中存在两个基本结构：

## (1) Green函数（resolvent结构）

\[
G(z)=(z-H)^{-1}
\]

其中：

- \(z \in \mathbb{C}\)：谱参数（energy / frequency）
- 极点 \(z=\epsilon_n\)：本征能

---

## (2) 路径演化算符（braid / holonomy）

严格定义为：

\[
R(u)=\mathcal{T}\exp\left(-i\int_0^u H(s)ds\right)
\]

其中：

- \(u\)：路径参数（braid / adiabatic parameter）
- \(H(s)\)：参数依赖 Hamiltonian
- \(\mathcal{T}\)：路径有序算符

---

# 2. 路径演化的微分形式

路径演化满足：

\[
\frac{d}{du}R(u)=-iH(u)R(u)
\]

初值条件：

\[
R(0)=I
\]

---

# 3. 非对易性结构

## (1) Abelian情况

若满足：

\[
[H(u_1),H(u_2)]=0
\]

则：

\[
R(u)=\exp\left(-i\int_0^u H(s)ds\right)
\]

---

## (2) 非Abelian情况（一般情形）

若：

\[
[H(u_1),H(u_2)]\neq 0
\]

则必须使用：

\[
R(u)=\mathcal{T}\exp\left(-i\int_0^u H(s)ds\right)
\]

---

# 4. Dyson展开

路径演化可展开为：

\[
R(u)=I
-i\int_0^u ds_1 H(s_1)
\]

\[
+(-i)^2\int_0^u ds_1\int_0^{s_1}ds_2 H(s_1)H(s_2)+\cdots
\]

该结构体现：

> 非对易序列结构（ordered operator algebra）

---

# 5. Green函数的谱结构

对哈密顿量对角化：

\[
H|n\rangle=\epsilon_n|n\rangle
\]

Green函数可写为：

\[
G(z)=\sum_n \frac{|n\rangle\langle n|}{z-\epsilon_n}
\]

---

# 6. Schur补结构（P/Q分解）

Hilbert空间分解：

\[
\mathcal{H}=\mathcal{H}_P \oplus \mathcal{H}_Q
\]

哈密顿量块结构：

\[
H=\begin{pmatrix}
H_{PP} & H_{PQ} \\
H_{QP} & H_{QQ}
\end{pmatrix}
\]

Green函数满足：

\[
G_{PP}(z)=(z-H_{PP}-\Sigma(z))^{-1}
\]

其中：

\[
\Sigma(z)=H_{PQ}(z-H_{QQ})^{-1}H_{QP}
\]

---

# 7. 物理意义（Green函数侧）

## (1) 虚跃迁结构

\[
\Sigma(z): P \rightarrow Q \rightarrow P
\]

## (2) 物理效应

- 实部：能级重整化
- 虚部：寿命 / broadening

---

# 8. 路径演化的几何意义

R(u) 是一个：

\[
\boxed{\text{non-Abelian holonomy}}
\]

对应：

- connection transport
- braid group representation
- adiabatic evolution

---

# 9. Green函数与路径演化的对偶关系

两者通过以下关系连接：

\[
G(z)=i\int_0^\infty du\, e^{izu}R(u)
\]

等价于：

\[
G(z)=\mathcal{L}[R(u)]
\]

因此：

| 对象 | 类型 |
|------|------|
| \(z\) | 谱参数 |
| \(u\) | 路径参数 |
| \(G(z)\) | resolvent结构 |
| \(R(u)\) | path-ordered evolution |

---

# 10. 统一结构总结

整个体系可以分为两层：

## (1) 谱空间（Green函数）

\[
G(z)=(z-H)^{-1}
\]

→ 描述能级结构与虚跃迁

---

## (2) 路径空间（braid演化）

\[
R(u)=\mathcal{T}\exp\left(-i\int_0^u H(s)ds\right)
\]

→ 描述非对易几何输运

---

# 11. 核心结论

\[
\boxed{
\text{Green函数（谱） + 路径演化（几何） = 完整non-Abelian结构}
}
\]

其中：

- Green函数编码：bulk虚传播（Σ(z)）
- R(u)编码：braid / holonomy结构

---

