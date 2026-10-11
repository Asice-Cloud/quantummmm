# Model

---

## 一、出发点：有效模型

### 1.1 物理场景

研究 ABS 存在时Majorana 的 braiding 性质。系统包含 5 个 Majorana 模式：

- $\gamma_1,\gamma_2,\gamma_3$：MZM 模式（$\gamma_2$ 和 $\gamma_3$ 将被交换）
- $\gamma_a,\gamma_b$：量子点 ancilla 模式 

有效哈密顿量：

$$
H_{EM}(t) = iE_d\gamma_a\gamma_b + iE_1\gamma_1\gamma_2 + i|t_2|\gamma_a\gamma_2 - i|t_1|\gamma_b\gamma_1 - i|t_3|\gamma_a\gamma_3
$$

| 参数 | 含义 | 理想 MZM | ABS |
|---|---|---|---|
| $E_1$ | $\gamma_1$–$\gamma_2$ 杂化能 | 0 | $\neq 0$ |
| $t_1$ | $\gamma_1$–ancilla 耦合 | 0 | $\neq 0$ |
| $t_2,t_3$ | 门控编织耦合 | 时变 | 时变 |
| $E_d$ | 量子点能级 | 时变 | 时变 |

### 1.2 三段 braiding 协议

按原文第二节，每段时长为 $\tau$，三个步骤完成一次 $\gamma_2,\gamma_3$ 交换，耗时 $3\tau$；将三段协议重复两遍完成两次交换，总时间为 $6\tau$。论文的目标重叠为 $|\langle\psi_1^+(0)|\psi_1^-(6\tau)\rangle|$。

论文 Fig. 1(d) 横轴的 $100/\mathrm{meV}$ 是 $\tau$ 的单位（取 $\hbar=1$），并不表示 $\tau$ 在扫描中固定。本项目当前推导保留 $\tau$ 的符号，把它作为与门控幅度 $t_c$ 一样的可控设计量；此前数值核对取 $\tau=100/\mathrm{meV}$，此时一次交换为 $300/\mathrm{meV}$，两次交换为 $600/\mathrm{meV}$。一般情况下分别为 $3\tau$ 和 $6\tau$；$12\tau$ 对应四次交换。

与 Fig. 1(d) 比较时，按图注取 $E_1=0.01\,\mathrm{meV}$；正文的 $0.001\,\mathrm{meV}$ 按已确认的笔误处理。

每段使用局部时间 $s\in[0,\tau]$，门控函数为 $f_\pm(s)=\frac{1\pm\cos(\pi s/\tau)}{2}$：

```
Step 1 (0→τ):   G1 关 → t₂ 打开, E_d → 0         γ₂ 移入量子点
Step 2 (τ→2τ):  G2 关, G1 开 → t₃ 打开, t₂ 关闭  γ₃ 接管, γ₂ 退出
Step 3 (2τ→3τ): G2 开 → t₃ 关闭, E_d 恢复        γ₃ 退出，控制回到初始设置
                                    ─────────────
                                    结果: γ₂↔γ₃ 交换
```

控制在 $3\tau$ 后重复，量子态继续从第一轮的终态演化；第二轮不重新置初态。

---

## 二、为什么是 so(5)？

### 2.1 双线性生成元的李代数闭包

定义 $X_{ij}=i\gamma_i\gamma_j$（共 $C_2^5=10$ 个）。它们满足标准 $so(5)$ 对易关系：

$$[X_{ij}, X_{kl}] = 2i(\delta_{jk}X_{il} - \delta_{ik}X_{jl} - \delta_{jl}X_{ik} + \delta_{il}X_{jk})$$

### 2.2 论文哈密顿量在 so(5) 基底中的投影

固定五维索引顺序为 $(1,2,3,a,b)$，在矩阵下标中以 $4=a,5=b$。按升序指标定义 $X_{ij}=i\gamma_i\gamma_j$，则 $H_{EM}$ 恰好是 5 个生成元的线性组合：

$$H_{EM} = E_1 X_{12} - |t_2| X_{24} + |t_1| X_{15} + |t_3| X_{34} + E_d X_{45}$$

三个耦合项的符号来自原始乘积的顺序：$\gamma_a\gamma_2=-\gamma_2\gamma_a$、$\gamma_b\gamma_1=-\gamma_1\gamma_b$、$\gamma_a\gamma_3=-\gamma_3\gamma_a$。这里的幅度表示相应时刻门控后的耦合。

**三段协议本质上是在 10 维 so(5) 李代数里的分段时变演化，不预设二能级投影。**

## so(5)表示

Cl(5) Gamma 矩阵

5 个 $2\times2$ 四元数矩阵：

$$\Gamma_1=\begin{pmatrix}0&1\\1&0\end{pmatrix},\;
\Gamma_2=\begin{pmatrix}0&-\mathbf i\\ \mathbf i&0\end{pmatrix},\;
\Gamma_3=\begin{pmatrix}0&-\mathbf j\\ \mathbf j&0\end{pmatrix},\;
\Gamma_4=\begin{pmatrix}0&-\mathbf k\\ \mathbf k&0\end{pmatrix},\;
\Gamma_5=\begin{pmatrix}1&0\\0&-1\end{pmatrix}$$

满足 $\{\Gamma_i,\Gamma_j\}=2\delta_{ij}$。



 旋量生成元：$\Sigma_{ij} = \frac14[\Gamma_i, \Gamma_j]$

以 $\Sigma_{12}$ 为例显式计算：

$$\Gamma_1\Gamma_2 = \begin{pmatrix}0&1\\1&0\end{pmatrix}\begin{pmatrix}0&-\mathbf i\\\mathbf i&0\end{pmatrix}
= \begin{pmatrix}\mathbf i&0\\0&-\mathbf i\end{pmatrix}$$

$$\Gamma_2\Gamma_1 = \begin{pmatrix}0&-\mathbf i\\\mathbf i&0\end{pmatrix}\begin{pmatrix}0&1\\1&0\end{pmatrix}
= \begin{pmatrix}-\mathbf i&0\\0&\mathbf i\end{pmatrix}$$

$$[\Gamma_1,\Gamma_2] = \Gamma_1\Gamma_2 - \Gamma_2\Gamma_1 = \begin{pmatrix}2\mathbf i&0\\0&-2\mathbf i\end{pmatrix}$$

$$\Sigma_{12} = \frac14[\Gamma_1,\Gamma_2] = \begin{pmatrix}\mathbf i/2&0\\0&-\mathbf i/2\end{pmatrix}$$

同理可算得全部 5 个活跃生成元的分块形式：

| $\Sigma_{12}(E_1)$ | $\Sigma_{24}(-|t_2|)$ | $\Sigma_{15}(+|t_1|)$ | $\Sigma_{34}(+|t_3|)$ | $\Sigma_{45}(E_d)$ |
|---|---|---|---|---|---|
| $\begin{pmatrix}\mathbf i/2&0\\0&-\mathbf i/2\end{pmatrix}$ | $\begin{pmatrix}\mathbf j/2&0\\0&\mathbf j/2\end{pmatrix}$ | $\begin{pmatrix}0&-1/2\\1/2&0\end{pmatrix}$ | $\begin{pmatrix}-\mathbf i/2&0\\0&-\mathbf i/2\end{pmatrix}$ | $\begin{pmatrix}0&\mathbf k/2\\\mathbf k/2&0\end{pmatrix}$ |

### 旋量矩阵与五维演化生成元的对应

将四元数表示成 $2\times2$ 复矩阵后，上述 $\Gamma_i$、$\Sigma_{ij}$ 是 $4\times4$ 复旋量矩阵。标量虚数 $i$ 与四元数单位 $\mathbf i,\mathbf j,\mathbf k$ 要区分。因为 $\Gamma_i\Gamma_j=2\Sigma_{ij}$，哈密顿量的旋量表示为

$$
H_{\rm spin}=2i\left(E_1\Sigma_{12}-|t_2|\Sigma_{24}
+|t_1|\Sigma_{15}+|t_3|\Sigma_{34}+E_d\Sigma_{45}\right).
$$

对应的五维实生成元是

$$
K=2\left(E_1L_{12}-|t_2|L_{24}+|t_1|L_{15}
+|t_3|L_{34}+E_dL_{45}\right),\qquad
L_{ij}=e_i e_j^T-e_j e_i^T.
$$

这正是当前 $\partial_t\mathcal B=[K,\mathcal B]$ 使用的生成元。原矩阵表的三个耦合标签已按升序指标修正；矩阵本身及现有 PDE 程序的符号通过了独立核对。完整推导和检验见 [有效哈密顿量核对](fixed_protocol_pde/effective_hamiltonian_audit.md)。
