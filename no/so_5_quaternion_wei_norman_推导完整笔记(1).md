# so(5) Quaternion Wei–Norman 分解模型（完整推导笔记）

> 说明：本文整理为一套自洽的 Wei–Norman 参数化推导框架，用于 quaternion 表示下的 so(5) 动力学分解。结构按照 Lie 代数闭合 → 指数分解 → 参数微分方程 → 与物理哈密顿对应的标准流程组织。

---

# 1. 目标问题

考虑一个动力学演化算符

\[
U(t) \in SO(5)
\]

满足

\[
\dot U(t) = H(t) U(t), \quad H(t) \in \mathfrak{so}(5)
\]

目标是进行 Wei–Norman 分解：

\[
U(t) = \prod_{k=1}^{10} \exp\big(x_k(t) X_k\big)
\]

其中 \(\{X_k\}\) 为 so(5) 的生成元（10维），并选取 quaternion-compatible basis。

---

# 2. so(5) 与 quaternion 分解结构

so(5) 可分解为：

\[
\mathfrak{so}(5) \simeq \mathfrak{so}(4) \oplus \mathbb{R}^4
\]

进一步：

\[
\mathfrak{so}(4) \simeq \mathfrak{su}(2)_L \oplus \mathfrak{su}(2)_R
\]

因此可选 basis：

- 左 quaternion：\(\vec{L} = (L_x, L_y, L_z)\)
- 右 quaternion：\(\vec{R} = (R_x, R_y, R_z)\)
- 混合生成元：\(\{Q_a\}_{a=1}^4\)

总计 10 个生成元。

---

# 3. Wei–Norman 分解 ansatz

选取有序分解：

\[
U(t) = e^{x_1 Q_1} e^{x_2 Q_2} e^{x_3 Q_3} e^{x_4 Q_4}
        e^{y_1 L_1} e^{y_2 L_2} e^{y_3 L_3}
        e^{z_1 R_1} e^{z_2 R_2} e^{z_3 R_3}
\]

定义参数向量：

\[
\theta(t) = (x_a, y_i, z_i)
\]

---

# 4. Wei–Norman 核心公式

对时间导数：

\[
\dot U U^{-1} = H(t)
\]

展开：

\[
\dot U U^{-1} = \sum_i \dot\theta_i \; \mathrm{Ad}_{e^{\theta_{<i}}}(X_i)
\]

其中：

\[
\mathrm{Ad}_g(X) = g X g^{-1}
\]

---

# 5. so(5) Lie 代数结构常数

定义：

\[
[X_a, X_b] = f_{ab}^{\ \, c} X_c
\]

关键结构：

### (1) quaternion SU(2) 子代数

\[
[L_i, L_j] = \epsilon_{ijk} L_k
\]

\[
[R_i, R_j] = \epsilon_{ijk} R_k
\]

\[
[L_i, R_j] = 0
\]

---

### (2) 混合生成元

\[
[L_i, Q_a] = (\sigma_i)_{ab} Q_b
\]

\[
[R_i, Q_a] = (\tilde\sigma_i)_{ab} Q_b
\]

其中 \(\sigma_i, \tilde\sigma_i\) 是 quaternion 左右作用矩阵。

---

# 6. BCH 展开与 Ad 联合结构

Wei–Norman 关键步骤：

\[
e^A B e^{-A} = B + [A,B] + \frac{1}{2}[A,[A,B]] + \cdots
\]

因此：

\[
\mathrm{Ad}_{e^{x_i X_i}}(X_j) = e^{x_i \mathrm{ad}_{X_i}} X_j
\]

其中：

\[
\mathrm{ad}_X(Y) = [X,Y]
\]

---

# 7. 投影方程（核心步骤）

将：

\[
H(t) = \sum_k h_k(t) X_k
\]

与：

\[
\dot U U^{-1}
\]

逐基展开，得到非线性 ODE：

\[
\dot \theta_i = F_i(\theta_1, ..., \theta_{i-1}; h(t))
\]

结构为**下三角系统**：

- 每个参数只依赖前面变量
- 保证可逐层积分

---

# 8. quaternion 表示简化

定义 quaternion basis：

\[
q = a_0 + a_1 i + a_2 j + a_3 k
\]

so(4) 映射：

\[
SO(4) \simeq (SU(2)_L \times SU(2)_R)/\mathbb{Z}_2
\]

作用：

\[
q \mapsto u_L q u_R^{-1}
\]

对应生成元：

\[
L_i q = \frac{1}{2} \sigma_i q,
\quad
R_i q = \frac{1}{2} q \sigma_i
\]

---

# 9. so(5) 中 extra 维度（关键结构）

额外 4 个生成元 \(Q_a\)：

- 连接 quaternion space 与 extra dimension
- 形成 spinor-like representation

结构：

\[
[Q_a, Q_b] = L_{ab} + R_{ab}
\]

\[
[L_i, Q_a] = (\sigma_i Q)_a
\]

---

# 10. Wei–Norman 最终方程结构

最终系统写为：

\[
\dot x_a = h_a(t) + \sum_{b<a} C_{ab}(x) h_b(t)
\]

\[
\dot y_i = h^{(L)}_i(t) + f_i(x,y)
\]

\[
\dot z_i = h^{(R)}_i(t) + g_i(x,y,z)
\]

---

# 11. 物理解释

该分解对应：

- quaternion rotation（SU(2)_L × SU(2)_R）
- extra Majorana / ABS coupling sector
- so(5) embedding of extended fermionic Hilbert space

---

# 12. 结论

Wei–Norman 方法将 so(5) 动力学转化为：

- 可积的非线性 triangular ODE system
- quaternion × mixed generator decomposition
- 适用于 topological / Majorana / ABS coupling models

---

# 附录：推荐下一步

如果继续推进，可以做：

1. 显式写出 10×10 adjoint representation matrix
2. 用 PRB111 的具体 Hamiltonian 投影 h_k(t)
3. 推导 ABS coupling 的 effective gauge field
4. 数值解 Wei–Norman ODE system

