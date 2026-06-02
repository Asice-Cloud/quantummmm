### 各段的李代数闭包

在 10 个 so(5) 生成元 $X_1,\dots,X_{10}$ 中：

| $X_1$               | $X_2$               | $X_3$               | $X_4$               | $X_5$               | $X_6$               | $X_7$               | $X_8$               | $X_9$               | $X_{10}$            |
| ------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- | ------------------- |
| $i\gamma_1\gamma_2$ | $i\gamma_1\gamma_3$ | $i\gamma_2\gamma_3$ | $i\gamma_1\gamma_a$ | $i\gamma_a\gamma_2$ | $i\gamma_a\gamma_3$ | $i\gamma_1\gamma_b$ | $i\gamma_2\gamma_b$ | $i\gamma_3\gamma_b$ | $i\gamma_a\gamma_b$ |

**Step 1**：活跃 $\{X_1, X_5, X_7, X_{10}\}$

通过对易子：
$$[X_1,X_5]\to X_4,\quad [X_1,X_7]\to X_8,\quad [X_5,X_{10}]\to X_8,\quad [X_7,X_{10}]\to X_4$$

闭包：$\{X_1, X_4, X_5, X_7, X_8, X_{10}\}$，维数 **6**。
这是 $\{\gamma_1,\gamma_2,\gamma_a,\gamma_b\}$ 上的 $so(4)\cong su(2)\oplus su(2)$。
$\gamma_3$ 为旁观者spectator。

**Step 2**：活跃 $\{X_1, X_5, X_6, X_7\}$

对易子链：
$$[X_1,X_5]\to X_4,\;[X_1,X_7]\to X_8,\;[X_5,X_6]\to X_3,\;[X_1,X_6]\to X_2,\;[X_6,X_7]\to X_9,\dots$$

闭包：全部 10 个生成元，维数 **10**——完整的 $so(5)$。commutant 平凡（仅有标量），
无法做固定块对角化。

**Step 3**：活跃 $\{X_1, X_6, X_{10}\}$

$X_1$ 与 $X_6,X_{10}$ 对易（$\gamma_1,\gamma_2$ 与 $\gamma_3,\gamma_a,\gamma_b$ 无重合指标）。
$\{X_6,X_9,X_{10}\}$ 形成封闭 su(2) 三元组：
$$[X_6,X_{10}]\to X_9,\quad [X_9,X_6]\to X_{10},\quad [X_{10},X_9]\to X_6$$

闭包：$\{X_1, X_6, X_9, X_{10}\}$，维数 **4**。
结构：$u(1)_{X_1}\oplus su(2)_{\{X_6,X_9,X_{10}\}}$。

---

| 段     | 活跃生成元           | 闭包       | 维数 | 代数结构                       |
| ------ | -------------------- | ---------- | ---- | ------------------------------ |
| Step 1 | $X_1,X_5,X_7,X_{10}$ | $+X_4,X_8$ | 6    | $so(4)\cong su(2)\oplus su(2)$ |
| Step 2 | $X_1,X_5,X_6,X_7$    | 全部 10 个 | 10   | 完整 $so(5)\cong sp(2)$        |
| Step 3 | $X_1,X_6,X_{10}$     | $+X_9$     | 4    | $u(1)\oplus su(2)$             |

> **只有 Step 3 的 $su(2)$ 子块能写成单位四元数闭式。**

### 符号

| 符号                      | 类型                  | 含义                                                         |
| ------------------------- | --------------------- | ------------------------------------------------------------ |
| $\gamma_1,\dots,\gamma_5$ | Majorana 算符         | 五个 Majorana 费米子，$\{\gamma_i,\gamma_j\}=2\delta_{ij}$   |
| $\Gamma_1,\dots,\Gamma_5$ | $2\times2$ 四元数矩阵 | `Cl(5)` Gamma 矩阵，$\{\Gamma_i,\Gamma_j\}=2\delta_{ij}$     |
| $\Sigma_{ij}$             | $2\times2$ 四元数矩阵 | `so(5)` 生成元 = $\frac14[\Gamma_i,\Gamma_j]$，10 个，张成 $\mathfrak{sp}(2)$ |
| $K(t)$                    | $2\times2$ 四元数矩阵 | 旋量表示的演化生成元，$K(t)=\sum h_{ij}(t)\Sigma_{ij}\in\mathfrak{sp}(2)$ |
| $A(t)$                    | 四元数                | $K$ 的左上 $1\times1$ 块，$A\in\operatorname{Im}\mathbb H$（纯虚四元数） |
| $B(t)$                    | 四元数                | $K$ 的右上 $1\times1$ 块，$B\in\mathbb H$（一般四元数）      |
| $C(t)$                    | 四元数                | $K$ 的左下 $1\times1$ 块，$C\approx -\bar B$                 |
| $D(t)$                    | 四元数                | $K$ 的右下 $1\times1$ 块，$D\in\operatorname{Im}\mathbb H$   |
| $U(t)$                    | $2\times2$ 四元数矩阵 | 旋量演化算符，$\dot U=KU$，$U(0)=I$，$U\in Sp(2)$            |
| $X(t)$                    | 四元数                | $U$ 的左上 $1\times1$ 块                                     |
| $Y(t)$                    | 四元数                | $U$ 的右上 $1\times1$ 块                                     |
| $Z(t)$                    | 四元数                | $U$ 的左下 $1\times1$ 块                                     |
| $W(t)$                    | 四元数                | $U$ 的右下 $1\times1$ 块                                     |
| $q(t)$                    | 四元数                | **Riccati 变量**：$q=ZX^{-1}$，4 实分量，$q(0)=0$            |
| $R(t)$                    | $5\times5$ 实矩阵     | Majorana 演化矩阵，$\gamma_i(t)=R_{ij}\gamma_j(0)$，$R\in SO(5)$ |
| $R_{123}$                 | $3\times3$ 实矩阵     | $R$ 在 $\{\gamma_1,\gamma_2,\gamma_3\}$ 上的子块             |
| $E_1$                     | 实数                  | $\gamma_1$–$\gamma_2$ 杂化能（ABS 特征参数）                 |
| $t_1,t_2,t_3$             | 实数（时变）          | $\gamma_i$ 与 ancilla 的门控耦合强度                         |
| $E_d$                     | 实数（时变）          | 量子点能级                                                   |
| $\tau$                    | 实数                  | 每段 braid 步的时长                                          |

核心关系链：
$$K=\begin{pmatrix}A&B\\C&D\end{pmatrix},\quad
U=\begin{pmatrix}X&Y\\Z&W\end{pmatrix},\quad
\dot U=KU,\quad
q=ZX^{-1},\quad
\dot q=C+Dq-qA-qBq.$$