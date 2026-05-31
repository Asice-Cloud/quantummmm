# 从 five-Majorana full-space 到 `so(5) ≅ sp(2)` 的四元数完整推导

这份文档只做一件事：

- 不先投影、不过早假设子空间；
- 先从论文 PRB111 的 five-Majorana 有效模型出发；
- 严格推出 `so(5)` 的 commutator 闭包；
- 再给出 `sp(2)` 的具体四元数矩阵表示；
- 最后讨论这套表示如何和论文的三段门控协议对应。

换句话说，这是一份 **full-space** 的推导，不是 reduced model。

---

## 1. 起点：论文的 five-Majorana 有效哈密顿量

论文的有效模型写成

$$
H_{EM}(t)
=
 iE_d(t)\gamma_a\gamma_b
+ iE_1\gamma_1\gamma_2
+ i|t_2(t)|\gamma_a\gamma_2
- i|t_1(t)|\gamma_b\gamma_1
- i|t_3(t)|\gamma_a\gamma_3.
$$

这里有五个 Majorana 算符：

$$
\Gamma(t)=
\begin{pmatrix}
\gamma_1(t)\\
\gamma_2(t)\\
\gamma_3(t)\\
\gamma_a(t)\\
\gamma_b(t)
\end{pmatrix},
\qquad
\gamma_i^\dagger=\gamma_i,
\qquad
\{\gamma_i,\gamma_j\}=2\delta_{ij}.
$$

我们的任务不是把它先压成二能级，而是直接研究这五个 Majorana 在完整李代数中的闭包。

---

## 2. Majorana 双线性生成元与 `so(5)` 闭包

定义双线性生成元

$$
X_{ij}:= i\gamma_i\gamma_j,\qquad i<j.
$$

这 5 个 Majorana 一共给出

$$
\binom52=10
$$

个独立双线性生成元。

### 2.1 直接计算 commutator

利用 Majorana 反对易关系，可得

$$
[\gamma_i\gamma_j,\gamma_k\gamma_l]
=2\big(
\delta_{jk}\gamma_i\gamma_l
-\delta_{ik}\gamma_j\gamma_l
-\delta_{jl}\gamma_i\gamma_k
+\delta_{il}\gamma_j\gamma_k
\big).
$$

乘上 `i` 因子后，得到生成元的李代数结构：

$$
[X_{ij},X_{kl}]
=
2i\big(
\delta_{jk}X_{il}
-\delta_{ik}X_{jl}
-\delta_{jl}X_{ik}
+\delta_{il}X_{jk}
\big).
$$

这就是标准的 `so(5)` 结构常数公式。

### 2.2 一个具体基底

为了和论文变量一一对应，取下面这组基底：

$$
X_1=X_{12}=i\gamma_1\gamma_2,
\qquad
X_2=X_{13}=i\gamma_1\gamma_3,
\qquad
X_3=X_{23}=i\gamma_2\gamma_3,
$$

$$
X_4=X_{1a}=i\gamma_1\gamma_a,
\qquad
X_5=X_{a2}=i\gamma_a\gamma_2,
\qquad
X_6=X_{a3}=i\gamma_a\gamma_3,
$$

$$
X_7=X_{1b}=i\gamma_1\gamma_b,
\qquad
X_8=X_{2b}=i\gamma_2\gamma_b,
\qquad
X_9=X_{3b}=i\gamma_3\gamma_b,
$$

$$
X_{10}=X_{ab}=i\gamma_a\gamma_b.
$$

这 10 个生成元构成 `so(5)` 的一组实基底。

---

## 3. 论文哈密顿量在 `so(5)` 基底中的写法

把 `H_{EM}(t)` 直接投到上述基底，得到

$$
H_{EM}(t)
= E_1 X_1 + |t_2(t)|X_5 - |t_1(t)|X_7 - |t_3(t)|X_6 + E_d(t)X_{10}.
$$

这里没有做任何 reduced projection，也没有先假设 `SU(2)` 或 quaternion 子空间。

### 3.1 这意味着什么

这说明论文的控制变量

$$
E_1,\quad t_1(t),\quad t_2(t),\quad t_3(t),\quad E_d(t)
$$

只是 `so(5)` 生成元的时间依赖系数。也就是说，论文协议本质上就是在 `so(5)` 里的一个分段时变线性组合。

---

## 4. full-space 的 Majorana 运动方程

### 4.1 Heisenberg 方程

对于任意 Majorana 分量，有

$$
\dot\gamma_i(t)=i[H_{EM}(t),\gamma_i(t)].
$$

使用

$$
[i\gamma_i\gamma_j,\gamma_k]=2(\delta_{jk}\gamma_i-\delta_{ik}\gamma_j),
$$

可知五个 Majorana 始终闭合在线性空间里，因此可以写成

$$
\dot\Gamma(t)=\mathcal A(t)\Gamma(t),
$$

其中 `\mathcal A(t)` 是一个实反对称的 `5\times 5` 矩阵，故

$$
\mathcal A(t)\in so(5).
$$

### 4.2 具体的 `5\times 5` 反对称矩阵

取基底顺序 `(
\gamma_1,\gamma_2,\gamma_3,\gamma_a,\gamma_b
)`，则上面的 Hamiltonian 对应的矩阵为

$$
\mathcal A(t)=
\begin{pmatrix}
0 & 2E_1 & 0 & 0 & 2|t_1(t)| \\
-2E_1 & 0 & 0 & -2|t_2(t)| & 0 \\
0 & 0 & 0 & 2|t_3(t)| & 0 \\
0 & 2|t_2(t)| & -2|t_3(t)| & 0 & 2E_d(t) \\
-2|t_1(t)| & 0 & 0 & -2E_d(t) & 0
\end{pmatrix}.
$$

于是五个 Majorana 的演化就是

$$
\dot\Gamma(t)=\mathcal A(t)\Gamma(t).
$$

这就是“用 commutator 研究所有 Majorana”的最直接含义：**完全不先投影，只在 full `so(5)` 中做线性闭包。**

### 4.3 这一步和论文三段流程的关系

论文的三段门控只是在不同时间区间给出不同的 `\mathcal A(t)`：

- Step 1：`t_2(t)` 打开，`E_d(t)` 从 `E_0` 压到 0；
- Step 2：`t_3(t)` 打开，`t_2(t)` 关回去；
- Step 3：`t_3(t)` 关回去，`E_d(t)` 恢复到 `E_0`。

所以 full-space 的精确解永远是

$$
\Gamma(t)=\mathcal T\exp\left(\int_0^t\mathcal A(t')dt'\right)\Gamma(0).
$$

---

## 5. `so(5) ≅ sp(2)` 的四元数矩阵表示

现在转到四元数语言。这里要强调：我们不是换了物理体系，而是把同一个 `so(5)` 线性系统写成 `sp(2)` 的四元数矩阵。

### 5.1 四元数记号

记四元数单位为

$$
\mathbf i^2=\mathbf j^2=\mathbf k^2=\mathbf ijk=-1.
$$

设 `\mathbb H` 为四元数数域，`\operatorname{Im}\mathbb H` 表示纯虚四元数。

### 5.2 `sp(2)` 的标准定义

`sp(2)` 可写成所有 `2\times 2` 四元数反厄米矩阵：

$$
\mathfrak{sp}(2)
=
\left\{
\begin{pmatrix}
 u & q \\
 -\bar q & v
\end{pmatrix}
:
 u,v\in\operatorname{Im}\mathbb H,\ q\in\mathbb H
\right\}.
$$

这里 `u,v` 给出 6 个自由度，`q` 给出 4 个自由度，总共 `3+3+4=10` 个自由度，正好和 `so(5)` 一致。

这就是 `so(5)\cong sp(2)` 的具体四元数形式。

### 5.3 一组显式基底

定义

$$
D(u,v):=
\begin{pmatrix}
 u & 0 \\
 0 & v
\end{pmatrix},
\qquad
O(q):=
\begin{pmatrix}
 0 & q \\
 -\bar q & 0
\end{pmatrix},
$$

其中 `u,v\in\operatorname{Im}\mathbb H`，`q\in\mathbb H`。

取一组实基底：

$$
L_a:=\frac12 D(\mathbf e_a,0),
\qquad
R_a:=\frac12 D(0,\mathbf e_a),
\qquad a=1,2,3,
$$

$$
Q_0:=\frac12 O(1),
\qquad
Q_a:=\frac12 O(\mathbf e_a),
\qquad a=1,2,3.
$$

这里 $\mathbf e_1=\mathbf i,\mathbf e_2=\mathbf j,\mathbf e_3=\mathbf k$。

这 10 个基底正好是

- `L_a`：3 个左虚四元数方向；
- `R_a`：3 个右虚四元数方向；
- `Q_0,Q_a`：4 个混合方向。

### 5.4 commutator 关系

对任意纯虚四元数 `u,v` 和任意四元数 `q,q_1,q_2`，有

$$
[D(u,v),D(u',v')]=D([u,u'],[v,v']),
$$

$$
[D(u,v),O(q)]=O(uq-qv),
$$

$$
[O(q_1),O(q_2)]
=
D(q_2\bar q_1-q_1\bar q_2,\;\bar q_2 q_1-\bar q_1 q_2).
$$

因此 `sp(2)` 的 commutator 结构完整闭合，并且与 `so(5)` 同构。

### 5.5 这和 `so(5)` 的关系

`so(5)` 和 `sp(2)` 都是 10 维的实简单李代数，而且它们属于同一个复化类型：

$$
\mathfrak{so}(5,\mathbb C)\cong\mathfrak{sp}(2,\mathbb C).
$$

在实形式上，`so(5)\cong sp(2)`。

所以我们可以把同一个 five-Majorana 代数，一方面看成 `so(5)` 的双线性闭包，另一方面看成 `sp(2)` 的四元数矩阵。

---

## 6. Clifford / spinor 的四元数实现

如果想要更“矩阵化”的版本，可以给出一个 `Cl(5)` 的四元数表示。定义五个 `2\times 2` 四元数矩阵

$$
\Gamma_1=
\begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix},\qquad
\Gamma_2=
\begin{pmatrix}
0 & -\mathbf i \\
\mathbf i & 0
\end{pmatrix},
$$

$$
\Gamma_3=
\begin{pmatrix}
0 & -\mathbf j \\
\mathbf j & 0
\end{pmatrix},\qquad
\Gamma_4=
\begin{pmatrix}
0 & -\mathbf k \\
\mathbf k & 0
\end{pmatrix},
$$

$$
\Gamma_5=
\begin{pmatrix}
1 & 0 \\
0 & -1
\end{pmatrix}.
$$

直接可验证

$$
\{\Gamma_A,\Gamma_B\}=2\delta_{AB}\,\mathbb 1,\qquad A,B=1,\dots,5.
$$

于是双生成元

$$
\Sigma_{AB}:=\frac14[\Gamma_A,\Gamma_B]
$$

给出一个 `so(5)` 的具体四元数 spinor 表示。由于这些矩阵都是四元数反厄米的，它们实际上落在 `sp(2)` 里。

这说明：

- `so(5)` 的双线性闭包；
- `sp(2)` 的四元数矩阵；
- `Cl(5)` 的 spinor 表示；

其实是同一个结构的三种写法。

---

## 7. 从 commutator 角度研究所有 Majorana

### 7.1 operator-space 版本

如果把 Majorana 视为 operator-space 向量 `\Gamma(t)`，那么它的演化完全由

$$
\dot\Gamma(t)=\mathcal A(t)\Gamma(t)
$$

给出。这个矩阵 `\mathcal A(t)` 的条目就是哈密顿量中五个耦合参数的线性组合。

这意味着：

- 五个 Majorana 的行为没有被丢掉；
- 所有信息都在 `so(5)` 的 commutator 中；
- 不需要先把问题压成二能级。

### 7.2 state-space 版本

如果在 `sp(2)` 的 spinor 表示中讨论幺正演化，则有

$$
\dot U(t)=K(t)U(t),
\qquad
K(t)\in\mathfrak{sp}(2).
$$

这里 `K(t)` 就是 `-iH_{EM}(t)` 在四元数矩阵中的表示。

于是同一个问题有两种等价视角：

1. `so(5)` 作用在五个 Majorana 上；
2. `sp(2)` 作用在 spinor / quaternion 表示上。

这就是“从整体出发”研究所有 Majorana 的严格方式。

---

## 8. 和论文三段门控协议如何对应

把论文三段协议直接代入 `K(t)\in sp(2)`，就是让系数函数 `E_d(t), t_1(t), t_2(t), t_3(t), E_1` 在三段中分段变化。

### Step 1

- `t_2(t)` 打开；
- `E_d(t)` 从 `E_0` 降到 0；
- 其余项 `E_1,t_1(t),t_3(t)` 按论文规定保持或关闭。

### Step 2

- `t_3(t)` 打开；
- `t_2(t)` 退去；
- `E_d(t)` 保持在零附近。

### Step 3

- `t_3(t)` 关闭；
- `E_d(t)` 恢复到 `E_0`。

于是整个演化始终是

$$
U(t)=\mathcal T\exp\left(\int_0^t K(t')dt'\right)
\in Sp(2),
$$

或者等价地写成 `SO(5)` 上的路径有序旋转。

---

## 9. 什么时候四元数真的给出解析解

这是最关键的结论。

### 9.1 一般情形

如果 `K(t)` 在一个段内的不同时间不对易，那么不能指望一个单独的四元数闭式。必须保留

$$
\mathcal T\exp\left(\int K dt\right)
$$

或用 Wei–Norman 参数方程。

### 9.2 能真正简化的情形

只有当某一段经过 full-space 推导后，确实被证明闭合到一个真实的 `su(2)` 子代数，并且该段内的旋转方向固定时，才能写成单位四元数：

$$
q(t)=\cos\frac{\Theta(t)}{2}-\hat n\sin\frac{\Theta(t)}{2}.
$$

此时 `\hat n` 是固定方向，`\Theta(t)` 是累计角度。

这时四元数不只是表示，而是**解析闭式**。

### 9.3 对 PRB111 三段协议的判断

- **generic ABS 情形**：第二轮 commutator shell 已经把 `so(5)` 补满，因此整段协议不能直接压成一个单独的单位四元数。
- **`E_1=t_1=0` 的理想情形**：full protocol 降到 `so(4)`，各段还有可能进一步落到不同的 `su(2)` 小块；这时可以写出局部四元数旋转。
- **`E_1=0, t_1\neq 0` 的情形**：仍然不能先验地当成某个二能级旋转，是否简化要继续回到 full-space commutator 结构中验证。

---

## 10. 最终结论

1. **从整体出发是正确的**：先在 full `so(5)` 里用 commutator 研究五个 Majorana 的完整行为。
2. **`so(5) ≅ sp(2)` 是精确的**：它给出一个四元数矩阵的忠实表示，不是约化。
3. **四元数可以写得很具体**：`sp(2)` 就是 `2\times 2` 的四元数反厄米矩阵；`Cl(5)` 的 spinor 也能用四元数矩阵显式表示。
4. **解析解不自动出现**：只有当某一段被 full-space 证明落入真正的 `su(2)` 子块时，四元数才升级成闭式解。
5. **对论文的三段门控协议**：一般仍然是 `Spin(5)` 上的时间有序演化；四元数是精确表示工具，不是默认的求解器。

如果要继续往下做，下一步最自然的是把 `\mathcal A(t)` 的三个分段矩阵 `\mathcal A^{(1)}(t),\mathcal A^{(2)}(t),\mathcal A^{(3)}(t)` 明确写出，然后检查每一段是否真的闭合到某个 `su(2)` 子块。

---

## 11. 论文三段协议的 full-space 分段矩阵

为了和论文的门控步骤一致，记
 
$$
x:=\frac{\pi t}{\tau},
\qquad
f_-(x):=\frac{1-\cos x}{2},
\qquad
f_+(x):=\frac{1+\cos x}{2}.
$$

下面用 `t_c` 表示 `t_2` 和 `t_3` 打开的峰值，用 `t_{1c}` 表示 `G_1` 关联的 `t_1` 峰值。这样写只是为了把门控 schedule 说清楚；如果你想采用更保守的记号，也可以把 `t_{1c}` 继续写回 `t_1`。

### 11.1 Step 1

按论文的门控顺序，Step 1 中 `t_2` 打开、`E_d` 从 `E_0` 压到 0，同时 `G_1` 相关的通道也处于“打开”那一支。对应的分段函数取为

$$
|t_2^{(1)}(t)|=f_-(x)\,t_c,
\qquad
|t_1^{(1)}(t)|=f_-(x)\,t_{1c},
\qquad
|t_3^{(1)}(t)|=0,
\qquad
E_d^{(1)}(t)=f_+(x)\,E_0.
$$

于是 `\dot\Gamma(t)=\mathcal A^{(1)}(t)\Gamma(t)` 中的 `5\times 5` 反对称矩阵为

$$
\mathcal A^{(1)}(t)=
\begin{pmatrix}
0 & 2E_1 & 0 & 0 & 2|t_1^{(1)}(t)| \\
-2E_1 & 0 & 0 & -2|t_2^{(1)}(t)| & 0 \\
0 & 0 & 0 & 0 & 0 \\
0 & 2|t_2^{(1)}(t)| & 0 & 0 & 2E_d^{(1)}(t) \\
-2|t_1^{(1)}(t)| & 0 & 0 & -2E_d^{(1)}(t) & 0
\end{pmatrix}.
$$

### 11.2 Step 2

Step 2 中，`t_3` 接管，`t_2` 退去，`E_d` 处于零附近，同时 `G_1` 处于“关闭”那一支，所以 `t_1` 回到退去的分支。对应分段函数取为

$$
|t_2^{(2)}(t)|=f_+(x)\,t_c,
\qquad
|t_1^{(2)}(t)|=f_+(x)\,t_{1c},
\qquad
|t_3^{(2)}(t)|=f_-(x)\,t_c,
\qquad
E_d^{(2)}(t)=0.
$$

于是

$$
\mathcal A^{(2)}(t)=
\begin{pmatrix}
0 & 2E_1 & 0 & 0 & 2|t_1^{(2)}(t)| \\
-2E_1 & 0 & 0 & -2|t_2^{(2)}(t)| & 0 \\
0 & 0 & 0 & 2|t_3^{(2)}(t)| & 0 \\
0 & 2|t_2^{(2)}(t)| & -2|t_3^{(2)}(t)| & 0 & 0 \\
-2|t_1^{(2)}(t)| & 0 & 0 & 0 & 0
\end{pmatrix}.
$$

### 11.3 Step 3

Step 3 中，`t_3` 再次退去，`E_d` 恢复到 `E_0`。按最直接的论文读法，`G_1` 已经回到关闭支，因此 `t_1` 不再作为独立打开通道保留。对应分段函数取为

$$
|t_2^{(3)}(t)|=0,
\qquad
|t_1^{(3)}(t)|=0,
\qquad
|t_3^{(3)}(t)|=f_+(x)\,t_c,
\qquad
E_d^{(3)}(t)=f_-(x)\,E_0.
$$

于是

$$
\mathcal A^{(3)}(t)=
\begin{pmatrix}
0 & 2E_1 & 0 & 0 & 0 \\
-2E_1 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 2|t_3^{(3)}(t)| & 0 \\
0 & 0 & -2|t_3^{(3)}(t)| & 0 & 2E_d^{(3)}(t) \\
0 & 0 & 0 & -2E_d^{(3)}(t) & 0
\end{pmatrix}.
$$

### 11.4 这三段矩阵说明了什么

这三段的共同点是：它们都严格作用在五个 Majorana 生成的 full-space 上，没有预设任何 reduced qubit。

不同点是：

- Step 1 里 `\gamma_1,\gamma_2,\gamma_a,\gamma_b` 的链路最完整；
- Step 2 里 `\gamma_3` 进入主链路；
- Step 3 里 `\gamma_3` 和 `\gamma_a,\gamma_b` 的恢复项主导返回过程。

所以，如果要从解析角度推进，最自然的下一步不是直接把它压成单位四元数，而是先检查每一段矩阵是否进一步掉入某个真实的 `su(2)` 子块。

---

## 12. 对应的 `Sp(2)` 四元数矩阵写法

现在把上面的 `\mathcal A^{(s)}(t)` 通过固定的 `so(5)\cong sp(2)` 同构，写成四元数矩阵生成元。

### 12.1 一般形式

对每一段 `s=1,2,3`，记对应的 `sp(2)` 生成元为

$$
K^{(s)}(t):=-iH^{(s)}(t)\in\mathfrak{sp}(2).
$$

因此可以写成

$$
K^{(s)}(t)=
\begin{pmatrix}
u^{(s)}(t) & q^{(s)}(t) \\
-\overline{q^{(s)}(t)} & v^{(s)}(t)
\end{pmatrix},
$$

其中 `u^{(s)}(t),v^{(s)}(t)\in\operatorname{Im}\mathbb H`，`q^{(s)}(t)\in\mathbb H`。

### 12.2 用基底展开

在上一节定义的 `sp(2)` 基底 `\{L_a,R_a,Q_0,Q_a\}` 下，可写成

$$
K^{(s)}(t)
=
\sum_{a=1}^3 u_a^{(s)}(t)L_a
+ \sum_{a=1}^3 v_a^{(s)}(t)R_a
+ q_0^{(s)}(t)Q_0
+ \sum_{a=1}^3 q_a^{(s)}(t)Q_a.
$$

这里的 10 个实系数 `u_a^{(s)},v_a^{(s)},q_0^{(s)},q_a^{(s)}` 都是 `E_1, t_1^{(s)}(t), t_2^{(s)}(t), t_3^{(s)}(t), E_d^{(s)}(t)` 的线性组合；一旦你选定上面 `so(5)\cong sp(2)` 的固定同构，这种线性关系就唯一确定。

### 12.3 三段门控的 `Sp(2)` 读法

于是论文三段协议在 `Sp(2)` 里就是

$$
U^{(s)}(t)=\mathcal T\exp\left(\int_0^t K^{(s)}(t')dt'\right),
\qquad s=1,2,3.
$$

全程演化就是这三段的顺序乘积。

### 12.4 解析化的判据

这一步的意义在于：

- 如果某一段的 `K^{(s)}(t)` 在固定同构下进一步缩成一个单一的 `su(2)` 生成方向，那么它可以写成单位四元数旋转；
- 如果不行，那就必须保留 `Sp(2)` 的时间有序指数。

所以这份 full-space 推导现在已经把“论文的三段含时演化”精确地写成了 `SO(5)` / `Sp(2)` 的矩阵问题；下一步若要做解析简化，就只剩下逐段检查 commutator 闭包这一件事。

## 13. 三段矩阵的闭包检查

下面直接用第 2 节的 `so(5)` 基底

$$
X_1=X_{12},\quad X_2=X_{13},\quad X_3=X_{23},\quad X_4=X_{1a},\quad X_5=X_{a2},\quad X_6=X_{a3},\quad X_7=X_{1b},\quad X_8=X_{2b},\quad X_9=X_{3b},\quad X_{10}=X_{ab}.
$$

### 13.1 Step 1

Step 1 的 active generators 是

$$
\{X_1,X_5,X_7,X_{10}\}=
\{X_{12},X_{a2},X_{1b},X_{ab}\}.
$$

由第 2.1 节的 commutator 公式可见，闭包会生成

$$
[X_{12},X_{a2}]\propto X_{1a},\qquad [X_{12},X_{1b}]\propto X_{2b},
$$

$$
[X_{1a},X_{1b}]\propto X_{ab},\qquad [X_{1a},X_{a2}]\propto X_{12},
$$

$$
[X_{2b},X_{a2}]\propto X_{ab},\qquad [X_{2b},X_{1b}]\propto X_{12},
$$

以及

$$
[X_{1a},X_{ab}]\propto X_{1b},\qquad [X_{2b},X_{ab}]\propto X_{a2}.
$$

因此 Step 1 的最小闭包是六维子代数

$$
\operatorname{span}\{X_{12},X_{1a},X_{a2},X_{1b},X_{2b},X_{ab}\}\cong so(4)\cong su(2)\oplus su(2).
$$

这说明 Step 1 不是单一的 `su(2)` 块，所以一般不能直接压成一个单位四元数旋转。

### 13.2 Step 2

Step 2 的 active generators 是

$$
\{X_1,X_5,X_7,X_6\}=
\{X_{12},X_{a2},X_{1b},X_{a3}\}.
$$

此时 commutator 进一步给出

$$
[X_{12},X_{a2}]\propto X_{1a},\qquad [X_{12},X_{1b}]\propto X_{2b},
$$

$$
[X_{a2},X_{a3}]\propto X_{23},\qquad [X_{1a},X_{a3}]\propto X_{13},
$$

$$
[X_{23},X_{2b}]\propto X_{3b},\qquad [X_{1a},X_{1b}]\propto X_{ab}.
$$

连同原来的四个生成元，这已经把十个基底全部生成出来了，因此 Step 2 的闭包是整个

$$
so(5).
$$

所以 Step 2 绝不是 `su(2)` 子块，而是完整的 full-space 演化。

### 13.3 Step 3

Step 3 的 active generators 是

$$
\{X_1,X_6,X_{10}\}=
\{X_{12},X_{a3},X_{ab}\}.
$$

这里有一个关键点：

$$
[X_{12},X_{a3}]=[X_{12},X_{ab}]=0.
$$

而时间依赖的两项满足

$$
[X_{a3},X_{ab}]\propto X_{3b},\qquad
[X_{ab},X_{3b}]\propto X_{a3},\qquad
[X_{3b},X_{a3}]\propto X_{ab}.
$$

因此

$$
\operatorname{span}\{X_{a3},X_{3b},X_{ab}\}\cong su(2),
$$

而 `X_{12}` 只是一个与其对易的 `u(1)` 因子。也就是说，Step 3 的结构是

$$
su(2)\oplus u(1).
$$

这一步可以先把 `X_{12}` 的相位因子单独分离出来，剩余部分再用 `SU(2)` / 单位四元数参数化。

### 13.4 结论

三段协议的层级是：

- Step 1：`so(4)\cong su(2)\oplus su(2)`，不是单一 `su(2)`；
- Step 2：完整 `so(5)`，不能降成 `su(2)`；
- Step 3：`su(2)\oplus u(1)`，其中时间依赖的非平凡部分确实落在 `su(2)` 里。

所以，如果要谈单位四元数，唯一严格成立的入口是 Step 3，且还要先把那个与 `su(2)` 对易的 `u(1)` 相位因子剥离掉。

## 14. 实际解析路线

如果目标是尽可能得到解析解，而不是先验压缩模型，那么最稳妥的做法就是把 `Sp(2)` 作为主框架，单位四元数只作为满足闭包条件时的局部简化。

### 14.1 先保留 `Sp(2)` 的精确演化

全程先写成

$$
U(t)=\mathcal T\exp\left(\int_0^t K(t')dt'\right),
\qquad K(t)\in\mathfrak{sp}(2).
$$

这里不提前假设存在全局 `SU(2)` 结构，也不把三段协议强行压成单一单位四元数旋转。

### 14.2 再逐段做闭包检查

对每一段 `K^{(s)}(t)`，先检查它在固定同构下生成的最小李代数是什么：

- 若它真的闭到 `su(2)`，就可以写成单位四元数闭式；
- 若它只是 `su(2)\oplus u(1)`，则先把可对易的 `u(1)` 相位因子剥离，再处理 `su(2)` 部分；
- 若它仍然是更大的 `so(4)` 或完整 `so(5)`，就不要强行做单位四元数化。

### 14.3 不闭合的段继续留在 `Sp(2)`

对 Step 1 和 Step 2，这种做法尤其必要，因为它们并不天然落在单一 `su(2)` 上。此时可以继续保留 `Sp(2)` 的时间有序指数；如果后续要进一步解析化，再考虑 Wei–Norman 分解或 Magnus 展开。

### 14.4 最终策略

所以当前最合理的路线是：

1. 主体上保留 `Sp(2)` 的精确表示；
2. 分段检查闭包，只在真正的 `su(2)` 段上写单位四元数；
3. 对其余段保持 `Sp(2)`，必要时再用 Wei–Norman 或 Magnus 做进一步分析。
