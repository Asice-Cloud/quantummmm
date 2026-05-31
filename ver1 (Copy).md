# PRB111 (Phys. Rev. B 111, 205411)

# Majorana Effective Hamiltonian
# → so(5) Lie Closure
# → Wei–Norman Preparation

---

# 文档目的

本笔记目标：

严格从 PRB111 的 effective Majorana Hamiltonian 出发，

建立：

\[
H(t)
\rightarrow
\mathfrak{so}(5)
\rightarrow
\mathrm{Spin}(5)
\rightarrow
\text{Wei–Norman}
\]

的完整数学框架。

---

# 重要说明

本笔记仅保留已经严格推导的内容。

以下内容尚未使用：

- SU(2)投影
- Quaternion表示
- Hopf纤维化
- Berry连接
- Adiabatic近似

因为这些内容均未从论文直接推出。

---

# Part 1
# PRB111 Effective Hamiltonian

论文有效模型：

\[
H_{EM}(t)
=
iE_d(t)\gamma_a\gamma_b
+
iE_1\gamma_1\gamma_2
+
i t_2(t)\gamma_a\gamma_2
-
i t_1(t)\gamma_b\gamma_1
-
i t_3(t)\gamma_a\gamma_3.
\]

Majorana满足：

\[
\{\gamma_i,\gamma_j\}
=
2\delta_{ij}.
\]

---

# Part 2
# Majorana集合

定义顺序：

\[
(1,2,3,a,b).
\]

因此系统包含：

\[
\{
\gamma_1,
\gamma_2,
\gamma_3,
\gamma_a,
\gamma_b
\}.
\]

共有：

\[
N=5
\]

个Majorana。

---

# Part 3
# Majorana双线性生成元

定义：

\[
L_{ij}
=
i\gamma_i\gamma_j.
\]

满足：

\[
L_{ij}
=
-L_{ji}.
\]

因此独立生成元数目：

\[
\binom52=10.
\]

---

# Part 4
# so(5)基底

选取：

\[
B_1=L_{12}
\]

\[
B_2=L_{13}
\]

\[
B_3=L_{23}
\]

\[
B_4=L_{1a}
\]

\[
B_5=L_{2a}
\]

\[
B_6=L_{3a}
\]

\[
B_7=L_{1b}
\]

\[
B_8=L_{2b}
\]

\[
B_9=L_{3b}
\]

\[
B_{10}=L_{ab}.
\]

---

# Part 5
# Majorana Lie代数

从

\[
\{\gamma_i,\gamma_j\}=2\delta_{ij}
\]

严格推导得到：

\[
[L_{ij},L_{kl}]
=
2i
(
\delta_{jk}L_{il}
-\delta_{ik}L_{jl}
-\delta_{jl}L_{ik}
+\delta_{il}L_{jk}
).
\]

这是标准：

\[
\mathfrak{so}(5)
\]

结构常数公式。

---

# Part 6
# PRB111 Hamiltonian投影

根据基底定义：

\[
L_{a2}
=
-L_{2a}
=
-B_5
\]

\[
L_{a3}
=
-L_{3a}
=
-B_6
\]

\[
L_{b1}
=
-L_{1b}
=
-B_7
\]

\[
L_{ab}
=
B_{10}.
\]

因此：

\[
H(t)
=
E_1B_1
+
t_2(t)B_5
-
t_3(t)B_6
-
t_1(t)B_7
+
E_d(t)B_{10}.
\]

---

# Part 7
# 驱动向量 h(t)

写成：

\[
H(t)
=
\sum_{k=1}^{10}
h_k(t)B_k.
\]

得到：

\[
h(t)
=
(
E_1,
0,
0,
0,
t_2(t),
-t_3(t),
-t_1(t),
0,
0,
E_d(t)
).
\]

---

# Part 8
# Lie Closure

Hamiltonian初始生成元集合：

\[
S_0=
\{
B_1,
B_5,
B_6,
B_7,
B_{10}
\}.
\]

---

计算：

\[
[B_1,B_5]
=
2iB_4
\]

产生：

\[
B_4.
\]

---

计算：

\[
[B_1,B_7]
=
-2iB_8
\]

产生：

\[
B_8.
\]

---

计算：

\[
[B_6,B_4]
=
-2iB_2
\]

产生：

\[
B_2.
\]

---

计算：

\[
[B_2,B_1]
=
2iB_3
\]

产生：

\[
B_3.
\]

---

计算：

\[
[B_3,B_8]
=
-2iB_9
\]

产生：

\[
B_9.
\]

---

最终：

\[
\{
B_1,\ldots,B_{10}
\}
\]

全部出现。

因此：

\[
\operatorname{Lie}(S_0)
=
\mathfrak{so}(5).
\]

---

# Part 9
# 已完成的结构常数

so(3)部分：

\[
[B_1,B_2]
=
-2iB_3
\]

\[
[B_1,B_3]
=
2iB_2
\]

\[
[B_2,B_3]
=
-2iB_1.
\]

---

a-sector：

\[
[B_1,B_4]
=
-2iB_5
\]

\[
[B_1,B_5]
=
2iB_4
\]

\[
[B_2,B_4]
=
-2iB_6
\]

\[
[B_2,B_6]
=
2iB_4
\]

\[
[B_3,B_5]
=
-2iB_6
\]

\[
[B_3,B_6]
=
2iB_5.
\]

---

b-sector：

\[
[B_1,B_7]
=
-2iB_8
\]

\[
[B_1,B_8]
=
2iB_7
\]

\[
[B_2,B_7]
=
-2iB_9
\]

\[
[B_2,B_9]
=
2iB_7
\]

\[
[B_3,B_8]
=
-2iB_9
\]

\[
[B_3,B_9]
=
2iB_8.
\]

---

与 B10 的耦合：

\[
[B_{10},B_4]
=
2iB_7
\]

\[
[B_{10},B_5]
=
2iB_8
\]

\[
[B_{10},B_6]
=
2iB_9
\]

\[
[B_{10},B_7]
=
-2iB_4
\]

\[
[B_{10},B_8]
=
-2iB_5
\]

\[
[B_{10},B_9]
=
-2iB_6.
\]

---

# Part 10
# 论文Braiding Protocol

定义：

\[
s(x)
=
\frac{1-\cos(\pi x/\tau)}2.
\]

---

Step 1

\[
0<t<\tau
\]

\[
t_2=s\,t_c
\]

\[
t_1=s\,t_1
\]

\[
t_3=0
\]

\[
E_d=(1-s)E_0.
\]

---

Step 2

\[
\tau<t<2\tau
\]

\[
t_2=(1-s)t_c
\]

\[
t_1=(1-s)t_1
\]

\[
t_3=s\,t_c
\]

\[
E_d=0.
\]

---

Step 3

\[
2\tau<t<3\tau
\]

\[
t_2=0
\]

\[
t_1=0
\]

\[
t_3=(1-s)t_c
\]

\[
E_d=sE_0.
\]

---

一次交换完成。

论文研究：

\[
6\tau
\]

对应两次连续交换。

---

# Part 11
# Wei–Norman准备

选取顺序：

\[
U
=
e^{x_1B_1}
e^{x_2B_2}
e^{x_3B_3}
e^{x_4B_4}
e^{x_5B_5}
e^{x_6B_6}
e^{x_7B_7}
e^{x_8B_8}
e^{x_9B_9}
e^{x_{10}B_{10}}.
\]

---

满足：

\[
\dot U
=
H(t)U.
\]

---

需要计算：

\[
\operatorname{ad}_{B_i}
\]

以及：

\[
e^{x_i\operatorname{ad}_{B_i}}.
\]

---

# 当前状态

已完成：

✔ Hamiltonian投影

✔ Lie closure

✔ so(5)基底

✔ 大部分结构常数

下一步：

1. 完整结构常数表
2. 全部 ad(B_i) 矩阵
3. Wei–Norman矩阵 M(x)
4. 10维非线性ODE
5. Spin(5)/Quaternion表示