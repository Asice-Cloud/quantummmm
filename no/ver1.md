# PRB111 (Phys. Rev. B 111, 205411, 2025)
# Effective Majorana Hamiltonian 的 Wei–Norman 重写与解析解笔记

---

# 1. 研究目标

论文研究：

> ABS（Andreev Bound States）与有限重叠 Majorana Zero Modes (MZMs)
> 的非阿贝尔编织性质。

核心问题：

有限重叠会产生额外耦合

\[
E_1
\]

与

\[
t_1
\]

从而破坏理想 Majorana 编织。

---

# 2. 论文 Effective Hamiltonian

论文构造了如下 Effective Majorana 模型：

\[
H_{EM}(t)
=
iE_d \gamma_a \gamma_b
+
iE_1 \gamma_1 \gamma_2
+
i|t_2(t)|\gamma_a\gamma_2
-
i|t_1(t)|\gamma_b\gamma_1
-
i|t_3(t)|\gamma_a\gamma_3
\]

其中：

\[
\gamma_i^\dagger=\gamma_i
\]

为 Majorana 算符。

---

## 2.1 Majorana 集合

系统包含：

\[
\{
\gamma_1,
\gamma_2,
\gamma_3,
\gamma_a,
\gamma_b
\}
\]

共 5 个 Majorana。

---

## 2.2 各项物理意义

### QD 能级

\[
iE_d\gamma_a\gamma_b
\]

表示量子点内部耦合。

---

### 有限长度导致的 Majorana 重叠

\[
iE_1\gamma_1\gamma_2
\]

对应 Majorana splitting。

这是论文最重要参数之一。

---

### 编织控制项

\[
i t_2 \gamma_a\gamma_2
\]

\[
-i t_3 \gamma_a\gamma_3
\]

控制交换路径。

---

### 有限重叠产生的额外耦合

\[
-i t_1 \gamma_b\gamma_1
\]

这是导致非理想编织的重要来源。

---

# 3. Majorana 双线性代数

定义：

\[
X_{ij}
=
i\gamma_i\gamma_j
\]

共有：

\[
\binom{5}{2}=10
\]

个独立生成元。

因此：

\[
\mathfrak{so}(5)
\]

自然出现。

---

# 4. so(5) 基底选择

定义：

\[
X_1=i\gamma_1\gamma_2
\]

\[
X_2=i\gamma_2\gamma_3
\]

\[
X_3=i\gamma_3\gamma_1
\]

---

\[
X_4=i\gamma_a\gamma_1
\]

\[
X_5=i\gamma_a\gamma_2
\]

\[
X_6=i\gamma_a\gamma_3
\]

\[
X_7=i\gamma_b\gamma_1
\]

---

\[
X_8=i\gamma_a\gamma_b
\]

\[
X_9=i\gamma_b\gamma_2
\]

\[
X_{10}=i\gamma_b\gamma_3
\]

---

# 5. Hamiltonian 投影到 so(5)

写成：

\[
H(t)
=
\sum_{k=1}^{10}
h_k(t)X_k
\]

逐项比较：

---

## E1 项

\[
iE_1\gamma_1\gamma_2
=
E_1 X_1
\]

因此

\[
h_1=E_1
\]

---

## t2 项

\[
it_2\gamma_a\gamma_2
=
t_2 X_5
\]

因此

\[
h_5=t_2
\]

---

## t3 项

\[
-it_3\gamma_a\gamma_3
=
-t_3 X_6
\]

因此

\[
h_6=-t_3
\]

---

## t1 项

\[
-it_1\gamma_b\gamma_1
=
-t_1 X_7
\]

因此

\[
h_7=-t_1
\]

---

## Ed 项

\[
iE_d\gamma_a\gamma_b
=
E_d X_8
\]

因此

\[
h_8=E_d
\]

---

# 6. 最终 h_k

唯一非零：

\[
h_1=E_1
\]

\[
h_5=t_2
\]

\[
h_6=-t_3
\]

\[
h_7=-t_1
\]

\[
h_8=E_d
\]

其余：

\[
h_2=h_3=h_4=h_9=h_{10}=0
\]

---

# 7. Wei–Norman 分解

构造：

\[
U(t)
=
e^{x_1X_1}
e^{x_2X_2}
e^{x_3X_3}
e^{x_4X_4}
e^{x_5X_5}
e^{x_6X_6}
e^{x_7X_7}
e^{x_8X_8}
e^{x_9X_9}
e^{x_{10}X_{10}}
\]

---

满足：

\[
\dot U
=
H U
\]

---

# 8. Wei–Norman 方程

定义：

\[
\theta=
(x_1,\ldots,x_{10})
\]

则：

\[
\dot\theta_i
=
h_i
+
\sum_{j<i}
C_{ij}(\theta)h_j
\]

其中：

\[
C_{ij}
\]

由

\[
[X_i,X_j]
\]

决定。

---

# 9. SU(2) 编织子代数

观察：

\[
X_1=i\gamma_1\gamma_2
\]

\[
X_2=i\gamma_2\gamma_3
\]

\[
X_3=i\gamma_3\gamma_1
\]

满足：

\[
[X_1,X_2]
=
2X_3
\]

\[
[X_2,X_3]
=
2X_1
\]

\[
[X_3,X_1]
=
2X_2
\]

即：

\[
\mathfrak{su}(2)
\]

---

定义：

\[
J_x=X_2
\]

\[
J_y=X_3
\]

\[
J_z=X_1
\]

---

# 10. 理想 Majorana 编织

当：

\[
E_1=0
\]

\[
t_1=0
\]

时：

系统只剩交换路径。

对应：

\[
U_{braid}
=
e^{-i\frac{\pi}{2}J_x}
\]

---

这是理想 non-Abelian braid。

---

# 11. E1 导致的动态相位

加入：

\[
H_{E_1}
=
E_1J_z
\]

得到：

\[
U_{E_1}
=
\exp
\left(
-iJ_z
\int E_1(t)dt
\right)
\]

定义：

\[
\phi_E
=
\int E_1(t)dt
\]

---

# 12. t1 导致的附加旋转

将：

\[
X_7=i\gamma_b\gamma_1
\]

投影到逻辑比特空间后，

可视为额外旋转轴。

记作：

\[
J_y
\]

方向扰动。

于是：

\[
H_{t_1}
=
t_1J_y
\]

---

# 13. 有效 SU(2) 哈密顿量

得到：

\[
H_{eff}
=
E_1J_z
+
t_1J_y
\]

---

写成向量形式：

\[
H_{eff}
=
\mathbf B\cdot\mathbf J
\]

其中：

\[
\mathbf B
=
(0,t_1,E_1)
\]

---

# 14. 有效场强

定义：

\[
\Omega
=
\sqrt{E_1^2+t_1^2}
\]

---

单位方向：

\[
\hat n
=
\frac{(0,t_1,E_1)}
{\Omega}
\]

---

# 15. 解析时间演化

若 E1 与 t1 近似常数：

\[
U(t)
=
\exp
\left[
-i
\Omega t
\,
(\hat n\cdot J)
\right]
\]

---

# 16. 总编织算符

因此：

\[
U
=
\exp
\left(
-i\frac{\pi}{2}J_x
\right)
\,
\exp
\left[
-i
\int
(E_1J_z+t_1J_y)
dt
\right]
\]

---

# 17. 动态相位

定义：

\[
\Theta_{dyn}
=
\int_0^{6\tau}
\sqrt{E_1^2(t)+t_1^2(t)}
\,dt
\]

---

# 18. Fidelity 振荡

编织保真度：

\[
F(\tau)
=
|
\langle\psi_f|
U(\tau)
|\psi_i\rangle
|^2
\]

---

近似得到：

\[
F(\tau)
\sim
\sin^2
\Theta_{dyn}(\tau)
\]

---

因此振荡周期：

\[
T
=
\frac{\pi}
{\sqrt{E_1^2+t_1^2}}
\]

---

# 19. PRB111 核心物理结论

有限重叠 Majorana：

\[
E_1\neq0
\]

导致动态相位。

---

额外耦合：

\[
t_1\neq0
\]

导致旋转轴偏转。

---

二者共同产生：

\[
\Theta_{dyn}
=
\int
\sqrt{E_1^2+t_1^2}
dt
\]

---

从而导致：

\[
F(\tau)
\]

随编织时间振荡。

---

# 20. 最终总结

论文中的 ABS 与有限重叠 MZM

均可理解为：

\[
SU(2)
\]

理想编织

+

\[
E_1
\]

动态相位

+

\[
t_1
\]

附加旋转

的组合系统。

其解析时间演化可写为：

\[
\boxed{
U
=
e^{-i\frac{\pi}{2}J_x}
e^{-i\int(E_1J_z+t_1J_y)dt}
}
\]

而编织偏离理想结果的根源为：

\[
\boxed{
\Theta_{dyn}
=
\int
\sqrt{E_1^2+t_1^2}
dt
}
\]

这正是 PRB111 数值结果中振荡行为的来源。