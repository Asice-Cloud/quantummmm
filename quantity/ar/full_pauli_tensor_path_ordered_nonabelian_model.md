# 完整 Pauli 张量积路径演化 non-Abelian 可计算模型



---



\[
\mathrm{End}(V\otimes V)
\]

的标准基为：

\[
\boxed{
\{\sigma^\alpha\otimes\sigma^\beta\}
}
\]

其中：

\[
\alpha,\beta\in\{0,x,y,z\}
\]

因此任意二元Hamiltonian都可写为：

\[
\boxed{
h(u)
=
\sum_{\alpha,\beta}
h_{\alpha\beta}(u)
\sigma^\alpha\otimes\sigma^\beta
}
\]

其中：

- \(u\)：路径参数 / braid parameter
- \(h_{\alpha\beta}(u)\)：路径依赖耦合

---

# 2. 路径演化算符

定义局域路径演化：

\[
\boxed{
R(u)
=
\mathcal T
\exp\left(
-i\int_0^u h(s)ds
\right)
}
\]

其中：

\[
\mathcal T
\]

表示路径序排序。

由于：

\[
[h(s_1),h(s_2)]\neq0
\]

因此必须保留路径序结构。

---

# 3. Dyson展开

定义：

\[
R(u)=I+R^{(1)}+R^{(2)}+R^{(3)}+\cdots
\]

---

# 4. 一阶项

由Dyson展开：

\[
\boxed{
R^{(1)}
=
-i\int_0^u ds_1\,h(s_1)
}
\]

代入Pauli展开：

\[
R^{(1)}
=
-i\int_0^u ds_1
\sum_{\alpha,\beta}
h_{\alpha\beta}(s_1)
\sigma^\alpha\otimes\sigma^\beta
\]

交换求和与积分：

\[
\boxed{
R^{(1)}
=
- i
\sum_{\alpha,\beta}
\left(
\int_0^u h_{\alpha\beta}(s_1)ds_1
\right)
\sigma^\alpha\otimes\sigma^\beta
}
\]

---

# 5. 二阶项

Dyson展开二阶项：

\[
\boxed{
R^{(2)}
=
(-i)^2
\int_0^u ds_1
\int_0^{s_1} ds_2
h(s_1)h(s_2)
}
\]

代入Pauli展开：

\[
h(s_1)
=
\sum_{\alpha,\beta}
h_{\alpha\beta}(s_1)
\sigma^\alpha\otimes\sigma^\beta
\]

\[
h(s_2)
=
\sum_{\mu,\nu}
h_{\mu\nu}(s_2)
\sigma^\mu\otimes\sigma^\nu
\]

因此：

\[
R^{(2)}
=
(-i)^2
\int ds_1ds_2
\sum_{\alpha,\beta}
\sum_{\mu,\nu}

h_{\alpha\beta}(s_1)
h_{\mu\nu}(s_2)

(\sigma^\alpha\sigma^\mu)
\otimes
(\sigma^\beta\sigma^\nu)
\]

于是：

\[
\boxed{
R^{(2)}
=
-
\sum_{\alpha,\beta,\mu,\nu}
\left(
\int_0^u ds_1
\int_0^{s_1}ds_2
h_{\alpha\beta}(s_1)
h_{\mu\nu}(s_2)
\right)

(\sigma^\alpha\sigma^\mu)
\otimes
(\sigma^\beta\sigma^\nu)
}
\]

---

# 6. 三阶项

Dyson展开：

\[
\boxed{
R^{(3)}
=
(-i)^3
\int_0^u ds_1
\int_0^{s_1} ds_2
\int_0^{s_2} ds_3

h(s_1)h(s_2)h(s_3)
}
\]

代入Pauli展开：

\[
R^{(3)}
=
(-i)^3
\int ds_1ds_2ds_3

\sum_{a,b,c}

h_a(s_1)
h_b(s_2)
h_c(s_3)

P_aP_bP_c
\]

其中：

\[
P_a=\sigma^{\alpha_a}\otimes\sigma^{\beta_a}
\]

---

# 7. 三体嵌入

定义三体空间：

\[
V^{\otimes3}
\]

定义局域嵌入Hamiltonian：

\[
\boxed{
h_{12}(u)
=
\sum_{\alpha,\beta}
h_{\alpha\beta}(u)
\sigma^\alpha\otimes\sigma^\beta\otimes I
}
\]

\[
\boxed{
h_{23}(u)
=
\sum_{\mu,\nu}
h_{\mu\nu}(u)
I\otimes\sigma^\mu\otimes\sigma^\nu
}
\]

对应演化算符：

\[
\boxed{
R_{12}(u)
=
\mathcal T
\exp\left(
-i\int_0^u h_{12}(s)ds
\right)
}
\]

\[
\boxed{
R_{23}(u)
=
\mathcal T
\exp\left(
-i\int_0^u h_{23}(s)ds
\right)
}
\]

---

# 8. Yang–Baxter deviation

定义：

\[
\boxed{
\Delta(u,v)
=
R_{12}(u)
R_{23}(u+v)
R_{12}(v)
-
R_{23}(v)
R_{12}(u+v)
R_{23}(u)
}
\]

若：

\[
\Delta=0
\]

则Yang–Baxter关系成立。

---

# 9. 最低阶展开

将Dyson展开代入：

\[
R_{12}=I+A_{12}+B_{12}+\cdots
\]

\[
R_{23}=I+A_{23}+B_{23}+\cdots
\]

其中：

\[
A_{12}
=
-i\int h_{12}(s)ds
\]

\[
A_{23}
=
-i\int h_{23}(s)ds
\]

---

# 10. 一阶项

YBE deviation一阶项：

\[
A_{12}+A_{23}+A_{12}
-
(A_{23}+A_{12}+A_{23})
\]

完全抵消。

因此：

\[
\boxed{
\Delta^{(1)}=0
}
\]

---

# 11. 二阶项

二阶项包括：

\[
A_{12}A_{23}
\]

\[
A_{23}A_{12}
\]

等。

全部整理后仍完全抵消。

因此：

\[
\boxed{
\Delta^{(2)}=0
}
\]

---

# 12. 三阶项

最低非平凡项来自：

\[
A_{12}A_{23}A_{12}
-
A_{23}A_{12}A_{23}
\]

整理得到：

\[
\boxed{
\Delta
\sim
\int ds_1ds_2
[h_{12}(s_1),h_{23}(s_2)]
}
\]

这是non-Abelian性的最低阶来源。

---

# 13. 计算交换子

代入：

\[
h_{12}(s_1)
=
\sum_{\alpha,\beta}
h_{\alpha\beta}(s_1)
\sigma^\alpha\otimes\sigma^\beta\otimes I
\]

\[
h_{23}(s_2)
=
\sum_{\mu,\nu}
h_{\mu\nu}(s_2)
I\otimes\sigma^\mu\otimes\sigma^\nu
\]

于是：

\[
[h_{12}(s_1),h_{23}(s_2)]
=
\sum
h_{\alpha\beta}(s_1)
h_{\mu\nu}(s_2)
[
\sigma^\alpha\otimes\sigma^\beta\otimes I,
I\otimes\sigma^\mu\otimes\sigma^\nu
]
\]

由于只有第二个site重叠：

\[
\boxed{
[h_{12}(s_1),h_{23}(s_2)]
=
\sum
h_{\alpha\beta}(s_1)
h_{\mu\nu}(s_2)
\sigma^\alpha
\otimes
[\sigma^\beta,\sigma^\mu]
\otimes
\sigma^\nu
}
\]

---

# 14. Pauli交换关系

使用：

\[
[\sigma^a,\sigma^b]
=
2i\epsilon_{abc}\sigma^c
\]

因此：

\[
\boxed{
[h_{12}(s_1),h_{23}(s_2)]
=
2i
\sum
h_{\alpha\beta}(s_1)
h_{\mu\nu}(s_2)
\epsilon_{\beta\mu\gamma}
\sigma^\alpha
\otimes
\sigma^\gamma
\otimes
\sigma^\nu
}
\]

于是：

\[
\boxed{
\Delta
\sim
2i
\int ds_1ds_2
\sum
h_{\alpha\beta}(s_1)
h_{\mu\nu}(s_2)
\epsilon_{\beta\mu\gamma}
\sigma^\alpha
\otimes
\sigma^\gamma
\otimes
\sigma^\nu
}
\]

---

# 15. non-Abelian量化

定义Frobenius non-Abelian measure：

\[
\boxed{
\mathcal N
=
\sqrt{
\mathrm{Tr}(\Delta^\dagger\Delta)
}
}
\]

利用Pauli正交关系：

\[
\boxed{
\mathrm{Tr}(P_aP_b)
=
2^N\delta_{ab}
}
\]

得到：

\[
\boxed{
\mathcal N
\sim
\sqrt{
\int ds_1ds_2
\sum
|h_{\alpha\beta}(s_1)|^2
|h_{\mu\nu}(s_2)|^2
\epsilon_{\beta\mu\gamma}^2
}
}
\]

1) Dyson 第三阶的精确形式（保留有序积分常数）：
\[
R^{(3)}(u)=(-i)^3\int_0^u ds_1\int_0^{s_1} ds_2\int_0^{s_2} ds_3\;h(s_1)h(s_2)h(s_3).
\]

2) YBE 偏差的最低非平凡阶为三阶，精确三阶差为三重时间有序积分的差：
\[
\begin{aligned}
\Delta^{(3)}={}&(-i)^3\Big[\int_0^u ds_1\int_0^{u+v} dt\int_0^v ds_3\\
&\quad\;h_{12}(s_1)h_{23}(t)h_{12}(s_3)
-\int_0^v ds_1'\int_0^{u+v} dt'\int_0^u ds_3'\h_{23}(s_1')h_{12}(t')h_{23}(s_3')\Big].
\end{aligned}
\]

3) 把包含重疊站的項重寫為交換子形式时，会出现主导的交换子项
\[
[h_{12}(s),h_{23}(t)]=2i\sum_{\alpha,\beta,\mu,\nu,\gamma}h_{\alpha\beta}(s)h_{\mu\nu}(t)\epsilon_{\beta\mu\gamma}\;\sigma^\alpha\otimes\sigma^\gamma\otimes\sigma^\nu.
\]
代入前因子 $(-i)^3$ 与上述 $2i$，主导代数常数为
\[
(-i)^3\cdot 2i = -2.
\]
因此交换子贡献在 Δ 中带有明确的数值常数 $-2$，时间依赖由由三重有序积分合併後得到的核 $K(s,t)$（由 θ 函数組成）決定。

4) 对 Frobenius 范数的更精确缩放：主导贡献满足
\[
\|\Delta^{(3)}\|_F^2\sim 4\cdot 2^3\sum_{\alpha,\beta,\mu,\nu,\gamma}\Big|\int\!\!\int K(s,t)\,h_{\alpha\beta}(s)h_{\mu\nu}(t)\,ds\,dt\Big|^2+\cdots
\]
其中 $4$ 来自 $(-2)^2$，$2^3$ 是三个位点的迹因子。于是量级为 $\|\Delta^{(3)}\|_F\sim O(J^3)$，而 $\mathrm{Tr}(\Delta^\dagger\Delta)\sim O(J^6)$。

5) 小结：文中原先把非阿贝尔度量写成两项平方和的近似（类似 $J^4$ 量级）是粗略估计；严格 Dyson 展开显示最低非平凡阶是三阶（$J^3$），需在文中注明这一点并保留时间核与系数以便数值比较。

建议在后续数值验证中直接用上面给出的 $\Delta^{(3)}$ 三重积分形式或在数值模拟中直接构造有限尺寸的 $R$ 矩阵并计算 $\Delta$ 与 $\|\Delta\|_F$，以验证系数与尺度。
}
\]

由于：

\[
\epsilon_{\beta\mu\gamma}\neq0
\iff
\beta\neq\mu
\]

因此最终：

\[
\boxed{
\mathcal N
\sim
\sqrt{
\int ds_1ds_2
\sum_{\beta\neq\mu}
|h_{\alpha\beta}(s_1)|^2
|h_{\mu\nu}(s_2)|^2
}
}
\]

---

# 16. 物理意义

因此：

\[
\boxed{
\text{non-Abelian性来自路径上的Pauli通道不对易累计}
}
\]

即：

\[
\boxed{
\sigma^x,
\sigma^y,
\sigma^z
\text{之间的竞争}
}
\]

产生了braid deformation。

---

# 17. 特殊情况

## 17.1 单Pauli通道

若：

\[
h(u)=J(u)\sigma^x\otimes\sigma^x
\]

则：

\[
[h_{12},h_{23}]=0
\]

因此：

\[
\boxed{
\mathcal N=0
}
\]

系统为Abelian。

---

## 17.2 多Pauli竞争

若：

\[
h(u)
=
J_x(u)\sigma^x\otimes\sigma^x
+
J_y(u)\sigma^y\otimes\sigma^y
\]

则：

\[
[\sigma^x,\sigma^y]\neq0
\]

因此：

\[
\boxed{
\mathcal N>0
}
\]

系统出现non-Abelian braid structure。

---

# 18. 与Schur complement联系

定义Hamiltonian分块：

\[
H=
\begin{pmatrix}
H_{PP}&H_{PQ}\\
H_{QP}&H_{QQ}
\end{pmatrix}
\]

Green函数：

\[
G(z)=(z-H)^{-1}
\]

Schur complement：

\[
\boxed{
\Sigma(z)
=
H_{PQ}(z-H_{QQ})^{-1}H_{QP}
}
\]

有效Hamiltonian：

\[
\boxed{
H_{eff}=H_{PP}+\Sigma(z)
}
\]

因此：

\[
h_{\alpha\beta}(u)
\rightarrow
h_{\alpha\beta}^{eff}(u)
}
\]

最终：

\[
\boxed{
ABS leakage
\Rightarrow
Pauli-channel mixing
\Rightarrow
non-Abelian deformation
}
\]

---

# 19. 最终统一结构

整个理论最终闭合为：

\[
\boxed{
h_{\alpha\beta}(u)
\Rightarrow
R(u)
\Rightarrow
[h_{12},h_{23}]
\Rightarrow
\Delta_{YB}
\Rightarrow
\mathcal N
}
\]

并且：

\[
\boxed{
\mathcal N
\sim
\sqrt{
\int ds_1ds_2
\sum_{\beta\neq\mu}
|h_{\alpha\beta}(s_1)|^2
|h_{\mu\nu}(s_2)|^2
}
}
\]

---

# 20. 最终物理图像

整个理论核心为：

\[
\boxed{
\text{局域Pauli张量积Hamiltonian}
\Rightarrow
\text{路径演化}
\Rightarrow
\text{Yang–Baxter braid structure}
\Rightarrow
\text{non-Abelian统计}
}
\]

而：

\[
\boxed{
\Sigma(z)
=
\text{ABS引入的额外Pauli mixing}
}
\]

最终导致：

\[
\boxed{
\text{braid deformation}
}
\]

即：

\[
\boxed{
	ext{spectral topology}
\leftrightarrow
	ext{operator algebra geometry}


## 22. Exact third-order calibration note

The third-order kernel formula above gives the correct operator structure and parameter dependence, but the comparison against the exact numerical YBE deviation shows that an overall calibration constant is still required for absolute normalization.

Define the exact third-order Frobenius norm and the kernel estimate by

$$
\mathcal N_{\mathrm{exact}}^{(3)}=\|\Delta^{(3)}_{\mathrm{exact}}\|_F,
\qquad
\mathcal N_{\mathrm{kernel}}^{(3)}=
\left[
4\cdot 2^3
\sum_{\alpha,\beta,\mu,\nu,\gamma}
\left|\int\!\!\int K(s,t)
 h_{\alpha\beta}(s)h_{\mu\nu}(t)\,ds\,dt\right|^2
\right]^{1/2}.
$$

Then the numerical comparison is summarized by

$$
\mathcal N_{\mathrm{exact}}^{(3)} \approx C_{\mathrm{cal}}\,\mathcal N_{\mathrm{kernel}}^{(3)}.
$$

So the formula is correct as a **third-order shape-preserving estimator**, while the final absolute value must be fixed by calibration in the chosen projection convention.

In the current numerical implementation used in this repository, the fitted calibration constant is approximately

$$
C_{\mathrm{cal}} \approx 5.09\times 10^2.
$$

This is a convention-dependent conversion factor between the kernel estimator and the exact third-order Frobenius norm.



# 修正与精确常数因子（补充说明）

$$
\mathcal N
\sim
\sqrt{
\int ds_1ds_2
\sum
|h_{\alpha\beta}(s_1)|^2
|h_{\mu\nu}(s_2)|^2
\epsilon_{\beta\mu\gamma}^2
}
$$



注：上面的简化式为启发性估计。下面给出更精确的低阶分析并指出被简化的常数因子与时间核。

1) Dyson 第三阶的精确形式（保留有序积分常数）：

$$
R^{(3)}(u)=(-i)^3\int_0^u ds_1\int_0^{s_1} ds_2\int_0^{s_2} ds_3\;h(s_1)h(s_2)h(s_3).
$$

2) YBE 偏差的最低非平凡阶为三阶，精确三阶差为三重时间有序积分的差：

$$
\begin{aligned}
\Delta^{(3)}={}&(-i)^3\Big[\int_0^u ds_1\int_0^{u+v} dt\int_0^v ds_3\\
&\quad\;h_{12}(s_1)h_{23}(t)h_{12}(s_3)
-\int_0^v ds_1'\int_0^{u+v} dt'\int_0^u ds_3'\;h_{23}(s_1')h_{12}(t')h_{23}(s_3')\Big].
\end{aligned}
$$

3) 把包含重疊站的項重寫為交換子形式时，会出现主导的交换子项

$$
[h_{12}(s),h_{23}(t)]=2i\sum_{\alpha,\beta,\mu,\nu,\gamma}h_{\alpha\beta}(s)h_{\mu\nu}(t)\epsilon_{\beta\mu\gamma}\;\sigma^\alpha\otimes\sigma^\gamma\otimes\sigma^\nu.
$$

代入前因子 $(-i)^3$ 与上述 $2i$，主导代数常数为
$$
(-i)^3\cdot 2i = -2.
$$
因此交换子贡献在 Δ 中带有明确的数值常数 $-2$，时间依赖由由三重有序积分合併後得到的核 $K(s,t)$（由 θ 函数組成）決定。

4) 对 Frobenius 范数的更精确缩放：主导贡献满足

$$
\|\Delta^{(3)}\|_F^2\sim 4\cdot 2^3\sum_{\alpha,\beta,\mu,\nu,\gamma}\Big|\int\!\!\int K(s,t)\,h_{\alpha\beta}(s)h_{\mu\nu}(t)\,ds\,dt\Big|^2+\cdots
$$

其中 $4$ 来自 $(-2)^2$，$2^3$ 是三个位点的迹因子。于是量级为 $\|\Delta^{(3)}\|_F\sim O(J^3)$，而 $\mathrm{Tr}(\Delta^\dagger\Delta)\sim O(J^6)$。

5) 小结：文中原先把非阿贝尔度量写成两项平方和的近似（类似 $J^4$ 量级）是粗略估计；严格 Dyson 展开显示最低非平凡阶是三阶（$J^3$），需在文中注明这一点并保留时间核与系数以便数值比较。

建议在后续数值验证中直接用上面给出的 $\Delta^{(3)}$ 三重积分形式或在数值模拟中直接构造有限尺寸的 $R$ 矩阵并计算 $\Delta$ 与 $\|\Delta\|_F$，以验证系数与尺度。

---

## 显式时间核 $K(s,t)$ 的 θ-函数表示

为了把三重时间有序积分整理成双重积分乘上交换子形式，可以把
\(\Delta^{(3)}\) 表示为
$$
\Delta^{(3)} = (-i)^3\int_0^{u} ds\int_0^{u+v} dt\;K(s,t)\,[h_{12}(s),\,h_{23}(t)],
$$
其中时间核 $K(s,t)$ 来自把第三个积分（在两项中分别为对 $s_3\in[0,v]$ 和 $s_1'\in[0,u]$ 的积分）用 Heaviside 步函数（记为 $\Theta$）表示相对次序后合并得到。令 $\Theta(x)$ 为 Heaviside 函数（$\Theta(x)=1$ 若 $x>0$，$\Theta(0)=1/2$ 可选），经代数化简得：

$$
\begin{aligned}
K(s,t)={}&\;\int_0^v d s'\big[\Theta(s-t)+\Theta(t-s')\big]\n+\;-\int_0^u d s'\big[\Theta(t-s)+\Theta(s'-t)\big]\\
={}&\;v\,\Theta(s-t)+\min(t,v) - u\,\Theta(t-s)-\max(0,u-t).
\end{aligned}
$$

解释：

- 第一项来源于 $A_{12}(u)A_{23}(u+v)A_{12}(v)$ 中固定\(s,t\) 后对第三个变量 $s'\in[0,v]$ 的积分，分解为两部分：一部分在 $s>t$ 时对 $s'$ 的总权重为 $v$（对应 $\Theta(s-t)$），另一部分为 $\int_0^v\Theta(t-s')ds'=\min(t,v)$。
- 第二项来源于交换的三重积分 $A_{23}(v)A_{12}(u+v)A_{23}(u)$，对 $s'\in[0,u]$ 的积分给出 $u\Theta(t-s)+\max(0,u-t)$ 项，作为减项出现。

因此 $K(s,t)$ 给出了明确的时间依赖核，可以用于精确评估
\(\Delta^{(3)}\) 中时间秩序对非阿贝尔强度的貢獻。在数值实现中，直接用该 $K(s,t)$ 计算双重积分通常比显式构造三重乘积更高效。