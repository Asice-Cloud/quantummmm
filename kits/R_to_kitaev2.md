### 从常数 YBE 解到指数形式的 R，再到 Kitaev 链（草稿）

这份笔记是对 [R_to_Kitaev.md](R_to_Kitaev.md) 的一个「指数形式」版本：

1. 先从满足常数 Yang–Baxter 方程的 SU(2) 不变两比特算符出发，回顾其线性形式
	$$
	R = a I + b\,\sigma^x\otimes\sigma^x + c\,\sigma^y\otimes\sigma^y + d\,\sigma^z\otimes\sigma^z,
	$$
	并把它改写为某个两体厄米算符的指数
	$$
	R = e^{iK},\qquad K=\sum_{\alpha,\beta}c_{\alpha\beta}\,\sigma^\alpha\otimes\sigma^\beta.
	$$
2. 然后把 $K$ 当作局域哈密顿密度，嵌入到整条链 $H=\sum_j K_{j,j+1}$，用 Jordan–Wigner 映射到费米子，得到 Kitaev‑型链参数 $t,\Delta,\mu$ 与 $c_{\alpha\beta}$ 的关系。


---

#### 1. 常数 YBE 解的线性形式回顾

在 [R_to_Kitaev.md](R_to_Kitaev.md) 里，我们在 SU(2) 共变性
$$
(U\otimes U)R(U^\dagger\otimes U^\dagger)=R,\qquad U\in SU(2)
$$
和常数 Yang–Baxter 方程的约束下，把 $R$ 写成对称、反对称投影算符 $P_3, P_1$ 的线性组合
$$
R = a P_3 + b P_1,
$$
这等价于在 Pauli 张量积基底下写成
$$
R = a I + b\,\sigma^x\otimes\sigma^x + c\,\sigma^y\otimes\sigma^y + d\,\sigma^z\otimes\sigma^z.
$$
这里 $a,b,c,d\in\mathbb C$ 满足一组来自 YBE 的多项式约束（原文给出的 12 条方程）。

从物理上，这只是说明：在 SU(2) 不变、作用在两个自旋 $1/2$ 的情形下，满足 YBE 的常数解 $R$ 总可以写成 $I$ 与三个各向同性自旋耦合项 $\sigma^\alpha\otimes\sigma^\alpha$ 的线性组合。


#### 2. 把线性形式重新参数化为指数形式

现在希望把同一个 $R$ 写成指数
$$
R = e^{iK},\qquad K=\sum_{\alpha,\beta} c_{\alpha\beta}\,\sigma^\alpha\otimes\sigma^\beta.
$$
为保持和原来的 SU(2) 不变情形平行，这里先考虑「对角」的情形：
$$
K = J_x\,\sigma^x\otimes\sigma^x + J_y\,\sigma^y\otimes\sigma^y + J_z\,\sigma^z\otimes\sigma^z,
$$
其中 $J_x,J_y,J_z\in\mathbb R$，因此 $K$ 是厄米的，两比特上的 $R=e^{iK}$ 是幺正算符。

一个关键的代数事实：
$$
\sigma^x\otimes\sigma^x,\ \sigma^y\otimes\sigma^y,\ \sigma^z\otimes\sigma^z
$$
两两对易，且都满足 $(\sigma^\alpha\otimes\sigma^\alpha)^2=I$。于是指数可以因式分解：
$$
\begin{aligned}
R &= e^{iJ_x\sigma^x\otimes\sigma^x}\ e^{iJ_y\sigma^y\otimes\sigma^y}\ e^{iJ_z\sigma^z\otimes\sigma^z},\\
e^{iJ_x\sigma^x\otimes\sigma^x} &= \cos J_x\,I + i\sin J_x\,\sigma^x\otimes\sigma^x,\\
e^{iJ_y\sigma^y\otimes\sigma^y} &= \cos J_y\,I + i\sin J_y\,\sigma^y\otimes\sigma^y,\\
e^{iJ_z\sigma^z\otimes\sigma^z} &= \cos J_z\,I + i\sin J_z\,\sigma^z\otimes\sigma^z.
\end{aligned}
$$

把这三个因子相乘后，$R$ 仍然落在四维线性子空间
$$
\operatorname{span}\{I,\sigma^x\otimes\sigma^x,\sigma^y\otimes\sigma^y,\sigma^z\otimes\sigma^z\}
$$
里，因此可以重新写成
$$
R = a(J) I + b(J)\,\sigma^x\otimes\sigma^x + c(J)\,\sigma^y\otimes\sigma^y + d(J)\,\sigma^z\otimes\sigma^z.
$$
也就是说：

- 线性参数 $(a,b,c,d)$ 和指数参数 $(J_x,J_y,J_z)$ 之间存在一个一一对应（在一定范围内），
- YBE 的代数约束原来是对 $(a,b,c,d)$ 的多项式条件，现在等价于对 $(J_x,J_y,J_z)$ 的某些三角函数关系。

为了后面具体比较，下面把 $a(J),b(J),c(J),d(J)$ 的三角函数表达式写出来。记
$$
A=\sigma^x\otimes\sigma^x,\quad B=\sigma^y\otimes\sigma^y,\quad C=\sigma^z\otimes\sigma^z,
$$
它们满足
$$
A^2=B^2=C^2=I,\qquad AB=BA=-C,\ BC=CB=-A,\ CA=AC=-B.
$$
写
$$
c_x=\cos J_x,\ s_x=\sin J_x,\ \text{同理 }c_y,s_y,c_z,s_z,
$$
则有
$$
e^{iJ_xA}=c_x I + i s_x A,\quad e^{iJ_yB}=c_y I + i s_y B,\quad e^{iJ_zC}=c_z I + i s_z C.
$$
直接按 $R=e^{iJ_xA}e^{iJ_yB}e^{iJ_zC}$ 展开并用上面的乘法关系化简，可得
$$
R = a(J) I + b(J)A + c(J)B + d(J)C,
$$
其中
$$
\begin{aligned}
a(J) &= c_x c_y c_z + i\, s_x s_y s_z,\\
b(J) &= i\,c_z s_x c_y + s_z c_x s_y,\\
c(J) &= i\,c_z c_x s_y + s_z s_x c_y,\\
d(J) &= c_z s_x s_y + i\,s_z c_x c_y.
\end{aligned}
$$
这给出了从指数参数 $(J_x,J_y,J_z)$ 到线性参数 $(a,b,c,d)$ 的显式三角映射；反过来若给定满足 YBE 的 $(a,b,c,d)$，原则上可以通过对上式求逆（例如先对 $R$ 对角化，取对数）得到 $(J_x,J_y,J_z)$，不过解析表达式会比较复杂，这里不展开。

在后面做 Jordan–Wigner 时，只需要知道「$R$ 仍然是 $I+\sigma^x\sigma^x+\sigma^y\sigma^y+\sigma^z\sigma^z$ 的线性组合」，于是原文中关于 $t,\Delta,\mu$ 的推导可以原封不动地搬过来，只是把 $(a,b,c,d)$ 换成 $(a(J),b(J),c(J),d(J))$。


#### 3. 选择生成元 K 作为哈密顿密度

指数形式还有一个更自然的物理读法：与其把 $H$ 直接取成 $\sum_j R_{j,j+1}$，不如把 **指数里的生成元** $K$ 解释为局域哈密顿密度：
$$
H = \sum_{j=1}^{L-1} K_{j,j+1},\qquad K_{j,j+1}=I^{\otimes(j-1)}\otimes K\otimes I^{\otimes(L-j-1)}.
$$
这样做有两个好处：

1. $K$ 本身就是「两体自旋耦合」的标准形式，易于做 Jordan–Wigner；
2. 若以后要从时间演化 $U(\tau)=e^{-iH\tau}$ 里构造 braid 算符 $U=e^{(\pi/4)\gamma_a\gamma_b}$，$H$ 正好是一个标准的二次/四次 Majorana 哈密顿量。

因此，下面的 Jordan–Wigner 推广部分都以 $K$ 为起点。


---

#### 4. 把 K 写成升降算符形式

仍然从对角的 SU(2) 形式出发：
$$
K = J_x\,\sigma^x_i\sigma^x_{i+1} + J_y\,\sigma^y_i\sigma^y_{i+1} + J_z\,\sigma^z_i\sigma^z_{i+1}.
$$
引入升降算符
$$
\sigma^\pm = \tfrac12(\sigma^x\pm i\sigma^y),\qquad \sigma^x=\sigma^++\sigma^-,\ \sigma^y=\frac{\sigma^+-\sigma^-}{i}.
$$
直接代数运算给出（在相邻的 $i,i+1$ 上）
$$
\begin{aligned}
\sigma^x_i\sigma^x_{i+1}+\sigma^y_i\sigma^y_{i+1} &= 2\big(\sigma^+_i\sigma^-_{i+1}+\sigma^-_i\sigma^+_{i+1}\big),\\
\sigma^x_i\sigma^x_{i+1}-\sigma^y_i\sigma^y_{i+1} &= 2\big(\sigma^+_i\sigma^+_{i+1}+\sigma^-_i\sigma^-_{i+1}\big).
\end{aligned}
$$
从而 $K$ 可写成
$$
\begin{aligned}
K_{i,i+1} &= J_x\,\sigma^x_i\sigma^x_{i+1} + J_y\,\sigma^y_i\sigma^y_{i+1} + J_z\,\sigma^z_i\sigma^z_{i+1}\\
&= (J_x+J_y)\big(\sigma^+_i\sigma^-_{i+1}+\sigma^-_i\sigma^+_{i+1}\big) \\
&\quad + (J_x-J_y)\big(\sigma^+_i\sigma^+_{i+1}+\sigma^-_i\sigma^-_{i+1}\big) + J_z\,\sigma^z_i\sigma^z_{i+1}.
\end{aligned}
$$
这一步和 [R_to_Kitaev.md](R_to_Kitaev.md) 中对 $(b,c,d)$ 的处理完全平行，只是把 $(b,c,d)$ 换成了 $(J_x,J_y,J_z)$。


#### 5. Jordan–Wigner 映射：从 K 到整条链的费米子形式

Jordan–Wigner 映射为
$$
\sigma^+_j=c_j^{\dagger}e^{i\pi\sum_{k<j}n_k},\qquad
\sigma^-_j=c_j e^{i\pi\sum_{k<j}n_k},\qquad
\sigma^z_j=2n_j-1,
$$
其中 $n_j=c_j^\dagger c_j$ 是占据数算符。

对最近邻 $i,i+1$，串算符在乘积中抵消，给出熟悉的结果：
$$
\begin{aligned}
\sigma^+_i\sigma^-_{i+1}+\sigma^-_i\sigma^+_{i+1}&\mapsto c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i,\\
\sigma^+_i\sigma^+_{i+1}+\sigma^-_i\sigma^-_{i+1}&\mapsto c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i,\\
\sigma^z_i\sigma^z_{i+1}&\mapsto (2n_i-1)(2n_{i+1}-1)\\
&=4n_in_{i+1}-2(n_i+n_{i+1})+1.
\end{aligned}
$$
代入 $K_{i,i+1}$：
$$
\begin{aligned}
K_{i,i+1} &= (J_x+J_y)\big(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i\big) \\
&\quad + (J_x-J_y)\big(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i\big)\\
&\quad + J_z\big(4n_in_{i+1}-2(n_i+n_{i+1})+1\big).
\end{aligned}
$$

把 $K$ 作为哈密顿密度，整条链的哈密顿量为
$$
H = \sum_{i=1}^{L-1}K_{i,i+1}.
$$
于是二次费米子部分可以识别为 Kitaev 链的跳跃与配对：
$$
H_{\text{quadratic}} = \sum_i\Big[\,t\,(c_i^{\dagger}c_{i+1}+h.c.)+\Delta\,(c_i^{\dagger}c_{i+1}^{\dagger}+h.c.)\Big] + (\text{密度线性项}),
$$
其中（在最简单的归一化下一种选择是）
$$
t = J_x+J_y,\qquad \Delta = J_x-J_y.
$$
而 $J_z$ 产生的
$$
J_z\big[4n_in_{i+1}-2(n_i+n_{i+1})+1\big]
$$
一方面给出最近邻相互作用 $4J_z\,n_in_{i+1}$，另一方面其线性密度部分
$$
-2J_z(n_i+n_{i+1})
$$
在链上求和后对应一个有效化学势 $\mu\propto 4J_z$（加上可吸收的常数偏移）。

总结为：
$$
t \propto (J_x+J_y),\qquad \Delta\propto (J_x-J_y),\qquad \mu\text{ 的线性部分来自 }J_z\text{（再加上重整化）}.
$$


#### 6. 把「指数形式」和原来的 YBE 约束联系起来

从 YBE 的角度看，原来对 $(a,b,c,d)$ 的代数约束
$$
F_k(a,b,c,d)=0,\quad k=1,\dots,12
$$
现在可以通过
$$
R=e^{iK(J_x,J_y,J_z)}\quad\Longrightarrow\quad (a,b,c,d)=(a(J),b(J),c(J),d(J))
$$
重写成对 $(J_x,J_y,J_z)$ 的一组三角关系
$$
F_k\big(a(J),b(J),c(J),d(J)\big)=0.
$$
典型情形包括：

- XXX 点：$J_x=J_y=J_z$，对应 $R\propto P$（交换算符）那条可积曲线；
- XY/自由 Kitaev 类：$J_z=0$，仅存在 $\sigma^x\sigma^x,\sigma^y\sigma^y$，对应无相互作用的 Kitaev 链；
- Ising/XXZ 相互作用：$J_x=J_y=0$ 或 $J_x=J_y\ne J_z$ 等，产生含四费米子相互作用的链。

从 Kitaev 链的角度看：

- 指数形式 $R=e^{iK}$ 并没有改变「Jordan–Wigner 后得到 Kitaev 链」的结构，只是把原来线性参数 $(b,c,d)$ 换成了某些角度 $(J_x,J_y,J_z)$ 的三角函数组合；
- YBE 的可积性条件在 $(J_x,J_y,J_z)$ 空间里选出一些特殊的曲线/面，对应 XXX、XY、Ising 等熟悉的可积族；
- 在这些子流形上，Kitaev 链的参数 $(t,\Delta,\mu)$ 也就随着 $(J_x,J_y,J_z)$ 沿着可积流动。


> 备注：上面对 $t,\Delta,\mu$ 的识别都忽略了整体常数以及可能的符号约定（例如有些文献把 Kitaev 链写成 $-t\sum c^\dagger c_{j+1}$），这些都可以通过简单的归一化或整体负号来调整，不改变拓扑相结构。

