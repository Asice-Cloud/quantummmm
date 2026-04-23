为了和普遍写法(如$U=exp(\frac{\pi}{4}\gamma_a \gamma_b)$保持同样的形式， 可以把它改写为某个两体厄米算符的指数， 这也是更通用的写法：
$$
R = e^{iH_P},\qquad H_P\in\mathfrak{su}(4) \\
\qquad
Hp =\sum_{\mu,\mu} a_{\mu\mu}\,\sigma^\mu_j\otimes\sigma^\mu_{j+1}
$$
然后把 $H_p$ 当作局域哈密顿密度，嵌入到整条链 $H=\sum_j H_p^{j,j+1}$.

下面我们先针对对角的简化的模型进行jordan wigner变换，取
$$
Hp =\sum_{\mu,\mu} a_{\mu\mu}\,\sigma^\mu_j\otimes\sigma^\mu_{j+1}
$$
其中
$$
\sigma^x\otimes\sigma^x,\ \sigma^y\otimes\sigma^y,\ \sigma^z\otimes\sigma^z
$$

两两对易，且都满足 $(\sigma^\alpha\otimes\sigma^\alpha)^2=I$。


#### Jordan Wigner

记升降算符 $\sigma^\pm=\tfrac{1}{2}(\sigma^x\pm i\sigma^y)$，于是等价写法：
$$
\sigma^x=\sigma^++\sigma^-,\qquad \sigma^y=\frac{\sigma^+-\sigma^-}{i}.
$$
直接计算得（对相邻站点 $i,i+1$）：
$$
\begin{aligned}

\sigma^x_i\sigma^x_{i+1}+\sigma^y_i\sigma^y_{i+1}&=2\big(\sigma^+_i\sigma^-_{i+1}+\sigma^-_i\sigma^+_{i+1}\big),\\

\sigma^x_i\sigma^x_{i+1}-\sigma^y_i\sigma^y_{i+1}&=2\big(\sigma^+_i\sigma^+_{i+1}+\sigma^-_i\sigma^-_{i+1}\big).

\end{aligned}
$$
Jordan Wigner transformation:
$$
\sigma^+_j=c_j^{\dagger}e^{i\pi\sum_{k<j}n_k},\qquad \sigma^-_j=c_j e^{i\pi\sum_{k<j}n_k},\qquad \sigma^z_j=2n_j-1.
$$
计算最近邻的交换类项（例如 $\sigma^+_i\sigma^-_{i+1}$）
$$
\begin{aligned}
\sigma^+_i\sigma^-_{i+1}
&= c_i^{\dagger} e^{i\pi\sum_{k<i}n_k}\; c_{i+1} e^{i\pi\sum_{k<i+1}n_k} \\
&= c_i^{\dagger} c_{i+1}\; e^{i\pi\sum_{k<i}n_k}\,e^{i\pi\sum_{k<i+1}n_k}.
\end{aligned}
$$
注意到 $\sum_{k<i+1}n_k=\sum_{k<i}n_k + n_i$，==>
$$
e^{i\pi\sum_{k<i}n_k}\,e^{i\pi\sum_{k<i+1}n_k}=e^{i\pi(2\sum_{k<i}n_k + n_i)}=e^{i\pi n_i},
$$
因为 $e^{2i\pi\sum_{k<i}n_k}=1$。 ==>
$$
\sigma^+_i\sigma^-_{i+1} = c_i^{\dagger} c_{i+1}\; e^{i\pi n_i}.
$$
对于$e^{iπ n_i}=(-1)^{n_i}$, 占据数基$ |… n_i,n_{i+1} …⟩$ 上有$ n_i,n_{i+1}∈{0,1}$,  $c_i^† c_{i+1} $只有在$ n_i=0$ 且$ n_{i+1}=1$ 时产生非零结果而在这情形下$ e^{iπ n_i}=1$；其它情形下要么两边同时为 0，要么因占据限制而抵消。配对项 $c_i^† c_{i+1}^† $的情形同理（只有在$ n_i=n_{i+1}=0$ 时非零，此时$ e^{iπ n_i}=1$）。

于是对最近邻 $i,i+1$，串算符 $e^{i\pi\sum_{k<i}n_k}$ 在乘积中相互抵消，因此有
$$
\begin{aligned}
\sigma^+_i\sigma^-_{i+1}+\sigma^-_i\sigma^+_{i+1}&\mapsto c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i,\\
\sigma^+_i\sigma^+_{i+1}+\sigma^-_i\sigma^-_{i+1}&\mapsto c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i,\\
\sigma^z_i\sigma^z_{i+1}&\mapsto (2n_i-1)(2n_{i+1}-1)=4n_in_{i+1}-2(n_i+n_{i+1})+1.
\end{aligned}
$$


对应的局部 R‑矩阵为
$$
R_{j,j+1} = e^{iH_P^{(j,j+1)}}.
$$
总哈密顿量取为
$$
H = \sum_{j=1}^{N-1} H_P^{(j,j+1)}.\\
Hp =\sum_{\mu,\mu} a_{\mu\mu}\,\sigma^\mu_j\otimes\sigma^\mu_{j+1}
$$
代入近邻项的哈密顿量 $H_p^{i,i+1}$：
$$
\begin{aligned}
H_p^{i,i+1} &= J_x\,\sigma^x_i\sigma^x_{i+1} + J_y\,\sigma^y_i\sigma^y_{i+1} + J_z\,\sigma^z_i\sigma^z_{i+1}\\
&= (J_x+J_y)\big(\sigma^+_i\sigma^-_{i+1}+\sigma^-_i\sigma^+_{i+1}\big) \\
&\quad + (J_x-J_y)\big(\sigma^+_i\sigma^+_{i+1}+\sigma^-_i\sigma^-_{i+1}\big) + J_z\,\sigma^z_i\sigma^z_{i+1}.
\end{aligned}
$$
再带入$\sigma \rightarrow c$
$$
\begin{aligned}
H_p^{i,i+1} &= (J_x+J_y)\big(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i\big) \\
&\quad + (J_x-J_y)\big(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i\big)\\
&\quad + J_z\big(4n_in_{i+1}-2(n_i+n_{i+1})+1\big).
\end{aligned}
$$
整条链的哈密顿量为 $H = \sum_{i=1}^{L-1}H_{i,i+1}.$

对照
$$
H= -\mu\sum_j\Big(n_j-\tfrac12\Big) - t\sum_j\big(c_j^{\dagger}c_{j+1}+h.c.\big) + \Delta\sum_j\big(c_j c_{j+1}+h.c.\big)+\cdots
$$
所有来自 $c_{xx}\sigma^x_j\sigma^x_{j+1}$ 与 $c_{yy}\sigma^y_j\sigma^y_{j+1}$ 的贡献，在 JW 之后都是严格二次的费米子算符：
$$
H_{\text{quadratic}} = \sum_i\Big[\,t\,(c_i^{\dagger}c_{i+1}+h.c.)+\Delta\,(c_i^{\dagger}c_{i+1}^{\dagger}+h.c.)\Big] + (\text{密度线性项}),
$$
因此并入 $H_{\mathrm{quad}}$，决定 $t,\Delta$ 等参数. 而 $J_z$ 产生的
$$
J_z\big[4n_in_{i+1}-2(n_i+n_{i+1})+1\big]
$$
一方面给出最近邻相互作用 $4J_z\,n_in_{i+1}$，另一方面其线性密度部分$-2J_z(n_i+n_{i+1})$

所以 $c_{zz}\,\sigma^z_j\sigma^z_{j+1}$ 贡献可以自然写成为：
$$
H_{\mathrm{int}}\supset 4c_{zz}\sum_j n_jn_{j+1},\qquad H_{\mathrm{gauge}}\supset -2c_{zz}\sum_j(n_j+n_{j+1})+\text{const}.
$$


上面自然的出现三类项：

- $H_{\mathrm{quad}}$ 包括所有**二次费米子项**（即所有形如 $c^\dagger c$, $c^\dagger c^\dagger$ 及其厄米共轭的项），因此化学势项 $\mu\sum_j n_j$ 属于 $H_{\mathrm{quad}}$；这些项决定能带、谱隙和零模式结构。

- $ H_{\mathrm{int}}$ 包括所有**四次及以上**的多体作用项（如 $n_jn_{j+1}$ 等）。

- $H_{\mathrm{gauge}}$ 包括纯常数项、纯局域常数能量位移或仅修正的项。

因此总哈密顿量自然分解为
$$
H = H_{\mathrm{quad}} + H_{\mathrm{int}} + H_{\mathrm{gauge}},
$$



#### 非对角形式

如果把生成元写成更一般的形式
$$
K = \sum_{\mu,\nu\in\{0,x,y,z\}} c_{\mu\nu}\,\sigma^\mu\otimes\sigma^\nu,
\qquad R=e^{iK},
$$
- $\sigma^x_i\sigma^x_{i+1}$, $\sigma^y_i\sigma^y_{i+1}$：如前，组合为 hopping 和 pairing 的二次項，字符串因子在最近邻乘积中抵消。

- 交叉项（例如） $\sigma^x_i\sigma^y_{i+1}$：展开为升降算符后经过 JW，可得到带复系数的二次组合（hopping/pairing）：
	$$
	\sigma^x_i\sigma^y_{i+1}\xrightarrow{\mathrm{JW}} \tfrac{1}{i}\big(c_i^\dagger c_{i+1}^\dagger - c_i^\dagger c_{i+1} + c_i c_{i+1}^\dagger - c_i c_{i+1}\big),
	$$
	因此,此类交叉项仍然并入 $H_{\mathrm{quad}}$。

- 含 $z$ 的混合项（例） $\sigma^x_i\sigma^z_{i+1}$：写成
$$
	\sigma^x_i\sigma^z_{i+1}=(\sigma^+_i+\sigma^-_i)(2n_{i+1}-1),
$$

JW 后留下串因子：
$$
	\sigma^x_i\sigma^z_{i+1}\xrightarrow{\mathrm{JW}} e^{i\pi\sum_{k<i}n_k}\,(c_i^\dagger+c_i)\,(2n_{i+1}-1),
$$

这个算符在费米子表述中含有非局域的串算符，不能被直接写成纯二次 BdG 形式（通常会产生高阶或非局域效应）。

- $\sigma^z_i\sigma^z_{i+1}$：如前
$$
	\sigma^z_i\sigma^z_{i+1}\xrightarrow{\mathrm{JW}} 4n_in_{i+1}-2(n_i+n_{i+1})+1,
$$
给出四费米子、线性密度与常数三类贡献。

**这一部分完整结果放在附录**，因为后面需要用到的是以上的哈密顿形式可以被自然分成三部分这个结论。

### 附录：JW 映射到 Majorana 双线性（显式对照）

在第 1、2 站点上定义 Majorana 算符
$$
\gamma_1=c_1+c_1^{\dagger},\qquad \gamma_2=i(c_1-c_1^{\dagger}),\qquad
\gamma_3=c_2+c_2^{\dagger},\qquad \gamma_4=i(c_2-c_2^{\dagger}).
$$
对常见的最近邻二体基元，Jordan–Wigner 后再写成 Majorana 双线性，有如下显式映射（按本文约定的因子与编号）：
$$
\begin{aligned}
\sigma^x\otimes\sigma^x &\mapsto \tfrac{1}{2}\,i\big(\gamma_1\gamma_4 + \gamma_2\gamma_3\big),\\
\sigma^y\otimes\sigma^y &\mapsto \tfrac{1}{2}\,i\big(-\gamma_1\gamma_4 + \gamma_2\gamma_3\big),\\
\sigma^x\otimes\sigma^y &\mapsto \tfrac{1}{2}\,i\big(\gamma_1\gamma_4 - \gamma_2\gamma_3\big),\\
\sigma^y\otimes\sigma^x &\mapsto \tfrac{1}{2}\,i\big(\gamma_1\gamma_4 + \gamma_2\gamma_3\big),\\
\sigma^z\otimes\sigma^z &\mapsto 4n_1n_2 - 2(n_1+n_2) + 1.
\end{aligned}
$$

因此，当 $H_P$ 仅由 XX, YY, XY, YX 四类二次基元主导时，$H_P$ 的二次部分可写为（去掉常数项）：
$$
H_P^{(2)} = i\sum_{a<b} A_{ab}\,\gamma_a\gamma_b
= \tfrac{i}{2}\Big[ (c_{xx}+c_{xy}+c_{yx}+c_{yy})\,\gamma_1\gamma_4 + (c_{xx}-c_{xy}+c_{yx}-c_{yy})\,\gamma_2\gamma_3\Big],
$$
或按线性代数写成 $A_{ab}$ 关于 $c_{\mu\nu}$ 的线性组合（不同项的系数由上面映射直接读出）。

（注：这里采用了在映射中常见的 $1/2$ 因子与 $i$ 前置的习惯，方便将 Pauli 基元直接对应到 $i\gamma_a\gamma_b$ 形式；含 $\sigma^z\sigma^z$ 的项会产生四次项与线性密度项，因此不属于纯双线性子空间。）


---------------------

例子：选一个简单的情形：
$$
H_P\big|_{\langle j,j+1\rangle}
 = J_x\,\sigma^x_j\sigma^x_{j+1}
 + J_y\,\sigma^y_j\sigma^y_{j+1}
 + J_z\,\sigma^z_j\sigma^z_{j+1}
 + h\,(\sigma^z_j+\sigma^z_{j+1})
 + \varepsilon\,I,
$$
也就是只含最常见的 XX/YY/ZZ 最近邻耦合和一个 on-site 磁场（$h$）以及常数项（$\varepsilon$）。对开放链
$$
H = \sum_{j=1}^{L-1} H_P\big|_{\langle j,j+1\rangle}
$$
做 Jordan–Wigner 映射后，利用前面的公式可以一步得到
$$
\begin{aligned}
H = \sum_{j=1}^{L-1}\Big[(J_x+J_y)\big(c_j^\dagger c_{j+1}+c_{j+1}^\dagger c_j\big)
+ (J_x-J_y)\big(c_j^\dagger c_{j+1}^\dagger+c_{j+1}c_j\big) \\
+ 4J_z\,n_jn_{j+1} + (-2J_z+2h)\,(n_j+n_{j+1})+ (J_z-2h+\varepsilon)\Big].
\end{aligned}
$$
把它按前面的三部分重写，就是
$$
\begin{aligned}
H_{\mathrm{quad}} &= \sum_{j=1}^{L-1}\Big[ t\,(c_j^\dagger c_{j+1}+c_{j+1}^\dagger c_j)
 + \Delta\,(c_j^\dagger c_{j+1}^\dagger+c_{j+1}c_j)\Big] - \mu\sum_j n_j,\\
H_{\mathrm{int}} &= U\sum_{j=1}^{L-1} n_jn_{j+1},\\
H_{\mathrm{gauge}} &= (L − 1) (J_z − 2h + \epsilon).
\end{aligned}
$$
其中参数 $t,\Delta,U,\mu$ 都是 $J_x,J_y,J_z,h,\varepsilon$ 的线性组合；例如在上面的具体例子中可以取
$$
t=J_x+J_y,\qquad \Delta=J_x-J_y,\qquad U=4J_z,
$$
而 $\mu$ 则由 $J_z,h$ 的线性组合给出（精确表达式需要把 $\sum_j(n_j+n_{j+1})$ 在整条链上展开并单独处理端点，但对体相和拓扑结构来说只是一种化学势的重新定义）。



#### 完整的 JW 展开和哈密顿分解

下面把最近邻生成元在通用系数 $c_{\mu\nu}$ 下的 JW 展开。记 $n_j=c_j^\dagger c_j$，串算符 $S_j=\exp(i\pi\sum_{k<j}n_k)$，$\sigma^\pm=\tfrac12(\sigma^x\pm i\sigma^y)$。对键 $\langle i,i+1\rangle$：

局域基元的 JW 映射（最近邻）
- $I\otimes I\mapsto 1$
- $\sigma^x_i\otimes I\mapsto S_i(c_i^\dagger+c_i)$
- $\sigma^y_i\otimes I\mapsto S_i(-i c_i^\dagger + i c_i)$
- $\sigma^z_i\otimes I\mapsto 2n_i-1$
- 同理作用在 $i+1$ 上的单体：$I\otimes\sigma^x_{i+1}\mapsto S_{i+1}(c_{i+1}^\dagger+c_{i+1})$, 等等。

两体项（最近邻）
- $\sigma^x_i\sigma^x_{i+1}\mapsto c_i^\dagger c_{i+1}^\dagger + c_i^\dagger c_{i+1} + c_i c_{i+1}^\dagger + c_i c_{i+1}$
- $\sigma^y_i\sigma^y_{i+1}\mapsto -c_i^\dagger c_{i+1}^\dagger + c_i^\dagger c_{i+1} + c_i c_{i+1}^\dagger - c_i c_{i+1}$
- $\sigma^x_i\sigma^y_{i+1}\mapsto -i c_i^\dagger c_{i+1}^\dagger + i c_i^\dagger c_{i+1} - i c_i c_{i+1}^\dagger + i c_i c_{i+1}$
- $\sigma^y_i\sigma^x_{i+1}\mapsto -i c_i^\dagger c_{i+1}^\dagger - i c_i^\dagger c_{i+1} + i c_i c_{i+1}^\dagger + i c_i c_{i+1}$
- $\sigma^z_i\sigma^z_{i+1}\mapsto 4n_i n_{i+1}-2(n_i+n_{i+1})+1$
- 含单侧 $z$ 的混合项：
    - $\sigma^x_i\sigma^z_{i+1}=(\sigma^+_i+\sigma^-_i)(2n_{i+1}-1)\mapsto S_i(c_i^\dagger+c_i)(2n_{i+1}-1)$
    - $\sigma^y_i\sigma^z_{i+1}\mapsto S_i(-i c_i^\dagger+i c_i)(2n_{i+1}-1)$
    - $\sigma^z_i\sigma^x_{i+1}\mapsto (2n_i-1)S_{i+1}(c_{i+1}^\dagger+c_{i+1})$
    - $\sigma^z_i\sigma^y_{i+1}\mapsto (2n_i-1)S_{i+1}(-i c_{i+1}^\dagger+i c_{i+1})$

于是键 $\langle i,i+1\rangle$ 的 JW 表达为线性叠加：
$$
\begin{aligned}
K^{(i,i+1)}_{JW} =
&\; c_{00} \\
&+ c_{x0} S_i(c_i^\dagger+c_i) + c_{y0} S_i(-i c_i^\dagger+i c_i) + c_{z0}(2n_i-1) \\
&+ c_{0x} S_{i+1}(c_{i+1}^\dagger+c_{i+1}) + c_{0y} S_{i+1}(-i c_{i+1}^\dagger+i c_{i+1}) + c_{0z}(2n_{i+1}-1) \\
&+ c_{xx}(c_i^\dagger c_{i+1}^\dagger + c_i^\dagger c_{i+1} + c_i c_{i+1}^\dagger + c_i c_{i+1}) \\
&+ c_{yy}(-c_i^\dagger c_{i+1}^\dagger + c_i^\dagger c_{i+1} + c_i c_{i+1}^\dagger - c_i c_{i+1}) \\
&+ c_{xy}(-i c_i^\dagger c_{i+1}^\dagger + i c_i^\dagger c_{i+1} - i c_i c_{i+1}^\dagger + i c_i c_{i+1}) \\
&+ c_{yx}(-i c_i^\dagger c_{i+1}^\dagger - i c_i^\dagger c_{i+1} + i c_i c_{i+1}^\dagger + i c_i c_{i+1}) \\
&+ c_{xz} S_i(c_i^\dagger+c_i)(2n_{i+1}-1) + c_{yz} S_i(-i c_i^\dagger+i c_i)(2n_{i+1}-1) \\
&+ c_{zx} (2n_i-1)S_{i+1}(c_{i+1}^\dagger+c_{i+1}) + c_{zy} (2n_i-1)S_{i+1}(-i c_{i+1}^\dagger+i c_{i+1}) \\
&+ c_{zz}(4n_i n_{i+1}-2(n_i+n_{i+1})+1).
\end{aligned}
$$

全链哈密顿量为 $H_{JW}=\sum_{i=1}^{L-1} K^{(i,i+1)}_{JW}$。

三部分（按算符次数/性质）自然分解：
- $H_{\mathrm{quad}}$：所有不含 $S_j$ 的二次项（hopping/pairing）与线性密度项，系数为 $c_{xx},c_{yy},c_{xy},c_{yx},c_{z0},c_{0z},\dots$ 的线性组合；。
- $H_{\mathrm{int}}$：所有纯四次项，主要来自 $c_{zz}\,4n_i n_{i+1}$，以及多项联合可能产生的高阶项。
- $H_{\mathrm{string}}$（或 $H_{\mathrm{nonlocal}}$）：所有含串前缀 $S_j$ 的项（如含单侧 $z$ 的混合项及单体 x/y），在费米子表示上是非局域的，需要微扰、平均场或数值方法处理（这是ai的建议，暂时未尝试）。
- $H_{\mathrm{gauge}}$：常数项与可吸收的能量零点（如 $c_{00}$、$c_{zz}$ 的常数贡献等）。


需要给出
$c_{\mu\nu}$ 到 Kitaev‑链参数的线性映射矩阵 $M$, 对应$t,\Delta,\mu ...$
\
---

下面把完整的线性映射与可解性判定追加在此：

约定 16 分量系数向量（列向量）：
$$
C = [c_{xx},\; c_{yy},\; c_{xy},\; c_{yx},\; c_{zz},\; c_{z0},\; c_{0z},\; c_{x0},\; c_{0x},\; c_{y0},\; c_{0y},\; c_{xz},\; c_{zx},\; c_{yz},\; c_{zy},\; c_{00}]^T.
$$

目标物理参数向量：
$$
p = [t,\; \Delta,\; U,\; \mu,\; C_{\mathrm{perbond}}]^T,
$$
其中 $t$ 为最近邻 hopping，$\Delta$ 为配对幅度，$U$ 为最近邻相互作用强度，$\mu$ 为体化学势（哈密顿写作 $-\mu\sum_j n_j$），$C_{\mathrm{perbond}}$ 為每键常数能量偏移。

显式线性关系 $p = M\cdot C$（按列顺序与上面 $C$ 一致）。为了便于直接代入，这里给出显式分量形式以及对应的 5×16 符号矩阵：

逐分量公式：
- $t = c_{xx} + c_{yy} + i\, ( c_{xy} - c_{yx} ),$
- $\Delta = c_{xx} - c_{yy} - i\, ( c_{xy} + c_{yx} ),$
- $U = 4\,c_{zz},$
- $\mu = 4\,c_{zz} - 2\,( c_{z0} + c_{0z} ),$
- $C_{\mathrm{perbond}} = c_{zz} - ( c_{z0} + c_{0z} ) + c_{00}.$

对应矩阵（行顺序 $[t,\Delta,U,\mu,C_{\mathrm{perbond}}]$，列顺序如 $C$）：
$$
M = \begin{pmatrix}
 1 &  1 &  i & -i & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 1 & -1 & -i & -i & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 &  0 &  0 &  0 & 4 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 &  0 &  0 &  0 & 4 & -2 & -2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 &  0 &  0 &  0 & 1 & -1 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}.
$$

注：上面包含复系数以允许一般的 XY/混合项；若所有 $c_{\mu\nu}$ 为实并满足有关的厄米性约束（例如 $c_{xy}=-c_{yx}$ 等），则 $t,\Delta$ 将是实数或成对共轭，哈密顿仍为厄米。

串算符（非局域项）触发分量：
- 下列任何分量非零会在 JW 映射后引入字符串前缀 $S_j=\exp(i\pi\sum_{k<j} n_k)$，从而产生非局域项：
	- 单体 x/y 分量： $c_{x0},\; c_{y0},\; c_{0x},\; c_{0y}$；
	- 含单侧 z 的混合项： $c_{xz},\; c_{yz},\; c_{zx},\; c_{zy}$。

可解性快速判定规则：
- 自由二次（Kitaev/XY，可 Fourier + Bogoliubov 严格对角化）：当并且仅当所有字符串分量都为 0 且 $c_{zz}=0$，即
$$
c_{x0}=c_{0x}=c_{y0}=c_{0y}=c_{xz}=c_{zx}=c_{yz}=c_{zy}=0,\quad c_{zz}=0.
$$

- XXZ（Bethe‑ansatz 可积，非自由但可积）：若所有字符串分量 = 0 且 $c_{zz}\neq0$，模型等价于 XXZ/Ising‑like 链，可用 Bethe‑ansatz 或已知可积方法处理；参数为 $U=4c_{zz}$，并且 $t,\Delta$ 由上式给出。

- 一般不可积/非局域：若任一字符串分量非零，则 JW 后出现非局域串项，通常不属于上述两类，需要数值、平均场或逐项近似分析；特殊可解子情形仍可能存在但需逐一检查。

自动化建议：我可以把这一套判定写成小脚本（Python），接受任意 16 个 $c_{\mu\nu}$ 并输出 $p$ = $(t,\Delta,U,\mu,C_{\mathrm{perbond}})$ 与可解性类别（并列出触发串项的分量）。如需我现在添加到仓库，请回复“添加脚本”。

---

（已追加完毕）

---

测试用最小修改参数（可直接被工具解析）：

c_xx = 1.1
c_yy = 1.0
c_xy = 0.0
c_yx = 0.0
c_zz = 0.0
c_z0 = 0.0
c_0z = 0.0
c_00 = 0.0

（将上述 c_{μν} 代入映射可得到 Δ = 0.1, t = 2.1, μ = 0）
