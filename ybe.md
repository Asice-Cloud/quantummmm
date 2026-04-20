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



下面给出常用的子集 $c_{\mu\nu}$ 到 Kitaev‑链参数的线性映射矩阵 $M$。为简洁起见，我们在这里只列出直接产生对称实的跳跃/配对、密度线性项、最近邻相互作用与键常数的分量；含串前缀 $S_j$ 的单侧 $x,y$ 项（例如 $c_{x0},c_{0x},c_{y0},c_{0y}$）以及产生手征/虚系数的交叉项（例如 $c_{xy},c_{yx}$ 对实对称 $t,\Delta$ 的贡献会出现在虚部或反对称型项）并不计入下表的五个输出量，但在文中已经按项列出并归类为 $H_{\mathrm{string}}$ 或可产生手征跳跃的项。

定义列向量的 $c$ 分量顺序为
$$
c = [c_{xx},\; c_{yy},\; c_{xy},\; c_{yx},\; c_{zz},\; c_{z0},\; c_{0z},\; c_{00}]^T.
$$
我们把目标参数排列为
$$
p = [t,\;\Delta,\;\mu_{\text{site}},\; U,\; E_{\text{bond}}]^T,
$$
其中：
- $t$：最近邻实对称跳跃系数（乘以 $c_i^\dagger c_{i+1}+\mathrm{h.c.}$）；
- $\Delta$：最近邻实配对系数（乘以 $c_i^\dagger c_{i+1}^\dagger+\mathrm{h.c.}$）；
- $\mu_{\text{site}}$：每一格点前的线性密度系数（即 $\sum_j\mu_{\text{site}}\,n_j$ 的系数）；
- $U$：最近邻密度‑密度相互作用系数（乘以 $n_i n_{i+1}$）；
- $E_{\text{bond}}$：每一条键的常数能量位移（键常数，方便把能量写成 $\sum_{\langle i,i+1\rangle} E_{\text{bond}}$）。

在上面约定下，线性映射 $p=M\,c$ 的显式矩阵为：
$$
M=\begin{pmatrix}
1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
1 & -1& 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & -4& 2 & 2 & 0 \\
0 & 0 & 0 & 0 & 4 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & -1& -1& 1
\end{pmatrix}.
$$

注解（来自上文逐项展开的系数）
- 行 1 ($t$)：$t=c_{xx}+c_{yy}$（对称实跳跃）；交叉项 $c_{xy},c_{yx}$ 产生反对称/虚部跳跃（手征分量），不计入此标量 $t$。
- 行 2 ($\Delta$)：$\Delta=c_{xx}-c_{yy}$（实配对）；交叉项可生成虚部配对成分。
- 行 3 ($\mu_{\text{site}}$)：每格点的线性密度系数为 $\mu_{\text{site}}=2(c_{z0}+c_{0z})-4c_{zz}$（从 $\sigma^z\otimes I, I\otimes\sigma^z$ 和 $\sigma^z\sigma^z$ 的展开得到；边界效应在开链时需单独处理）。
- 行 4 ($U$)：最近邻密度交互作用 $U=4c_{zz}$ 来自 $\sigma^z_i\sigma^z_{i+1}\mapsto 4n_i n_{i+1}+\cdots$。
- 行 5 ($E_{\text{bond}}$)：每键常数 $E_{\text{bond}}=c_{00}-c_{z0}-c_{0z}+c_{zz}$（来自 $I\otimes I$、单体 $z$ 的常数项与 $\sigma^z\sigma^z$ 的常数项）。