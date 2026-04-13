## Solution of Yang Baxter Equation and Integrable Sub manifold

### Content

1. **YBE的可行解及辫子群**
2. 可积子流形
3. **时间演化braid -> 空间变换**

### 1.YBE的可行解及辫群关系





Quantum dynamical ybe (QDYBE)的定义$^{[1]}$为
$$
R: V\otimes V \rightarrow V \otimes V\\
(R⊗id)∘(id⊗R)∘(R⊗id)=(id⊗R)∘(R⊗id)∘(id⊗R)
$$
这里面给出了的一个general form是形如：
$$
R = \sum_{\alpha}E_{\alpha \alpha} \otimes E_{\alpha \alpha} + \sum_{\alpha \ne \beta}E_{\alpha \alpha} \otimes E_{\beta \beta} + \sum_{\alpha \ne \beta}E_{\alpha \beta} \otimes E_{\beta \alpha}
$$
借助这个简单的例子，下面将R写成了
$$
R=\sum_{\mu,\nu\in\{0,x,y,z\}} a_{\mu\nu}\,\sigma^\mu_j\otimes\sigma^\nu_{j+1}
$$
的形式，相比之前的四元数模型，是将控制门形式改写成了现在的相互作用形式。



但是QDYBE的约束都非常复杂，所以这里先尝试了R是在SU(2)作用下不变，即$(U \otimes U)R(U^{\dagger}\otimes U^{\dagger}) = R, U\in SU(2)$, 或者表示成李代数形式 $[R,\Delta J] = 0, \forall J \in su(2)$ , 

由此得出了$R = a P_3 + b P_1$的形式，其中$P_3 = \frac{1+P}{2}, P_1 = \frac{1-P}{2} , P:=P(\alpha \otimes \beta) = \beta \otimes \alpha$ ， 即对应这SU(2)群的分解

$2 \times 2 = 3+1$的形式，$P_i$则表示了投影到对称、反对称部分， 因为投影算符的性质，这个R是满足YBE的。

我们下面也会用到这个简化的对角形式：
$$
R=a\,I + b\,\sigma^x\sigma^x + c\,\sigma^y\sigma^y + d\,\sigma^z\sigma^z.
$$


#### 满足YBE的R作为算符

格点体系的全 Hilbert 空间为单点态空间 $V$ 的张量积 $V^{\otimes L}$。任何作用在两点上的算符$R:\;V\otimes V\to V\otimes V$
可唯一嵌入到全链为局域算符$^{[2]}$：
$$
R_{i,i+1}=I^{\otimes(i-1)}\otimes R\otimes I^{\otimes(L-i-1)}.
$$
作用在第 i, i+1 位用张量积表示,是局域两体算符。

如果$R$ 满足 Yang–Baxter 方程，则这些局域嵌入满足局域 YBE，从而给出 braid‑group 的表示和可积传输矩阵的局部生成元——这正是把代数生成元 $b_i$ 映为 $\mathrm{id}^{\otimes(i-1)}\otimes R\otimes\mathrm{id}^{\otimes(n-i-1)}$ 的代数依据$^{[2]}$。满足辫群关系的参数(a,b,c,d)见附录.

R满足YBE后，即满足了短程算符，下面简单证明该形式也可以满足长程算符：

假设 \(\{i,j\}\cap\{k,l\}=\varnothing\)。任取一对基元项
$$
T_{ij}^{\alpha\beta} = \sigma_i^\alpha\sigma_j^\beta,\qquad

T_{kl}^{\gamma\delta} = \sigma_k^\gamma\sigma_l^\delta,
$$
它们只在四个两两不同的 site 上非平凡。由于不同 site 上的泡利对易，
$$
T_{ij}^{\alpha\beta} T_{kl}^{\gamma\delta}

= \sigma_i^\alpha\sigma_j^\beta\sigma_k^\gamma\sigma_l^\delta

= \sigma_k^\gamma\sigma_l^\delta\sigma_i^\alpha\sigma_j^\beta

= T_{kl}^{\gamma\delta} T_{ij}^{\alpha\beta}.
$$


于是对所有 \(\alpha,\beta,\gamma,\delta\) 都有
$$
[T_{ij}^{\alpha\beta},T_{kl}^{\gamma\delta}] = 0.
$$


由于对易子的双线性，这立刻推出
$$
[R_{ij},R_{kl}] = 0,\qquad

\{i,j\}\cap\{k,l\}=\varnothing,
$$
短程上会出现 $\sigma_i^x\sigma_i^z = i\sigma_i^y,\qquad\sigma_i^z\sigma_i^x = -i\sigma_i^y$ ,因此是三元关系




#### 改写形式

为了和普遍写法(如$U=exp(\frac{\pi}{4}\gamma_a \gamma_b)$保持同样的形式， 可以把它改写为某个两体厄米算符的指数， 这也是更通用的写法：
$$
R = e^{iH_P},\qquad H_P\in\mathfrak{su}(4) \\
\qquad
H_P = \sum_{i,j\in\{0,x,y,z\}} c_{ij}\,\sigma_i\otimes\sigma_j,\quad (\sigma_0 = I)
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



### 2.可积子流形

使用配置空间来构造辫子操作，先定义如下：

- $\Sigma$：底层二维流形，描述**样品整体的拓扑类型和边界条件**。
    - 平面 / 盘：单个有限样品、有物理边界；
    - 环面 $T^2$：周期边界条件（无边缘），如无限大晶体 + 周期边界...；
    - 高 genus 曲面：等效于在系统中引入若干的柄，相当于多个非平凡环路，增加全局拓扑简并度。

- $X=\{x_i\}$：缺陷 / 涡旋 / 端点的集合，描述**哪些位置存在局域拓扑缺陷，从而绑定 Majorana 零模或任意子**。
    - p+ip 超导中的涡旋核心；
    - Kitaev 蜂窝模型中的 $\pi$-flux plaquette（vison）；
    - 一维拓扑超导链的端点，等等。



例如在 2D p+ip 或 Kitaev 蜂窝模型中，引入若干缺陷（涡旋、vison 等），记其在底空间 $\Sigma$ 上的位置为
$$
X = (x_1,\dots,x_n)\in \mathcal C = \{X\mid x_i\in\Sigma,\ x_i\neq x_j\}.
$$
配置空间 \(\mathcal C\)的坐标为 \(\lambda=(\lambda^1,\lambda^2,\dots)\)。在每个 \(\lambda\) 上，选定一个有限维“低能子空间” \(\mathcal H_0(\lambda)\subset\mathcal H\)，由一组正交归一的本征态
$$
\{|\psi_a(\lambda)\rangle\}_{a=1}^{N_0}
$$
张成（例如零模/基态简并空间）。$\mathcal H_0(X) = \mathrm{span}\{\text{零能模态/简并基态}\}.$

这些子空间拼接成一个 Hilbert 丛
$$
\pi:\mathcal E\to\mathcal C,\qquad \pi^{-1}(X)=\mathcal H_0(X).
$$

定义投影算符
$$
P(\lambda) = \sum_{a=1}^{N_0} |\psi_a(\lambda)\rangle\langle\psi_a(\lambda)|.
$$

定义Berry 联络与曲率为
$$
A_{ij}=i\langle\psi_i(x)|\,\mathrm d\psi_j(x)\rangle
$$
Lie 代数值：
$$
\mathcal A = \sum_{\mu} A_\mu(\lambda)\,d\lambda^\mu,\quad A_\mu(\lambda)\in\mathfrak g,
$$
其中 \(\mathfrak g\) 在我们的 Majorana/BdG 场景下通常是 \(so(2n)\) 或其子代数。曲率为
$$
F = d\mathcal A + \mathcal A\wedge\mathcal A.\\
F_{ij}=dA_{ij}+A_{ij}\wedge A_{ij}, \\
$$

使用投影算符写出来就是：
$$
F_{\mu\nu}(X) = P[\partial_\mu P,\partial_\nu P]P\\
=\partial_\mu A_{\nu}-\partial_\nu A_{\mu}+[A_{\mu},A_{\nu}]
$$
(这里思路源自于Non-Abelian Stokes Theorem$^{[3]}$,比原文中多了dA的形式,具体推导在附录，因为这个构造复杂这里不展开了)

在规范变换
$$
|\psi_a(\lambda)\rangle \mapsto \sum_b |\psi_b(\lambda)\rangle\,U_{ba}(\lambda),
\quad U(\lambda)\in U(N_0)
$$
下，联络按
$$
A_\mu \mapsto U^{-1}A_\mu U + i\,U^{-1}\partial_\mu U
$$
变换，标准的规范联络变换律。



沿着闭合路径 $\gamma:[0,1]\to\mathcal C$的平行移动（holonomy）给出拓扑演化算符：
$$
U_R[\gamma]=\mathcal P\exp\Big(-\int_\gamma A\Big),\\
= \mathcal P\exp\Bigl(-\int_0^1 A_\mu(\gamma(s))\,\dot X^\mu(s)\,ds\Bigr)
$$
其中 $\mathcal P$ 表示路径有序指数,是在零模子空间上的幺正算符。下面的部分会说明这就是MCG群的生成元（Dehn twist）



> 若曲率 $F=0$ 或仅取中心值，则 $U[\gamma]$ 只依赖于路径的同伦类，[$\gamma$] 对应于 braid group / mapping class group 的一个元素，从而在零模子空间上给出一个拓扑不变量（例如 Ising 的 half twist / Dehn twist 矩阵）。简而言之就是，$F=0$ => holonomy 给出辫子群表示。

这一理念可以由非阿贝尔stokes解释，也就是下文重点使用的方法：

联络的曲率 $F$ 在考虑的区域上严格为零（平直），则沿闭合曲线的平行搬运（holonomy）只依赖于曲线的同伦类；记联络为 $A$，闭合曲线 $\gamma$ 的 holonomy 写为
$$
U[\gamma]=\mathcal P\exp\Big(\oint_\gamma A\Big).
$$
利用非阿贝尔 Stokes（面次序指数），可以把它写成所围曲面 $S$ 上曲率的面次序指数：
$$
U[\gamma]=\mathcal P_S\exp\Big(\iint_S F\Big).
$$
若 $F\equiv0$，右侧对任意两个以同一同伦类为边界的填充面均相同，故 $U[\gamma]$ 仅取决于 $[\gamma]$（同伦类）。在参数空间或带标记点的配置空间情形，基本群 $\pi_1$ 即为编织群（braid group）；在带洞的曲面情形，相关的是映射类群（mapping class group）。因此平直联络直接给出这些群在零模子空间上的表示：每个同伦类对应一个幺正算符，并且该算符由拓扑资料确定，成为拓扑不变量。

若曲率不是严格为零但在该子空间上为标量倍的单位算子（即 $F(x)=\lambda(x)I$，称为中心化曲率），则面次序指数退化为一个标量因子乘以单位算符，holonomy 仍由同伦类决定，但只给出群的投影表示（projective representation），这在带整体相位不可观测时很常见。

直观举例：在 Ising‑型或任意子系统中，若在有效基态子空间上联接近似平直，则绕粒子交换的闭环给出的 $U[\gamma]$ 正是 braid group 的矩阵表示；在曲面上对基圈做 Dehn twist，得到的矩阵即映射类群的表示（half‑twist / Dehn twist 矩阵），这些矩阵在联接平直或曲率中心化时只由拓扑数据决定。

- 交换两个 $\sigma$ 任何子（按融合基底 ${1,\psi}$）的 Braid 矩阵为对角形式：

$$
B_{\sigma\sigma}=\mathrm{diag}(R^{\sigma\sigma}_1,\;R^{\sigma\sigma}_\psi)=\mathrm{diag}(e^{-i\pi/8},\;e^{3i\pi/8}).
$$

​    这里 $R^{\sigma\sigma}_a$ 是交换两个 $\sigma$ 后若它们融合到通道 $a\in\{1,\psi\}$ 得到的相因子。

- Dehn twist / half‑twist（绕单个任何子做 $2\pi$ 自旋）由拓扑自旋 $\theta_a$ 给出：例如

$$
\theta_\sigma=e^{i\pi/8},\qquad \theta_\psi=-1,\qquad \theta_1=1.
$$

​    因此在对应的零模子空间上，绕包含某个任何子的圈做 Dehn twist 会乘以相应的 $\theta$（若作用在局域标记的子空间上则为标量因子），这就是映射类群MCG元素在该子空间上的具体实现。





(下面的思路未经数值模拟纯脑补)

引入 \(R=e^{iH_P}\) 之后，自然的下一步想法可以是围绕某个“理想可积点” \(H_P^{(0)}\) 做小扰动(我们前面写成了三个部分分解)：
$$
H_P = H_P^{(0)} + \varepsilon V,
$$
其中 \(H_P^{(0)}\) 满足量子 YBE（或至少经典 YBE），\(V\) 是一般的 su(4) 扰动。

- 经典 YBE 提供了一个对 \(V\) 的线性化条件：哪些方向 \(V\) 保持 \(O(\varepsilon^2)\) 曲率为零，哪些方向立即产生非零曲率；
- 对应“哪些扰动只改变有效 (t,\Delta,\mu) 而不破坏拓扑/平坦结构，哪些扰动会引入相互作用/非对易结构，从而拉高实现同一个 Dehn twist 所需的电路复杂度”。

### 3.时间演化braid$\rightarrow$空间变换



> **命题** 在一个带谱隙的 2+1 维拓扑相中，给定任意 braid word $\beta\in B_n$，存在一族 Hamiltonian‑level 的纯空间操作（由 Dehn twist / 几何畸变 / 局域 $R=e^{iH_P}$ 门重构给出），使得在零模/编码子空间上，其非阿贝尔 Berry holonomy 与沿时间编织任意子的 holonomy 等价（在 SU(2) 内共轭，误差由曲率/YBE 偏差、非绝对慢走与离散化精度控制）。

即：沿着连续路径演化的H(λ),例如通过绝热/ quasi‑adiabatic 演化实现，现在想要证明这些Braid的操作可以用空间的变换来表示。



下面依次在三个层次上论证：

1. 纯拓扑：$B_n \cong \mathrm{MCG}(D_n,\partial D)$，且 $\mathrm{MCG}$ 由 Dehn twist 生成；
2. R/Berry：在 $R = exp(iH\_P)$​ 框架中构造空间操作并用 Berry 几何控制误差。
3. ~~TQFT：世界线 braid 表示与 mapping class 表示一致~~(这一部分是查的，未经实际推导验证，放在附录了)



#### 拓扑证明：辫群与 Dehn twist 生成的 mapping class group 等价

首先是MCG（Mapping Class Group）：一张给定拓扑表面 Σ（例如带 n 个标记点的圆盘 D_n）上的映射类群，记作 MCG(Σ)，定义为把保持边界和标记点的自同胚按同伦/同胚同类（isotopy）分掉后的群。直观上它是“把表面通过切‑扭‑粘变形”的同伦类集合. 用具体公式描述如下：

记 $D\subset\mathbb R^2$ 为闭圆盘，$\{p_1,\dots,p_n\}\subset\mathrm{int}(D)$ 为 n 个标记点。我们考虑保持边界点集 $\partial D$ 不动的同胚的同伦类：
$$
\mathrm{MCG}(D_n,\partial D)
 := \pi_0\bigl(\{f:D\to D\mid f\text{为同胚},\ f(\partial D)=\partial D,\ f(\{p_i\})=\{p_i\}\}\bigr).
$$

另一方面，n 股辫群 $B_n$ 定义为 n 个点在区间 $[0,1]$ 上的编织同伦类。



**定理 1.1.** 有自然同构
$$
\Phi: B_n \simeq  \mathrm{MCG}(D_n,\partial D).
$$

可以这样理解：

- 从一个 braid worldline 配置出发，把时间方向压缩，将其投影到空间圆盘上，得到一个保持边界的同胚的同伦类；
- 反过来，从一个 mapping class 出发，可构造等价的“点在时间中运动”的轨迹；
- 这两个构造互为逆，且群结构（连接、复合）对应。

因此，拓扑上任意 braid word $\beta\in B_n$ 都等价于某个 mapping class $[f]\in \mathrm{MCG}(D_n,\partial D)$。



接下来考虑 mapping class group 的生成元。

**定义 1.2（Dehn twist）.** 对一条简单闭曲线 $c\subset D_n$，沿着 $c$ 取一个窄带邻域，将其切开，在一侧沿曲线方向旋转 $2\pi$ 后再粘回，得到新的曲面同胚。其同伦类记为 $T_c\in\mathrm{MCG}(D_n,\partial D)$，称为沿 $c$ 的 Dehn twist。



可以写出形式如下： 设一个环面 $A = S^1 \times [0,1]$, 映射$T: A \rightarrow A$
$$
T(\theta,t) = (\theta+2\pi,t)
$$
Dehn twist写为：
$$
T: (\alpha,\beta) \rightarrow (\alpha+\beta,\beta)
$$

Half twist:=
$$
\sqrt T := (θ , t ) \rightarrow (θ + π + πt, t ).
$$
以下图片来自论文$^{4}$ 

![image-20260413031213485](/home/asice-cloud/.config/Typora/typora-user-images/image-20260413031213485.png)



**定理 1.3（Dehn–Lickorish 生成定理的球面多穿孔版本）.** 对于多穿孔圆盘（或等价的球面带标记点），存在有限族简单闭曲线 $\{c_j\}$，使得
$$
\mathrm{MCG}(D_n,\partial D)
 = \langle T_{c_j}^{\pm1}\mid j\in J\rangle
$$
即**任意 mapping class 都可写成有限个 Dehn twist 的乘积**。换言之，Dehn twist就是MCG生成元

需要补充的是：Dehn twist并不对应辫子群的生成元(或者是我们之前构造的R)，与之对应的是half twist.这一点可以由图片看出

1. **局域交换/绕行生成元（half twist 类型）**：
    $$
    X_{ab}(\theta)=\frac{\theta}{2}\,\gamma_a\gamma_b,\qquad U_{ab}(\theta)=e^{X_{ab}(\theta)}.
    $$

    - 在 Ising 任意子理论中，$\theta=\pi/2$ 对应标准的 $\sigma$‑$\sigma$ 交换（up to phase）。
    - 这类幺正正是“half twist”的 Majorana 表示：把两个端点/缺陷在配置空间中交换一次。 

2. **沿某条非平凡曲线的 Dehn twist 生成元**：

    - 把 Dehn twist $T_\gamma$ 看成沿某闭合曲线 $\gamma$ 绕行一圈的 mapping class 元；

    - 在 Majorana 表示中，$T_\gamma$ 的作用往往可以写成若干 $X_{ab}$ 的组合：
        $$
        \rho(T_\gamma)=\exp\Big(\frac12\sum_{a<b}\Theta^{(\gamma)}_{ab}\,\gamma_a\gamma_b\Big).
        $$

    - 这里 $\Theta^{(\gamma)}$ 是沿配置空间路径的联络积分：
        $$
        \Theta^{(\gamma)}_{ab}=\oint_\gamma\Omega_{ab}(\lambda)\,d\lambda.
        $$



**同构 $\Phi$ 与生成定理合并即得**：

> **推论 1.4（纯拓扑版主命题）.** 对任意 braid word $\beta\in B_n$，存在简单闭曲线 $c_1,\dots,c_m$ 及符号 $\varepsilon_k=\pm1$，使得
> $$
> \Phi(\beta) = [f] = T_{c_1}^{\varepsilon_1}\cdots T_{c_m}^{\varepsilon_m}\in\mathrm{MCG}(D_n,\partial D).
> $$
> 换言之，每一个 braid 都等价于有限个 Dehn twist 的组合，而每个 Dehn twist 是一个**纯空间几何操作**。

这一步纯拓扑论证,仅说明“时间中的编织”和“空间中的扭转”在 mapping class group 意义下是同一个对象的不同代表。





#### $R = exp(iH\_P)$ / Berry ：从抽象表示到具体 Hamiltonian 空间操作

借助TQFT的理论：给定一个 2+1 维拓扑量子场论，在一个带 n 个穿孔的空间截面上，其 Hilbert 空间记为
$$
\mathcal H_0(D_n; \{\sigma_i\}),
$$
其中 $\sigma_i$ 是附着在各穿孔上的 anyon 类型（例如 Ising TQFT 中的 $\sigma$）。



我们现在在具体的晶格拓扑相，实现 TQFT 的 Hilbert 空间与幺正表示。假设：

1. 有一个局域自旋/费米体系，其低能有效理论实现了某个拓扑序，对应的零模/编码子空间为 $\mathcal H_0$；

2. 在每条键 $\langle ij\rangle$ 上有一个两体生成元
$$
     H_P^{(ij)} = \sum_{a,b} c_{ab}^{(ij)}\,\sigma_a\otimes\sigma_b,\qquad R_{ij}=e^{iH_P^{(ij)}};
$$

3. 总哈密顿量 $H(\lambda)$ 由这些 $H_P^{(ij)}(\lambda)$ 组成，参数 $\lambda$ 取值在某个带谱隙的参数空间 $\Omega$ 内。

这里的哈密顿量有两层物理意义：

- “上层的空间”：$H(\lambda)$ 只是一个具体晶格模型，用来实现给定的 2+1 维 TQFT。它的零模/编码子空间 
$\mathcal H_0(\lambda)$ 携带了前面讨论的 braid/MCG 幺正表示；
- “具体形式”：每个局域块 $H_P^{(ij)}$ 描述的是微观上的配对、跳跃、规范约束等局域物理。我们在前面它写成
$$
H_P^{(ij)} \simeq H_{\mathrm{quad}}^{(ij)} + H_{\mathrm{int}}^{(ij)} + H_{\mathrm{gauge}}^{(ij)},
$$

这里采用这种自然分解，$H_{\mathrm{quad}}$ 选出 Kitaev/p+ip/honeycomb 那个拓扑通道，$H_{\mathrm{int}}, H_{\mathrm{gauge}}$ 则是偏离理想可积/YBE 点的局域扰动，它们用来定义 $\varepsilon_{\mathrm{YBE}}$、复杂度曲率等。



前面的第二部分已经说明：

- 每个 $H_P^{(ij)}$ 在 JW/Majorana 映射后自然分解为:
$H_{\mathrm{quad}} + H_{\mathrm{int}} + H_{\mathrm{gauge}}$，前者给出 BdG 拓扑结构，后者描述偏离可积/平坦的扰动；

- YBE / classical YBE 在 $H_P$ 空间中刻画出可积/近平坦子流形，使 Berry 曲率 $F$ 小，保证 holonomy 对路径同伦“刚性”。

对参数空间 $\Omega$ 中的路径 $\gamma$，其 Berry holonomy 定义为
$$
U[\gamma] = \mathcal P\exp\Bigl(-\int_\gamma A\Bigr),\qquad
 A_\mu(X) = iP\,\partial_\mu P\,P,\qquad
$$
（$P$ 为投影到 $\mathcal H_0(\lambda)$）。在配置空间的语境下，$\gamma$ 可以理解为：

- A. 随时间移动缺陷/序参量纹理的路径（worldline 图像）；
- B. 在空间上重写配对图 / Kekulé 畸变等空间操作。



##### 3.0命题

这里将要证明以下不等式：(系数吸收进了$\lesssim$)
$$
\bigl\|U_{\mathrm{spatial}} - e^{i\phi}U_{\mathrm{top}}\bigr\|
\;\lesssim\; \varepsilon_{\mathrm{YBE}} + \kappa\,\mathcal A_{\gamma,\gamma'} + \varepsilon_{\mathrm{adiab}} + \varepsilon_{\mathrm{Trotter}},
$$
其中：

- $U_{\mathrm{top}}$ 对应抽象 TQFT 中的 $\rho_{\mathrm{WL}}(\beta)$ 或 $\rho_{\mathrm{MCG}}([f])$；
- $U_{\mathrm{spatial}}$ 是由有限轮局域 $R=e^{iH_P}$ 门、配对图重写或 Kekulé/Dirac 畸变构成的空间操作；
- $\varepsilon_{\mathrm{YBE}}$ 测量所选 $H_P$ 偏离可积/平坦子流形的程度；
- $\kappa$ 为 Berry 曲率上界，$\mathcal A_{\gamma,\gamma'}$ 为两条同伦路径间“填充面积”；
- $\varepsilon_{\mathrm{adiab}}$ 来自时间演化不完全 adiabatic；
- $\varepsilon_{\mathrm{Trotter}}$ 来自用离散 R 门近似连续生成元的误差。



这一不等式的证明骨架依赖于四个引理：

1. **引理 A（平坦/近平坦联络下的路径同伦稳定性）**：在 $\|F\|\le\kappa$ 区域内，同伦路径给出 holonomy 之差 $\sim\kappa\times$ 填充面积；
2. **引理 B（adiabatic 时间演化与 Berry holonomy 等价）**：在带谱隙的情形，真实时间演化与 Berry holonomy 只差 $\varepsilon_{\mathrm{adiab}}$；
3. **引理 C（F,R‑word 与 $e^{iH_P}$ 嵌入一致性）**：在理想可积点上，$e^{iH_P^{(0)}}$ 在 $\mathcal H_0$ 上等同于 TQFT 的 R‑符号，小扰动下拓扑表示稳定；
4. **引理 D（空间重构电路逼近 quasi‑adiabatic 演化）**：quasi‑adiabatic 生成元可由有限轮局域 R 门逼近，误差 $\varepsilon_{\mathrm{Trotter}}$ 可控。



下面将上述四个引理的内容和主不等式的拼接过程展开说明。

##### 3.1 设定与记号

令参数/配置空间为光滑流形 $\Omega$，对每个 $X\in\Omega$ 有带谱隙的哈密顿量 $H(X)$，其零模/编码子空间投影为
$$
P(X): \mathcal H \to \mathcal H_0(X),\qquad P^2=P=P^\dagger.
$$
在 $\Omega$ 上定义 Berry 联络与曲率
$$
 A_\mu(X) = iP\,\partial_\mu P\,P,\qquad
F_{\mu\nu}(X) = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu]
                         = P[\partial_\mu P,\partial_\nu P]P.
$$
对任意光滑路径 $\gamma:[0,1]\to\Omega$，Berry holonomy 定义为
$$
U[\gamma] = \mathcal P\exp\Bigl(-\int_0^1 A_\mu(X(t))\,\dot X^\mu(t)\,\mathrm dt\Bigr).
$$
我们假设在考虑的区域内存在统一界
$$
\|F_{\mu\nu}(X)\| \le \kappa,\qquad \forall X\in \Omega.
$$

设 $\gamma$ 是“理想世界线”路径（anyon 真实编织或 TQFT 抽象路径），$\gamma'$ 是通过空间重构/Kekulé/Dirac 畸变等实现的“空间操作路径”。二者端点相同，且在保持谱隙的区域内同伦。记 $S$ 为两条路径围成的有向曲面，其面积为 $\mathcal A_{\gamma,\gamma'}$。

Wilson 环：
$$
U[\partial S] = \mathcal P\exp\Bigl(-\oint_{\partial S}A\Bigl).
$$
$U[γ] $是任意路径的路径有序指数；$U[∂S] $​是当 γ 为某面 S 的边界时的特殊情况。若有两条端点相同的路径 γ 和 γ'，拼成的闭合路径就是 ∂S，且$U[∂S]=U[γ′] U[γ]^{−1}$

因此若 $U[∂S]=I$（例如曲率 F 在该面上为 0），则 $U[γ′]=U[γ]$，即 holonomy 仅依照同伦类定。



##### 3.2 引理 A：近平坦联络下的路径同伦稳定性

**引理 A.** 若 $\|F\|\le\kappa$，则存在常数 $C_A$，使得对任意同伦的两条路径 $\gamma,\gamma'$，有
$$
\bigl\|U[\gamma']-U[\gamma]\bigr\| \le C_A\,\kappa\,\mathcal A_{\gamma,\gamma'}.
$$

*证明* 这是非阿贝尔 Stokes 公式$^{[3]}$的一个定量化版本。

- 由 Berry 联络的 Wilson 线定义和 path‑ordering，可将 $U[\gamma']U[\gamma]^{-1}$ 写成沿着由 $\gamma$、$\gamma'$ 和若干“连接线”组成的闭合回路 $\partial S$ 的 Wilson 环；

非阿贝尔 Stokes 定理在此可写作更具体的表面序形式：取基点 $x_0\in S$ 并定义把 $x_0$ 平行传输到 $x$ 的 Wilson 线
$$
U(x,x_0)=\mathcal P\exp\Bigl(-\int_{x_0}^x A\Bigr).
$$
那么曲率 $F=dA+A\wedge A$ 的表面序表示为
$$
U[\partial S]=\mathcal P_S\exp\Bigl(-\iint_S U(x,x_0)^{-1}F(x)U(x,x_0)\,dS_x\Bigr),
$$
其中 $\mathcal P_S$ 表示对面元按某种固定顺序進行排列（离散化为小面元后按该顺序相乘并取极限）。

离散化的直观近似为：把 $S$ 划分为小面元 $\{\Delta S_k\}$（中心点 $x_k$），选参考路径 $p_k$ 把基点连到 $x_k$，则
$$
U[\partial S]\approx\prod_k^{\mathcal P_S}\exp\big(-U(p_k)^{-1}F(x_k)U(p_k)\,\Delta S_k\big).
$$
对表面有序指数做 Dyson 级数展开（令 $\tilde F(x)=U(x,x_0)^{-1}F(x)U(x,x_0)$），前两项为
$$
U[\partial S]=I -\iint_S \tilde F(x)\,dS_x - \iint_{S_>}\tilde F(x_1)\tilde F(x_2)\,dS_{x_1}dS_{x_2} + O(\mathrm{Area}^3),
$$
其中 $S_>$ 表示按表面序的有序二重积分；高阶项包含非对易的嵌套积分与对易子。



取算符范数并假设 $\sup_{x\in S}\|F(x)\|\le\kappa$，可得主阶估计
$$
\|U[\partial S]-I\| \le C\,\kappa\,\mathrm{Area}(S) + O(\kappa^2\mathrm{Area}(S)^2),
$$
因此在小面积或曲率有界情況下有
$$
\|U[\partial S]-I\| \lesssim C \cdot \kappa\,\mathrm{Area}(S).
$$
由于 $U[\partial S]=U[\gamma']U[\gamma]^{-1}$，因为$U[\gamma]$是酉的，
$$
U[γ'] − U[γ] = (U[∂S] U[γ]) − U[γ] = (U[∂S] − I) U[γ]\\
\\
\|U[\gamma']-U[\gamma]\| = \|U[\partial S]-I\| \lesssim C_A  \kappa\,\mathcal A_{\gamma,\gamma'},
$$
从而导出引理 A 中所需的不等式

形象的说明，是在曲率有界的近平坦联络下，同伦的两条路径给出的 holonomy 之差被“曲率上界 × 填充面积”控制，从而确保在 YBE/近平坦子流形附近，空间变形路径与世界线路径的 Berry 相位/矩阵几乎相同。



##### 3.3 引理 B：adiabatic 时间演化与 Berry holonomy

**引理 B.** 设 $H(X(t))$ 是一条光滑、谱隙 $\Delta>0$ 统一有界的时间参数路径，对应的时间演化算符为 $U_{\mathrm{dyn}}(t)$。若演化足够慢（总演化时间 $T$ 足够大），则存在常数 $C_B$ 使得
$$
\bigl\|U_{\mathrm{dyn}}(T)P(X(0)) - U[\gamma]P(X(0))\bigr\| \le C_B\,\varepsilon_{\mathrm{adiab}},
$$
其中 $\varepsilon_{\mathrm{adiab}} \sim \mathcal O(\|\dot H\|/\Delta^2)$ 之类的 adiabatic 小参数。

这是 adiabatic 定理的标准形式（Kato、Avron–Seiler–Yaffe 等）：

- 利用瞬时本征态展开，将时间演化分解为 dynamical 相位与几何相位；
- adiabatic 极限下，跃迁到激发态的幅度受 $\|\dot H\|/\Delta^2$ 控制；
- 限制到初态在基态/零模子空间的情形，激发态贡献可忽略，剩下的几何部分正是 Berry holonomy，误差给出上式中的 $\varepsilon_{\mathrm{adiab}}$。

因此，真实的绝热时间演化在编码子空间上与 Berry holonomy 等价，只差一个可控小量。



##### 3.4 引理 C：F,R‑数据与 $e^{iH_P}$ 嵌入的一致性

**引理 C.** 在理想可积点 $H_P^{(0)}$ 上，对每一条键 $(ij)$，存在 $H_P^{(0,ij)}$ 使得
$$
R_{ij}^{(0)} := e^{iH_P^{(0,ij)}}\big|_{\mathcal H_0}
$$
在编码子空间 $\mathcal H_0$ 上实现 TQFT 给出的 R‑矩阵；同时，若 $H_P^{(ij)}=H_P^{(0,ij)}+V^{(ij)}$ 的扰动满足 $\|V^{(ij)}\|\le\epsilon$ 且沿整条路径谱隙保持，则对应的 Berry holonomy 与理想 TQFT 表示相差至多
$$
\varepsilon_{\mathrm{YBE}} \lesssim C_C\,\epsilon
$$
同样可以使用非阿贝尔stokes定理证明。证明太长放在附录



##### 3.5 引理 D：空间重构电路逼近 quasi‑adiabatic 演化

**引理 D.** 设存在一个 quasi‑adiabatic 生成元 $K(t)$，它是由局域算符指数衰减加权积分得到的（Hastings–Wen 方案），使得
$$
U_{\mathrm{QA}} = \mathcal T\exp\Bigl(-i\int_0^T K(t)\,\mathrm dt\Bigr)
$$
在编码子空间上实现沿路径 $\gamma'$ 的 adiabatic/Berry 演化。则对任意给定精度 $\delta>0$，可以构造有限轮局域 R‑门电路 $U_{\mathrm{spatial}}$，使得
$$
\bigl\|U_{\mathrm{spatial}}-U_{\mathrm{QA}}\bigr\| \le \varepsilon_{\mathrm{Trotter}} \le \delta.
$$

*证明思路. *

用一个平滑的滤波核把投影导数卷进去，得到 quasi‑adiabatic 生成元 K(t)。直观上它是“平滑版”的瞬时几何生成器，用来把基态子空间沿路径搬运，而不会把能量推到激发态

- quasi‑adiabatic 生成元 $K(t)$ 是局域和的积分，满足 Lieb–Robinson 速度界，因而可在空间上分块、在时间上离散；
- 用 Trotter–Suzuki 分解将时间序列上的指数拆分为各个局域块的指数乘积；
- 每一块指数 $e^{-i h_\alpha \Delta t}$ 可由有限轮局域 $R=e^{iH_P}$ 门近似实现（通过门集完备性或直接识别）；
- 误差由标准 Trotter 估计控制，随时间步长 $\Delta t$ 和分块截断尺度衰减，从而给出 $\varepsilon_{\mathrm{Trotter}}$。

这一步把“连续的 quasi‑adiabatic 演化”真正压缩为“有限深度、有限范围的空间重构电路”，具体形式可以是配对图重写、Kekulé 模式拨动、branch cut/涡旋线段的重新排布等。



##### 3.6 主不等式的拼接

最后说明如何由上述四个引理得到主不等式。分几步进行：

1. **理想拓扑门 vs 理想 Berry holonomy.**

- 在 TQFT 层，有
         $$U_{\mathrm{top}} \simeq \rho_{\mathrm{WL}}(\beta) = \rho_{\mathrm{MCG}}(\Phi(\beta)), $$
         其中 $\simeq$ 忽略整体相位；

- 在理想可积点，由引理 C 和 B，可在某条“世界线路径” $\gamma$ 上构造理想哈密顿量 $H^{(0)}(X)$，其 adiabatic/Berry 演化 $U[\gamma]$ 在编码子空间上等同于 $U_{\mathrm{top}}$，至多差一个 $\varepsilon_{\mathrm{YBE}}+\varepsilon_{\mathrm{adiab}}$ 的小量。

2. **世界线路径 vs 空间路径.**

- 在参数/配置空间中构造对应的空间操作路径 $\gamma'$（用 Dehn twist、配对图重写、Kekulé/Dirac 纹理变形等实现），保证 $\gamma$ 与 $\gamma'$ 在保持谱隙的区域内同伦；

- 由引理 A，二者 Berry holonomy 差异满足
    $$
         \bigl\|U[\gamma']-U[\gamma]\bigr\| \lesssim \kappa\,\mathcal A_{\gamma,\gamma'}.
    $$

3. **理想 Berry vs quasi‑adiabatic/电路实现.**

- 由引理 B 和 quasi‑adiabatic 构造，可用 $U_{\mathrm{QA}}$ 逼近 $U[\gamma']$，误差吸收到 $\varepsilon_{\mathrm{adiab}}$；

- 再由引理 D，用有限轮 R‑门电路 $U_{\mathrm{spatial}}$ 逼近 $U_{\mathrm{QA}}$，误差记为 $\varepsilon_{\mathrm{Trotter}}$。

4. **三角不等式合并误差.** 将上述三步的偏差用三角不等式串联：
    $$
    \begin{aligned}
         \bigl\|U_{\mathrm{spatial}}-e^{i\phi}U_{\mathrm{top}}\bigr\|
         &\le \bigl\|U_{\mathrm{spatial}}-U_{\mathrm{QA}}\bigr\|
                         + \bigl\|U_{\mathrm{QA}}-U[\gamma']\bigr\|
                         + \bigl\|U[\gamma']-U[\gamma]\bigr\|\\
         &\quad + \bigl\|U[\gamma]-e^{i\phi}U_{\mathrm{top}}\bigr\|.
        \end{aligned}
    $$
    

    下面给出每一项的可取上界（常数为模型/几何常数，已在正文各处定义）：

    Trotter / 电路误差（引理 D）：令 quasi‑adiabatic 生成元被分解为局域项且每项范数上界为 $k_*$，总演化时间 $T$，Trotter 步数 $m$，空间截断半径 $R_k$，则存在常数 $C_D,C'_D,\mu>0$ 使得
    $$
            \bigl\|U_{\mathrm{spatial}}-U_{\mathrm{QA}}\bigr\| \le C_D\frac{T\,k_*^2}{m} + C'_D e^{-\mu R_k}.
    $$
    

    - 绝热 / quasi‑adiabatic 误差（引理 B）：令 $\Delta$ 为沿路径的最小谱隙，记 $L_0$ 为基态能量相关项尺度，$L_V$ 为扰动导数尺度，则存在常数 $C_B,C''_B$ 使得

    $$
            \bigl\|U_{\mathrm{QA}}-U[\gamma']\bigr\| \le C_B\frac{L_0+L_V}{\Delta^2} + C''_B\frac{\epsilon_{\mathrm{tot}}}{\Delta}.
    $$

    

    - 曲率 / 同伦面积项（引理 A，Stokes 精确界）：令曲面 $S$ 为 $\gamma,\gamma'$ 所围曲面、面面积为 $\mathcal A_{\gamma,\gamma'}$，并设 $\|F\|\le\kappa$，则存在几何常数 $C_A,C_A'>0$ 使得

    $$
            \bigl\|U[\gamma']-U[\gamma]\bigr\| \le C_A\,\kappa\,\mathcal A_{\gamma,\gamma'}\;\exp\Big(C_A'\,\kappa\,\mathcal A_{\gamma,\gamma'}\Big).
    $$

    

    在小曲率/小面积近似（$\kappa\,\mathcal A_{\gamma,\gamma'}\ll1$）下，可线性化为 $\lesssim C_A\,\kappa\,\mathcal A_{\gamma,\gamma'}$。

    

    - YBE / R‑嵌入偏差（引理 C）：令路径长度（或门数）为 $L$，扰动范数尺度 $\epsilon_{\mathrm{tot}}=\sum_{\alpha}\epsilon_{\alpha}$，则存在常数 $C_C,C'_C$ 使得

    $$
            \bigl\|U[\gamma]-e^{i\phi}U_{\mathrm{top}}\bigr\| \le C_C\,L\,\frac{\epsilon_{\mathrm{tot}}}{\Delta} + C'_C\frac{\epsilon_{\mathrm{tot}}^2}{\Delta^2}.
    $$

    

    把上述不等式代入三角不等式并合并常数，可得到一个带显式常数的主不等式：
    $$
    \begin{aligned}
        \|U_{\mathrm{spatial}}-e^{i\phi}U_{\mathrm{top}}\| &\le C_D\frac{T\,k_*^2}{m} + C'_D e^{-\mu R_k} \\
        &\quad + C_B\frac{L_0+L_V}{\Delta^2} + C''_B\frac{\epsilon_{\mathrm{tot}}}{\Delta} \\
        &\quad + C_A\,\kappa\,\mathcal A_{\gamma,\gamma'}\,e^{C_A'\,\kappa\,\mathcal A_{\gamma,\gamma'}} \\
        &\quad + C_C\,L\,\frac{\epsilon_{\mathrm{tot}}}{\Delta} + C'_C\frac{\epsilon_{\mathrm{tot}}^2}{\Delta^2}.
        \end{aligned}
    $$
    

    注：如果在应用情形中满足 $\kappa\,\mathcal A_{\gamma,\gamma'}\ll1$ 且 $\epsilon_{\mathrm{tot}}/\Delta\ll1$，可忽略指数项和二阶项，主不等式退化为文中先前的简洁形式（把所有模型常数合并为“≲”符号）。这就得到开头给出的主不等式形式。



##### 3.7总结

把纯拓扑层的 $B_n\cong\mathrm{MCG}(D_n)$、TQFT 层的 $\rho_{\mathrm{WL}}=\rho_{\mathrm{MCG}}\circ\Phi$ 与上述 R/Berry 层的误差估计结合，即得本文命题，用形象一点描述：

> 任意 braid $\beta$ 在 TQFT 中对应的拓扑门 $U_{\mathrm{top}}$ 都可以在具体晶格拓扑相中通过一族**几何/空间操作**实现：包括 Dehn twist、branch cut/Kekulé/Dirac 涡旋畸变、配对图重写以及局域 $R = exp(iH\_P)$ 门电路。只要所选 $H_P$ 方向处在 YBE/近平坦子流形附近，且操作过程保持谱隙不闭合，这些空间操作在零模/编码子空间上与“时间中的 braid”给出同一个 SU(2) 共轭类的幺正，误差由上式中的四项可控。



### reference

[1]*The Dynamical Yang-Baxter Equation, Representation Theory, and Quantum Integrable Systems*

[2]*nLab https://ncatlab.org/nlab/show/Yang-Baxter+equation* 

[3]Non-Abelian Stokes Theorem . December 1995  .DOI:10.1142/9789812831323_0017

[4]Instantaneous braids and Dehn twists in topologically ordered states.  Guanyu Zhu , Ali Lavasani, and Maissam Barkeshli. DOI: 10.1103/PhysRevB.102.075105



最后一章的定理的所有出处：

- **辫群、mapping class group 与 Dehn twist 生成**  
    - J. S. Birman, *Braids, Links, and Mapping Class Groups*, Annals of Mathematics Studies 82, Princeton Univ. Press (1974).  
    - B. Farb and D. Margalit, *A Primer on Mapping Class Groups*, Princeton Univ. Press (2011)

- **Reshetikhin–Turaev 型 TQFT 与世界线/MCG 表示**  
    - N. Reshetikhin and V. G. Turaev, *Invariants of 3-manifolds via link polynomials and quantum groups*, Invent. Math. 103, 547–597 (1991).  
    - V. G. Turaev, *Quantum Invariants of Knots and 3-Manifolds*, de Gruyter (1994/2010).  
    - B. Bakalov and A. Kirillov Jr., *Lectures on Tensor Categories and Modular Functors*, AMS (2001)，关于模范畴、F/R‑符号与 MCG 表示的统一框架。

- **Berry 相、非阿贝尔 Berry 联络与几何解释**  
    - B. Simon, *Holonomy, the quantum adiabatic theorem, and Berry's phase*, Phys. Rev. Lett. 51, 2167–2170 (1983).  
    - F. Wilczek and A. Zee, *Appearance of Gauge Structure in Simple Dynamical Systems*, Phys. Rev. Lett. 52, 2111–2114 (1984)，给出非阿贝尔 Berry 联络与曲率的标准形式。

- **绝热定理与编码子空间上的 adiabatic 演化**  
    - T. Kato, *On the adiabatic theorem of quantum mechanics*, J. Phys. Soc. Jpn. 5, 435–439 (1950).  
    - J. E. Avron, R. Seiler and L. G. Yaffe, *Adiabatic theorems and applications to the quantum Hall effect*, Commun. Math. Phys. 110, 33–49 (1987).  
    - S. Teufel, *Adiabatic Perturbation Theory in Quantum Dynamics*, Lecture Notes in Mathematics 1821, Springer (2003)，系统讨论了有谱隙情形下的绝热定理和误差估计。

- **拓扑序稳定性与局域扰动下的鲁棒性**  
    - M. B. Hastings and X.-G. Wen, *Quasi-adiabatic continuation of quantum states: The stability of topological ground-state degeneracy and emergent gauge invariance*, Phys. Rev. B 72, 045141 (2005).  
    
    - S. Bravyi, M. B. Hastings and S. Michalakis, *Topological quantum order: stability under local perturbations*, J. Math. Phys. 51, 093512 (2010).  

    - S. Michalakis and J. P. Zwolak, *Stability of frustration-free Hamiltonians*, Commun. Math. Phys. 322, 277–302 (2013)。
    
        上面是引理 C 中“沿可积/YBE 子流形的小扰动不破坏拓扑表示”基础。
    
- **Lieb–Robinson 界与 quasi‑adiabatic 连续性、有限深电路逼近**  
    - E. H. Lieb and D. W. Robinson, *The finite group velocity of quantum spin systems*, Commun. Math. Phys. 28, 251–257 (1972)——Lieb–Robinson 界最初形式。  
    - M. B. Hastings, *Lieb-Schultz-Mattis in higher dimensions*, Phys. Rev. B 69, 104431 (2004)；以及与 X.-G. Wen 合作的 2005 年工作，上述 quasi‑adiabatic 生成元 $K(t)$ 的构造即源于此。  
    - B. Nachtergaele and R. Sims, *Lieb-Robinson bounds in quantum many-body physics*, in *Entropy and the Quantum*, Contemp. Math. 529, AMS (2010)，对 Lieb–Robinson 界及其在 quasi‑adiabatic 连续性中的应用做了综述。  
        这些是引理 D 中“将连续演化压缩为有限深局域电路”的步骤。

### Appendix



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



#### 可积子流形部分补充

**(a) 从平行移动到曲率 $F=dA+A\wedge A$**

在配置空间的某个局部坐标片中，用坐标 $x^\mu$ 描述点。零模子空间上的 Berry 平行移动算符 $U(t)$ 满足一阶微分方程
$$
\frac{\mathrm d}{\mathrm dt}U(t)
=-A_\mu\big(x(t)\big)\,\dot x^\mu(t)\,U(t),
$$
其中 $A=A_\mu(x)\,\mathrm d x^\mu$ 是 $\mathfrak u(N)$‑值 1‑形式。形式解写作路径有序指数
$$
U[\gamma]=\mathcal P\exp\Big(-\int_\gamma A\Big).
$$

取二维坐标 $(x,y)$，考虑以 $x_0$ 为左下角、边长 $\delta x,\delta y$ 的无穷小矩形回路 $\square$。沿每条边的平行移动近似为
$$
U_x(\delta x)\approx\mathbf1-A_x(x_0)\delta x,\qquad
U_y(\delta y)\approx\mathbf1-A_y(x_0)\delta y,
$$
其中 $A_x,A_y$ 是 $A$ 在 $x,y$ 方向的分量。总 holonomy 是四条边有序乘积
$$
U_{\square}
\approx U_x(\delta x)\,U_y(\delta y)\,U_x(-\delta x)\,U_y(-\delta y).
$$
把 $A_x,A_y$ 在 $x_0$ 附近展开，保留到 $\mathcal O(\delta x\,\delta y)$，可以得到标准结果：
$$
U_{\square}
\approx \mathbf1 - F_{xy}(x_0)\,\delta x\,\delta y,
$$
其中
$$
F_{xy}=\partial_xA_y-\partial_yA_x+[A_x,A_y].
$$
用微分形式记号写成
$$
F=dA+A\wedge A.
$$
也就是说，$F$ 正是“无穷小回路 holonomy 偏离单位的首要项”。这个定义在规范变换下满足 $F\to gFg^{-1}$ 的协变性，是刻画联络局部曲率的自然对象。

**(b) $F=0$ 时 Berry 联络是纯规范，holonomy 给出辫子群表示**

设 $U\subset\mathcal C$ 是一个单连通开集，若在 $U$ 上曲率恒为零：
$$
F|_U\equiv 0,
$$
则有标准结论：存在光滑幺正矩阵 $g:U\to U(N)$，使得
$$
A=g^{-1}\,\mathrm d g,
$$
即 $A$ 是一个纯规范联络。

证明思路（略述）：选定参考点 $x_0\in U$，对任一点 $x\in U$ 选一条路径 $\gamma:x_0\to x$，令
$$
g(x):=U[\gamma].
$$
因为 $F=0$，从 Ambrose–Singer 定理可知不同路径给出的 $U[\gamma]$ 只依赖端点（同伦相关路径 holonomy 相同），从而 $g(x)$ 定义良好，且直接满足 $A=g^{-1}\mathrm d g$。

在规范变换 $|\psi\rangle\mapsto g|\psi\rangle$ 下，联络变为
$$
A' = gAg^{-1}+g\,\mathrm d g^{-1}=0.
$$
因此在这个规范下，沿任意起终点相同且**可缩回一点**的局部回路 $\gamma\subset U$，有
$$
U'[\gamma]=\mathcal P\exp\Big(-\int_\gamma A'\Big)=\mathbf1.
$$
也就是说：在 $F=0$ 的区域里，所有局部小回路的 holonomy 都是单位元，所有非平庸的 holonomy 只能来自配置空间的非平凡基本群 $\pi_1(\mathcal C)$ 的元素。

在我们的设置中，$n$ 个缺陷的配置空间 $\mathcal C$ 的基本群是辫子群 $B_n$，于是平坦联络 $A$ 的 holonomy 给出了一个群表示
$$
\rho: B_n\cong\pi_1(\mathcal C)\longrightarrow U\big(\mathcal H_{\text{zero modes}}\big),
\qquad
[\gamma]\longmapsto U_R[\gamma].
$$
特别地，对辫子生成元 $\sigma_i$，记
$$
B_i:=\rho(\sigma_i),
$$
则因为 $\rho$ 是群同态，必有
$$
B_iB_{i+1}B_i=B_{i+1}B_iB_{i+1},
\qquad
[B_i,B_j]=0\ \ (|i-j|\ge2),
$$
这就是辫子关系在零模 Hilbert 空间中的算符形式。



####  TQFT 层证明：世界线 braid 表示与 mapping class 表示一致

现在假设给定一个 2+1 维拓扑量子场论，在一个带 n 个穿孔的空间截面上，其 Hilbert 空间记为
$$
\mathcal H_0(D_n; \{\sigma_i\}),
$$
其中 $\sigma_i$ 是附着在各穿孔上的 anyon 类型（例如 Ising TQFT 中的 $\sigma$）。

在这样的 TQFT 中，存在两种等价的幺正表示构造：

1. **世界线表示**：

    - 把 braid $\beta\in B_n$ 画成 2+1 维时空中的 worldlines；

    - 按照 F,R‑符号的 Reshetikhin–Turaev 规则计算振幅，得到
        $$
        \rho_{\mathrm{WL}}(\beta)\in U(\mathcal H_0).
        $$

2. **mapping class 表示**：

    - 把 mapping class $[f]\in\mathrm{MCG}(D_n,\partial D)$ 看作“空间手术”操作：切‑扭‑粘；

    - 由 3D TQFT 的公理，该操作诱导
        $$
        \rho_{\mathrm{MCG}}([f])\in U(\mathcal H_0).
        $$

标准 TQFT 结构（见 Reshetikhin–Turaev、Walker 等）告诉我们：

**定理 a.1（表示一致性）.** 在上述构造下，有
$$
\rho_{\mathrm{WL}}(\beta) = \rho_{\mathrm{MCG}}(\Phi(\beta))
$$
至多差一个整体 U(1) 相位。这意味着：

- “时间中沿 worldlines 编织 anyon” 与
- “在空间截面上做对应的 mapping class（Dehn twist word）”

在 Hilbert 空间上给出的是**同一个幺正表示**。

结合推论 1.4 即得：

> **推论 a.2（TQFT 版主命题）.** 在任意 2+1 维 TQFT 中，对任意 braid word $\beta\in B_n$，存在简单闭曲线 $c_k$ 及符号 $\varepsilon_k$，使得在 $\mathcal H_0$ 上有
> $$
> \rho_{\mathrm{WL}}(\beta) = \rho_{\mathrm{MCG}}\Bigl(\prod_k T_{c_k}^{\varepsilon_k}\Bigr)
> $$
> （至多差一整体相位）。
>
> 也就是说，在 TQFT 的层面 **“时间中的 braid” 与若干“空间 Dehn twist” 的组合本来就是同一类幺正操作。**



#### 引理C证明

用非阿贝尔 Stokes 定理把两个闭路对应的 holonomy 差异用曲率的面积积分控制，并把曲率的范数上界以哈密顿量的谱隙与微扰分量的范数显式量化。

假设和定义
- 参数空间点记为 $\lambda\in\Lambda$，沿参数路径变动的哈密顿为
  $$H(\lambda)=H_0(\lambda)+V(\lambda).$$
- 基态编码投影为 $P(\lambda)$，理想情形的投影记为 $P_0(\lambda)$（对应 $H_0$）。记 $\Delta_0:=\inf_{\lambda}\mathrm{gap}(H_0(\lambda))>0$。假设
  $$\sup_{\lambda}\|V(\lambda)\|=:\epsilon_{\mathrm{tot}}<\Delta_0/2.$$
- 把微扰细分为三类局域分量：
  $$V(\lambda)=H_{\mathrm{int}}(\lambda)+H_{\mathrm{string}}(\lambda)+H_{\mathrm{gauge}}(\lambda),$$
  并定义各自范数上界 $\epsilon_{\alpha}:=\sup_{\lambda}\|H_{\alpha}(\lambda)\|$，使 $\epsilon_{\mathrm{tot}}=\sum_\alpha\epsilon_{\alpha}$。
- 假设 $H(\lambda)$ 在参数上 $C^1$（足够光滑），且对参数导数有界：
  $$
  M_1:=\sup_{\lambda,a}\|\partial_a H(\lambda)\|<\infty.
  $$

命题:
  在上述假设下，对于任何两条同边界的闭路所张成曲面 $S$，相应 holonomy 差异有上界
$$
  \|U_{\mathrm{pert}}(\partial S)-U_0(\partial S)\| \le C\,A_S\,\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}\,\frac{M_1^2}{\Delta_0^4} + R(\epsilon_{\mathrm{tot}},\Delta_0),
$$
其中 $A_S$ 为曲面面积，$C$ 与格点几何有关且与各分量的局域/支撑半径通过常数因子耦合，余项 $R(\cdot)$ 量级 $O(\epsilon_{\mathrm{tot}}^2/\Delta_0^6)$。

证明:

(1) 投影的解析表示与基本界定。
对每个 $\lambda$，在复平面选取闭曲线 $\Gamma$ 包围基态分离谱且与最近谱点距离至少 $\Delta:=\Delta_0-\|V\|\ge\Delta_0/2$。投影有留数表示
$$
  P(\lambda)=-\frac{1}{2\pi i}\oint_{\Gamma}R(z;\lambda)\,dz,\\ R(z;\lambda)=(H(\lambda)-z)^{-1}.
$$
  因此对参数微分
$$
  \partial_a P = -\frac{1}{2\pi i}\oint_{\Gamma} R(z)\,(\partial_a H)\,R(z)\,dz.
$$
  对任意 $z\in\Gamma$ 有 $\|R(z)\|\le 1/\mathrm{dist}(z,\mathrm{spec}(H))\le 2/\Delta_0$，故
$$
  \|\partial_a P\| \le C_1\frac{\|\partial_a H\|}{\Delta_0^2}\le C_1\frac{M_1}{\Delta_0^2},
$$
 其中 $C_1=|\Gamma|/(2\pi)\cdot 4$ 为几何常数（可选取相同的 $\Gamma$ 以统一常数）。

(2) 曲率的表示与范数界。
曲率定义为 $F_{ab}=P[\partial_a P,\partial_b P]P$。因此
$$
  \|F_{ab}\|\le 2\,\|\partial_a P\|\,\|\partial_b P\| \le C_2\frac{\|\partial_a H\|\,\|\partial_b H\|}{\Delta_0^4}\le C_2\frac{M_1^2}{\Delta_0^4}.
$$
所以曲率在曲面上的面积积分满足
$$
\iint_S\|F\| \le C_2\,A_S\,\frac{M_1^2}{\Delta_0^4}.
$$
(3) 投影的微扰展开与 $\delta P$ 的界。

使用解析函数微扰论或留数差分，得到
$$
\delta P:=P-P_0 = -\frac{1}{2\pi i}\oint_{\Gamma}\big(R(z;H)-R(z;H_0)\big)dz
$$

$$
R(z)-R_0(z)=-R_0(z)V R(z)=-\sum_{n\ge1}(-1)^{n-1}R_0(z)\big(VR_0(z)\big)^n,
$$

当对围道 $\Gamma$ 上的点满足 $\sup_{z\in\Gamma}\|R_0(z)V\|<1$ 时，上述 Neumann 级数收敛（在我们的假设 $\|V\|<\Delta_0/2$ 下通常成立，因为 $\|R_0(z)\|\le 2/\Delta_0$）。代入投影的留数表达式并取范数得
$$
\|\delta P\|=\Big\| -\frac{1}{2\pi i}\oint_{\Gamma}\big(R-R_0\big)dz\Big\| \le \frac{|\Gamma|}{2\pi}\sup_{z\in\Gamma}\sum_{n\ge1}\|R_0(z)\|^{n+1}\,\|V\|^n.
$$
当记 $q:=\sup_{z\in\Gamma}\|R_0(z)\|\,\|V\|<1$ 时，上述级数求和并简化得
$$
\|\delta P\| \le \frac{|\Gamma|}{2\pi}\frac{\|R_0\|^2\,\|V\|}{1-\|R_0\|\,\|V\|} \le \frac{|\Gamma|}{2\pi}\frac{(2/\Delta_0)^2\,\|V\|}{1-2\|V\|/\Delta_0}=:C_3\frac{\|V\|}{\Delta_0}.
$$


所以在常见的小微扰情形 $\|V\|\ll\Delta_0$ 下，可进一步用简化常数写作 $\|\delta P\|\lesssim C_3'\,\|V\|/\Delta_0$。

对导数差量的估计，从导数的留数表达式出发：
$$
\partial_a P - \partial_a P_0 = -\frac{1}{2\pi i}\oint_{\Gamma} \big(R(\partial_a H)R - R_0(\partial_a H_0)R_0\big)\,dz.
$$
把被积子重写为三项之和：
$$
\begin{aligned}

R(\partial_a H)R - R_0(\partial_a H_0)R_0 &= (R-R_0)(\partial_a H)R 

 + R_0(\partial_a H - \partial_a H_0)R 

 + R_0(\partial_a H_0)(R-R_0).

\end{aligned}
$$
对每一项使用范数不等式和下面的 resolvent 差界可得：

- 第一与第三项（含 $R-R_0$）：利用 $\|R-R_0\|\le \|R_0\|^2\,\|V\|/(1-\|R_0\|\,\|V\|)$ 以及 $\|R\|\le \|R_0\|/(1-\|R_0\|\,\|V\|)$，得到典型贡献

$$
\Big\|\frac{1}{2\pi i}\oint_{\Gamma}(R-R_0)(\partial_a H)R\,dz\Big\| \lesssim \frac{|\Gamma|}{2\pi}\frac{\|R_0\|^2\,\|V\|}{1-\|R_0\|\,\|V\|}\cdot\|\partial_a H\|\cdot\|R\| 

  \lesssim C'\frac{\|\partial_a H\|\,\|V\|}{\Delta_0^3(1-2\|V\|/\Delta_0)^2}.
$$

- 第二项（含 $\partial_a H-\partial_a H_0=\partial_a V$）：直接用 $\|R_0\|\le 2/\Delta_0$ 和 $\|R\|\le 2/\Delta_0$，得到

$$
\Big\|\frac{1}{2\pi i}\oint_{\Gamma}R_0(\partial_a V)R\,dz\Big\| \le \frac{|\Gamma|}{2\pi}\|R_0\|\,\|\partial_a V\|\,\|R\| \le C''\frac{\|\partial_a V\|}{\Delta_0^2}.
$$

综上合并主要项，并在 $\|V\|/\Delta_0$ 小的前提下忽略更高阶项，得到常用的简化界
$$
\|\partial_a P - \partial_a P_0\| \le C_4\Big(\frac{\|\partial_a H\|\,\|V\|}{\Delta_0^3} + \frac{\|\partial_a V\|}{\Delta_0^2}\Big) + O\Big(\frac{\|V\|^2\,\|\partial_a H\|}{\Delta_0^4}\Big).
$$
其中 $\|\partial_a V\|$ 可進一步以各分量 $H_{\alpha}$ 的界表示（例如 $\|\partial_a V\|\le\sum_{\alpha}\|\partial_a H_{\alpha}\|$），從而把导数差的界写成各微扰成分的贡献之和。



我们还是用系数统一表示为以下形式：
$$
\|\delta P\| \le C_3\frac{\|V\|}{\Delta_0}=C_3\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}.
$$
$$
  \|\partial_a P - \partial_a P_0\| \le C_4\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}\frac{\|\partial_a H\|}{\Delta_0^2} + O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^3}\Big).
$$

(4) 曲率差的展开与分量依赖。
 将 $P=P_0+\delta P$ 代入 $F=P[\partial_a P,\partial_b P]P$，展开差量得到由 $\delta P$ 与 $\partial P-\partial P_0$ 控制的若干项。

$(P=P_0+\delta P)，(\partial P=\partial P_0+\delta(\partial P))$代入并得到若干项（用三角不等式分项上界）
$$
  \|\widetilde F - F_0\| \le C_5\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}\frac{M_1^2}{\Delta_0^4} + O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^6}\Big).
$$
  若需把 $\epsilon_{\mathrm{tot}}$ 拆成各分量:
$$
\epsilon_{\mathrm{tot}}=\sum_{\alpha}\epsilon_{\alpha},\qquad \|\partial_a H\|\lesssim \|\partial_a H_0\| + \sum_{\alpha}C_{\alpha}\,\epsilon_{\alpha},
$$
  其中 $C_{\alpha}$ 包含支撑半径与几何权重（对长串项 $C_{\mathrm{string}}$ 较大）。因此
  $$\|\widetilde F - F_0\| \lesssim \sum_{\alpha}C'_{\alpha}\frac{\epsilon_{\alpha}}{\Delta_0}\frac{M_1^2}{\Delta_0^4}+O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^6}\Big).$$

(5) 曲率差到 holonomy 差的转换（Stokes）
非阿贝尔 Stokes 给出($\mathcal S$表示表面的序)
$$
U(\partial S)=\mathcal S\exp\Big(\iint_S W(s,t)F(s,t)W(s,t)^{-1}dS\Big).
$$
把曲面算子写为
$$
X:=\iint_S W(s,t)F(s,t)W(s,t)^{-1}\,dS,

\qquad X_0:=\iint_S W(s,t)F_0(s,t)W(s,t)^{-1}\,dS,
$$


其中共轭因子 $W(s,t)$ 为酉，不改变范数。由定义有
$$
\|X-X_0\| \le \iint_S\|F-F_0\|,\qquad \|X\|\le \iint_S\|F\|,\;\; \|X_0\|\le \iint_S\|F_0\|.
$$
利用 Duhamel 表示或算符指數的基本积分恒等，可以得到严格的差分界：
$$
e^{X}-e^{X_0}=\int_0^1 e^{(1-s)X_0}(X-X_0)e^{sX}\,ds
$$
从而
$$
\|e^{X}-e^{X_0}\| \le \|X-X_0\|\,e^{\|X\|+\|X_0\|}.
$$
因此
$$
\|U_{\mathrm{pert}}(\partial S)-U_0(\partial S)\| \le \Big(\iint_S\|\widetilde F-F_0\|\Big)\,\exp\Big(\iint_S\|F\|+\iint_S\|F_0\|\Big).
$$


代入  (2)–(4) 中得到的曲率与曲率差的上界（把 $\iint_S\|\widetilde F-F_0\|$ 用 Step 4 的界代替，並把 $\iint_S\|F\|,\iint_S\|F_0\|$ 用 (2)的界代替），得到带指数因子的更精确界：
$$
\|U_{\mathrm{pert}}-U_0\| \le \Big(C_6\,A_S\,\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}\,\frac{M_1^2}{\Delta_0^4}\Big)\exp\Big(C_7\,A_S\,\frac{M_1^2}{\Delta_0^4}\Big) + O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^6}\Big).
$$
在常用的小曲率情形（即 $A_S\,M_1^2/\Delta_0^4\ll1$ 且 $\epsilon_{\mathrm{tot}}/\Delta_0\ll1$）下，指數項可展開為 $1+o(1)$，主阶项退化为线性形式，得到下面所写的近似表达（常数 $C$ 可由 $C_6,C_7,C_5,C_2$ 组合得到）。

 

对于小范数场，指数映射线性化可得
$$
\|U_{\mathrm{pert}}(\partial S)-U_0(\partial S)\| \le \iint_S\|\widetilde F-F_0\|\,+\,O\Big(\big(\iint_S\|F\|\big)^2\Big).
$$
 代入 (4)的界，得到命题不等式（常数 $C$ 可用 $C_5,C_2$ 与几何因子组合表示）：
$$
  \|U_{\mathrm{pert}}(\partial S)-U_0(\partial S)\| \le C\,A_S\,\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}\,\frac{M_1^2}{\Delta_0^4} + O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^6}\Big).
$$



常数与改善方向（注解）

- 上述降冪 $\Delta_0^{-5}$（整体看似）来自对 $\partial P$ 的双重 resolvent 估计；更精细的交换子重排或使用分块对角化（Schrieffer–Wolff）可在某些情形下把分母冪次降低到 $\Delta_0^{-2}$ 或 $\Delta_0^{-3}$，但需额外假设（例如 $\partial_a H$ 与 $[H_0,V]$ 的结构性小）。
- 若要数值对接，本证明中出现的常数 $C,C_1,\dots$ 与几何权重 $C_{\alpha}$ 应由具体模型（格点拓扑、算符支撑半径）评估。建议把各分量的范数 $\epsilon_{\alpha}$ 由 `kits/compute_jw_mapping.py` 输出，并用本文不等式估算期待的 $\varepsilon_{\mathrm{YBE}}$。

结论：在假设 $\epsilon_{\mathrm{tot}}<\Delta_0/2$ 且 $H(\lambda)$ 平滑的前提下，Stokes‑based 的严格推导把 holonomy / YBE 偏差以曲率差的面积积分表示，并可用谱隙与各微扰分量范数给出如上明确上界，从而把主文中抽象的 $\epsilon_{\mathrm{tot}}$ 具体化为对应分量的数值可验证界。



#### 引理 D 

从准绝热生成元的构造出发，严格推导出生成元分解为局域项、单项范数上界 $k_*$ 与采用 $m$ 步 Trotter 展开时的误差上界，同时估计因空间截断引入的指数小尾项。

用一个平滑的滤波核把投影导数卷进去，得到 quasi‑adiabatic 生成元 K(t)。直观上它是“平滑版”的瞬时几何生成器，用来把基态子空间沿路径搬运，而不会把能量推到激发态

假设与符号：
- 投影族 $P(t)$（$t\in[0,T]$）光滑，对应哈密顿 $H(t)$ 在整个区间有谱隙 $\Delta:=\inf_t\mathrm{gap}(H(t))>0$。

 - 定义 $\dot P:=\partial_t P$，并采用光滑滤波核 $W(s)$（快速衰减且奇偶性适当）构造 quasi‑adiabatic 生成元
	$$
	K(t)=i\int_{-\infty}^{\infty} W(s) e^{iH(t)s}\,\dot P(t)\,e^{-iH(t)s}\,ds.
	$$
	

 常用选择使得 $\hat W(\omega)=\int e^{i\omega s}W(s)ds$ 在 $|\omega|<\Delta/2$ 附近等于 $-1/\omega$，以产生把基态子空间沿时间演化平移的生成元。



**Lieb–Robinson 引理**：在有局域相互作用的量子格点体系里，“信息”和局域扰动不会瞬间传到远处，而是以有限速度传播。

形式上，对两个分别作用在距离 d 的局域算符 A, B，有常数$ C, μ, v$（Lieb–Robinson 速度）使得随时间 t 的对易子满足

$∥[A(t),B]∥≤C ∥A∥∥B∥ e^{−μ(d−v∣t∣)}.$

当 v∣t∣<d*v*∣*t*∣<*d*,在光锥之外时，对易子被指数压制——远处几乎不受影响。结果：算符在有限时间内保持“准局域性”，这正是把连续生成元截断为局域项并只产生指数小尾项的理论基础。



1：生成元的几何与局域性（准局域性）

1.1. 由于 $\dot P$ 为局域算符（由 $H$ 的局域性与有限导数保证），且 $e^{iH s} O e^{-iH s}$ 随 $s$ 的增长以 Lieb–Robinson 光速扩散，卷积核 $W(s)$ 的快速衰减保证了积分产生的一个准局域算符。更精确地，对任何局域算符支撑半径为 $r$ 的 $O$，存在常数 $c,\mu,v_{LR}$，使得
$$
    \|[e^{iHs} O e^{-iHs}, X]\|\le c\,\|O\|\,\|X\|\,e^{-\mu(d-r-v_{LR}|s|))}
$$
（$d$ 为算符支撑间距）。把此不等式与 $W(s)$ 卷积可以导出：若截断 $K(t)$ 的支撑到半径 $R_k$（即只保留与基点距离小于 $R_k$ 的项），截断误差为
$$
    \|K(t)-K_{R_k}(t)\| \le C' e^{-\mu' R_k},
$$


其中 $\mu'$ 与 $\mu$、$W$ 的衰减常数以及 $v_{LR}$ 有关。此即所谓空间截断的指数小尾项（Lieb–Robinson tail）。

2：生成元每项的范数上界 $k_*$

2.1. 从 $K(t)$ 的表达式及留数傅里叶分析，可把 $K(t)$ 展为局域项的和：
$$
K(t)=\sum_{\alpha} k_{\alpha}(t),
$$
每个 $k_{\alpha}$ 支撑在直径 $\lesssim R_k$ 的小区域内。对单项范数上界的常规估计给出
$$
\|k_{\alpha}(t)\| \le C_{\alpha}\int |W(s)|\,\|e^{iHs}\dot P e^{-iHs}\|\,ds \lesssim C_{\alpha}'\frac{\|\dot H\|}{\Delta},
$$
​    其中使用了 $\dot P$ 与 $\dot H$ 的线性关系（通过投影导数的解析表示，见下），以及 $W$ 的频域特性导致一个 $1/\Delta$ 因子。于是定义
$$
k_*:=\sup_{t,\alpha}\|k_{\alpha}(t)\| \lesssim C_{\mathrm{qa}}\frac{\|\dot H\|}{\Delta}.
$$
​    常数 $C_{\mathrm{qa}}$ 含有 $W$ 的 L1 范数与局域性几何因子。

3：Trotter 时间离散化误差

3.1. 把时间分为 $m$ 段，每段长度 $\delta t=T/m$，令 $K_j:=K(t_j)$（取中点或左端点均可），考虑近似
$$
\mathcal T e^{-i\int_0^T K(t)dt} \approx \prod_{j=1}^m e^{-iK_j\delta t}.
$$
​    使用基础的时间切分误差分析或 Duhamel 展开与交换子估计，可得到二阶型（或一阶）误差上界：
$$
\Big\|\mathcal T e^{-i\int_0^T K(t)dt} - \prod_{j=1}^m e^{-iK_j\delta t}\Big\| \le C_T\frac{T}{m}\sup_{j,j'}\|[K_j,K_{j'}]\| + O\Big(\frac{T^2}{m^2}\max\|K\|^3\Big).
$$
​    交换子的范数可用局域分解与 $k_*$ 上界估计：若 $K_j=\sum_\alpha k_{\alpha}^j$，则
$$
\|[K_j,K_{j'}]\| \le \sum_{\alpha,\beta}\|[k_{\alpha}^j,k_{\beta}^{j'}]\| \lesssim C''\,k_*^2
$$
​    （非重叠支撑对 commutator 为零，重叠项数由局域结构给出常数因子）。因此
$$
\Big\|\mathcal T e^{-i\int_0^T K(t)dt} - \prod_{j=1}^m e^{-iK_j\delta t}\Big\| \le C_D\frac{T k_*^2}{m} + O\Big(\frac{T^2 k_*^3}{m^2}\Big).
$$
4：把局域指数分解为基本门（R 门）与累积误差

4.1. 每个局域指数 $e^{-i k_{\alpha}\delta t}$ 由有限数目的 $R$ 类门（本研究中 $R=e^{iH_P^{(ij)}}$）近似实现，近似误差可由 Trotter‑Suzuki 或 Lie–Trotter 展开控制；若对每个局域项采用常数轮次的门分解，则总门数与局域项数成正比，累积误差为每项门分解误差乘以项数，故不改变时间离散误差的主要参数依赖（$m,k_*,T$）。

5：空间截断误差（Lieb–Robinson tail）

5.1. 如步骤 1 所述，将每个 $k_{\alpha}$ 截断到半径 $R_k$ 的局域算符会产生尾项，其范数界为
$$
    \|k_{\alpha}-k_{\alpha}^{(R_k)}\| \le C_3 e^{-\mu R_k}.
$$


​    把所有项加总得到总截断尾项
$$
    \|K-K^{(R_k)}\| \le C_4 e^{-\mu R_k}
$$
​    并可把该尾项对时间演化的影响用 Dyson / Duhamel 展开证明为同阶指数小量，故在最终误差中以加项的形式出现。

总结（引理 D 结论）

综上，存在常数 $C_D,C'_D,\mu>0$ 使得
$$
\varepsilon_{\mathrm{Trotter}}:=\Big\|\mathcal T e^{-i\int_0^T K(t)dt} - \prod_{j=1}^m \prod_{\alpha} e^{-i k_{\alpha}^j\delta t}\Big\| \le C_D\frac{T k_*^2}{m} + C'_D e^{-\mu R_k} + O\Big(\frac{T^2 k_*^3}{m^2}\Big),
$$
​    且
$$
k_*\lesssim C_{\mathrm{qa}}\frac{\|\dot H\|}{\Delta},
$$
