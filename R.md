## Solution of Yang Baxter Equation and Integrable Sub manifold

### Content

1. **YBE的可行解及辫子群**
2. **配置空间、Dehn twist 和 half twist**
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
但是QDYBE的约束都非常复杂，所以这里先尝试了R是在SU(2)作用下不变，即$(U \otimes U)R(U^{\dagger}\otimes U^{\dagger}) = R, U\in SU(2)$, 或者表示成李代数形式 $[R,\Delta J] = 0, \forall J \in su(2)$ , 

由此得出了$R = a P_3 + b P_1$的形式，其中$P_3 = \frac{1+P}{2}, P_1 = \frac{1-P}{2} , P:=P(\alpha \otimes \beta) = \beta \otimes \alpha$ ， 即对应这SU(2)群的分解

$2 \times 2 = 3+1$的形式，$P_i$则表示了投影到对称、反对称部分， 因为投影算符的性质，这个R是满足YBE的。



借助这个简单的例子，下面将R写成了
$$
R=a\,I + b\,\sigma^x\sigma^x + c\,\sigma^y\sigma^y + d\,\sigma^z\sigma^z.
$$
的形式，将张量积$\sigma \otimes \sigma$简写做了$\sigma \sigma$， 然后考虑满足QDYBE$(a,b,c,d)$的约束(在最后面部分给出）。 相比之前的四元数模型，是将控制门形式改写成了现在的相互作用形式。



#### 推广到算符

格点体系的全 Hilbert 空间为单点态空间 $V$ 的张量积 $V^{\otimes L}$。任何作用在两点上的算符$R:\;V\otimes V\to V\otimes V$
可唯一嵌入到全链为局域算符：
$$
R_{i,i+1}=I^{\otimes(i-1)}\otimes R\otimes I^{\otimes(L-i-1)}.
$$
作用在第 i, i+1 位用张量积表示,是局域两体算符。

如果$R$ 满足 Yang–Baxter 方程，则这些局域嵌入满足局域 YBE，从而给出 braid‑group 的表示和可积传输矩阵的局部生成元——这正是把代数生成元 $b_i$ 映为 $\mathrm{id}^{\otimes(i-1)}\otimes R\otimes\mathrm{id}^{\otimes(n-i-1)}$ 的代数依据$^{[2]}$。满足辫群关系的参数(a,b,c,d)见附录.

R满足YBE后，即满足了短程算符，下面简单证明长程算符的形式为：
$$
R_{ij}= aI + b\,\sigma_i^x\sigma_j^x + c\,\sigma_i^y\sigma_j^y + d\,\sigma_i^z\sigma_j^z
$$
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

这里有
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
H_P = \sum_{\mu,\nu\in\{0,x,y,z\}} c_{\mu\nu}\,\sigma^\mu_j\otimes\sigma^\nu_{j+1}
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

### 2.配置空间、Dehn twist 、half twist







### reference

[1]*The Dynamical Yang-Baxter Equation, Representation Theory, and Quantum Integrable Systems*
[2]*nLab https://ncatlab.org/nlab/show/Yang-Baxter+equation* 



### Appendix

(这里由程序解出12个约束后的结果)

 Yang–Baxter 给出的代数约束（对复参数可分实、虚部分）：
$$
a d (b-c)=0,\\
b c (a-d)=0,\\
a b c - a b d - a c d + b c d = 0.
$$

