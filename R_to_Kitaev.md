### 从 Yang Baxter 到 Kitaev 链

#### YBE

*这一部分参考了文献*

*The Dynamical Yang-Baxter Equation, Representation Theory, and Quantum Integrable Systems*
*以及nLab https://ncatlab.org/nlab/show/Yang-Baxter+equation* 



Quantum dynamical ybe (QDYBE)的定义为

$R: V\otimes V \rightarrow V \otimes V$

$(R⊗id)∘(id⊗R)∘(R⊗id)=(id⊗R)∘(R⊗id)∘(id⊗R)$



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
的形式，将张量积$\sigma \otimes \sigma$简写做了$\sigma \sigma$， 然后考虑满足QDYBE$(a,b,c,d)$的约束(在最后面部分给出）。 相比之前的四元数模型，是将控制门形式改写成了现在的相互作用形式。然后根据

![image-20260325191530677](/home/asice-cloud/.config/Typora/typora-user-images/image-20260325191530677.png)

可以将R进行推广。



#### 从 $R$ 推广为 $R_{i,i+1}$ 

格点体系的全 Hilbert 空间为单点态空间 $V$ 的张量积 $V^{\otimes L}$。任何作用在两点上的算符$R:\;V\otimes V\to V\otimes V$
可唯一嵌入到全链为局域算符：
$$
R_{i,i+1}=I^{\otimes(i-1)}\otimes R\otimes I^{\otimes(L-i-1)}.
$$
 作用在第 i, i+1 位用张量积表示,是局域两体算符。

如果$R$ 满足 Yang–Baxter 方程，则这些局域嵌入满足局域 YBE，从而给出 braid‑group 的表示和可积传输矩阵的局部生成元——这正是把代数生成元 $b_i$ 映为 $\mathrm{id}^{\otimes(i-1)}\otimes R\otimes\mathrm{id}^{\otimes(n-i-1)}$ 的代数依据。

关于考虑使用 Jordan–Wigner ：自旋算符映射到费米子时通常出现串算子 $e^{i\pi\sum_{k<j}n_k}$,在链上求和得总哈密顿 H = Σ_j R_{j,j+1}。Jordan–Wigner 的串因子在最近邻乘积中抵消，但对于最近邻两点 $i,i+1$，这些串因子在算符乘积中互相抵消，因此 $R_{i,i+1}$ 映射后仍是局域的二体/四体费米子算符（没有产生非局域长串）。这保证了把自旋层面的 $R_{i,i+1}$ 直接翻译为 Kitaev局域项的可行性。



#### 从 $R_{i,i+1}$ 到 $t,\,\Delta,\,\mu$ 

这一部分是在 Jordan–Wigner映射与链上求和约定下，证明对于

![{\displaystyle H=-\mu \sum _{j=1}^{N}\left(c_{j}^{\dagger }c_{j}-{\frac {1}{2}}\right)+\sum _{j=1}^{N-1}\left[-t\left(c_{j+1}^{\dagger }c_{j}+c_{j}^{\dagger }c_{j+1}\right)+|\Delta |\left(c_{j+1}^{\dagger }c_{j}^{\dagger }+c_{j}c_{j+1}\right)\right]}](https://wikimedia.org/api/rest_v1/media/math/render/svg/cd0eab38d1481b89a6fc627707733bede18873c8)

在现在的模型
$$
R_{i,i+1}=a\,I + b\,\sigma^x_i\sigma^x_{i+1} + c\,\sigma^y_i\sigma^y_{i+1} + d\,\sigma^z_i\sigma^z_{i+1}.
$$
下，得到：
$$
t\propto(b+c),\qquad \Delta\propto(b-c),\qquad \mu\text{ 的线性部分由 }a,d\text{ 给出（可被重整化）}.
$$
并给出常数因子在常见归一化下的具体值。



**一、算符**
$$
R_{i,i+1}=a\,I + b\,\sigma^x_i\sigma^x_{i+1} + c\,\sigma^y_i\sigma^y_{i+1} + d\,\sigma^z_i\sigma^z_{i+1}.
$$

***符号说明***

 $\sigma^\alpha$（$\alpha=x,y,z$）是每个格点上的 Pauli 矩阵，作用在单点的二维态空间 $V$ 上。（在 $|\uparrow\rangle=(1,0)^T,\;|\downarrow\rangle=(0,1)^T$ 基下）：
$$
\sigma^x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
\sigma^y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad
\sigma^z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$
定义升降算符
$$
\sigma^\pm=\tfrac12(\sigma^x\pm i\sigma^y),\qquad
\sigma^+=\begin{pmatrix}0&1\\0&0\end{pmatrix},\;\sigma^- =\begin{pmatrix}0&0\\1&0\end{pmatrix}.
$$


局域标示 $\sigma^\alpha_i$ 表示这个矩阵嵌入到全链 $V^{\otimes L}$ 的第 $i$ 位（即 $I^{\otimes(i-1)}\otimes\sigma^\alpha\otimes I^{\otimes(L-i)}$,因此项 $\sigma_i^a \sigma_{i+1}^a$ 实际是两个格点上的张量积算符 $\sigma^a \otimes \sigma^a$（嵌入到全链上）。）。它们满足代数关系
$$
\{\sigma^\alpha,\sigma^\beta\}=2\delta_{\alpha\beta}I,\qquad [\sigma^\alpha,\sigma^\beta]=2i\epsilon_{\alpha\beta\gamma}\sigma^\gamma.
$$




在 Jordan–Wigner 映射中用到的关系为：
$$
\sigma^+_j\mapsto c_j^{\dagger}e^{i\pi\sum_{k<j}n_k},\qquad \sigma^-_j\mapsto c_j e^{i\pi\sum_{k<j}n_k},\qquad \sigma^z_j\mapsto 2n_j-1.
$$

（最近邻 $i,i+1$，串算符在乘积中抵消，使局域自旋两体算符在费米子表示下仍为局域项。）



**二、用升降算符展开（代数恒等）**
记 $\sigma^\pm=\tfrac{1}{2}(\sigma^x\pm i\sigma^y)$，于是等价写法：
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



**三、Jordan–Wigner 映射**
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

对于$e^{iπ n_i}=(-1)^{n_i}$, 占据数基$ |… n_i,n_{i+1} …⟩$ 上有$ n_i,n_{i+1}∈{0,1}$,  $c_i^† c_{i+1} $只有在$ n_i=0$ 且$ n_{i+1}=1 时产生非零结果$而在这情形下$ e^{iπ n_i}=1$；其它情形下要么两边同时为 0，要么因占据限制而抵消。配对项 $c_i^† c_{i+1}^† $的情形同理（只有在$ n_i=n_{i+1}=0$ 时非零，此时$ e^{iπ n_i}=1$）。



于是对最近邻 $i,i+1$，串算符 $e^{i\pi\sum_{k<i}n_k}$ 在乘积中相互抵消，因此有
$$
\begin{aligned}
\sigma^+_i\sigma^-_{i+1}+\sigma^-_i\sigma^+_{i+1}&\mapsto c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i,\\
\sigma^+_i\sigma^+_{i+1}+\sigma^-_i\sigma^-_{i+1}&\mapsto c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i,\\
\sigma^z_i\sigma^z_{i+1}&\mapsto (2n_i-1)(2n_{i+1}-1)=4n_in_{i+1}-2(n_i+n_{i+1})+1.
\end{aligned}
$$



**四、逐项代入并识别二次项, 代入上述等式到 $R_{i,i+1}$：**
$$
\begin{aligned}
R_{i,i+1} &= a\,I + b\,\sigma^x_i\sigma^x_{i+1} + c\,\sigma^y_i\sigma^y_{i+1} + d\,\sigma^z_i\sigma^z_{i+1} \\
&= a\,I + (b+c)\big(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i\big) + (b-c)\big(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i\big)\\
&\quad + d\big(4n_in_{i+1}-2(n_i+n_{i+1})+1\big).
\end{aligned}
$$


在单个键（bond）处，二次费米子项的系数即如上所示：
- hopping (单键处)：系数 $(b+c)$ * $(c_i^{\dagger}c_{i+1}+h.c.)$；
- pairing (单键处)：系数 $(b-c)$ * $(c_i^{\dagger}c_{i+1}^{\dagger}+h.c.)$。

因此局域识别给出 $t\propto(b+c)$、$\Delta\propto(b-c)$（比例常数由你对总体哈密顿量的归一化决定）。



**五、链上求和与化学势项的精确系数**
取常见的做法把哈密顿量定义为所有键的和：
$$
H=\sum_{j=1}^{L-1} R_{j,j+1}.
$$
(开链，周期边界类似但要处理双计数)

把第四项中关于 $d$ 的线性密度部分写出：单键产生 $-2d(n_j+n_{j+1})$。把所有键求和：对于体内站点（非边界），$n_j$ 出现在键 $(j-1,j)$ 和 $(j,j+1)$ 的线性项中，各贡献 $-2d$，合计为 $-4d$。因此链上线性项为
$$
H_{\rm lin}^{(d)} = -4d\sum_{j=1}^L n_j + \text{(边界处的修正)}.
$$
同时，常数项与 $a$ 的贡献为每个键贡献 $a$；求和后产生整体能级位移 $a\times(L-1)$，该常数可以通过重定义参考能量或并入 $\mu$ 的常数偏移中吸收。

若采用kitaev chain惯例写法
$$H= -\mu\sum_j\Big(n_j-\tfrac12\Big) - t\sum_j\big(c_j^{\dagger}c_{j+1}+h.c.\big) + \Delta\sum_j\big(c_j c_{j+1}+h.c.\big)+\cdots,$$
则将上面链上求和结果与该形式比较，可选择直接取归一化使得
$$t = b+c,\qquad \Delta = b-c,$$
并识别
$$
\mu = 4d + \mu_0,\qquad \mu_0\text{ 来自 }a\text{（常数偏移，可吸收）}.
$$
（若在定义中对 $t$ 或 $\Delta$ 加了额外的负号或 $1/2$ 因子，上述等号右边应乘以相应的归一化常数。）



**六、关于 $d$ 项的“可重整化”说明**
四费米子项 $4d\,n_in_{i+1}$ 是相互作用项，不属于二次自由理论；但是在低能或平均场近似中，此项会带来两个重要效应：

- 平均场分解 $n_in_{i+1}\to n_i\langle n_{i+1}\rangle + \langle n_i\rangle n_{i+1}-\langle n_i\rangle\langle n_{i+1}\rangle$ 直接生成线性项，等效于修改 $\mu$（即重整化 $\mu$）；
- 在严格的重整化群分析中，四费米子在低维（1D）下可能是边缘/微扰量，会把系统流到不同的低能固定点，从而间接改变有效的二次项系数（包括对 $\mu,t,\Delta$ 的重整化）。

补充推导：

1)  $\sigma^z_i\sigma^z_{i+1}$ 的展开：
$$
\sigma^z_i\sigma^z_{i+1}=(2n_i-1)(2n_{i+1}-1)=4n_i n_{i+1}-2(n_i+n_{i+1})+1.
$$

2) 对二次密度项作平均场分解：
$$
n_i n_{i+1}\approx n_i\langle n_{i+1}\rangle + \langle n_i\rangle n_{i+1} - \langle n_i\rangle\langle n_{i+1}\rangle.
$$

3) 把上式代入单键贡献 $d(4n_i n_{i+1}-2(n_i+n_{i+1})+1)$ 得：
$$
d\big[4(n_i\langle n_{i+1}\rangle + \langle n_i\rangle n_{i+1}-\langle n_i\rangle\langle n_{i+1}\rangle)-2(n_i+n_{i+1})+1\big].
$$

4) 若系統平移不变（$\langle n_i\rangle=\bar n$），合并所有键对同一处的贡献，得到线性项系数：
$$
H^{(d),\mathrm{MF}}_{\rm lin}=4d(2\bar n-1)\sum_j n_j + {\rm const.}
$$
5) 相当于重整了化学势
$$
\mu_{\rm eff}=\mu-4d(2\bar n-1),
$$
剩余的涨落项 $(n_i-\bar n)(n_{i+1}-\bar n)$ 表示未被平均场捕获的相互作用，可能在 1D 中导致重要效应，需用更严格的约束

**结论**

- 直接把 $H=\sum_j R_{j,j+1}$ 作为 Kitaev‑型哈密顿量时，可取 $t=b+c,\,\Delta=b-c$，并且链上 $d$ 的线性贡献给出 $\mu\supset 4d$，$a$ 只产生可吸收的常数偏移。
- 若把 $t,\Delta$ 标准化为带负号的常见写法（例如 $-t\sum c_j^{\dagger}c_{j+1}$），只需在等号前乘上所用的总体负号（即 $t_{paper}=- (b+c)$ 等）。

最终得到：
$$
t \propto (b+c),\qquad \Delta \propto (b-c),\qquad \mu\text{ 由 }a,d\text{ 的线性项给出（可重整化）}.
$$

总结就是，张量积形式体现算符是局域的两体算符；b,c 控制二次自由部分（hopping/pairing → 对应邻位 Majorana 双线性），d 同时生成化学势修正与最近邻相互作用（四费米）



**示例(ai 生成)**

- 若要得到自由（非相互作用）Kitaev 链，需要消去四费米子项，即要求
    $d=0.$ 
    这时 $R_{i,i+1}$ 仅是二次费米子算符（加上常数项），可直接对应 $t,\Delta,\mu$。
- 若 $d\neq0$ 则得到含相互作用的扩展模型（例如 XXZ‑类的相互作用），不再是自由马约拉那链。

例 1：取 $a=0,d=0,b=1,c=1$，则
$$t\propto 2,\quad \Delta\propto 0,$$
即纯跳跃自由费米子链。对应 Kitaev 链的配对为 0。

例 2：取 $a=0,d=0,b=1,c=-1$，则
$$t\propto 0,\quad \Delta\propto 2,$$
即纯配对项（最大化 p‑wave 配对）。



### **附录：YBE 代数约束与典型物理情形**

(这里由程序解出12个约束后的结果)

 Yang–Baxter 给出的代数约束（对复参数可分实、虚部分）：
$$
a d (b-c)=0,\\
b c (a-d)=0,\\
a b c - a b d - a c d + b c d = 0.
$$


下表列出这些约束下若干常见物理情形及其意义：

| 物理情形 | 代数约束/参数关系 | 物理含义与备注 |
|---:|---|---|
| XXX（SU(2) 不变） | $b=c=d=t$ → 约束化为 $t^2(a-t)=0$ | 若 $t\neq0$，需 $a=t$，即 $R\propto P$（置换/交换算符）。这是 XXX 点（可 Baxter 化得到 $R(u)=uI+P$）。 |
| XXZ（平面各向同性） | $b=c\neq d$：代数约束通常要求若 $b\neq0$ 则 $a=d$，进一步退化到 $a=b=c=d$（常数情形下难以得到非平凡的常数 R） | 常见的 XXZ 可积族为带谱参数的情形（Baxter 化，双曲/三角函数），但作为常数 4 参数截面通常受限。 |
| XY / 自由 Kitaev‑样族（无相互作用） | $d=0$ 且常见解为 $a=0,d=0$（即 $R=b\,\sigma^x\sigma^x + c\,\sigma^y\sigma^y$ 为一族解） | 通过 Jordan–Wigner 映射得到二次费米子：跳跃 $t\propto(b+c)$，配对 $\Delta\propto(b-c)$。这是 free‑fermion（Kitaev/XY）类模型，适合马约拉那对角化。 |
| 纯跳跃或纯配对（单项） | 仅 b 或仅 c 非零（其余为 0） | 简单独立耦合，往往满足 YBE（方程大多退化为 0），对应纯 hopping 或纯 pairing 的局域项。 |
| 只有 $\sigma^z\sigma^z$（Ising‑型） | $b=c=0$，任意 $a,d$ 满足约束 | 对应 Ising/XXZ 中的相互作用项（四费米子项在 JW 后出现）。可包含化学势偏移与最近邻相互作用。 |
| 平凡/退化情形 | 任意两个或更多参数为 0，如 $a=0,d=0$、或 $b=0,c=0$ 等 | 多数代数方程自动满足，给出容易处理的可积或可对角化子空间。

**小结**：代数约束表明可积（满足常数 YBE）的四参数 $R$ 是受限的——非平凡的 SU(2)‑不变点（XXX）与某些退化/单项耦合属于解集；而一般同时包含跳跃与配对且含相互作用（任意 $b,c,d$）的自由 Kitaev 链并不总是出现在常数 R 的解集中（但可在带谱参数族或通过其他映射下出现）。

- 带谱参数的情形（Baxter 化）：常数 R 解集受限，许多重要可积模型（如 XXZ）并不是常数 R 的非平凡解，而是通过引入谱参数 $u$ 得到 $R(u)$（Baxter 化）。在这种情形下，参数间的代数关系由解析函数（双曲/三角）调节，物理上对应可调的各向异性或耦合强度。

    
  
  

```
num_equations=12
Eq1: 4*re(b)*re(a*d) - 4*re(c)*re(a*d) - 4*im(b)*im(a*d) + 4*im(c)*im(a*d)
Eq2: 4*re(b)*im(a*d) - 4*re(c)*im(a*d) + 4*re(a*d)*im(b) - 4*re(a*d)*im(c)
Eq3: -4*re(b)*re(a*d) + 4*re(c)*re(a*d) + 4*im(b)*im(a*d) - 4*im(c)*im(a*d)
Eq4: -4*re(b)*im(a*d) + 4*re(c)*im(a*d) - 4*re(a*d)*im(b) + 4*re(a*d)*im(c)
Eq5: -4*re(a)*re(b*c) + 4*re(d)*re(b*c) + 4*im(a)*im(b*c) - 4*im(d)*im(b*c)
Eq6: -4*re(a)*im(b*c) + 4*re(d)*im(b*c) - 4*re(b*c)*im(a) + 4*re(b*c)*im(d)
Eq7: -4*re(a*b*c) + 4*re(a*b*d) + 4*re(a*c*d) - 4*re(b*c*d)
Eq8: -4*im(a*b*c) + 4*im(a*b*d) + 4*im(a*c*d) - 4*im(b*c*d)
Eq9: 4*re(a)*re(b*c) - 4*re(d)*re(b*c) - 4*im(a)*im(b*c) + 4*im(d)*im(b*c)
Eq10: 4*re(a)*im(b*c) - 4*re(d)*im(b*c) + 4*re(b*c)*im(a) - 4*re(b*c)*im(d)
Eq11: 4*re(a*b*c) - 4*re(a*b*d) - 4*re(a*c*d) + 4*re(b*c*d)
Eq12: 4*im(a*b*c) - 4*im(a*b*d) - 4*im(a*c*d) + 4*im(b*c*d)
```



### 附录：简单验证

主要是验证微观参数是否按线性关系映射到有效 Kitaev 链，即$ t ∝ (b+c)、Δ ∝ (b−c)、μ ≈ 4d + μ0$，以及模型是否能重现 Kitaev 链的拓扑性质与 Majorana 零模。

使用模型中的$(a,b,c,d,\mu_0)$去构建了 4x4 BdG single-particle Hamiltonian ,然后尝试验证能否复现这些极限情况。

1.k‑space能隙随化学势 μ 的最小正能量曲线，在 |μ|≈2|t| 间隙闭合

![image-20260326152420884](/home/asice-cloud/.config/Typora/typora-user-images/image-20260326152420884.png)

2.以及验证拓扑不变量（winding 在拓扑区取非零值，在相界变为 0；Pfaffian 符号应在相界翻转，两者一致则验证良好）

![image-20260326152706049](/home/asice-cloud/.config/Typora/typora-user-images/image-20260326152706049.png)

3. 绘制测得的有效参数 t_eff 对 (b+c) 以及 Δ_eff 对 (b−c) 的散点并画线性拟合。JSON 包含拟合斜率/截距/R²

![image-20260326153029818](/home/asice-cloud/.config/Typora/typora-user-images/image-20260326153029818.png)

![image-20260326153035767](/home/asice-cloud/.config/Typora/typora-user-images/image-20260326153035767.png)



然后其实本来想结合一下文献，看看能不能用现有模型通过一些变换然后和文献匹配， 昨天使用了Non-Abelian statistics of Majorana zero modes in the presence of an Andreev bound state 这个学长的文献，然后用现在的模型是可以推导出类似的结论，就是减小杂化或缩短交换时间等可降低动力学相位 Θ的累计，提高与理想拓扑交换的重合。但是我这边没法恢复回来，暂时没有复现可恢复高保真的效果。所以这一部分现在没有整理完，之后再换种思路试试怎么做变换。

