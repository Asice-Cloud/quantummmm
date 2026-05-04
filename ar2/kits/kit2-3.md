### 从 1D R 到远程算符与 2D 编织：系统整理

本节把前面分散在 R_to_Kitaev、kit2、kit2-2 中的内容系统整理，回答两个核心问题：

1. 从 1D 上最近邻的
	 $$R_{i,i+1}=aI+b\,\sigma^x_i\sigma^x_{i+1}+c\,\sigma^y_i\sigma^y_{i+1}+d\,\sigma^z_i\sigma^z_{i+1}$$
	 出发，在自旋→费米→Majorana+Z2 的映射下，**远程算符 $R_{ij}$（$j>i+1$）的精确形式是什么**？
2. 在 2D/Z2 框架和量子门语言下，**如何用局域 R 构造“端点 i,j 之间的远程作用”和“编织算符”**？

整体结构：

- 1D：最近邻 R 的 Jordan–Wigner 映射与长程自旋算符的 JW 串形式；
- Majorana+Z2：把 JW 串重写成“端点 Majorana × 路径上的 Z2 链变量乘积”；
- 门序列：用最近邻交换和最近邻 R 拼出只作用在端点 i,j 的复合幺正门；
- 2D：在二维图上用路径 $\gamma$ 推广上述结构，并与编织/世界线的几何图像对应起来。


#### 5.1 记号约定与总体关系

- 自旋 Hilbert 空间中的两体算符 $R_{ij}^{(\text{spin})}$：
	$$
	R_{ij}^{(\text{spin})}=aI+b\,\sigma^x_i\sigma^x_j+c\,\sigma^y_i\sigma^y_j+d\,\sigma^z_i\sigma^z_j,
	$$
	这是最初在自旋基底上定义的两体算符，对任意一对格点 $(i,j)$ 都可以这样写（最近邻 $j=i+1$ 只是其中一个特例）。

- 费米/Majorana 展开下的 $R_{ij}$：用 Jordan–Wigner（1D）或四 Majorana + Z2 规范场（2D）映射，把 $R_{ij}^{(\text{spin})}$ 展开为
$$
	R_{ij} = C_{ij}\,I + D_{ij}^{(\text{dens})} + Q_{ij}^{(2)}(c,\gamma)\,\Xi_{ij},
$$
其中：
- $C_{ij}\,I$ 是常数项；
- $D_{ij}^{(\text{dens})}$ 收集所有“密度/四费米”项，如 $n_in_j,\ n_i,\ n_j$ 等；
- $Q_{ij}^{(2)}(c,\gamma)$ 是只含端点自由度的二次算符（在费米语言里是 $c_i,c_j$ 的 hopping/pairing 组合，在 Majorana 语言里是 $i c_i c_j$ 或其线性组合）；
- $\Xi_{ij}$ 是从 i 到 j 的 Z2 串：在 1D 中是 JW 串 $e^{i\pi\sum_{k=i}^{j-1}n_k}$，在 2D Z2 语言中是路径 $\gamma$ 上链变量的乘积 $u_\gamma$。

- 字符串算符 $B_{ij}(\gamma)$：专指上式中“端点二次项 × 串”的那一块，
	$$
	B_{ij}(\gamma)=u_\gamma\,(i c_i c_j),
	$$
	它可以理解为
	$$
	Q_{ij}^{(2)}(c,\gamma)\,\Xi_{ij}\ \text{在 Majorana+Z2 规范选择下的一个标准代表},
	$$
	一般不包含 $C_{ij}I$ 和 $D_{ij}^{(\text{dens})}$。在 1D 上，$u_\gamma$ 与 JW 串 $e^{i\pi\sum_{k=i}^{j-1}n_k}$ 等价；在 2D 上则显式记录路径/拓扑信息。

- 门意义的 $R^{(\text{gate})}_{ij}$：用最近邻 R 和 SWAP 串成的、在代数上只作用在因子 i,j 上的复合幺正门（见第 3 节公式）。$R_{ij}$、$B_{ij}(\gamma)$ 是“哈密顿量/算符层面”的对象，而 $R^{(\text{gate})}_{ij}$ 是“量子门/演化层面”的对象。

- 编织算符 $U_\gamma$：沿路径 $\gamma$ 乘起的一串局域 Majorana 交换门（见第 4 节），在零模子空间上给出任何子的编织表示。它和 $B_{ij}(\gamma)$ 一样，都依赖路径 $\gamma$，但 $B_{ij}(\gamma)$ 是“静态的字符串算符”，$U_\gamma$ 则是“由时间演化得到的编织算符”。

#### 5.2 1D：自旋→费米下的 $R_{ij}$ 具体形式

对最近邻 $j=i+1$，在 1.1 节已经给出
$$
\begin{aligned}
R_{i,i+1}
&=(b+c)\,(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i)
 +(b-c)\,(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i)\\
&\quad +d\,(4n_in_{i+1}-2(n_i+n_{i+1})+1)+aI,
\end{aligned}
$$
可以识别出：
- 常数项：$aI+d\,I$ 中的一部分；
- 密度/四费米项：$4d\,n_in_{i+1}-2d(n_i+n_{i+1})$；
- 端点二次项：
$$
	Q_{i,i+1}^{(2)}(c)= (b+c)(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i)
	+(b-c)(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i),
$$
对应 Kitaev 链里的 hopping 和 pairing。

对远程 $j>i+1$，由 Jordan–Wigner 映射
$$
\sigma^+_k=c_k^{\dagger}e^{i\pi\sum_{m<k}n_m},\qquad
\sigma^-_k=c_k e^{i\pi\sum_{m<k}n_m}
$$
直接计算得到
$$
\sigma^+_i\sigma^-_j = c_i^{\dagger} c_j\; e^{i\pi\sum_{k=i}^{j-1}n_k},\quad
\sigma^-_i\sigma^+_j = c_i c_j^{\dagger}\; e^{i\pi\sum_{k=i}^{j-1}n_k},
$$
以及
$$
\sigma^+_i\sigma^+_j = c_i^{\dagger} c_j^{\dagger}\; e^{i\pi\sum_{k=i}^{j-1}n_k},\quad
\sigma^-_i\sigma^-_j = c_i c_j\; e^{i\pi\sum_{k=i}^{j-1}n_k}.
$$
再用
$$
\sigma^x=\sigma^++\sigma^-,\qquad \sigma^y=\tfrac{1}{i}(\sigma^+-\sigma^-),
$$
展开 $\sigma^x_i\sigma^x_j,\sigma^y_i\sigma^y_j$。先在自旋升降算符层面整理：
$$
\begin{aligned}
\sigma^x_i\sigma^x_j
&=(\sigma^+_i+\sigma^-_i)(\sigma^+_j+\sigma^-_j)\\
&=\sigma^+_i\sigma^+_j+\sigma^+_i\sigma^-_j+\sigma^-_i\sigma^+_j+\sigma^-_i\sigma^-_j,\\[4pt]
\sigma^y_i\sigma^y_j
&=-\,(\sigma^+_i-\sigma^-_i)(\sigma^+_j-\sigma^-_j)\\
&=-\sigma^+_i\sigma^+_j+\sigma^+_i\sigma^-_j+\sigma^-_i\sigma^+_j-\sigma^-_i\sigma^-_j.
\end{aligned}
$$
因此可以写成
$$
\sigma^x_i\sigma^x_j+\sigma^y_i\sigma^y_j
=2(\sigma^+_i\sigma^-_j+\sigma^-_i\sigma^+_j),\\
\sigma^x_i\sigma^x_j-\sigma^y_i\sigma^y_j
=2(\sigma^+_i\sigma^+_j+\sigma^-_i\sigma^-_j).
$$
把前面的 JW 结果代入，可得四个升降算符组合都带有同一个串 $e^{i\pi\sum_{k=i}^{j-1}n_k}$，从而
$$
\begin{aligned}
b\,\sigma^x_i\sigma^x_j+c\,\sigma^y_i\sigma^y_j
&=\tfrac{b+c}{2}(\sigma^x_i\sigma^x_j+\sigma^y_i\sigma^y_j)
 +\tfrac{b-c}{2}(\sigma^x_i\sigma^x_j-\sigma^y_i\sigma^y_j)\\
&= (b+c)(\sigma^+_i\sigma^-_j+\sigma^-_i\sigma^+_j)
 +(b-c)(\sigma^+_i\sigma^+_j+\sigma^-_i\sigma^-_j)\\
&=\Big[(b+c)(\text{hopping 二次项})+(b-c)(\text{pairing 二次项})\Big] e^{i\pi\sum_{k=i}^{j-1}n_k},
\end{aligned}
$$
其中“hopping 二次项”由 $c_i^{\dagger}c_j, c_j^{\dagger}c_i$ 组成，“pairing 二次项”由 $c_i^{\dagger}c_j^{\dagger}, c_j c_i$ 组成。更具体地，可以取
$$
Q_{ij}^{(2)}(c)=(b+c)(c_i^{\dagger}c_j+c_j^{\dagger}c_i)+(b-c)(c_i^{\dagger}c_j^{\dagger}+c_j c_i),
$$
而所有由于重排产生的符号可以吸收到 $C_{ij},D_{ij}^{(\text{dens})}$ 或规范选择中。再加上 $aI$ 与 $d\,\sigma^z_i\sigma^z_j$ 产生的常数和密度/四费米项，就得到远程 $R_{ij}^{(\text{spin})}$ 可写为
$$
R_{ij}=C_{ij}I+D_{ij}^{(\text{dens})}+Q_{ij}^{(2)}(c)\,e^{i\pi\sum_{k=i}^{j-1}n_k},\qquad j>i+1,
$$
其中 $Q_{ij}^{(2)}(c)$ 是由 $c_i^{\dagger}c_j,c_j^{\dagger}c_i,c_i^{\dagger}c_j^{\dagger},c_j c_i$ 线性组合成的二次算符。这就把“完整的远程 $R_{ij}$”拆成了前面总结中的 $C_{ij},D_{ij}^{(\text{dens})},Q_{ij}^{(2)}(c),\Xi_{ij}$ 四部分，其中
$$
\Xi_{ij}=e^{i\pi\sum_{k=i}^{j-1}n_k}
$$
就是 1D 中的 Z2 串。

从这个角度看，1D 中的“字符串算符部分”可以写成
$$
R_{ij}^{(\text{string-part})}=Q_{ij}^{(2)}(c)\,e^{i\pi\sum_{k=i}^{j-1}n_k},
$$
而 $C_{ij}I+D_{ij}^{(\text{dens})}$ 则对拓扑/编织性质通常不起决定性作用（只是改变能量或引入相互作用）。

从完全类似的计算可以看出：即便在最近邻情形（$j=i+1$，此时 JW 串退化为 1），相邻的二次项 $Q_{i,i+1}^{(2)}(c)$ 与 $Q_{i+1,i+2}^{(2)}(c)$ 也是**一般不对易**的。例如在最简单的 hopping 模型中，取
$$
H_{12}=c_1^{\dagger}c_2+c_2^{\dagger}c_1,\qquad
H_{23}=c_2^{\dagger}c_3+c_3^{\dagger}c_2,
$$
直接计算可得
$$
[H_{12},H_{23}] = c_1^{\dagger}c_3 - c_3^{\dagger}c_1\neq 0,
$$
说明“共享端点”的最近邻二次项在费米语言中同样是不对易的。这与我们在 Majorana 语言中得到的 $[\gamma_1\gamma_2,\gamma_2\gamma_3]\neq 0$ 完全一致，只是用不同基底表述而已。

#### 5.3 2D：Majorana+Z2 语言下的 $B_{ij}(\gamma)$ 与 $R_{ij}$

在 2D Z2+Majorana 框架中，每条最近邻 a 型键 $(p,q)$ 上有
$$
\sigma^a_p\sigma^a_q = u^{(a)}_{pq}(i c_p c_q),
$$
若某一类键上的耦合常数来自于上面的 $R_{pq}^{(\text{spin})}$ 中的 $b,c$ 项，则其 Majorana 部分正是 $i c_p c_q$ 的线性组合。对一条从 i 到 j 的路径 $\gamma$，定义
$$
u_\gamma=\prod_{(k,l)\in\gamma}u_{kl},
$$
则我们抽象出远程 Majorana 键
$$
B_{ij}(\gamma)=u_\gamma(i c_i c_j),
$$
它在 1D 上对应于
$$
Q_{ij}^{(2)}(c)\,e^{i\pi\sum_{k=i}^{j-1}n_k}
$$
这一部分；在 2D 上则是“端点 Majorana × Z2 串”的自然推广，用来抓住 $R_{ij}$ 在拓扑/零模子空间中最重要的那一块。

从这个意义上，我们可以把 Majorana+Z2 语言下的 $R_{ij}$ 写成
$$
R_{ij}=C_{ij}I+D_{ij}^{(\text{dens})}+\lambda_{ij}\,B_{ij}(\gamma)+\cdots,
$$
其中 $\lambda_{ij}$ 由 $a,b,c,d$ 决定，“$\cdots$” 代表可能存在的其他高阶或不同类型的端点耦合，而 $B_{ij}(\gamma)$ 就是那一类最简单的“字符串 × 端点 Majorana”结构。后文在讨论远程交换和编织时，主要关心的正是这一部分。

#### 5.4 几种远程对象的角色对比

在上述具体表达式基础上，可以总结：

1) **1D 自旋→费米：带 JW 串的远程 $R_{ij}$**  
完整的 $R_{ij}$ 在 1D JW 语言中总可以写成
$$
R_{ij}=C_{ij}I+D_{ij}^{(\text{dens})}+Q_{ij}^{(2)}(c)\,e^{i\pi\sum_{k=i}^{j-1}n_k},
$$
其中最后一项就是“远程字符串部分”。

2) **Majorana+Z2：路径字符串算符 $B_{ij}(\gamma)$（$R_{ij}$ 的二次部分）**  
在引入链路 Z2 变量 $u_{kl}$ 后，对从 i 到 j 的路径 $\gamma$ 定义
$$
u_\gamma = \prod_{(k,l)\in\gamma} u_{kl},
$$
则对应的长程 Majorana 键写作
$$
B_{ij}(\gamma) = u_\gamma\,(i c_i c_j).
$$
这在 1D 上与“$R_{ij}$ 的 JW 串部分”等价，即 $R_{ij}$ 被展开到费米/Majorana 语言后，其端点二次项乘以 JW 串的那一部分；在 2D 上则自然编码了路径/拓扑依赖性，是本工作中讨论远程算符与编织时的标准形式。

3) **门序列意义的远程两体门 $R^{(\text{gate})}_{ij}$**  
若底层硬件/模型只允许最近邻 R 或 SWAP，远程 i,j 上的两体门可以通过交换序列构造：
$$
S_{i\to j}=S_{j-1,j} S_{j-2,j-1}\cdots S_{i,i+1},
$$
$$
R^{(\text{gate})}_{ij}=S_{i\to j}\,R_{j-1,j}\,S_{i\to j}^{\dagger}.
$$
这一定义在算符代数上等价于“只作用在第 i、第 j 个因子的两体变换”，但在物理实现上是由一串最近邻门组合而成的复合幺正。

4) **2D 上沿路径的编织算符 $U_\gamma$**  
在 2D Z2+Majorana 框架中，沿路径 $\gamma$ 的“局域交换”可以选取为局域 Majorana 双线性上的演化
$$
U_{kl}^{(\text{braid})} \approx \exp\Big(\pm\frac{\pi}{4}\,\gamma_k\gamma_l\Big),
$$
则整条路径的编织算符为
$$
U_\gamma=\prod_{(k,l)\in\gamma} U_{kl}^{(\text{braid})}.
$$
在零模简并子空间中，$U_\gamma$ 与字符串算符 $B_{ij}(\gamma)$ 一起，描述了“沿路径搬运/绕圈”对态空间的非平庸作用，两者是本工作中讨论远程算符与 2D 编织的两个互补视角。

用
$$
\sigma^x=\sigma^++\sigma^-,\qquad \sigma^y=\frac{\sigma^+-\sigma^-}{i}
$$

展开，可见远程 $\sigma^x_i\sigma^x_j,\sigma^y_i\sigma^y_j$ 的结构和最近邻时完全类似，只是每一项都被 **同一个串** 乘上：

$$
\sigma^x_i\sigma^x_j\sim(\text{hopping 二次项})\times e^{i\pi\sum_{k=i}^{j-1}n_k},
$$
$$
\sigma^y_i\sigma^y_j\sim(\text{pairing 二次项})\times e^{i\pi\sum_{k=i}^{j-1}n_k}.
$$

因此，对于 1D 链上远距 $i<j$，自旋层面的
$$
R_{ij}=aI+b\,\sigma^x_i\sigma^x_j+c\,\sigma^y_i\sigma^y_j+d\,\sigma^z_i\sigma^z_j
$$
在费米表示下具有“端点二次算符 × 中间 JW 串”的统一结构：
$$
R_{ij} \sim (\text{端点的 hopping/ pairing / 密度})\; e^{i\pi\sum_{k=i}^{j-1}n_k}.
$$

这给出了 1D 上远程算符的精确形式：
- **端点**：和最近邻情况一样，是 $c_i,c_j$ 的二次组合（跳跃/配对/密度）；
- **路径因子**：从 i 到 j 的 JW 串 $e^{i\pi\sum_{k=i}^{j-1}n_k}$，相当于一条 Z2 路径上的字符串。


 

### 2. Majorana + Z2 视角：把串重写成链变量乘积

在 2D/更一般格子上，更方便的语言是 Majorana + Z2 规范场（kit2 中已详细给出）：

- 在每个格点 j 上引入四个 Majorana：$b^x_j,b^y_j,b^z_j,c_j$；
- 映射自旋算符为
	$$
	\sigma^a_j = i\,b^a_j c_j,\qquad a=x,y,z;
	$$
- 对一条类型为 a 的键 $(i,j)$ 定义链路变量
	$$
	u^{(a)}_{ij}=-i b^a_i b^a_j,\qquad (u^{(a)}_{ij})^2=1;
	$$
- 则键项可以写成
	$$
	\sigma^a_i\sigma^a_j = u^{(a)}_{ij}\,(i c_i c_j).
	$$

在最近邻情形下，这正好把原来的自旋两体项变成“Z2 链变量 × 端点 Majorana 的二次项”，哈密顿量可以写为
$$
H=\frac i4\sum_{ij} A_{ij}(u)\,c_i c_j,
$$
其中 $A_{ij}(u)$ 由 $J^{(a)}_{ij}u^{(a)}_{ij}$ 组成，是实反对称矩阵。

要把 1D 中的 JW 串重写成 Z2 链变量的乘积，可以：

- 选一条从 i 到 j 的路径 $\gamma$（在 1D 上就是 $(i,i+1,\dots,j)$）；
- 在路径每条最近邻边 $(k,l)$ 上都有对应的 Z2 链变量 $u_{kl}$；
- 定义沿路径的乘积
	$$
	u_\gamma = \prod_{(k,l)\in\gamma} u_{kl}.
	$$

在合适的规范选择和映射下，可以把 JW 串 $e^{i\pi\sum_{k=i}^{j-1}n_k}$ 重新解释为这样一个 $u_\gamma$（或者与之只差一个规范变换/整体因子）。于是 1D 上的远程算符可以用与 2D 同样的“端点 + 路径串”形式来写：
$$
B_{ij}(\gamma) = u_\gamma\,(i c_i c_j),
$$
其中 $\gamma$ 是从 i 到 j 的路径；在 1D 上只有一条最短路径，在 2D 上则有成族的路径（不同同伦类给出物理上不同的作用，与编织有关）。

这就是 kit2 中反复出现的结构：
- **局域最近邻**：$\sigma^a_i\sigma^a_j=u^{(a)}_{ij}(i c_i c_j)$；
- **远程键/字符串算符**：选择路径 $\gamma$，沿路径把最近邻链变量乘起来得到 $u_\gamma$，再与端点 Majorana 二次项相乘。


 

### 3. 门序列意义上的远程 $R_{ij}$：用最近邻交换实现

上一节讨论的是“算符/哈密顿量层面”的远程字符串算符。若从“操作（量子门）”角度出发，只允许最近邻 R 或 SWAP，也可以构造一个在自旋 Hilbert 空间里只作用在端点 i,j 的复合幺正门 $R^{(\text{gate})}_{ij}$。

设链的 Hilbert 空间为 $V^{\otimes L}$，抽象 R 作用在 $V\otimes V$ 上。最近邻嵌入为
$$
R_{k,k+1} = I^{\otimes(k-1)}\otimes R\otimes I^{\otimes(L-k-1)}.
$$

令 $S_{k,k+1}$ 为相邻两个格点的交换门（SWAP），在张量积上对应置换第 k 和 k+1 个因子。对 $i<j$，定义将第 i 个格点搬到 j 位置的交换序列
$$
S_{i\to j}=S_{j-1,j} S_{j-2,j-1}\cdots S_{i,i+1}.
$$
则复合门
$$
R^{(\text{gate})}_{ij}=S_{i\to j}\,R_{j-1,j}\,S_{i\to j}^{\dagger}
$$
在算符代数层面等价于一个“只作用在第 i、第 j 个因子上的两体算符”，但它是由一串最近邻操作拼出来的 gate，而不是哈密顿量中的单一局域项。这就是 1D 模型里在**操作意义**上实现远程 R 的自然方式。

这一构造可以在 2D 图上完全照搬：
- 选一条 2D 最近邻路径 $\gamma:i=i_0\to i_1\to\cdots\to i_n=j$；
- 把交换门 S 和最近邻 $R_{i_k i_{k+1}}$ 沿路径有序作用；
- 得到的复合幺正 $R^{(\text{gate})}_{ij}(\gamma)$ 只改变端点附近自旋自由度，其余格点的态通过交换“让路”。


 

### 4. 2D 上的路径、远程算符与编织

在 2D + Z2 + Majorana 框架中，上述两种远程构造自然合并到“路径 $\gamma$ 的几何信息”里：

- **Hamiltonian/算符意义**：对给定路径 $\gamma$，定义字符串算符
	$$
	B_{ij}(\gamma) = u_\gamma\,(i c_i c_j),\qquad u_\gamma = \prod_{(k,l)\in\gamma} u_{kl}.
	$$
	不同路径（尤其是绕行与不绕行其他拓扑缺陷）的 $u_\gamma$ 在零模子空间中给出不同的作用，对应世界线的不同绕缠数。

- **量子门/演化意义**：在每条边 $(k,l)\in\gamma$ 上，通过选择合适的 R 参数与演化时间，可以实现局域的 Majorana 交换门
	$$
	U_{kl}^{(\text{braid})} \approx \exp\Big(\pm\frac{\pi}{4}\,\gamma_k\gamma_l\Big),
	$$
	整条路径的编织算符为有序乘积
	$$
	U_\gamma = \prod_{(k,l)\in\gamma} U_{kl}^{(\text{braid})},
	$$
	在零模简并子空间中实现任何子的编织。不同同伦类的路径 $\gamma$ 给出互不相同的 $U_\gamma$，这就是非阿贝尔统计的来源。

综上：

- 在 1D 上，远程自旋算符不可避免地带有 JW 串 $e^{i\pi\sum_{k=i}^{j-1}n_k}$；
- 在 Majorana+Z2 语言中，这个串自然重写为路径上的 Z2 链变量乘积 $u_\gamma$，从而得到“端点 Majorana × Z2 字符串”的长程算符；
- 在操作层面，可以用最近邻交换和最近邻 R 串成等效的远程 gate $R^{(\text{gate})}_{ij}$；
- 在 2D 上，这些结构与任意子编织的路径 $\gamma$ 一一对应：远程算符和编织算符都变成“沿路径收集的 Z2/相位信息 × 端点上的局域作用”。

这完成了从 1D 最近邻 R→远程算符→2D Z2+Majorana→编织的系统整理，可作为 kit2 与 kit2-2 的补充说明与“远程算符”推导的集中展示。


 

### 5. 小结：几种标准的远程算符形式

为了便于在后文直接引用，这里把前面几种得到“远程作用”的典型形式集中列出。

先统一约定记号：

- 自旋层面的 $R_{ij}$：指
$$
	R_{ij}^{(	ext{spin})}=aI+b\,\sigma^x_i\sigma^x_j+c\,\sigma^y_i\sigma^y_j+d\,\sigma^z_i\sigma^z_j,
$$
这是最初在自旋 Hilbert 空间上定义的两体算符（无论 $j=i+1$ 还是 $j>i+1$）。

- 费米/Majorana 展开下的 $R_{ij}$：指把 $R_{ij}^{(\text{spin})}$ 用 JW 或四 Majorana 映射展开后得到的算符，其中一般包含：
	- 端点的二次项（$c_i,c_j$ 或 $\gamma$ 的 hopping/ pairing 型耦合），
	- 密度和四费米项，
	- 再乘以一条从 i 到 j 的 Z2 串（1D 中是 JW 串，2D 中是 $u_\gamma$）。

- 字符串算符 $B_{ij}(\gamma)$：专指“端点 Majorana 的二次项 × 路径上的 Z2 串”这一块，
$$
	B_{ij}(\gamma)=u_\gamma\,(i c_i c_j),
$$
它等于 $R_{ij}$ 在零模/拓扑子空间中最重要的那一部分，但一般不包括 $R_{ij}$ 里的恒等和密度项。

- 门意义的 $R^{(\text{gate})}_{ij}$：用最近邻 R 和 SWAP 串成的、在代数上只作用在因子 i,j 上的复合幺正门。

- 编织算符 $U_\gamma$：沿路径 $\gamma$ 乘起的一串局域 Majorana 交换门，在零模子空间上给出任何子的编织表示。

1) **1D 自旋→费米：带 JW 串的远程 $R_{ij}$**  
对 $i<j$，自旋层面的

$$
R_{ij}=aI+b\,\sigma^x_i\sigma^x_j+c\,\sigma^y_i\sigma^y_j+d\,\sigma^z_i\sigma^z_j
$$

在 Jordan–Wigner 映射下具有“端点二次算符 × 中间串”的结构：
$$
R_{ij} \sim (\text{端点的 hopping/pairing/密度})\; e^{i\pi\sum_{k=i}^{j-1}n_k}.
$$
这里 $e^{i\pi\sum_{k=i}^{j-1}n_k}$ 是从 i 到 j 的 JW 串，是 1D 上最直接的远程算符因子。

2) **Majorana+Z2：路径字符串算符 $B_{ij}(\gamma)$（$R_{ij}$ 的二次部分）**  
在引入链路 Z2 变量 $u_{kl}$ 后，对从 i 到 j 的路径 $\gamma$ 定义
$$
u_\gamma = \prod_{(k,l)\in\gamma} u_{kl},
$$
则对应的长程 Majorana 键写作
$$
B_{ij}(\gamma) = u_\gamma\,(i c_i c_j).
$$
这在 1D 上与“$R_{ij}$ 的 JW 串部分”等价，即 $R_{ij}$ 被展开到费米/Majorana 语言后，其端点二次项乘以 JW 串的那一部分；在 2D 上则自然编码了路径/拓扑依赖性，是本工作中讨论远程算符与编织时的标准形式。

3) **门序列意义的远程两体门 $R^{(\text{gate})}_{ij}$**  
若底层硬件/模型只允许最近邻 R 或 SWAP，远程 i,j 上的两体门可以通过交换序列构造：
$$
S_{i\to j}=S_{j-1,j} S_{j-2,j-1}\cdots S_{i,i+1},
$$
$$
R^{(\text{gate})}_{ij}=S_{i\to j}\,R_{j-1,j}\,S_{i\to j}^{\dagger}.
$$
这一定义在算符代数上等价于“只作用在第 i、第 j 个因子的两体变换”，但在物理实现上是由一串最近邻门组合而成的复合幺正。

4) **2D 上沿路径的编织算符 $U_\gamma$**  
在 2D Z2+Majorana 框架中，沿路径 $\gamma$ 的“局域交换”可以选取为局域 Majorana 双线性上的演化
$$
U_{kl}^{(\text{braid})} \approx \exp\Big(\pm\frac{\pi}{4}\,\gamma_k\gamma_l\Big),
$$
则整条路径的编织算符为
$$
U_\gamma=\prod_{(k,l)\in\gamma} U_{kl}^{(\text{braid})}.
$$
在零模简并子空间中，$U_\gamma$ 与字符串算符 $B_{ij}(\gamma)$ 一起，描述了“沿路径搬运/绕圈”对态空间的非平庸作用，两者是本工作中讨论远程算符与 2D 编织的两个互补视角。


 

### 6. 代数验证与辫子群条件下的 $a,b,c,d$

本节对上面远程算符和编织门的两个关键性质给出**严格代数证明**，并在此基础上说明：在“远程交换可交换”的辫子群条件下，对 $a,b,c,d$ 没有额外约束；真正约束 $a,b,c,d$ 的是局域 YBE/近邻辫子关系。

#### 6.1 远程字符串算符 $B_{ij}(\gamma)$ 的对易性

**命题 1.** 设 $B_{ij}(\gamma)=u_\gamma(i c_i c_j)$，其中 $u_\gamma=\prod_{(k,l)\in\gamma}u_{kl}$，$u_{kl}$ 两两对易且 $u_{kl}^2=1$，$c_m$ 是满足
$$
\{c_m,c_n\}=2\delta_{mn}
$$
的 Majorana 算符。若两条路径 $\gamma,\gamma'$ 的端点集合与边集合完全不相交，则
$$
[B_{ij}(\gamma),B_{kl}(\gamma')]=0.
$$

**证明.** 在假设下，$\{i,j\}\cap\{k,l\}=\varnothing$，且 $u_\gamma,u_{\gamma'}$ 是由互不相交边上的 $u_{mn}$ 乘积组成。于是：

1. Z2 部分：所有 $u_{mn}$ 彼此对易，且与所有 $c_p$ 对易（这是 Kitaev 模型中 $u_{mn}$ 的定义性质），因此
$$
[u_\gamma,u_{\gamma'}]=0,\qquad [u_\gamma, i c_k c_l]=[u_{\gamma'}, i c_i c_j]=0.
$$

2. Majorana 端点部分：记 $X=i c_i c_j$，$Y=i c_k c_l$。利用 $c_m^2=1$ 及 $c_m c_n=-c_n c_m$（$m\neq n$），直接计算
$$
XY = -c_i c_j c_k c_l,\qquad YX = -c_k c_l c_i c_j.
$$
由于四个指标两两不同，任意交换一次会引入一个负号，总共需要偶数次交换（4 次）把 $c_k c_l c_i c_j$ 排列成 $c_i c_j c_k c_l$，所以
$$
c_k c_l c_i c_j = (+1)\,c_i c_j c_k c_l.
$$
于是 $XY=YX$，即 $[X,Y]=0$，从而
$$
[i c_i c_j, i c_k c_l]=0.
$$

综上，由 $B_{ij}(\gamma)=u_\gamma X$、$B_{kl}(\gamma')=u_{\gamma'}Y$ 与上述对易关系立刻得到
$$
[B_{ij}(\gamma),B_{kl}(\gamma')]=0.
$$
命题得证。$\square$

此命题即“远程字符串算符可交换”的严格表达。它完全基于支撑不相交和 Z2 变量、Majorana 算符的代数关系，与 $a,b,c,d$ 的具体数值无关。

**补充：近邻字符串算符的非对易性（与上式相对）**  
上面的命题只适用于两条路径及端点**完全不相交**的情形。若两条字符串算符共享一个端点（例如 $B_{12}(\gamma_1)$ 与 $B_{23}(\gamma_2)$），则一般**不对易**。设
$$
B_{12}=u_{12}(i c_1 c_2),\qquad B_{23}=u_{23}(i c_2 c_3),
$$
其中 $u_{12},u_{23}$ 与所有 $c_j$ 对易，且彼此对易。则
$$
[B_{12},B_{23}] = u_{12}u_{23}\,[i c_1 c_2, i c_2 c_3].
$$
令 $\gamma_1=c_1,\gamma_2=c_2,\gamma_3=c_3$，上式右边正是
$$
u_{12}u_{23}\,[\gamma_1\gamma_2,\gamma_2\gamma_3],
$$
而在 6.2 节中我们已经显式计算出（见后文式 (6.2.")）
$$
[\gamma_1\gamma_2,\gamma_2\gamma_3]=2\gamma_1\gamma_3\neq 0.
$$
因此
$$
[B_{12},B_{23}]\neq 0,
$$
这说明：
- 对应**不相交**路径/端点，$B$ 算符可交换（命题 1）；
- 对应**共享端点**的近邻路径，$B$ 算符本身就是不对易的，其非对易性完全来源于端点 Majorana 双线性的非对易性，而 Z2 串因子只是给出一个对易的整体因子。

#### 6.2 Majorana 交换门的辫子关系

现在考虑 Majorana 交换门
$$
U_{ab}=\exp\Big(\theta\,\gamma_a\gamma_b\Big),
$$
其中 $\theta=\pm\frac{\pi}{4}$，$\gamma_m$ 为 Majorana 算符：$\{\gamma_m,\gamma_n\}=2\delta_{mn}$。设有三个位点 1,2,3，对应三个 Majorana $\gamma_1,\gamma_2,\gamma_3$，定义
$$
U_1:=U_{12}=\exp\Big(\theta\,\gamma_1\gamma_2\Big),\qquad
U_2:=U_{23}=\exp\Big(\theta\,\gamma_2\gamma_3\Big).
$$

我们证明：

**命题 2.** 取 $\theta=\pm\tfrac{\pi}{4}$ 时，有
$$
U_1 U_2 U_1 = U_2 U_1 U_2,
$$
即 $U_1,U_2$ 给出辫子群 $B_3$ 的一个表示。

证明分两步：先求共轭作用，再比较两侧的作用是否一致。

**(i) $U_{ab}$ 的共轭作用是旋转。**  
记 $K_{ab}:=\gamma_a\gamma_b$。注意 $K_{ab}^\dagger=-K_{ab}$，且 $K_{ab}^2=-1$。对 $\gamma_a$ 有
$$
\begin{aligned}
e^{\theta K_{ab}}\gamma_a e^{-\theta K_{ab}}
&= \gamma_a + \theta[K_{ab},\gamma_a] + \frac{\theta^2}{2!}[K_{ab},[K_{ab},\gamma_a]]+\cdots.
\end{aligned}
$$
由反对易关系可算得
$$
[K_{ab},\gamma_a]=[\gamma_a\gamma_b,\gamma_a]=-2\gamma_b,\\
[K_{ab},\gamma_b]=[\gamma_a\gamma_b,\gamma_b]=2\gamma_a,
$$
以及二重对易
$$
[K_{ab},[K_{ab},\gamma_a]]=[K_{ab},-2\gamma_b]=-4\gamma_a,
$$
继续下去得到一个二维线性子空间 $\mathrm{span}\{\gamma_a,\gamma_b\}$ 中的旋转。整理幂级数可得
$$
e^{\theta K_{ab}}\gamma_a e^{-\theta K_{ab}}
= \gamma_a\cos(2\theta)+\gamma_b\sin(2\theta),\\
e^{\theta K_{ab}}\gamma_b e^{-\theta K_{ab}}
= \gamma_b\cos(2\theta)-\gamma_a\sin(2\theta),
$$
而对与 $a,b$ 都不同的 $\gamma_c$ 有 $[K_{ab},\gamma_c]=0$，故
$$
e^{\theta K_{ab}}\gamma_c e^{-\theta K_{ab}}=\gamma_c,\qquad c\neq a,b.
$$

取 $\theta=\pm\tfrac{\pi}{4}$ 时，$2\theta=\pm\tfrac{\pi}{2}$，于是
$$
\cos(2\theta)=0,\quad \sin(2\theta)=\pm1.
$$
从而
$$
U_{ab}\gamma_a U_{ab}^\dagger = \pm\gamma_b,\\
U_{ab}\gamma_b U_{ab}^\dagger = \mp\gamma_a,\\
U_{ab}\gamma_c U_{ab}^\dagger = \gamma_c\ (c\neq a,b).
$$
这表明 $U_{ab}$ 在 $(\gamma_a,\gamma_b)$ 平面内实现一个 $\pm\tfrac{\pi}{2}$ 的旋转，对其他 $\gamma_c$ 不作用。

**(ii) $U_1,U_2$ 满足辫子关系。**  
只需比较 $U_1 U_2 U_1$ 与 $U_2 U_1 U_2$ 在基 $\{\gamma_1,\gamma_2,\gamma_3\}$ 上的共轭作用是否一致。利用上式：

- 先看 $U_1 U_2 U_1$ 作用下的 $\gamma_1$：
$$
U_1:\ \gamma_1\mapsto\pm\gamma_2;\\
U_2:\ \gamma_2\mapsto\pm\gamma_3;\\
U_1:\ \gamma_3\mapsto\gamma_3\ (\gamma_3 \text{ 与 }\gamma_1,\gamma_2 \text{ 无关}).
$$
整体上 $\gamma_1$ 被送到 $\pm\gamma_3$（符号是可整体吸收的相位）。

- 再看 $U_2 U_1 U_2$ 作用下的 $\gamma_1$：
$$
U_2:\ \gamma_1\mapsto\gamma_1;\\
U_1:\ \gamma_1\mapsto\pm\gamma_2;\\
U_2:\ \gamma_2\mapsto\pm\gamma_3.
$$
同样得到 $\gamma_1\mapsto\pm\gamma_3$。

对 $\gamma_2,\gamma_3$ 做同样计算，可验证 $U_1 U_2 U_1$ 与 $U_2 U_1 U_2$ 在三维实向量空间 $\mathrm{span}\{\gamma_1,\gamma_2,\gamma_3\}$ 上诱导的线性变换完全一致（仅差一个整体相位）。由于 Majorana 算符生成的算符代数是忠实表示，这意味着两侧算符本身也只可能相差一个总体 U(1) 相位，而在辫子群表示中允许这样的整体相位，因此可写为
$$
U_1 U_2 U_1 = e^{i\phi} U_2 U_1 U_2.
$$
通过直接比较两侧在 Fock 空间中一组规范基矢上的作用（或要求表示是单值的）可以进一步固定 $\phi=0$，从而得到严格的辫子关系
$$
U_1 U_2 U_1 = U_2 U_1 U_2.
$$
命题得证。$\square$

当两对 Majorana 完全不相交时（如 $U_{12}$ 与 $U_{34}$），由 (i) 可知它们在各自二维子空间内旋转且生成元对易，因而
$$
[\gamma_1\gamma_2,\gamma_3\gamma_4]=0 \;\Rightarrow\; U_{12}U_{34}=U_{34}U_{12},
$$
对应辫子群中 $|i-j|>1$ 时生成元对易的“远程可交换”条件。


相反地，当两对共享一个端点（例如 $\gamma_1\gamma_2$ 与 $\gamma_2\gamma_3$）时，生成元是**不对易**的。设
$$
A = \gamma_1\gamma_2,\qquad B = \gamma_2\gamma_3,
$$
直接计算有
$$
AB = \gamma_1\gamma_2\gamma_2\gamma_3 = \gamma_1\gamma_3,
$$
而
$$
\begin{aligned}
BA &= \gamma_2\gamma_3\gamma_1\gamma_2
	= \gamma_2(-\gamma_1\gamma_3)\gamma_2\\
	&= -\gamma_2\gamma_1\gamma_3\gamma_2
	= -(-\gamma_1\gamma_2)\gamma_3\gamma_2\\
	&= \gamma_1\gamma_2\gamma_3\gamma_2
	= \gamma_1\gamma_2(-\gamma_2\gamma_3)\\
	&= -\gamma_1\gamma_2\gamma_2\gamma_3
	= -\gamma_1\gamma_3.
\end{aligned}
$$
于是
$$
[A,B]=AB-BA=2\gamma_1\gamma_3\neq 0,
$$
这表明近邻生成元 $\gamma_1\gamma_2$ 与 $\gamma_2\gamma_3$ 不对易，从而对应的 $U_{12},U_{23}$ 也不对易，这正是局域辫子关系非平凡的来源。
