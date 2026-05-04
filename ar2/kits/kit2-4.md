### 从 1D R 到 2D 模型与编织：完整推导

本节在 `R_to_Kitaev`、`kit2`、`kit2-2`、`kit2-3` 的基础上，系统、一步一步推导：

1. 把 1D 上的常数两体算符
	$$
	R_{ij}^{(\text{spin})}=aI+b\,\sigma^x_i\sigma^x_j+c\,\sigma^y_i\sigma^y_j+d\,\sigma^z_i\sigma^z_j
	$$
	推广到一个 **二维格子上的局域哈密顿量**；
2. 在 Majorana+Z2 语言下写出 2D 哈密顿量的具体形式，并说明它如何“逐边”继承 1D 时 R 与 $(t,\Delta,\mu)$ 的对应关系；
3. 在这个 2D 模型内部，严格推导**远程算符**和**编织算符**的形式，并说明它们与 1D JW 串和 1D 远程 R\_{ij} 的关系。

整个思路保持和 1D 完全平行：

- 1D：点集是 $\{1,2,\dots,L\}$，边集只含 $(i,i+1)$，在每条边上放同一个两体算符 $R^{(\text{spin})}$，得到 1D 模型;
- 2D：点集是 2D 格子上的所有格点，边集是二维图上的所有最近邻边 $(i,j)$，同样在每条最近邻边上放同一个 $R^{(\text{spin})}$，得到 2D 模型；
- 远程：在 1D 上，远程 R\_{ij} 会自动带上 JW 串；在 2D 上，远程算符写成“端点 Majorana × 路径上 Z2 链变量乘积”，并在零模子空间中实现编织。

下面依次推导。

 

### 1. 把 1D 的 R 嵌入到 2D 哈密顿量：图论视角

先用图的语言重新表述 1D 的情形：

- 1D 链对应一条路径图 $\mathcal{G}_{1D}=(V_{1D},E_{1D})$，其中
	- $V_{1D}=\{1,2,\dots,L\}$；
	- $E_{1D}=\{(i,i+1)\mid i=1,\dots,L-1\}$；
- 1D 自旋哈密顿量可以写为
$$
	H_{1D}^{(\text{spin})} = \sum_{(i,j)\in E_{1D}} R_{ij}^{(\text{spin})},
	$$
其中每条边上都放同一个两体算符 $R^{(\text{spin})}$，只是支撑位置不同（嵌入到 Hilbert 空间的不同因子上）。

**推广到 2D 的关键点在于：**

- 不是“在 1D 上先推导一个远程 R\_{ij}，再把这个 R\_{ij} 拿去当 2D 的哈密顿量”，
- 而是**直接在二维图的每一条最近邻边上放同一个 $R^{(\text{spin})}$**，把 1D 的“线状图”换成二维图即可。

具体地，对一个给定的二维图 $\mathcal{G}_{2D}=(V_{2D},E_{2D})$（可以是方格、蜂窝等）：

- 每个格点 $i\in V_{2D}$ 搭一个自旋 Hilbert 空间 $\mathbb{C}^2$；
- 对每条最近邻边 $(i,j)\in E_{2D}$（例如方格上的水平/竖直边，蜂窝上的三种键），嵌入
$$
	R^{(\text{spin})}_{ij}
	=I^{\otimes (i-1)}\otimes R^{(\text{spin})}\otimes I^{\otimes (N-2-i)},
$$
即只作用在自旋 i 和 j 上，其余为恒等；

- 2D 自旋哈密顿量定义为
$$
	H_{2D}^{(\text{spin})}
	= \sum_{(i,j)\in E_{2D}} R_{ij}^{(\text{spin})}.
$$

此时，在每条边上的局域作用与 1D 情形完全相同，只是边的集合 $E$ 不再是“线”，而是 2D 图。

**与 1D 映射的关系：**

- 1D 中，我们把 $R_{i,i+1}^{(\text{spin})}$ Jordan–Wigner 到费米语言，识别出最近邻的 hopping、pairing、密度项，从而得到 $(t,\Delta,\mu)$；
- 在 2D 中，只是把这种“沿边的局域两体作用”从 1D 的边集 $E_{1D}$ 推广到 2D 的边集 $E_{2D}$，因此**每条边上的 $(t,\Delta,\mu)$ 与 1D 映射完全相同**，区别只在于动量空间从一维变成二维、动量函数 $\epsilon(\mathbf{k}),\Delta(\mathbf{k})$ 的结构更复杂。

为了更清楚，我们在下一节把 2D 中的 Majorana+Z2 形式写出来，再回头对比 1D 的映射。

 

### 2. 2D Majorana+Z2 映射：逐边继承 1D 的 R→$(t,\Delta,\mu)$ 对应

在 `kit2.md` 中我们已经给出 2D Kitaev 蜂窝模型的标准 Majorana+Z2 表述，这里用更接近本节记号的语言再组织一遍。

#### 2.1 四 Majorana 表示与链路变量

在每个格点 $j\in V_{2D}$ 引入四个 Majorana 算符 $b^x_j,b^y_j,b^z_j,c_j$，满足
$$
\{b^\alpha_j,b^\beta_k\}=2\delta_{jk}\delta^{\alpha\beta},\quad
\{c_j,c_k\}=2\delta_{jk},\quad
\{b^\alpha_j,c_k\}=0.
$$
自旋算符表示为
$$
\sigma^a_j = i\,b^a_j c_j,\qquad a=x,y,z.
$$

对每条 a 型键 $(i,j)\in E_{2D}^{(a)}$（例如蜂窝三种键 x,y,z）定义 Z2 链路变量
$$
u^{(a)}_{ij}=-i b^a_i b^a_j,\qquad (u^{(a)}_{ij})^2=1,
$$
则有著名的等式
$$
\sigma^a_i\sigma^a_j = u^{(a)}_{ij}\,(i c_i c_j).
$$

这样，任何包含 $\sigma^a_i\sigma^a_j$ 的两体项都可以写成“**Z2 链变量 × 端点 Majorana 二次项**”。

#### 2.2 把 $R_{ij}^{(\text{spin})}$ 改写成 Majorana+Z2 形式

对任意一条最近邻边 $(i,j)\in E_{2D}$，自旋空间上的两体算符是
$$
R_{ij}^{(\text{spin})}=aI+b\,\sigma^x_i\sigma^x_j+c\,\sigma^y_i\sigma^y_j+d\,\sigma^z_i\sigma^z_j.
$$
用上面的映射，一般可以写成
$$
\begin{aligned}
R_{ij}^{(\text{spin})}
&= aI 
  + b\,u^{(x)}_{ij}(i c_i c_j)
  + c\,u^{(y)}_{ij}(i c_i c_j)
  + d\,u^{(z)}_{ij}(i c_i c_j)\\
&= aI + J_{ij}(u)\,(i c_i c_j),
\end{aligned}
$$
其中
$$
J_{ij}(u) := b\,u^{(x)}_{ij} + c\,u^{(y)}_{ij} + d\,u^{(z)}_{ij}
$$
是一个取值为 $\mathbb{R}$ 的函数（在给定 Z2 背景 $u^{(a)}_{ij}=\pm1$ 下就是一个实数系数）。

于是 2D 自旋哈密顿量
$$
H_{2D}^{(\text{spin})} = \sum_{(i,j)\in E_{2D}} R_{ij}^{(\text{spin})}
$$
在 Majorana+Z2 语言下可以写成
$$
H_{2D}[c,u]
= a|E_{2D}|\,I + \frac{i}{4}\sum_{i,j} A_{ij}(u)\,c_i c_j,
$$
其中矩阵元
$$
A_{ij}(u) := 2J_{ij}(u)
$$
在 $(i,j)\in E_{2D}$ 时非零，其他为 0（并取反对称扩展 $A_{ji}=-A_{ij}$）。

这就是标准的“静态 Z2 背景上的自由 Majorana 哈密顿量”。

**与 1D R→$(t,\Delta,\mu)$ 的关系：**

- 在 1D 中，我们已经把最近邻 $R_{i,i+1}^{(\text{spin})}$ 用 JW 化成
	$$
	R_{i,i+1}
	=(b+c)\,(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i)
	 +(b-c)\,(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i)
	 + (\text{密度和常数项}),
	$$
	由此读出
	$$
	t\propto (b+c),\qquad \Delta\propto (b-c),
	$$
	并由 $a,d$ 的线性项匹配化学势 $\mu$。
- 在 2D 中，对于每条边 $(i,j)$，自旋两体项还是那四项，只是边的几何方向不同（方格/蜂窝）。因此**每条边上的二次部分系数和 1D 完全一样**，只是动量空间的色散关系从一维变成
	$$
	\epsilon(\mathbf{k}) \sim \sum_{(i,j)\in E_{2D}} (b+c)\,e^{i\mathbf{k}\cdot(\mathbf{r}_i-\mathbf{r}_j)},\\
	\Delta(\mathbf{k}) \sim \sum_{(i,j)\in E_{2D}} (b-c)\,e^{i\mathbf{k}\cdot(\mathbf{r}_i-\mathbf{r}_j)},
	$$
	从而给出二维拓扑超导的 BdG 形式（见 `kit2.md` 中的蜂窝格例子）。

总之：**从 1D 到 2D，R 与单边上的 $(t,\Delta,\mu)$ 映射规则并没有变，变的是图的几何和动量空间的维数。**

#### 2.3 显式的 2D BdG 形式：用 $a,b,c,d$ 写出 $t,\Delta,\mu$

为了把“R→2D 哈密顿量”的关系写得完全显式，我们回到费米算符形式。先看一条最近邻边 $(i,j)$ 上的两体算符：
$$
R_{ij}^{(\text{spin})}=aI+b\,\sigma^x_i\sigma^x_j+c\,\sigma^y_i\sigma^y_j+d\,\sigma^z_i\sigma^z_j.
$$
用升降算符
$$
\sigma^x=\sigma^+ + \sigma^-,\qquad
\sigma^y=-i\sigma^+ + i\sigma^-,
$$
可得
$$
\begin{aligned}
\sigma^x_i\sigma^x_j
&=\sigma^+_i\sigma^-_j+\sigma^-_i\sigma^+_j
 +\sigma^+_i\sigma^+_j+\sigma^-_i\sigma^-_j,\\
\sigma^y_i\sigma^y_j
&=\sigma^+_i\sigma^-_j+\sigma^-_i\sigma^+_j
 -\sigma^+_i\sigma^+_j-\sigma^-_i\sigma^-_j.
\end{aligned}
$$
于是
$$
\begin{aligned}
b\,\sigma^x_i\sigma^x_j+c\,\sigma^y_i\sigma^y_j
&=(b+c)\,(\sigma^+_i\sigma^-_j+\sigma^-_i\sigma^+_j)\\
&\quad +(b-c)\,(\sigma^+_i\sigma^+_j+\sigma^-_i\sigma^-_j).
\end{aligned}
$$
再用 1D 情形中的 JW 映射（对最近邻边无串）：
$$
\sigma^+_i\sigma^-_j\to c_i^\dagger c_j,\quad
\sigma^-_i\sigma^+_j\to c_j^\dagger c_i,\quad
\sigma^+_i\sigma^+_j\to c_i^\dagger c_j^\dagger,\quad
\sigma^-_i\sigma^-_j\to c_j c_i,
$$
以及
$$
\sigma^z_k=2n_k-1,\qquad
\sigma^z_i\sigma^z_j=(2n_i-1)(2n_j-1)
=4n_in_j-2(n_i+n_j)+1,
$$
得到单边上的费米哈密顿量形式
$$
\boxed{
\begin{aligned}
R_{ij}^{(\text{ferm})}
&= aI + d\,(4n_in_j-2(n_i+n_j)+1)\\
&\quad +(b+c)\,(c_i^\dagger c_j+c_j^\dagger c_i)\\
&\quad +(b-c)\,(c_i^\dagger c_j^\dagger+c_j c_i).
\end{aligned}}
$$

从这条边出发，2D 上的总费米哈密顿量为
$$
H_{2D}^{(\text{ferm})}
=\sum_{(i,j)\in E_{2D}} R_{ij}^{(\text{ferm})}.
$$
于是系数与参数 $a,b,c,d$ 的关系为：

- 最近邻 **hopping**：
	$$t_{ij}=b+c,$$
	沿每一条边 $(i,j)$ 上的 $c_i^\dagger c_j+c_j^\dagger c_i$；
- 最近邻 **配对**：
	$$\Delta_{ij}=b-c,$$
	沿每一条边上的 $c_i^\dagger c_j^\dagger+c_j c_i$；
- **密度-密度相互作用**：
	$$U_{ij}=4d,$$
	对应 $4d\,n_in_j$；
- **化学势/单体密度项**：
	每条边贡献 $-2d(n_i+n_j)$，所以对某个格点 $i$，若配位数为 $z_i$，则总的化学势系数为
	$$
	\mu_i = 2d\,z_i + \mu_0(a),
	$$
	其中 $\mu_0(a)$ 来自整体常数项 $aI+dI$ 的重新整理（通常可以吸收到能量基准中，只留下与 $n_i$ 成正比的部分）。

在 **自由费米子/YBE 自由点** 上，可以选择满足特定约束的 $(a,b,c,d)$ 使得四体项 $4d\,n_in_j$ 在物理子空间中被投影掉或抵消，这样 $H_{2D}^{(\text{ferm})}$ 就退化为标准的二维 Kitaev 型 BdG 哈密顿量：
$$
H_{2D}^{(\text{BdG})}
=\sum_{(i,j)\in E_{2D}}\Big[t_{ij}(c_i^\dagger c_j+c_j^\dagger c_i)
+\Delta_{ij}(c_i^\dagger c_j^\dagger+c_j c_i)\Big]
-\sum_i \mu_i n_i+\text{常数},
$$
其中 $t_{ij}=b+c,\ \Delta_{ij}=b-c,\ \mu_i$ 由 $d$ 和配位数 $z_i$ 决定，这就是 2D 情形下 “R→$(t,\Delta,\mu)$” 的**显式参数对应关系**。

 

### 3. 2D 上的远程算符：从 1D JW 串到 2D 路径字符串

#### 3.1 1D 回顾：远程 R\_{ij} = 端点二次项 × JW 串

在 `kit2-3.md` 中，我们已经严格推导了 1D 上远程 R\_{ij} 的形式：对 $i<j$，
$$
R_{ij}
=C_{ij}I+D_{ij}^{(\text{dens})}
 +Q_{ij}^{(2)}(c)\,e^{i\pi\sum_{k=i}^{j-1}n_k},
$$
其中
$$
Q_{ij}^{(2)}(c)=(b+c)(c_i^{\dagger}c_j+c_j^{\dagger}c_i)+(b-c)(c_i^{\dagger}c_j^{\dagger}+c_j c_i)
$$
是端点二次算符，而
$$
e^{i\pi\sum_{k=i}^{j-1}n_k}
$$
是从 i 到 j 的 JW 串。于是“字符串部分”
$$
R_{ij}^{(\text{string-part})}=Q_{ij}^{(2)}(c)\,e^{i\pi\sum_{k=i}^{j-1}n_k}
$$
和我们后来在 Majorana+Z2 语言中定义的
$$
B_{ij}(\gamma)=u_\gamma(i c_i c_j)
$$
完全平行：它们都表示“端点二次算符 × 路径上的 Z2 串”。

#### 3.2 2D：路径 $\gamma$ 与字符串算符 $B_{ij}(\gamma)$

在 2D Z2+Majorana 框架中，给定图 $\mathcal{G}_{2D}=(V_{2D},E_{2D})$ 和一个从 i 到 j 的路径
$$
\gamma:\ i=i_0 \to i_1 \to \cdots \to i_n=j,\qquad (i_k,i_{k+1})\in E_{2D},
$$
我们定义路径上的 Z2 串
$$
u_\gamma=\prod_{k=0}^{n-1} u_{i_k i_{k+1}},
$$
其中每个 $u_{i_k i_{k+1}}$ 是那条边上的链路变量（x,y,z 三种键中相应的一种）。

然后定义 2D 上的远程 Majorana 键/字符串算符：
$$
B_{ij}(\gamma)=u_\gamma(i c_i c_j).
$$

这就是 1D 中 $Q_{ij}^{(2)}(c)\,e^{i\pi\sum_{k=i}^{j-1}n_k}$ 的二维推广：

- 端点部分：仍然是 Majorana 二次项 $i c_i c_j$，对应 1D 中的 hopping/pairing 组合；
- 路径部分：从一条直线上的 JW 串，变成二维图上的路径 Z2 串 $u_\gamma$；
- 几何信息：不同路径（特别是是否绕过某些 Z2 flux/缺陷）给出不同的 $u_\gamma$，因此在零模子空间中产生不同的相位或矩阵，这就是编织的几何来源。

#### 3.3 远程 R\_{ij} 与 2D 哈密顿量的关系

在 2D 模型中，基本哈密顿量只包含**最近邻边**上的 $R_{ij}^{(\text{spin})}$，而不像 1D 那样考虑全图上的“抽象 R\_{ij}$”。远程算符 $B_{ij}(\gamma)$ 是在给定 2D 模型上构造出来的**附加算符**，用来描述：

- 端点 i,j 之间的长程关联；
- 沿路径 $\gamma$ 搬运零模/任何子时，零模子空间中的有效作用。

在 Majorana+Z2 语言下，可以把 2D 上的 $R_{ij}$ 展开成
$$
R_{ij}=C_{ij}I+D_{ij}^{(\text{dens})}+\lambda_{ij} B_{ij}(\gamma)+\cdots,
$$
其中 $\lambda_{ij}$ 由 $(a,b,c,d)$ 决定，“$\cdots$” 代表其他类型的局域耦合；在拓扑有效理论里，最重要的就是 $\lambda_{ij} B_{ij}(\gamma)$ 这一类“字符串 × 端点 Majorana”结构。

 

### 4. 2D 上的编织：从 $B_{ij}(\gamma)$ 到 $U_\gamma$

上面讨论的是“静态算符” $B_{ij}(\gamma)$。要讨论**编织**，需要看时间演化算符，即沿路径按顺序作用局域 Majorana 交换门。

#### 4.1 局域 Majorana 交换门

在 6.2 节中，我们已经证明：对任意一对 Majorana $\gamma_a,\gamma_b$，
$$
U_{ab}(\theta)=\exp\big(\theta\,\gamma_a\gamma_b\big)
$$
在 $\theta=\pm\tfrac{\pi}{4}$ 时实现一个“$\pm\tfrac{\pi}{2}$ 的旋转”，并给出了辫子关系
$$
U_{12}U_{23}U_{12}=U_{23}U_{12}U_{23},\qquad
U_{12}U_{34}=U_{34}U_{12}\ \text{(不相交时)}.
$$

在 2D Majorana+Z2 模型中，我们可以把**每条最近邻边**上的局域 Hamiltonian
$$
H_{ij}^{\text{(edge)}}(u) = J_{ij}(u)\,(i c_i c_j)
$$
视为生成元，取合适的演化时间 $\tau$，得到边上的局域交换门
$$
U_{ij}^{(\text{braid})}=\exp\Big(-i\tau H_{ij}^{\text{(edge)}}(u)\Big)
\approx \exp\Big(\pm\frac{\pi}{4}\,\gamma_i\gamma_j\Big),
$$
这里“$\approx$” 是在选取能级、约化到零模子空间、以及调参使得 $\tau J_{ij}$ 取到特定值的意义下成立（详情可参考 `R_to_Kitaev.md` 里从 R 构造单边 braid operator 的例子）。

#### 4.2 沿路径的编织算符 $U_\gamma$

给定一条二维路径
$$
\gamma:\ i=i_0\to i_1 \to\cdots\to i_n=j,
$$
在每条边 $(i_k,i_{k+1})$ 上施加一个局域交换门 $U_{i_k i_{k+1}}^{(\text{braid})}$，有序乘积得到
$$
U_\gamma=\prod_{k=0}^{n-1} U_{i_k i_{k+1}}^{(\text{braid})}.
$$

在零模简并子空间中，$U_\gamma$ 的作用正是“把一个任何子/零模沿路径 $\gamma$ 搬运、或绕某个缺陷/环路”的量子操作。由于局域门满足辫子关系：

- 当两条路径 $\gamma,\gamma'$ 及其端点不相交时，对应的 $U_\gamma,U_{\gamma'}$ 对易；
- 当路径只在一个点附近“交叉”或“相邻”时，$U_\gamma,U_{\gamma'}$ 满足辫子群的局域关系（类似 $U_{12},U_{23}$ 之间的关系）。

这与 6.1、6.2 中对 $B_{ij}(\gamma)$ 和 $U_{ab}$ 的代数验证完全一致，只是现在嵌入了一整个 2D 模型的几何背景。

#### 4.4 2D 上 R\_{ij} 编织后的显式形式：$R_{ij}^{(\text{braid})}(\gamma)$

上面 4.2 节是用抽象的 $U_\gamma$ 来描述编织。现在用 1D 中已经得到的 $Q_{ij}^{(2)}(c)$，把“编织后的 R\_{ij}$” 写成**显式含有 $a,b,c,d$ 的形式**。

在 1D 中我们已经证明：对 $i<j$，远程算符可以写成
$$
R_{ij}
=C_{ij}I+D_{ij}^{(\text{dens})}
 +Q_{ij}^{(2)}(c)\,e^{i\pi\sum_{k=i}^{j-1}n_k},
$$
其中
$$
Q_{ij}^{(2)}(c)
=(b+c)(c_i^\dagger c_j+c_j^\dagger c_i)
 +(b-c)(c_i^\dagger c_j^\dagger+c_j c_i).
$$

在 2D 中，给定一条从 i 到 j 的路径
$$
\gamma: i=i_0\to i_1\to\cdots\to i_n=j,
$$
定义路径 Z2 串 $u_\gamma$，远程 Majorana 字符串算符
$$
B_{ij}(\gamma)=u_\gamma(i c_i c_j).
$$
类比 1D 的结构，我们可以把 2D 上**沿路径 $\gamma$ 编织得到的远程 R\_{ij}$** 写成
$$
\boxed{
R_{ij}^{(\text{braid})}(\gamma)
= C_{ij}(a,d)I + D_{ij}^{(\text{dens})}(a,d)
 + Q_{ij}^{(2)}(c; b,c)\,u_\gamma,
}
$$
其中：

- $C_{ij}(a,d)$、$D_{ij}^{(\text{dens})}(a,d)$ 来自 $aI+d\,\sigma^z_i\sigma^z_j$ 展开后得到的常数项和密度/密度-密度项（与 2.3 节中单边的分析完全一致，只是现在 i,j 不要求是最近邻）；
- 
	$$
	Q_{ij}^{(2)}(c; b,c)
	=(b+c)(c_i^\dagger c_j+c_j^\dagger c_i)
	+(b-c)(c_i^\dagger c_j^\dagger+c_j c_i)
	$$
	是**端点二次算符部分**，系数只由 $(b,c)$ 决定，形式与 1D 完全一致；
- $u_\gamma$ 是二维路径上的 Z2 串，编码了编织路径绕过哪些 flux/缺陷。

如果希望把 $R_{ij}^{(\text{braid})}(\gamma)$ 直接写成 2D 自旋语言，可以利用
$$
\sigma^a_i\sigma^a_j = u^{(a)}_{ij}(i c_i c_j)
$$
的逆变换，把 $i c_i c_j$ 还原为“自旋双体 × Z2 链变量”的线性组合，得到类似
$$
R_{ij}^{(\text{braid})}(\gamma)
= a'I + b'\,\sigma^x_i W_x(\gamma)\sigma^x_j
 + c'\,\sigma^y_i W_y(\gamma)\sigma^y_j
 + d'\,\sigma^z_i W_z(\gamma)\sigma^z_j,
$$
其中 $b',c',d'$ 是 $a,b,c,d$ 的线性组合，而 $W_a(\gamma)$ 是沿路径 $\gamma$ 的相应 Wilson 线算符（用链路变量 $u^{(a)}$ 表示）。但在拓扑有效理论中，最核心的结构还是上式中显式写出的那一项：
$$
Q_{ij}^{(2)}(c; b,c)\,u_\gamma,
$$
也就是“由 $(b,c)$ 决定的端点二次算符 × 由路径 $\gamma$ 和 Z2 背景决定的字符串”。这就是 2D 上 R\_{ij} 在编织后的**具体算符形式**，而不仅仅是抽象的 $U_\gamma R_{ij}U_\gamma^{-1}$ 记号。

#### 4.3 与 2D 哈密顿量的对应

需要强调的是：

- 2D 哈密顿量 $H_{2D}[c,u]$ 本身**只包含最近邻边上的局域两体项**；
- 远程字符串算符 $B_{ij}(\gamma)$ 和编织算符 $U_\gamma$ 是在这个局域哈密顿量的基础上**额外构造出来的算符/演化**，用来描述拓扑激发、零模和任何子的长程性质；
- 在适当的低能有效理论和投影到零模子空间后，可以把 $H_{2D}$ 在这个子空间上的作用重写为**若干个 $B_{ij}(\gamma)$ 和 $U_\gamma$ 的组合**，从而建立“R→2D 哈密顿量→2D 编织”的完整链条。

 

### 5. 小结：2D 情况与 1D 的对应关系

- **局域哈密顿量层面**：
	- 1D：$H_{1D}^{(\text{spin})}=\sum_{(i,i+1)}R_{i,i+1}^{(\text{spin})}$；
	- 2D：$H_{2D}^{(\text{spin})}=\sum_{(i,j)\in E_{2D}}R_{ij}^{(\text{spin})}$；
	- Majorana+Z2 语言下，二者都变成“沿边的 Majorana 二次项”的求和，只是边集从一维变成二维。

- **R 与 $(t,\Delta,\mu)$ 的对应关系**：
	- 每条最近邻边上的 t、$\Delta$、$\mu$ 与 1D 完全一样，由 $(a,b,c,d)$ 决定；
	- 差别只在于 2D 的动量空间更复杂，能谱和拓扑不变量（Chern 数等）在 `kit2.md` 中已做分析。

- **远程算符层面**：
	- 1D：远程 R\_{ij} 的字符串部分是 $Q_{ij}^{(2)}(c)\,e^{i\pi\sum_{k=i}^{j-1}n_k}$；
	- 2D：远程字符串算符是 $B_{ij}(\gamma)=u_\gamma(i c_i c_j)$，其中 $u_\gamma$ 编码二维路径和绕行信息。

- **编织层面**：
	- 在 1D 上，路径只有一条直线，谈不上“绕圈”；
	- 在 2D 上，可以沿不同路径 $\gamma$ 进行局域交换，乘积得到编织算符 $U_\gamma$，其代数关系（远的可交换、近邻满足辫子关系）在 6.1,6.2 中已严格验证；
	- $B_{ij}(\gamma)$ 与 $U_\gamma$ 分别是“静态字符串算符”和“时间演化编织算符”两个互补视角。

因此，**2D 情况既保持了 1D 时“R→(t,\Delta,\mu)” 的逐边对应关系，又引入了新的几何自由度（路径/环路/拓扑），使得远程算符和编织真正成为二维拓扑结构的一部分，而不只是“远一点的局域作用”。**

