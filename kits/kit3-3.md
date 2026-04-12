## 3.3 由 R 参数构造的“上层空间”

这一节把前面讨论的想法正式整理成一个“参数空间在下、拓扑数据在上”的结构，用来研究 \((a,b,c,d)\) 的变化如何在更上层的空间里表现出来。

### 3.3.1 参数空间与拓扑相区

把满足 Yang–Baxter 方程的常数 R 写成
$$
R(a,b,c,d)
 \,=\, a\,I + b\,\sigma^x\sigma^x + c\,\sigma^y\sigma^y + d\,\sigma^z\sigma^z,
$$
相应地定义一个参数空间
$$
\mathcal M_R \subset \mathbb R^4,\qquad p=(a,b,c,d)\in\mathcal M_R.
$$

通过 $R(a,b,c,d)$ 与自旋链哈密顿量的对应，再经 Jordan–Wigner 变换，我们得到一族 1D/2D 的 BdG 哈密顿量
$$
p\mapsto H_{\mathrm{BdG}}(p),\qquad p\in\mathcal M_R.
$$

在 2D 情形下（例如方格/蜂窝格），可以在整个参数空间中标出“能隙不闭合”的区域，并用 Chern 数、零模数目以及对应的任意子类型来给出一个粗略的“拓扑相图”。形式上，我们把有隙区域分解为若干连通块
$$
\mathcal M_R^{(\alpha)}\subset\mathcal M_R,\qquad \alpha=\text{不同的拓扑相}.
$$
每个 $\mathcal M_R^{(\alpha)}$ 中，谱始终有隙，Chern 数与任意子范畴（例如 Ising / trivial / \dots）保持不变。

### 3.3.2 固定几何下的“R-上层空间”

选定一个底层几何数据：

- 一个二维流形 $\Sigma$（例如平面、环面等），
- 一组缺陷/涡旋位置 $X=\{x_i\}\subset\Sigma$，
- 一组背景 $\mathbb Z_2$ 规范场构型 $u=\{u_e\}$（在给定拓扑相区内，$u$ 可以视为固定的某一类自旋结构代表）。

对每一个参数点 $p\in\mathcal M_R$，我们借助前面 kit2、kit3-2 中的构造，得到相应的 Majorana+$\mathbb Z_2$ 模型和低能/零模子空间
$$
\mathcal H(p;\Sigma,X,u)\subset\mathcal H_{\mathrm{full}}(p;\Sigma,u).
$$

在格点描述中，每条边 $e=(ij)$ 上有一对 Majorana $\gamma_i,\gamma_j$，以及一个链变量 $u_e=\pm1$。我们定义局域算符

$$
U_e(u_e) \,=\, \exp\Big(\frac{\pi}{4}\,u_e\,\gamma_i\gamma_j\Big),
$$
并对一条路径 $\gamma$ 上的边做有序乘积
$$
U(\gamma,u) \,=\, \overrightarrow{\prod_{e\in\gamma}} U_e(u_e).
$$

把 $U(\gamma,u)$ 投影到拓扑子空间 $\mathcal H(p;\Sigma,X,u)$ 上，得到对路径同伦类 $[\gamma]\in\pi_1(\mathrm{Conf}_X(\Sigma))$ 的一个表示：

$$
\rho_p([\gamma]) \,=\, P_{\mathrm{top}}(p)\,U(\gamma,u)\,P_{\mathrm{top}}(p),
$$

其中 $P_{\mathrm{top}}(p)$ 是投影到零模/拓扑简并子空间的算符。

于是，对固定的底层几何 $(\Sigma,X,u)$，我们可以把“由 R 描述的上层空间”定义为一个随 $p$ 变化的对象族：

$$
\Big\{\mathcal H(p;\Sigma,X,u),\; \rho_p\Big\}_{p\in\mathcal M_R}.
$$

从拓扑场论的视角，这个族可以被理解为：

- 底：参数空间 $\mathcal M_R$；
- 纤维：在给定几何 $(\Sigma,X,u)$ 上的拓扑 Hilbert 空间 $\mathcal H(p;\Sigma,X,u)$；
- 结构：基本群/配置空间 $\pi_1(\mathrm{Conf}_X(\Sigma))$ 或 mapping class 群在该纤维上的表示 $\rho_p$。

也就是说，我们从单个 R 矩阵出发，不仅得到一个哈密顿量，还得到一个“以参数空间为基底、以拓扑 Hilbert 空间及其表示为纤维”的上层结构。

### 3.3.3 参数变化与“上层空间”的变化

现在来看 $(a,b,c,d)$ 的改变在这个结构上的体现。

1. **保持在同一拓扑相区 $\mathcal M_R^{(\alpha)}$ 内**

	只要谱的能隙不闭合，Chern 数和任意子类型（例如 Ising）保持不变。抽象地说，在同一相区内存在一族连续的幺正变换 $W(p\to p')$，把不同参数点的拓扑 Hilbert 空间与表示联系起来：
    $$
    W(p\to p')\,\mathcal H(p;\Sigma,X,u) \cong \mathcal H(p';\Sigma,X,u),
	$$

	且在这些幺正同构下，\(\rho_p\) 与 \(\rho_{p'}\) 是同构的表示。也就是说：

	- F、R、S、T 等范畴数据在整个 $\mathcal M_R^{(\alpha)}$ 内不变；
	- 变的只是“这些数据在具体格点 Majorana+$\mathbb Z_2$ 模型中的实现方式”。

	这种意义下，$\mathcal M_R^{(\alpha)}$ 上配备的是一个“平的局部系统”：沿着任意闭合的参数路径回到同一点，得到的拓扑 Hilbert 空间与表示只差一个整体的幺正基变换（没有新的拓扑相）。

2. **穿过相界：拓扑相的改变**

	当参数 $p$ 穿过能隙闭合的相界时，Chern 数、零模数目以及对应的任意子范畴会发生离散跳变。这时：

	- 纤维 $\mathcal H(p;\Sigma,X,u)$ 的维数和结构可能骤然改变；
	- 表示 $\rho_p$ 的同构类型也改变（例如从 Ising 范畴对应的表示变为 trivial 范畴对应的表示）。

	这可以被看作“上层空间在基底 $\mathcal M_R$ 的某个超曲面上发生了一次拓扑突变”。

3. **沿参数路径追踪上层结构**

	给定一条参数路径
	$$
	p(t)=(a(t),b(t),c(t),d(t)),\qquad t\in[0,1],
	$$
	我们可以在数值或解析上做如下事情：

	- 沿路径构造一族 BdG 哈密顿量 $H_{\mathrm{BdG}}(p(t))$，检查其能隙是否在整个区间内保持开启；
	- 对若干离散的 $t$ 值求出零模空间 $\mathcal H(p(t);\Sigma,X,u)$ 以及关键路径/辫子 $[\gamma]$ 的矩阵表示 $\rho_{p(t)}([\gamma])$；
	- 将这些矩阵与某个标准范畴（例如 Ising）的 F、R 等数据比较，识别在不同区段上的任意子类型是否相同；
	- 一旦发现 Chern 数或表示类型发生跳变，即可判定参数路径穿过了一个拓扑相界，意味着“上层空间的纤维类型”发生了变化。

从这个角度看，$(a,b,c,d)$ 的变化不仅仅是哈密顿量系数的变化，而是驱动了一个“以 R 为基底、以拓扑 Hilbert 空间和表示为纤维的上层空间”的形变与突变。接下来可以在具体例子（例如处于 Ising 区的某条参数路径）上，将上述一般框架具体化为可计算的模型和数值实验。

### 3.3.4 与标准 TQFT/MTC 纤维丛结构的对比

这里把上面定义的“R-上层空间”和常规 TQFT/MTC 中的纤维丛/局部系统结构做一个抽象对比，并精炼若干概念。

**(1) 标准 TQFT/MTC 的纤维丛视角**

在一个 $(2+1)$ 维的共形场论/TQFT（例如 Ising TQFT）中，通常有：

- 对每个带有打洞/缺陷的二维流形 $(\Sigma,X)$，给出一个拓扑 Hilbert 空间
	$$
	\mathcal H_{\mathrm{TQFT}}(\Sigma,X),
	$$
- 对 mapping class 群/编织群元素 $[\gamma]$，给出其在该 Hilbert 空间上的表示
	$$
	\rho_{\mathrm{TQFT}}([\gamma]) : \mathcal H_{\mathrm{TQFT}}(\Sigma,X)\to\mathcal H_{\mathrm{TQFT}}(\Sigma,X).
	$$

把 $(\Sigma,X)$ 看作“底”，$\mathcal H_{\mathrm{TQFT}}(\Sigma,X)$ 看作“纤维”，$\rho_{\mathrm{TQFT}}$ 就是在这些纤维上的一个平坦联络/局部系统。更具体地：

- 当改变缺陷配置 $X$ 的位置（在配置空间 $\mathrm{Conf}_n(\Sigma)$ 内移动）时，$\mathcal H_{\mathrm{TQFT}}$ 沿着配置空间形成一个局部系统；
- 当改变底流形的拓扑类型或做 mapping class 变换时，用 F、R、S、T 等数据在不同切割方式/不同流形之间粘合这些 Hilbert 空间。

在这个框架中，**既没有显式出现微观哈密顿量，也没有显式出现 R 矩阵**；所有结构都被抽象成范畴数据和几何数据。

**(2) “R-上层空间”的角色**

在本节的构造中：

- $\mathcal M_R$ 是 R 矩阵的参数空间（即“微观哈密顿量空间”的一个切片）；
- 对每个 $p\in\mathcal M_R$，我们都有一个“从微观出发构造出的 TQFT 型结构”
	$$
	\big(\mathcal H(p;\Sigma,X,u),\,\rho_p\big),
	$$
	它在适当的拓扑相区内与某个标准 TQFT/MTC（如 Ising）同构；
- “R-上层空间”恰好是在 $\mathcal M_R$ 上，把这些 TQFT 型结构纤维化起来的东西：
	$$
	p\mapsto \big(\mathcal H(p;\Sigma,X,u),\rho_p\big).
	$$

与标准 TQFT 的对比可以粗略地概括为：

- 标准 TQFT：**底是几何/缺陷配置空间**（例如 $(\Sigma,X)$ 的 moduli、$\mathrm{Conf}_n(\Sigma)$），纤维是拓扑 Hilbert 空间，联络由 F、R、S、T 等数据编码；
- 本文的 R-上层空间：**底是 R 的参数空间 $\mathcal M_R$**，纤维是“在某个固定几何 $(\Sigma,X,u)$ 上、由该 R 诱导出来的 TQFT 型结构”。

换句话说，我们在“几何/缺陷的上层”再加了一层：

- 最底层：$(\Sigma,X,u)$ 及其配置空间/基本群；
- 中间层：在固定 $p$ 下，由 Majorana+$\mathbb Z_2$ 模型给出的 TQFT 型结构 $(\mathcal H(p;\Sigma,X,u),\rho_p)$；
- 最上层：随 $p\in\mathcal M_R$ 变化的这整个 TQFT 型结构族。

**(3) 精炼若干定义与同构关系**

在拓扑相区 $\mathcal M_R^{(\alpha)}$ 内，如果我们知道这一相区对应的抽象 TQFT/MTC（例如 Ising），那么存在一族（可选到局部连续的）幺正变换
$$
W_{\alpha}(p): \mathcal H(p;\Sigma,X,u) \to \mathcal H_{\mathrm{TQFT}}^{(\alpha)}(\Sigma,X),\qquad p\in\mathcal M_R^{(\alpha)},
$$
使得在该相区内，对所有相关的 $[\gamma]$ 都有
$$
W_{\alpha}(p)\,\rho_p([\gamma])\,W_{\alpha}(p)^{-1}
 \,=\, \rho_{\mathrm{TQFT}}^{(\alpha)}([\gamma]).
$$

这里：

- $\mathcal H_{\mathrm{TQFT}}^{(\alpha)}(\Sigma,X)$ 是一旦给定相型 $\alpha$ 后就固定的拓扑 Hilbert 空间（与 $p$ 无关）；
- $\rho_{\mathrm{TQFT}}^{(\alpha)}$ 是该相型下的“标准” TQFT 表示（由 F、R、S、T 完全确定）。

这样，“R-上层空间”可以被理解为：

1. 在每个相区 $\mathcal M_R^{(\alpha)}$ 内，与一个固定的 TQFT/MTC 数据同构，只是坐标系随 $p$ 改变；
2. 在相界上，这种同构类型发生变化，对应于“更换了一个 TQFT/MTC”。

这就把“用 R 描述的上层空间”与常规定义的 TQFT/MTC 纤维丛结构放在了同一张图上：TQFT 负责描述“在给定相型里，上层空间的纤维是什么样”；而 $\mathcal M_R$ 和 R-上层空间，则负责描述“在微观参数变化时，我们在不同 TQFT 纤维之间如何移动和跳变”。

### 3.3.5 Ising 例子的三层结构示意图

下面用一个简单的示意图，把“R 参数空间 → 由 R 诱导的 TQFT 型结构 → 抽象 Ising TQFT”这三层关系画在一起（底层的几何/缺陷配置默认为固定的 $(\Sigma,X,u)$）。

$$
\begin{array}{ccc}
	ext{(顶层：微观参数)} & & \\
 p \in \mathcal M_R^{(\text{Ising})}
 & \xrightarrow{\;\text{R + Majorana+}\mathbb Z_2\;}&
 \big(\mathcal H(p;\Sigma,X,u),\,\rho_p\big)
 \\
 & & \quad\;\downarrow W_{\alpha}(p) \\
 & & \\
 & & \big(\mathcal H_{\mathrm{Ising}}(\Sigma,X),\,\rho_{\mathrm{Ising}}\big)
 \\
 & & \text{(中层：抽象 Ising TQFT 纤维)}
\\[6pt]
 & & \\
 & & \rho_{\mathrm{Ising}}:\;\pi_1(\mathrm{Conf}_X(\Sigma))\;\text{或}\;\mathrm{MCG}(\Sigma)\to \mathrm{U}(\mathcal H_{\mathrm{Ising}}) \\
 & & \text{(底层：几何/缺陷的配置空间与基本群)}
\end{array}
$$

其中：

- 顶层的 $p\in\mathcal M_R^{(\text{Ising})}$ 表示处在 Ising 拓扑相区内的一点（即某组 $(a,b,c,d)$）；
- 通过格点上的 R + Majorana+$\mathbb Z_2$ 构造，得到中间的“由 R 诱导的 TQFT 型结构” $(\mathcal H(p;\Sigma,X,u),\rho_p)$；
- 在该相区内，存在幺正基变换 $W_{\alpha}(p)$，把这一结构与标准的 Ising TQFT 纤维 $(\mathcal H_{\mathrm{Ising}}(\Sigma,X),\rho_{\mathrm{Ising}})$ 同构；
- 最底层是几何/缺陷配置 $(\Sigma,X)$ 及其配置空间/基本群（或 mapping class 群），Ising TQFT 给出了从这些群到 $\mathrm{U}(\mathcal H_{\mathrm{Ising}})$ 的表示 $\rho_{\mathrm{Ising}}$。

这样，在 Ising 具体例子中，“向下看”可以从 $p$ 沿着箭头，先到由 R 诱导的 Majorana+规范模型，再到抽象的 Ising TQFT；“向上看”则是从一个给定的 TQFT（Ising）出发，寻找哪些 R 参数区和微观模型可以实现它，以及在这些参数区之间如何通过形变或相变进行连接。

为了更直观地体现参数路径的作用，可以在顶层的 $\mathcal M_R^{(\text{Ising})}$ 里加入一条曲线
$$
\gamma_R:[0,1]\to \mathcal M_R^{(\text{Ising})},\qquad t\mapsto p(t)=(a(t),b(t),c(t),d(t)),
$$
例如保持 $a,d$ 不变，只在 $(b,c)$ 平面的一条线段上移动。对应地，这条路径在三层结构中的像可以示意为：

$$
\begin{array}{ccc}
\gamma_R(t) & & \\
 p(0) \longrightarrow p(t) \longrightarrow p(1)
 & \xrightarrow{\;\text{R + Majorana+}\mathbb Z_2\;}&
 \big(\mathcal H(p(0)),\rho_{p(0)}\big)
 \longrightarrow \big(\mathcal H(p(t)),\rho_{p(t)}\big)
 \longrightarrow \big(\mathcal H(p(1)),\rho_{p(1)}\big)
 \\
 & & \quad\;\downarrow W_{\alpha}(p(0)),\;W_{\alpha}(p(t)),\;W_{\alpha}(p(1)) \\
 & & \\
 & & \big(\mathcal H_{\mathrm{Ising}},\,\rho_{\mathrm{Ising}}\big)
\end{array}
$$

在整个路径 $\gamma_R$ 保持在 Ising 相区内且能隙不闭合的前提下：

- 底层的 Ising TQFT 纤维 $\mathcal H_{\mathrm{Ising}}$ 和表示 $\rho_{\mathrm{Ising}}$ 完全不变；
- 中间层的 $(\mathcal H(p(t)),\rho_{p(t)})$ 随 $t$ 连续变化，但在幺正变换 $W_{\alpha}(p(t))$ 下始终与同一个 $(\mathcal H_{\mathrm{Ising}},\rho_{\mathrm{Ising}})$ 同构；
- 顶层的 $p(t)$ 则记录了“同一种 Ising TQFT 在不同 R 参数下的具体实现方式”的变化。

如果把 $\gamma_R$ 延伸到穿出 Ising 区、跨越一个相界，那么上图的下方就会在某个 $t_c$ 处发生“纤维类型”的突变：$\mathcal H_{\mathrm{Ising}}$ 被另一种 TQFT 的 Hilbert 空间替换，$\rho_{\mathrm{Ising}}$ 被对应的新范畴表示替换，表示了一个真正的拓扑相变。




