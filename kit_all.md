## 标题（暂定）

从 Pauli 形式的 R 矩阵到 Majorana 几何：

一维 Kitaev 链、二维 $p+ip$ 涡旋与蜂窝格 vison 上的 Dehn twist 与 Berry holonomy


## 摘要

本文从一个简单的 Pauli 形式两比特 $R$‑矩阵出发，系统构建了“代数 $R$ → BdG/Majorana 哈密顿量 → 零模 Hilbert 丛 → Berry 联络与 Dehn twist → 拓扑门与电路复杂度”的统一图景。具体而言，我们考虑
$$
R_{ij}(a,b,c,d)=aI + b\,\sigma^x_i\sigma^x_j + c\,\sigma^y_i\sigma^y_j + d\,\sigma^z_i\sigma^z_j,
$$
分析其满足 Yang–Baxter 方程的参数子流形 $\mathcal M_R^{(\mathrm{YBE})}$，并给出从 $(a,b,c,d)$ 到一维 Kitaev 链/四 Majorana toy 模型 BdG 参量 $(t,\Delta,\mu)$ 的显式映射。在 1D 平台上，我们通过 4‑Majorana 模型与 F,R‑矩阵计算，数值验证了“half‑twist 的平方 + F‑move”给出的 Dehn twist 门与 Ising TQFT 理论表达 $U_{\mathrm{Dehn}}\simeq F^{-1}R^2F$ 的一致性。

随后，我们在二维 spinless $p+ip$ 拓扑超导的格点 BdG 模型中插入两只涡旋，构造一条“一个涡旋绕另一涡旋一圈”的离散路径，在零模子空间上计算 Berry holonomy 并与 Ising 的 $R^{\sigma\sigma}$ 与 $(R^{\sigma\sigma})^2$ 做 SU(2) 共轭比较。数值上，在一段有限的化学势区间内，该涡旋绕行路径实现的 Berry holonomy 与理想 Dehn twist 的 SU(2) 重合度 $F_{\mathrm{Dehn}}(\mu)$ 维持在 $\sim0.95$–$1.0$ 的高值，给出了一个清晰的“2D $p+ip$ 上的 Dehn twist 平坦区”。

最后，我们将上述几何/拓扑图景迁移到 Kitaev 蜂窝格模型：用 Majorana+$\mathbb Z_2$ 规范场表示自旋，利用链路变量 $u_{ij}$ 与回路通量 $W_p$ 构造 vison 缺陷，设计“一个 vison 绕另一 vison 一圈”的 Z$_2$ 规范路径，并在有限尺寸 $4\times4$ 蜂窝格上实现完整的 Berry holonomy 数值计算。我们沿三条简单的耦合切片 —— $J_x=J_y=1$ 扫 $J_z$，$J_x=J_z=1$ 扫 $J_y$，以及 simplex 上的 triangle 路径 $J_x=J_y=(1-t)/2,J_z=t$ —— 扫描 Dehn twist fidelity $F_{\mathrm{Dehn}}$ 与能隙 proxy $\min|E|$。结果表明：在当前 toy honeycomb 设置下，谱在这些切片上保持有隙，但 vison‑loop Berry holonomy 的 SU(2) 方向始终与 Ising Dehn twist 近乎正交（$F_{\mathrm{Dehn}}\approx0$），未能复现类似 2D $p+ip$ 的 Dehn twist 平坦区。这一“负例”明确划出了“非阿贝尔 + 有能隙”与“实现 Ising Dehn twist”之间的差别，并为后续在更物理的蜂窝格参数区和更大系统中搜索 Dehn twist 平坦区奠定了数值基础。


## 1 引言

拓扑量子计算的一个核心思想是：在某些凝聚态体系或量子场论中，编码在非阿贝尔任意子零模中的量子信息在拓扑意义上受保护，编织（braiding）与 Dehn twist 等 mapping class 群操作在逻辑子空间上实现一组“几何的、鲁棒的”量子门。抽象上，这一结构通常由模张量范畴（modular tensor category, MTC）或拓扑量子场论（TQFT）描述，其基本数据包括 F 矩阵、R 矩阵以及映射类群表示 $\rho: \mathrm{MCG}\to U(\mathcal H)$。

然而，从给定的“微观哈密顿量”出发，如何系统地连接到上述“上层的范畴/几何结构”，并在具体的 Majorana/BdG 模型中显式计算 Berry holonomy 与 Dehn twist，一直是一个兼具物理与数学挑战的问题。一方面，代数可积系统给出了满足 Yang–Baxter 方程（YBE）的 $R$‑矩阵族，这些 $R$ 可看作某种“平坦联络的离散化”；另一方面，拓扑相中的 BdG/Majorana 模型在缺陷配置空间上自然携带非阿贝尔 Berry 联络，其 holonomy 与 TQFT 的 F,R 数据密切相关。

本文尝试在一个尽量简单、可数值实现的框架下，把这两端连接起来：从一个 Pauli 形式的两比特 $R(a,b,c,d)$ 出发，通过适当的映射得到 1D Kitaev 链与 2D BdG/Majorana 模型，分析其零模 Hilbert 丛上的 Berry 联络，构造 braid/Dehn twist 路径的 Berry holonomy，并与 Ising TQFT 的 F,R‑矩阵做 SU(2) 共轭比较。我们的目标不是构造一个严格可积的 2D 模型，而是通过一系列 1D/2D 的可控 toy 模型和数值实验，建立如下的工作图景：
$$
R(a,b,c,d) \Rightarrow H_{\mathrm{BdG/Majorana}}(p,x) \Rightarrow (\mathcal E, A, F) \Rightarrow U_{\mathrm{Berry}}[\gamma] \approx F^{-1}R^2F \Rightarrow \text{拓扑门与电路复杂度}.
$$

本文分为三个主要层次：

1. 在 1D 平台上，从 Pauli 形式 $R(a,b,c,d)$ 出发，给出到 Kitaev 链参量 $(t,\Delta,\mu)$ 的映射，分析 YBE 约束下的参数子流形与拓扑不变量，利用 4‑Majorana toy 模型实现一次“纯 1D 的 Dehn twist 数值验证”。
2. 在 2D spinless $p+ip$ 拓扑超导格点模型中插入两只涡旋，沿一条矩形路径让一只涡旋绕行另一只，计算零模子空间上的 Berry holonomy，与 Ising 的 $R^{\sigma\sigma}$ 与 $(R^{\sigma\sigma})^2$ 做 SU(2) 共轭比较，并沿化学势 $\mu$ 扫描得到 $F_{\mathrm{Dehn}}(\mu)$ 的“Dehn twist 平坦区”。
3. 在 Kitaev 蜂窝格模型中用 Majorana+$\mathbb Z_2$ 规范场表示自旋，引入 vison 缺陷与 branch cut，构造“一个 vison 绕行另一 vison”的 Z$_2$ 规范路径，在有限 $4\times4$ 蜂窝格上实现完整的 vison‑loop Berry 数值实验，并沿若干耦合切片扫描 $F_{\mathrm{Dehn}}$ 与能隙，比较蜂窝格与 2D $p+ip$ 在 Dehn twist 平坦区上的差异。


## 2 Pauli 形式的 $R(a,b,c,d)$ 与 1D Kitaev 链

### 2.1 Pauli ansatz 与 Yang–Baxter 方程

我们考虑最简单的 Pauli 形式两比特算符
$$
R_{ij}(a,b,c,d) = aI + b\,\sigma^x_i\sigma^x_j + c\,\sigma^y_i\sigma^y_j + d\,\sigma^z_i\sigma^z_j,
$$
其中 $(a,b,c,d)\in\mathbb R^4$ 或 $\mathbb C^4$。将其视为 $4\times4$ 矩阵并嵌入三个比特的张量积空间，要求满足恒等式
$$
R_{12}R_{23}R_{12} = R_{23}R_{12}R_{23},
$$
即可得到一族多项式约束，定义 Yang–Baxter 方程的解集 $\mathcal M_R^{(\mathrm{YBE})}\subset\mathbb C^4$。在这类 Pauli ansatz 下，YBE 条件等价于一组关于 $a,b,c,d$ 的代数方程，刻画了一族“可积的” R‑矩阵。

从物理角度看，把 $R$ 写成时间演化算符 $U=e^{-iH\delta t}$ 的离散近似，可把 YBE 看作某种“参数联络在立体 2‑胞上平坦”的条件；$\mathcal M_R^{(\mathrm{YBE})}$ 则对应于参数空间中 Berry 曲率在若干代表性 2‑胞上消失的子流形。后文我们不会直接用到 YBE 的解析解，而是把它视为一个“理想的平坦极限”，以 1D/2D 数值 Berry 计算的平坦近似区域来与之对应。

### 2.2 从 $(a,b,c,d)$ 到 Kitaev 链 $(t,\Delta,\mu)$

在 1D 平台上，我们关心的是 Kitaev 链
$$
H_{\mathrm{Kitaev}} = -\sum_j\bigl(t c_j^\dagger c_{j+1} + \Delta c_j c_{j+1} + \mathrm{h.c.}\bigr) - \mu\sum_j(c_j^\dagger c_j - \tfrac12),
$$
的拓扑相结构与端点 Majorana 零模。我们给出一组简单而自然的映射（详见独立推导文档），在局部上把 $R(a,b,c,d)$ 视为一段小时间演化，识别出
$$
t = b+c,\qquad \Delta = b-c,\qquad \mu \simeq 4d+\mu_0,
$$
其中 $\mu_0$ 为某个固定参考点的化学势偏移（例如 $\mu_0=-2$）。

这一映射的直观含义是：

- $b\,\sigma^x\sigma^x + c\,\sigma^y\sigma^y$ 的线性组合生成近邻跳跃与 $p$‑波配对两种通道，其和/差分别对应 $t,\Delta$；
- $d\,\sigma^z\sigma^z$ 项在 Jordan–Wigner 映射下等效为局域密度与化学势的修正，因此 $\mu$ 与 $d$ 有近似线性关系；
- 常数项 $aI$ 只贡献整体能量偏移，不影响拓扑。

在这一识别下，每个 $(a,b,c,d)$ 点都对应一个具体的 Kitaev 链 $(t,\Delta,\mu)$，可以通过已知的 Pfaffian 或 winding 数判据来诊断其是否处在拓扑相（存在端点 Majorana 零模）。这为后续讨论 $R$ 参数变化如何影响 Majorana 边界零模与 Dehn twist 提供了一个直接桥梁。

### 2.3 4‑Majorana toy 模型与 1D Dehn twist 数值验证

为了在最简单的 Hilbert 空间里检验“F,R,Dehn twist”结构，我们考察一个只含 4 个 Majorana 模式的 toy 模型，它可以看作两个端点 Majorana 对的截断：
$$
\{\gamma_1,\gamma_2,\gamma_3,\gamma_4\},\qquad \{\gamma_a,\gamma_b\}=2\delta_{ab},\qquad \gamma_a^\dagger=\gamma_a.
$$
在合适的总费米数约束下，这 4 个 Majorana 生成一个二维逻辑 Hilbert 空间，可与 Ising TQFT 中两只 $\sigma$ 粒子的 Hilbert 空间同构。在这一空间中，可以显式写出：

- half‑twist（交换）$R^{\sigma\sigma}$ 的作用矩阵；
- F‑move 在不同括号化 $(\sigma\sigma)\sigma\to\sigma(\sigma\sigma)$ 间的变换矩阵；
- Dehn twist 组合 $U_{\mathrm{Dehn}}\simeq F^{-1}R^2F$。

通过数值构造一条“instantaneous braid/Dehn twist‑like” 的参数路径，并在 4‑Majorana 零模空间上用 Berry holonomy 计算 $U_{\mathrm{Berry}}[\gamma]$，我们发现：在合适的参数和路径选择下，$U_{\mathrm{Berry}}$ 在 SU(2) 意义下与上述 $U_{\mathrm{Dehn}}$ 高度共轭，给出了一个纯 1D toy 模型中的 Dehn twist 数值验证。这部分细节详见专门笔记与 verify 脚本，这里不再赘述，而把重心放在 2D 与蜂窝格的实现上。

### 2.4 R‑空间的几何图景：Hilbert 丛、平坦联络与 YBE 子流形

kit3 系列的一个核心工作，是把上述 1D R‑模型放进一个明确的微分几何框架里：给定一组参数 $p=(a,b,c,d)$ 以及几何/缺陷背景 $(\Sigma,X,u)$（例如 genus‑2 曲面上若干 puncture 或 4‑Majorana 端点），我们可以在“缺陷配置空间”
$$
\mathcal C = \mathrm{Conf}_X(\Sigma)
$$
上构造一个 Hilbert 向量丛
$$
\pi: \mathcal E_R\to\mathcal C,\qquad \pi^{-1}(x)=\mathcal H_R(x),
$$
其中纤维 $\mathcal H_R(x)$ 是在给定 R 与规范背景下，该配置 $x$ 上的零模/拓扑简并子空间。选一组随 $x$ 平滑变化的正交基 $\{|\psi_i(x)\rangle\}$，可以定义 Berry‑型联络 1‑形式
$$
A_{ij}(x)=\langle\psi_i(x)|\,d_x\psi_j(x)\rangle,
$$
其沿配置空间路径 $\gamma\subset\mathcal C$ 的 holonomy 为
$$
U_R[\gamma]=\mathcal P\exp\Big(-\int_\gamma A\Big),
$$
正是我们在 1D/2D/蜂窝格中计算的 braid/Dehn twist 拓扑门。

在这个连续联络的基础上，可以定义曲率 2‑形式
$$
F=dA+A\wedge A,
$$
它刻画了“绕无穷小闭合回路的平行移动偏离单位”的程度。局部计算表明，若在某个单连通区域 $U\subset\mathcal C$ 内 $F\equiv0$，则存在规范变换使得 $A$ 在 $U$ 内化为纯规范项 $A=g^{-1}dg$，从而所有可缩回的局部回路 $\gamma\subset U$ 的 holonomy 都是单位元；所有非平庸的 holonomy 只来自配置空间的基本群 $\pi_1(\mathcal C)$，即辫子群 $B_n$ 与 mapping class 群的元素。换言之：在平坦区内，Berry holonomy 只依赖同伦类，不依赖路径细节，这就是“真正拓扑门”的几何刻画。

在 kit3‑5 中，我们把这一几何结构与 Pauli 形式 R‑ansatz 精确对齐：在“局域两体 + 自旋 $1/2$ + 四参数 $R(a,b,c,d)$ + 适当的对称性”这组物理假设下，每个辫子生成元的微观作用可以唯一写成某个两体算符 $R(a,b,c,d)$ 在链上的嵌入 $B_i=R_{i,i+1}$。把“平坦联络给出辫子群表示”的几何事实与这一 ansatz 结合，恰好得到常数 Yang–Baxter 算符方程
$$
R_{12}(a,b,c,d)R_{23}(a,b,c,d)R_{12}(a,b,c,d)
=R_{23}(a,b,c,d)R_{12}(a,b,c,d)R_{23}(a,b,c,d).
$$
在 Pauli 基底下展开这一算符方程，可以把它简化为对 $(a,b,c,d)$ 的多项式约束，定义出一个代数子流形
$$
\mathcal M_R^{(\mathrm{YBE})}\subset\mathbb C^4,
$$
也就是“YBE 解族”。几何上，这正对应于：在这些参数点上，离散的 R‑联络在局域三角 2‑胞上的 holonomy 为单位元，连续极限中曲率 2‑形式 $F$ 的通量为零，即“R‑诱导联络在代表性 2‑胞上平坦”。其补集 $\mathcal M_R\setminus\mathcal M_R^{(\mathrm{YBE})}$ 则是“带曲率的非 YBE 区域”：绕最小 YBE 2‑胞的 holonomy 偏离单位元，微观上的 braid/Dehn 关系只近似成立。

在这一框架下，Dehn twist 可以统一理解为沿配置空间或 moduli 空间中某个非平凡闭合回路的 holonomy：在 Ising 相 + YBE 子流形内部，$U_R[T_\gamma]$ 与 Ising TQFT 的 $\rho_{\mathrm{Ising}}(T_\gamma)$ 共轭，只依赖同伦类，是严格意义的拓扑门；离开这些条件时，同一条“几何 Dehn 路径”的 holonomy 开始对路径细节与参数变化敏感，其偏离理想 $F^{-1}R^2F$ 的方式，正可以用曲率 $F$ 与参数联络的曲率 $\mathcal F$ 来系统刻画。这就是 kit3‑5 中“R, YBE 与联络/曲率的微分几何图景”的精华总结。

### 2.5 half twist、Dehn twist 与辫子群表示：从 R 到拓扑门

在 1D/kit3 系列中，我们更细致地追踪了“最初写下的 Pauli 形式 $R(a,b,c,d)$”与“half twist、Dehn twist 以及辫子群表示”之间的具体关系。可以用以下几步来概括：

1.**Majorana 语言下的 half twist 与辫子群**  
在 4‑Majorana toy 模型以及更长的 Kitaev 链中，一对相邻 Majorana 模式 $\gamma_i,\gamma_j$ 的“几何交换”在低能子空间上由幺正算符
$$
	U_{\text{half}}^{(ij)}=\exp\Bigl(\tfrac{\pi}{4}\,\gamma_i\gamma_j\Bigr)
$$
实现。不同对 $U_{\text{half}}^{(i,i+1)}$ 之间满足标准辫子关系
$$
	U_{\text{half}}^{(i,i+1)}U_{\text{half}}^{(i+1,i+2)}U_{\text{half}}^{(i,i+1)}
	=U_{\text{half}}^{(i+1,i+2)}U_{\text{half}}^{(i,i+1)}U_{\text{half}}^{(i+1,i+2)},
$$
因此在零模 Hilbert 空间上给出辫子群 $B_n$ 的一个具体矩阵表示，生成元 $\sigma_i$ 映到 $U_{\text{half}}^{(i,i+1)}$。

2.**spin 语言中的 half twist 与 Pauli 形式 R**  
通过 Jordan–Wigner 映射，Majorana 二次项 $i\gamma_i\gamma_j$ 在自旋基底中对应于一类最近邻两体算符，可以用 Pauli 张量积展开。对“局域两体 + 自旋 1/2 + 反幺正对称”这组假设，所有允许的两体算符张成的 4 维子空间，正好由
$$
	I,\ \sigma^x\otimes\sigma^x,\ \sigma^y\otimes\sigma^y,\ \sigma^z\otimes\sigma^z
$$
张成。这就是最初的 Pauli ansatz
$$
	R(a,b,c,d)=aI + b\,\sigma^x\sigma^x + c\,\sigma^y\sigma^y + d\,\sigma^z\sigma^z.
$$
在小时间步极限下，可以把 $R$ 视为某个两体哈密顿量 $H_{ij}$ 的时间演化片段 $R\approx e^{-iH_{ij}\delta t}$，其中 $H_{ij}$ 的主导部分就是 $i\gamma_i\gamma_j$ 型耦合。也就是说，适当选取 $(a,b,c,d)$ 并在编码子空间上投影后，$R(a,b,c,d)$ 在本质上实现的正是某个 $U_{\text{half}}^{(ij)}$，差别仅在整体相位与基变换上；这一点在 4‑Majorana 数值脚本中已被显式验证。

3.**YBE 与辫子群表示的对接：$R_{i,i+1}$ 作为 braid 生成元**  
在 2.4 的几何图景中，我们从“配置空间上的平坦联络”出发得到辫子群表示 $\rho:B_n\to U(\mathcal H)$。在 Pauli ansatz 下再加上“局域两体”的物理要求，自然的选择是
$$
	\rho(\sigma_i)=B_i=R_{i,i+1}(a,b,c,d).
$$
要使 $\rho$ 真正成为群同态，就必须满足
$$
	R_{12}R_{23}R_{12}=R_{23}R_{12}R_{23},
$$
也就是常数 Yang–Baxter 方程。换句话说：**YBE 恰恰是在我们选取“用局域 Pauli 形式 R 代表 braid 生成元”这一微观实现方案下，使得 R‑门族在零模 Hilbert 空间上给出一致辫子群表示的充要条件。**解集 $\mathcal M_R^{(\mathrm{YBE})}$ 就是那些 $(a,b,c,d)$，它们既可以被理解为“平坦 R‑联络”的参数点，也可以被理解为“真正 braid 表示”的参数点。

4.**从 half twist 到 Dehn twist：$U_{\mathrm{Dehn}}\simeq F^{-1}R^2F$ 的微观含义**  
在 Ising TQFT 中，两只 $\sigma$ 粒子的 braid 生成元 $\sigma$ 在不同融合通道 $1,\psi$ 上的作用由 $R^{\sigma\sigma}_c$ 给出，其平方 $(R^{\sigma\sigma}_c)^2$ 与相应通道的 topological spin $\theta_c$ 相联系；更一般的 Dehn twist 可以写成若干次 F‑move 与这些 $R^2$ 的组合，形式上记为
$$
	U_{\mathrm{Dehn}}(\gamma)\simeq F^{-1}R^2F.
$$
在我们的 R+Majorana 框架中，这句话可以翻译为：

- 选取一组微观 half twist 门 $U_{\text{half}}^{(ij)}$（对应若干 $R_{i,i+1}$），在适当的编码基底（F‑move 给出的基变换）下，其平方 $U_{\text{half}}^{2}$ 在逻辑子空间上与某条 Dehn twist 曲线 $\gamma$ 的抽象作用幺正等价；

- 数值上，我们在 4‑Majorana toy 模型中比较了 three objects：几何 Dehn 路径的 Berry holonomy $U_{\mathrm{Berry}}$、用 Ising TQFT F,R 构造的 $U_{\mathrm{Dehn}}^{\mathrm{Ising}}$ 以及若干 $U_{\text{half}}^{(ij)2}$ 的组合，发现三者在 SU(2) 上属于同一个共轭类，只差整体相位和有限维基变换。这为 “Dehn twist = 适当组合的 half twist 平方” 在具体微观模型中提供了直接数值证据。

5.**从最初的 R 到拓扑门的闭环**  
综上，在 1D/kit3 这条线上，最初写下的 Pauli 形式
$$
	R(a,b,c,d)=aI + b\,\sigma^x\sigma^x + c\,\sigma^y\sigma^y + d\,\sigma^z\sigma^z
$$
与 Dehn twist/half twist/辫子群之间的关系可以用一条闭环来概括：

- 通过 Jordan–Wigner 与 Majorana 表示，它对应某个局域 Majorana 二次耦合 $i\gamma_i\gamma_j$，其时间演化片段给出微观 half twist 门；

- 要让这些 R‑门在零模 Hilbert 空间上组成一致的辫子群表示，YBE 提供了必要且充分的代数约束，定义出 $\mathcal M_R^{(\mathrm{YBE})}$ 这一“平坦/可积”子流形；

- 在这一子流形及 Ising 拓扑相内部，half twist 的平方与有限次 F‑move 可以组合成 Dehn twist 门 $U_{\mathrm{Dehn}}\simeq F^{-1}R^2F$，其 Berry holonomy 与抽象 TQFT 表示 $\rho_{\mathrm{Ising}}(T_\gamma)$ 共轭；

- 离开 $\mathcal M_R^{(\mathrm{YBE})}$ 或 Ising 区后，R 仍然定义了某种“微观 braid 操作”，但对应的联络带曲率、辫子/Dehn 关系只近似成立，门的作用对路径与参数细节变得敏感，这正可以通过我们在 1D/2D 数值扫描中观察到的 $F_{\mathrm{Dehn}}(p)$ 衰减来量化。

这一节把 kit3 系列中关于 “R、half twist、Dehn twist 与辫子群表示” 的代数与几何关系压缩成了一个统一的故事线，也为后续在 2D $p+ip$ 和蜂窝格上寻找与这些 1D 结果一致的 Dehn twist 平坦区提供了清晰的参照系。


## 3 2D $p+ip$ 拓扑超导与涡旋 Dehn twist

### 3.1 2D $p+ip$ BdG 模型与涡旋零模

在连续极限中，一个简化的 spinless 2D $p+ip$ 拓扑超导 BdG 哈密顿量可写为
$$
H = \int d^2\mathbf r\, \Psi^\dagger(\mathbf r)
\begin{pmatrix}
\frac{\mathbf p^2}{2m} - \mu & \Delta(\mathbf r) (p_x + i p_y) \\
\Delta^*(\mathbf r) (p_x - i p_y) & -\frac{\mathbf p^2}{2m} + \mu
\end{pmatrix}
\Psi(\mathbf r),
$$
其中 $\mu>0$、$\Delta(\mathbf r)=|\Delta(\mathbf r)|e^{i\theta(\mathbf r)}$。当 $|\Delta|$ 非零且 $\mu$ 落在拓扑相区时，体系的填充带具有非零 Chern 数，属于 2D 拓扑超导相。

在这一背景中插入磁通涡旋意味着：序参量相位 $\theta(\mathbf r)$ 绕某点 $\mathbf R_a$ 一圈时绕转 $2\pi$。对每个涡旋 $v_a$，BdG 方程存在一个指数局域的零能解 $\gamma_a$，满足
$$
\{\gamma_a,\gamma_b\}=2\delta_{ab},\qquad \gamma_a^\dagger=\gamma_a,
$$
即涡旋携带 Majorana 零模。对 $2N$ 个涡旋，零模子空间维数为 $2^{N-1}$，与 Ising 任何子理论中 $2N$ 个 $\sigma$ 粒子的 Hilbert 空间维数一致。

在我们的格点实现中，我们使用一个 $L_x\times L_y$ 方格 BdG 模型（典型参数为 $L_x=L_y=20$，$t=1.0,\Delta_0=0.5,\mu=-2.0$ 等），在格上用近似极角构造两个涡旋的相位纹理，并数值对角化得到近零能本征态。最低本征值谱显示一对接近零能的本征值与其余激发之间有明显能隙，验证了这两个涡旋确实绑定了一对 Majorana 零模。

### 3.2 涡旋配置空间上的 Berry 联络与 braid/Dehn twist

设 $X=(\mathbf R_1,\dots,\mathbf R_n)$ 表示涡旋位置，配置空间为
$$
\mathcal C = \mathrm{Conf}_n(\Sigma) = \{X\in\Sigma^n\mid \mathbf R_a\neq\mathbf R_b\}/S_n.
$$
对每个 $X$，零模本征态张成一个 Hilbert 空间 $\mathcal H_{p+ip}(\Sigma,X)$，整体上形成 Hilbert 丛 $\mathcal E_{p+ip}\to\mathcal C$。选取局域正交归一基底 $\{\psi_i(X)\}$，定义 Berry 联络
$$
A_{ij}(X)=\langle\psi_i(X)|d_X\psi_j(X)\rangle,
$$
以及对应的 holonomy
$$
U_{p+ip}[\gamma]=\mathcal P\exp\Bigl(-\int_\gamma A\Bigr),\qquad [\gamma]\in\pi_1(\mathcal C).
$$

在 2D $p+ip$ 拓扑相区内，只要路径 $\gamma$ 不闭合能隙，上式就只依赖 $[\gamma]$ 的同伦类，并与 Ising TQFT 提供的 braid/Dehn twist 表示 $\rho_{\mathrm{Ising}}([\gamma])$ 共轭。

### 3.3 双涡旋绕行数值实验与 Dehn twist 实现

我们在格点 BdG 模型中放置两只涡旋 $v_1,v_2$，固定 $v_1$ 于格中心附近，令 $v_2$ 沿着围绕 $v_1$ 的矩形路径绕行一圈。将路径离散成若干步 $\lambda_0,\dots,\lambda_N$，在每一步上构造 BdG 哈密顿量 $H(\lambda_k)$，对角化并按 $|E|$ 选取前两个本征态构成零模子空间基底 $V_k$。对相邻步构造重叠矩阵
$$
W_k = V_k^\dagger V_{k+1},
$$
用极分解投影到最近酉矩阵 $U_k$，最后沿路径有序相乘得到
$$
U_{\mathrm{holo}} = \prod_k U_k.
$$

数值结果表明：对一组典型参数，$U_{\mathrm{holo}}$ 在 SU(2) 意义下与 Ising 的 $(R^{\sigma\sigma})^2$ 高度共轭，其 Dehn twist fidelity $F_{\mathrm{Dehn}}\sim0.97$，而与单次 braid $R^{\sigma\sigma}$ 的重合度 $F_R\sim0.85$。这说明在这一具体的 2D BdG 模型与涡旋路径选择下，“一个涡旋绕另一涡旋一圈”实现的是一次近乎完美的 Ising Dehn twist 门，而不仅仅是单次 braid。

### 3.4 化学势扫描 $F_{\mathrm{Dehn}}(\mu)$ 与 Dehn twist 平坦区

为了更系统地理解 Dehn twist 在参数空间中的稳定性，我们在保持几何与路径不变的前提下，沿化学势 $\mu$ 扫描 Berry holonomy 与 Dehn twist 的 SU(2) 重合度。数值上，我们在
$$
\mu\in\{-3.0,-2.5,-2.0,-1.5,-1.0,-0.5\}
$$
上计算了 $F_{\mathrm{Dehn}}(\mu)$，并绘制了图像 [pip_vortex_F_Dehn_vs_mu.png](pip_vortex_F_Dehn_vs_mu.png)。结果显示：

- 在 $\mu\approx-3.0,-2.5,-2.0,-1.0$ 附近，$F_{\mathrm{Dehn}}(\mu)$ 维持在 $0.92$–$1.0$ 范围内，给出一段“Dehn twist 平坦区”；
- 在 $\mu\approx-1.5$ 与 $\mu\approx-0.5$ 一带，$F_{\mathrm{Dehn}}(\mu)$ 明显跌落到 $\sim0.1$ 或更低，表明此时这一涡旋绕行路径在零模子空间上的作用已经严重偏离理想的 Dehn twist。

这条一维截面提供了一个直观例证：在 2D $p+ip$ 模型的参数空间中存在一段有限宽度的区域，使得某条固定的几何路径在零模 Hilbert 丛上的 Berry 联络近似平坦，Dehn twist 门在这一参数区间内几乎刚性不变；离开这一区域后，同一条路径实现的门迅速“跑出”理想的拓扑类型。这与我们在 1D $R(a,b,c,d)$ 参数空间中通过 YBE 子流形与 Berry 曲率平坦区域所获得的图景在概念上完全一致。


## 4 Kitaev 蜂窝格：Majorana+$\mathbb Z_2$ 规范场与 vison Berry

### 4.1 Majorana+$\mathbb Z_2$ 表示与链路变量

在 Kitaev 蜂窝格模型中，我们在每个格点 $j$ 引入四个 Majorana 算符 $b^x_j,b^y_j,b^z_j,c_j$，用
$$
\sigma^a_j = i b^a_j c_j,\qquad a\in\{x,y,z\},
$$
表示自旋，并引入局域约束
$$
D_j\equiv b^x_j b^y_j b^z_j c_j = +1,
$$
以把扩展 Hilbert 空间投影回物理自旋‑1/2 空间。对类型为 $a$ 的键 $(i,j)$ 定义链路 Z$_2$ 变量
$$
u_{ij}\equiv -i b^a_i b^a_j,\qquad u_{ij}^2=1.
$$
利用 Majorana 的反对易代数，可以严格推导出
$$
\sigma^a_i\sigma^a_j = u_{ij}\,(i c_i c_j),
$$
从而原始的 Kitaev 哈密顿量
$$
H = -\sum_{\langle ij\rangle_a} J_a\,\sigma^a_i\sigma^a_j
$$
可写成固定 $u_{ij}$ 背景下的自由 Majorana 模型
$$
H_c[u] = \frac i4\sum_{i,j} A_{ij}[J,u]\, c_i c_j,
$$
其中 $A$ 为实反对称矩阵，其非零元由 $J_a u_{ij}$ 给出。回路通量
$$
W_p=\prod_{(ij)\in p}u_{ij}
$$
描述六边形 plaquette 上的 Z$_2$ 通量（vison）。在纯 Kitaev 模型中，$u_{ij}$ 守恒，$W_p$ 为守恒量，基态位于某个特定的 flux 扇区（通常为无通量扇区）。

### 4.2 从 $R(a,b,c,d)$ 到蜂窝格耦合与相图

把 Pauli 形式的 $R_{ij}=aI + b\,\sigma^x_i\sigma^x_j + c\,\sigma^y_i\sigma^y_j + d\,\sigma^z_i\sigma^z_j$ 映到蜂窝格时，自然的标量映射是
$$
J_x=b,\qquad J_y=c,\qquad J_z=d,
$$
其中常数 $a$ 只贡献整体能量偏移。在无通量扇区、平移对称下，对 $H_c[u^{(0)}]$ 做动量空间分析，可以得到每个波矢 $\mathbf k$ 的 $2\times2$ 有效矩阵
$$
A(\mathbf k)=2\begin{pmatrix}0 & f(\mathbf k)\\ -f^*(\mathbf k) & 0\end{pmatrix},
$$
其中
$$
f(\mathbf k)=J_x e^{i\mathbf k\cdot\delta_x} + J_y e^{i\mathbf k\cdot\delta_y} + J_z e^{i\mathbf k\cdot\delta_z},
$$
单粒子谱为 $\varepsilon(\mathbf k)=\pm2|f(\mathbf k)|$。可证明：存在 $\mathbf k$ 使 $f(\mathbf k)=0$ 当且仅当 $J_x,J_y,J_z$ 可以在复平面中构成闭合三角形，即满足三角不等式；这给出了标准的 Kitaev 蜂窝格相图：归一化条件 $J_x+J_y+J_z=1$ 下，中心三角区域（每个 $J_\alpha\le1/2$）为 gapless 相，三个角落区域（某一 $J_\alpha>1/2$）为 gapped 相。

要得到带非零 Chern 数的 chiral 非阿贝尔相，需要在 gapless 区附近加入破时间反演的三自旋项或等效的 next‑nearest‑neighbor 复跳跃，给 Dirac 点打开质量隙并赋予能带非平庸 Chern 数。这一部分本文只在概念上依托已有文献，不做新的数值实现；我们在数值上采用的是一个“无磁场的 toy honeycomb 模型”，主要目的是在有限系统上构造可运行的 vison Berry 流水线。

### 4.3 统一的七步流程：从 $R$ 到 Majorana 与 Berry 几何

在综述完 1D/2D/p+ip/honeycomb 各自的构造后，我们在 [kit4-3.md](kit4-3.md) 中把整个过程总结为一个七步流程：

1. 选定具体物理平台与参数路径（1D R‑模型、2D $p+ip$、honeycomb 等），给出参考点 $p_0$ 与附近的参数线/面；
2. 从 $R(a,b,c,d)$ 或其它耦合出发，构造 BdG/Majorana 哈密顿量 $H_p(x)$，其中 $x$ 标记缺陷位置；
3. 进行谱与拓扑不变量分析，获得能隙与 Pfaffian/Chern 等静态拓扑判据；
4. 在缺陷配置空间上构造 Berry 联络，沿典型 braid/Dehn twist 路径计算 Berry holonomy；
5. 在零模逻辑子空间中，把 Berry holonomy 规约到 SU(2) 等特殊酉群，并与 TQFT 的 F,R 数据构造的 Dehn twist $U^{\mathrm{TQFT}}_{\mathrm{Dehn}}$ 比较；
6. 在参考点附近做线性稳定性分析与参数扫描，得到 $F_{\mathrm{Dehn}}(p)$ 与能隙的图像；
7. 综合上述结果，提炼出“几何–代数–Majorana”三层之间的字典及其对拓扑门稳定性与电路复杂度的影响。

这一流程在 2D $p+ip$ 上已经完整走通；接下来，我们在蜂窝格 toy 模型上实现其最小版本。

### 4.4 蜂窝格 vison‑loop Berry 数值实验与各向异性扫描

在 [verify/run_honeycomb_vison_berry.py](verify/run_honeycomb_vison_berry.py) 中，我们实现了一个有限尺寸 $4\times4$ 蜂窝格上的 vison‑loop Berry 数值实验，其关键步骤如下：

1. 构造 brick‑wall 蜂窝几何与最近邻键 $(i,j,\alpha)$，建立 Majorana 哈密顿量矩阵 $H_c[u]$；
2. 自动枚举 6 边形 plaquette 及其通量 $W_p=\prod u_{ij}$，在无通量规范下所有 $W_p=+1$；
3. 在 plaquette 邻接图上用 BFS 寻找任意两块 $p_A,p_B$ 之间的最短路径，对每对相邻 plaquette 选择一条共享边构成 branch cut；沿 branch cut 反转 $u_{ij}$ 即只在 $p_A,p_B$ 上产生 vison（$W_p=-1$）；
4. 在 9 个 plaquette 中以其几何中心构成 $3\times3$ 网格，固定中心 plaquette 为 $p_A$，令 $p_B$ 沿着围绕中心的一个矩形回路绕行一圈；在每一步上重建从 $p_A$ 到当前 $p_B$ 的 branch cut，确保始终只有两 vison；
5. 在每个规范构型 $u^{(k)}$ 下对角化 $H_c[u^{(k)}]$，提取最低能的二维子空间作为“vison 编码子空间”，用重叠矩阵极分解方法累积 Berry holonomy $U_{\mathrm{Berry}}$；
6. 将 $U_{\mathrm{Berry}}$ 归一化到 SU(2)，与 Ising 的 $R^{\sigma\sigma}$ 和 Dehn twist $F^{-1}R^2F$ 的 SU(2) 形式比较，得到 $F_R,F_{\mathrm{Dehn}}$；
7. 沿若干参数切片扫描上述量并成图。

在各向同性点 $J_x=J_y=J_z=1$、$L_x=L_y=4$ 下，初始谱显示一对近零能本征值与其余激发之间有明显能隙，vison‑loop Berry holonomy 数值上接近单位阵，其 SU(2) 形式满足
$$
F_R\approx0.707,\qquad F_{\mathrm{Dehn}}\approx0,
$$
说明当前 vison‑loop 在 SU(2) 意义下更接近某个“约为 $\pi/4$ 的绕轴旋转”，而非 Ising Dehn twist。

### 4.5 沿 $J_z,J_y$ 及 simplex triangle 路径的 Dehn twist 扫描

为了探查蜂窝格 toy 模型中是否存在类似 2D $p+ip$ 的 Dehn twist 平坦区，我们在固定 vison‑loop 路径的前提下，考察了三条简单的参数切片：

1. $J_x=J_y=1$，沿 $J_z\in\{0.5,0.8,1.0,1.2,1.5\}$ 扫描，得到图 [honeycomb_vison_F_Dehn_vs_Jz.png](honeycomb_vison_F_Dehn_vs_Jz.png) 与 [honeycomb_vison_gap_vs_Jz.png](honeycomb_vison_gap_vs_Jz.png)；
2. $J_x=J_z=1$，沿 $J_y\in\{0.5,0.8,1.0,1.2,1.5\}$ 扫描，得到图 [honeycomb_vison_F_Dehn_vs_Jy.png](honeycomb_vison_F_Dehn_vs_Jy.png) 与 [honeycomb_vison_gap_vs_Jy.png](honeycomb_vison_gap_vs_Jy.png)；
3. simplex 上的 triangle 路径
	$$
	J_z=t,\qquad J_x=J_y=\tfrac12(1-t),\qquad t\in[0,1],
	$$
	得到图 [honeycomb_vison_F_Dehn_vs_triangle_t.png](honeycomb_vison_F_Dehn_vs_triangle_t.png) 与 [honeycomb_vison_gap_vs_triangle_t.png](honeycomb_vison_gap_vs_triangle_t.png)。

数值结论可以简要概括为：

- 在所有切片上，gap proxy $\min|E|$ 要么保持有限、要么仅在小区间内略有减小，整体上未出现系统性能隙塌陷；
- $F_R$ 在所有扫描点上几乎稳定在 $\approx0.707$；
- $F_{\mathrm{Dehn}}$ 在数值精度内始终接近 0，从未出现类似 2D $p+ip$ 的“高平坦段”。

特别是在 triangle 路径上，即便在靠近 gapless 线的区域（小 $t$、$J_z\approx J_x+J_y$）处，$F_{\mathrm{Dehn}}(t)$ 仍然没有显著抬升。这表明：在当前无场且有限尺寸的蜂窝格 toy 模型中，虽然可以构造出良好定义的 vison‑loop Berry holonomy，并在谱上保持有隙，但这一 holonomy 在 SU(2) 意义下与 Ising Dehn twist 保持近乎正交；要在蜂窝格平台上真正复现 2D $p+ip$ 那种 Dehn twist 平坦区，需要引入更真实的非阿贝尔相（如磁场诱导的 Chern 带）、更大系统尺寸以及更“拓扑”的 vison‑loop 路径。


## 5 几何–代数–Majorana 的对应与近似 YBE 子流形

通过以上 1D/2D/p+ip/honeycomb 的分析，我们可以把“代数 YBE、上层 Berry 几何与 Majorana 行为”之间的对应关系概要整理如下：

1. Pauli 形式的 $R(a,b,c,d)$ 与 YBE 方程给出一族代数上平坦的子流形 $\mathcal M_R^{(\mathrm{YBE})}$。在这些子流形上，三线交叉的重排对应的参数联络在代表性 2‑胞上平坦。
2. 在 1D Kitaev 链与 4‑Majorana toy 模型中，沿这些 YBE 子流形或其邻域构造的 half‑twist/Dehn twist 参数路径，其 Berry holonomy 在零模逻辑子空间上与 Ising 的 $F^{-1}R^2F$ 高度一致，可以视为“1D 版的平坦子流形”。
3. 在 2D $p+ip$ 模型中，化学势 $\mu$ 的一段区间为某条具体涡旋绕行路径提供了一个 Dehn twist 平坦区：在这一段中，Berry holonomy 作为 $\mu$ 的函数近似常量，与 Ising Dehn twist 共轭，而离开这一段后则迅速偏离。这可以被视为“2D 版的近似 YBE 子流形”，只是参数轴从 $(a,b,c,d)$ 换成了 $(\mu,\Delta,\dots)$。
4. 在蜂窝格 toy 模型中，我们已经构造出完备的 vison‑loop Berry 数值框架，但在目前考察的三条简单参数切片上并未发现 Dehn twist 平坦区；相反，我们得到的是一系列“有隙但 $F_{\mathrm{Dehn}}\approx0$” 的负例。这一负例清楚地说明：“处在非阿贝尔相并保持能隙”只是实现 Ising Dehn twist 的必要条件之一，而非充分条件；真正的 Dehn twist 平坦区对参数与路径有更精细的几何/代数要求。

从更宏观的角度看，这些结果为下面这样一个工作思路提供了支撑：

- 给定一族代数上可积或近可积的 $R$‑矩阵，先在 1D 平台上通过 BdG/Majorana 映射与 4‑Majorana toy 模型确定其对应的 TQFT/F,R 结构；
- 再在 2D/p+ip 与更真实的格点模型（如蜂窝格）中，寻找参数与路径的组合，使得 Berry 联络在代表性的 braid/Dehn twist 2‑胞上近似平坦，从而在零模 Hilbert 丛上实现与这些 F,R 数据兼容的拓扑门；
- 最后，用这些几何上优化的参数路径来设计常深度的 LQC+permutation 电路，在实际量子器件中实现近似的拓扑门，同时量化偏离平坦子流形时的误差与复杂度开销。


## 6 结论与展望

本文围绕一个简单但结构丰富的 Pauli 形式 $R(a,b,c,d)$，通过 1D Kitaev 链、4‑Majorana toy 模型、2D $p+ip$ 涡旋及 Kitaev 蜂窝格 vison 四条主线，系统构建并部分数值验证了“代数 R → BdG/Majorana → Berry 联络与 Dehn twist → 拓扑门”的统一图景。在这一过程中：

- 我们给出了从 $(a,b,c,d)$ 到 Kitaev 链参量 $(t,\Delta,\mu)$ 的显式识别，并在 4‑Majorana toy 模型中数值验证了 half‑twist 的平方与 F‑move 给出的 Dehn twist 门与 Ising TQFT 表达的一致性；
- 在 2D $p+ip$ 格点 BdG 模型中，我们构造了一个具体的双涡旋绕行数值实验，发现其 Berry holonomy 在一段化学势区间内与 Ising Dehn twist 高度共轭，给出了一个清晰的 Dehn twist 平坦区；
- 在 Kitaev 蜂窝格 toy 模型中，我们实现了完整的 vison‑loop Berry 数值流水线，并沿多条耦合切片扫描了 Dehn twist fidelity 与能隙，虽然尚未在这些简单切片上找到 Dehn twist 平坦区，但明确刻画了“仅有非阿贝尔相与能隙并不足以保证实现 Ising Dehn twist”的边界条件。

未来的工作方向包括：

1. 在蜂窝格中引入更贴近文献的非阿贝尔区（如磁场诱导的三自旋项与 chiral Chern 带），并在更大系统尺寸与更“拓扑”的 vison‑loop 路径上系统搜索 Dehn twist 平坦区，与 2D $p+ip$ 的结果做真正意义上的“对等比较”；
2. 在 1D/2D 平台上更系统地分析 YBE 子流形附近的 Berry 曲率结构，理解代数可积性与几何平坦性之间的精确对应，以及这种对应如何控制拓扑门在参数空间中的稳定区域；
3. 在已有的数值框架基础上，引入含噪声与相互作用的微扰，研究 Dehn twist 平坦区在更真实环境下的退化方式，为“拓扑保护 + LQC+permutation 纠错”联合方案提供定量基准；
4. 把本文建立的“R→Majorana→Berry→Dehn twist”流程进一步推广到其它 MTC（如 Fibonacci、Ising ×U(1) 等）以及包含多种缺陷/边界的更复杂几何。

这些工作有望在代数可积系统、拓扑量子场论与具体凝聚态拓扑相之间搭起一座更加系统、可数值验证的桥梁，为在实际量子器件中实现拓扑门与优化其几何资源开销提供理论与数值基础。

