### 从 YBE anyon 编织到几何扭转：R = exp(iH\_P) 视角下的空间诱导拓扑

这一节把前面几份笔记（kit-exp, kit-new 等）中分散的内容，围绕一个清晰的问题重新组织：

> 如何把原本在 YBE / TQFT 中“抽象定义”的 anyon 编织（half twist, Dehn twist），通过 $R=e^{iH_P}$ 的指数表示和微扰，联系到具体晶格模型中**配置空间的几何扭转**，从而得到一种“由空间重构而非时间演化”产生的拓扑效应？

思路分三层：

1. 抽象层：YBE 与 anyon 编织；
2. 实现层：$R=e^{iH_P}$ 与微观哈密顿量；
3. 几何层：配置空间、Berry 曲率与“空间诱导”的拓扑。

下面依次梳理，并在每一层给出尽可能自洽的论证。

 

#### 1. 抽象层：从 YBE 到 anyon 编织

在抽象的可积模型 / TQFT 语言中，一个满足常数 YBE 的 R‑矩阵
$$
R: V\otimes V\to V\otimes V,
\qquad
R_{12}R_{13}R_{23}=R_{23}R_{13}R_{12},
$$
自然定义了 braid group 的一个表示：
$$
B_i := R_{i,i+1}: V^{\otimes n}\to V^{\otimes n},
$$
满足
$$
B_iB_{i+1}B_i = B_{i+1}B_iB_{i+1},\qquad B_iB_j=B_jB_i~~(|i-j|\ge 2).
$$
H=\sum_j H_{j,j+1}
在一个给定的 anyon 理论（如 Ising TQFT）中，选取某个类型 $\sigma$ 的 anyon，其两个 $\sigma$ 的融合空间
$$
V_{\sigma\sigma} = \mathrm{Hom}(\sigma\otimes\sigma,\cdot)
$$
是一个有限维向量空间。TQFT 给出 F、R 符号：

- R‑符号 $R^{\sigma\sigma}: V_{\sigma\sigma}\to V_{\sigma\sigma}$ 表示两个 $\sigma$ 交换（half twist）的幺正作用；
- F‑符号 $F^{\sigma\sigma\sigma}$ 则描述重新关联（不同的括号方式）时基的变换。

在这样的抽象理论中：

- half twist 对应单次交换 $\sigma\leftrightarrow\sigma$，由某个 R‑矩阵表示；
- Dehn twist 属于更大群（mapping class group）的元素，可以写成若干 F、R 的 word，例如在 Ising 理论中有
	$$
	U_{\mathrm{Dehn}}\sim F^{-1}(R^{\sigma\sigma})^2F,
	$$
	在适当的融合空间中给出一个 SU(2) 元素（剥去整体相位后）。

**这一层要点**：YBE 和 R‑矩阵起点是完全“抽象”的：它只告诉我们**某个 Hilbert 空间上的幺正表示**，并未涉及具体的空间几何或时间演化。我们要做的是把这个抽象的 R 与具体的晶格哈密顿量、缺陷配置空间联系起来。

 

#### 2. 实现层：R = exp(iH\_P) 与微观哈密顿量

我们在 kit-new 中采用了统一的指数表示：
$$
R = e^{iH_P},\qquad H_P\in\mathfrak{su}(4)
$$
（对两比特空间而言），其中
$$
H_P = \sum_{i,j\in\{0,x,y,z\}} c_{ij}\,\sigma_i\otimes\sigma_j,
$$
并把它放到具体的晶格模型上。

**2.1 H\_P 的自然分解**

对每一条键 $\langle j,j+1\rangle$ 上的 H\_P，在 Jordan–Wigner / Majorana 映射下，哈密顿量可自然拆分为
$$
H = H_{\mathrm{quad}} + H_{\mathrm{int}} + H_{\mathrm{gauge}},
$$
其中：

1. $H_{\mathrm{quad}}$：所有二次 Majorana 项 $i\gamma_a\gamma_b$ 的和，对应 1D/2D BdG / Kitaev 有效模型（给出 $t,\Delta,\mu,\dots$ 等参数）；
2. $H_{\mathrm{int}}$：四 Majorana 及以上的多体相互作用项，打破自由性和许多可积结构；
3. $H_{\mathrm{gauge}}$：只改变整体能量、全局 U(1) 相位或规范选择的项，对 Berry 曲率与拓扑数据本质上是“规范”或“平移”。

下面把这一“自然分解”在 1D JW 语言和 2D Majorana+Z2 语言中都推一遍，说明它并不是拍脑袋分三堆，而是按**费米子算符的次数和物理效果**自动出现的分解。

**(A) 1D：从 H\_P 到 JW 费米子，再到三部分**

从两比特生成元
$$
H_P = \sum_{\mu,\nu\in\{0,x,y,z\}} c_{\mu\nu}\,\sigma^\mu_j\otimes\sigma^\nu_{j+1}
$$
出发。按算符在 JW 之后的“费米子次数”来整理：

1. **$\sigma^x_j\sigma^x_{j+1},\;\sigma^y_j\sigma^y_{j+1}$ 型：产生 hopping/配对（二次项）**

在 [R_to_Kitaev.md](R_to_Kitaev.md) 已经推过
$$
\begin{aligned}
\sigma^x_j\sigma^x_{j+1}+\sigma^y_j\sigma^y_{j+1}&=2(\sigma^+_j\sigma^-_{j+1}+\sigma^-_j\sigma^+_{j+1}),\\
\sigma^x_j\sigma^x_{j+1}-\sigma^y_j\sigma^y_{j+1}&=2(\sigma^+_j\sigma^+_{j+1}+\sigma^-_j\sigma^-_{j+1}),
\end{aligned}
$$
而 JW 映射给出
$$
\sigma^+_j\sigma^-_{j+1}+\sigma^-_j\sigma^+_{j+1} \mapsto c_j^\dagger c_{j+1}+c_{j+1}^\dagger c_j,\\
\sigma^+_j\sigma^+_{j+1}+\sigma^-_j\sigma^-_{j+1} \mapsto c_j^\dagger c_{j+1}^\dagger+c_{j+1}c_j.
$$
因此，所有来自 $c_{xx}\sigma^x_j\sigma^x_{j+1}$ 与 $c_{yy}\sigma^y_j\sigma^y_{j+1}$ 的贡献，在 JW 之后都是**严格二次的费米子算符**（加上可能的常数因子）：
$$
\sim (c_j^\dagger c_{j+1}+h.c.) + (c_j^\dagger c_{j+1}^\dagger+h.c.),
$$
它们自然并入 $H_{\mathrm{quad}}$，决定 $t,\Delta$ 等参数。

2. **$\sigma^z_j\sigma^z_{j+1}$ 型：四费米子 + 线性密度 + 常数**

同样在 [R_to_Kitaev.md](R_to_Kitaev.md) 计算过
$$
\sigma^z_j\sigma^z_{j+1}=(2n_j-1)(2n_{j+1}-1)
 = 4n_jn_{j+1}-2(n_j+n_{j+1})+1.
$$
这里自然出现三类项：

- **四费米子项** $4n_jn_{j+1}$：这是严格的相互作用（四个 $c,c^\dagger$），属于 $H_{\mathrm{int}}$；
- **线性密度项** $-2(n_j+n_{j+1})$：在链上求和后给出化学势形如 $-\mu\sum_j n_j$，对拓扑分类来说只是在 $(t,\Delta)$ 空间上平移 $\mu$，可视为 $H_{\mathrm{quad}}$ 与 $H_{\mathrm{gauge}}$ 之间的“边界”，通常把其视作 $H_{\mathrm{gauge}}$ 一部分（调参用）；
- **常数项** $+1$：只给整体能量位移，完全属于 $H_{\mathrm{gauge}}$。

因此，所有 $c_{zz}\,\sigma^z_j\sigma^z_{j+1}$ 贡献可以**自动分裂**为：
$$
H_{\mathrm{int}}\supset 4c_{zz}\sum_j n_jn_{j+1},\qquad H_{\mathrm{gauge}}\supset -2c_{zz}\sum_j(n_j+n_{j+1})+\text{const}.
$$

3. **含 $I$ 或单体 $\sigma^z$ 的项：纯 on-site/gauge 调整**

H\_P 中若有
$$
c_{00}\,I\otimes I + c_{z0}\,\sigma^z_j\otimes I + c_{0z}\,I\otimes\sigma^z_{j+1},
$$
在 JW 下
$$
\sigma^z_j\mapsto 2n_j-1,
$$
因此这些都只给 on-site 密度项和整体常数：
$$
\sim \sum_j (\mu_j n_j) + \text{const}.
$$
它们**不产生新的 hopping/配对，也不产生四费米子相互作用**，因此完全可以收进 $H_{\mathrm{gauge}}$ 里——作为“重新定义化学势和零点能量”的自由度。

4. **含 $\sigma^x,\sigma^y$ 的单体：通常在对称/规范选择下被禁止或吸收**

例如 $c_{x0}\,\sigma^x_j\otimes I$ 在 JW 下会生成带长串的非局域项，一般在具有粒子数守恒、宇称对称或自旋旋转对称的物理模型中被禁止；在我们关心的 Kitaev/BdG 类模型里，通常直接令这些系数为零或通过局域基变换消去。因此它们要么不存在，要么被视作规范选择的一部分，同样并入 $H_{\mathrm{gauge}}$。

综合以上四类，可以看到：

- **所有“产生/湮灭一对邻近费米子”的项**自动落入 $H_{\mathrm{quad}}$（$t,\Delta$）；
- **所有“密度–密度 $(nn)$ 或更高阶”的项**自动落入 $H_{\mathrm{int}}$；
- **所有“只含 $I$ 或线性密度 $n_j$ 的项”**只改整体能量/化学势/规范选择，落入 $H_{\mathrm{gauge}}$。

这就是 1D 情形下“自然三分”的代数来源。

为了把上面的讨论写成一条**整链的显式 JW 后哈密顿量**，可以选一个简单但代表性的情形：只保留
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
H = &\sum_{j=1}^{L-1}\Big[ (J_x+J_y)\big(c_j^\dagger c_{j+1}+c_{j+1}^\dagger c_j\big)
 \
&\qquad\qquad + (J_x-J_y)\big(c_j^\dagger c_{j+1}^\dagger+c_{j+1}c_j\big) \\
&\qquad\qquad + 4J_z\,n_jn_{j+1} + (-2J_z+2h)\,(n_j+n_{j+1}) \\
&\qquad\qquad + (J_z-2h+\varepsilon)\Big].
\end{aligned}
$$
把它按前面的三部分重写，就是
$$
\begin{aligned}
H_{\mathrm{quad}} &= \sum_{j=1}^{L-1}\Big[ t\,(c_j^\dagger c_{j+1}+c_{j+1}^\dagger c_j)
 + \Delta\,(c_j^\dagger c_{j+1}^\dagger+c_{j+1}c_j)\Big] - \mu\sum_j n_j,\\
H_{\mathrm{int}} &= U\sum_{j=1}^{L-1} n_jn_{j+1},\\
H_{\mathrm{gauge}} &= \text{常数项（整体能量零点），以及边界上与 $n_j$ 相关的修正},
\end{aligned}
$$
其中参数 $t,\Delta,U,\mu$ 都是 $J_x,J_y,J_z,h,\varepsilon$ 的线性组合；例如在上面的具体例子中可以取
$$
t=J_x+J_y,\qquad \Delta=J_x-J_y,\qquad U=4J_z,
$$
而 $\mu$ 则由 $J_z,h$ 的线性组合给出（精确表达式需要把 $\sum_j(n_j+n_{j+1})$ 在整条链上展开并单独处理端点，但对体相和拓扑结构来说只是一种化学势的重新定义）。

约定（重要）：在本文档中我们按**算符次数（operator degree）**来分类哈密顿量的三部分：

$H_{\mathrm{quad}}$ 包括所有**二次费米子项**（即所有形如 $c^\dagger c$, $c^\dagger c^\dagger$ 及其厄米共轭的项），因此化学势项 $\mu\sum_j n_j$ 属于 H_{\mathrm{quad}}；这些项决定能带、谱隙和零模式结构。

$H_{\mathrm{int}}$ 包括所有**四次及以上**的多体作用项（如 $n_jn_{j+1}$ 等）。

$H_{\mathrm{gauge}}$ 包括纯常数项、纯局域常数能量位移或仅调整规范/基准的项。


**(B) 2D：在 Majorana+Z2 语言下的相同三分**

在 [kit2.md](kit2.md) 和 [kit2-exp.md](kit2-exp.md) 中，我们用四 Majorana + Z2 规范场映射了 2D 的 $
\sigma^a_i\sigma^a_j$：
$$
\sigma^a_i\sigma^a_j = u_{ij}^{(a)}\,(ic_ic_j),\qquad u_{ij}^{(a)}=-ib^a_ib^a_j,\quad a\in\{x,y,z\}.
$$
对一个纯 Kitaev‑型键哈密顿项
$$
K^{(a)}_{ij}=J_a\,\sigma^a_i\sigma^a_j
$$
有
$$
K^{(a)}_{ij}\mapsto J_a u_{ij}^{(a)}\,(ic_ic_j).
$$
因此：

- 所有 $J_a$ 控制的最近邻 Majorana 跳跃 $ic_ic_j$ 自然构成 $H_{\mathrm{quad}}$；
- 若在 H\_P 中加入 $\sigma^z_i\sigma^z_j$ 型的额外相互作用（超出纯 Kitaev 点），在 Majorana 语言下会产生 $(n_in_j)$ 或更高阶算符，自动落入 $H_{\mathrm{int}}$；
- 所有 on-site 项（如 $\sigma^z_i$、常数）只影响化学势或整体能量，归入 $H_{\mathrm{gauge}}$。

在 2D 上额外出现的是 Z2 链路变量 $u_{ij}$：

- **固定 $u_{ij}$ 背景**时，$H_{\mathrm{quad}}$ 仍然是自由 Majorana；不同 $\{u\}$ 扇区对应不同拓扑背景（flux 配置）；
- 若把 $u_{ij}$ 动力学化（加入翻转 $u_{ij}$ 的项），这些项本身通常也是四 Majorana 甚至更高阶的算符，例如 $u_{ij}ic_ic_j$、$W_p=\prod u_{ij}$，因此自然并入 $H_{\mathrm{int}}$ 或 $H_{\mathrm{gauge}}$（取决于是否只改变能量谱而不改零模结构）。

从这个角度看，**不论 1D 还是 2D，只要我们以“费米子次数 + 对拓扑数据的效应”作为分类标准，哈密顿量都会自动分裂成三类项**：

1. $H_{\mathrm{quad}}$：二次 Majorana，决定带结构、零模、bulk gap，与拓扑相结构和 braid/Dehn holonomy 直接相关；
2. $H_{\mathrm{int}}$：四次及以上，体现“相互作用/非自由”，在弱耦合下可视为对拓扑相的微扰，在强耦合时可能驱动相变；
3. $H_{\mathrm{gauge}}$：只改能量基准、化学势或规范选择，对拓扑不变量（Chern 数、topological spin 等）在相同相区内不产生本质影响。

这就是我们在 2.1 节宣称的“自然分三部分”的详细代数推导和物理含义。

在 4‑Majorana / 2‑site 极限点，我们选取一个特殊的
$$
H_P^{(0)}\quad\text{使得}\quad H_{\mathrm{quad}}^{(0)} = H_0 = \frac{i}{2}\gamma_2\gamma_3,\qquad H_{\mathrm{int}}^{(0)}\approx 0.
$$
此时时间演化算符
$$
U_0(\tau) = e^{-iH_0\tau} = \exp\Bigl(\frac{\pi}{4}\gamma_2\gamma_3\Bigr)
$$
恰好是 half twist 生成元。若我们选取 R 为
$$
R^{(0)} = e^{iH_P^{(0)}},
$$
则在适当的时间/尺度下 $R^{(0)}$ 在 2‑$\sigma$ 融合空间或 4‑Majorana 零模子空间中的作用，与抽象 TQFT 的 half twist 同构（差一个整体相位和基变换）。

**2.2 YBE 与“可积/近平坦子流形”**

在 $H_P$ 空间中，要求 R 满足常数 YBE 就等价于对 $H_P$ 的一组代数约束：小参数展开给出 classical YBE
$$
[H_{12},H_{13}] + [H_{12},H_{23}] + [H_{13},H_{23}] = 0,
$$
对应“线性化”的平坦联络条件，进一步要求完整的 operator YBE 成立，则得到真正的可积子流形。偏离这些子流形的扰动 $V$ 在几何上表现为 Berry 曲率 $F$ 变大，在数值上则表现为 Dehn twist 平台塌陷、braid fidelity 和 LQC 复杂度恶化。

**这一层结论**：我们已经把抽象的 R‑矩阵嵌入到了具体的晶格 H 中，且能够识别出“哪些 H\_P 给出了理想的 half twist / Dehn twist（可积/近平坦）”，哪些 H\_P 会引入强曲率和非可积效应。

 

#### 3. 几何层：配置空间、Berry 曲率与“空间诱导”的拓扑

接下来把视角从“时间演化”转到“配置空间几何”。

**3.1 缺陷配置空间与 Hilbert 丛**

在 2D p+ip 或 Kitaev 蜂窝模型中，引入若干缺陷（涡旋、vison 等），记其在底空间 $\Sigma$ 上的位置为
$$
X = (x_1,\dots,x_n)\in \mathcal C = \{X\mid x_i\in\Sigma,\ x_i\neq x_j\}.
$$
对每个 $X$ 有一个哈密顿量 $H(X)$ 和对应的低能子空间（零模/基态子空间）
$$
\mathcal H_0(X) = \mathrm{span}\{\text{零能模态/简并基态}\}.
$$
这些子空间拼接成一个 Hilbert 丛
$$
\pi:\mathcal E\to\mathcal C,\qquad \pi^{-1}(X)=\mathcal H_0(X).
$$

在局部选取一组正交归一基 $|\psi_a(X)\rangle$，定义投影 $P(X)=\sum_a |\psi_a(X)\rangle\langle\psi_a(X)|$。Berry 联络与曲率为
$$
 A_\mu(X) = iP\,\partial_\mu P\,P,\qquad
 F_{\mu\nu}(X) = P[\partial_\mu P,\partial_\nu P]P.
$$
给定一条配置空间中的闭合路径 $\gamma:[0,1]\to\mathcal C$，Berry holonomy
$$
U[\gamma] = \mathcal P\exp\Bigl(-\int_0^1 A_\mu(\gamma(s))\,\dot X^\mu(s)\,ds\Bigr)
$$
是在零模子空间上的幺正算符。

若曲率 $F=0$ 或仅取中心值，则 $U[\gamma]$ 只依赖于路径的同伦类，[$\gamma$] 对应于 braid group / mapping class group 的一个元素，从而在零模子空间上给出一个拓扑不变量（例如 Ising 的 half twist / Dehn twist 矩阵）。

**3.2 “时间中的运动” vs “空间中的重构”**

通常我们把 $\gamma$ 理解为“随时间缓慢移动缺陷”的路径：

- 参数 $s\in[0,1]$ 是 rescaled time；
- 系统哈密顿量 $H(X(s))$ 随时间变化，产生 adiabatic 演化；
- Berry holonomy 就是这段时间演化在零模子空间上的投影。

但几何上，$\gamma$ 只是 $\mathcal C$ 内的一条曲线 —— 我们同样可以把它理解成“空间上的重构方案”：

- 在 1D/2D 网络中，通过改变哪些 H\_P（或者等价地，打开/关闭哪些键、在哪些键施加 R），在空间上重写编码，使得**逻辑上等同**于“缺陷沿着 $\gamma$ 移动”；
- 如果我们以有限步和局域门来近似这一操作，每一步对应于在某个局部块上施加一个或若干 $e^{iH_P}$，这些操作在时间上可以高度并行（常数深度），而在配置空间上仍然实现了同一个同伦类的 $\gamma$。

**这一层的核心思想**：

> Berry holonomy 是配置空间路径的函数，而不是物理时间参数化的函数。只要我们在 Hilbert 丛 $\mathcal E$ 上实现了同一个同伦类的 $\gamma$，无论是通过“真实移动缺陷的时间演化”，还是通过“空间重构 + 局域 R 门的有限深序列”，拓扑上得到的 holonomy（在平坦/近平坦区域内）是相同的或仅差一个小扰动。

 

#### 4. “空间诱导拓扑”的准确性：从 R‑word 到几何扭转

现在我们把三层结构拼在一起，给出这一思路的较为精确的表述与论证框架。

**4.1 从抽象拓扑门到 R‑word**

给定一个抽象的拓扑门 $U_{\mathrm{top}}$（例如某个 braid word 或 Dehn twist），在 TQFT 层面我们可以把它写成 F、R 的有限 word：
$$
U_{\mathrm{top}} \simeq W(F,R^{\sigma\sigma})
$$
在相应的融合空间中定义。选择一个具体的 R‑矩阵解 $R^{\sigma\sigma}$，再通过
$$
R^{\sigma\sigma} = e^{iH_P^{(0)}}
$$
反推出对应的两体生成元 $H_P^{(0)}$。

在 $H_P$ 空间中，我们可以要求：

1. $H_P^{(0)}$ 尽量落在一个可积分或 classical YBE 残差很小的子流形上（近平坦层）；
2. 把 word 中的每个 R 替换为适当的 $e^{iH_P^{(0)}}$ 或其变体，形成一个在微观 Hilbert 空间中的有限串 $\widetilde{W}$。

**4.2 从 R‑word 到配置空间路径**

接下来，在具体模型（1D Kitaev 链、2D p+ip、蜂窝等）中布置缺陷和局域 H\_P：

1. 把 $H_P^{(0)}$ 放到适当的局域块上，确保局域 4‑Majorana / 2‑$\sigma$ 子空间的零模结构与 TQFT 融合空间同构；
2. 设计一个仅用有限轮局域门（每轮在空间上高度并行）的电路，实现 word $\widetilde{W}$ 在零模子空间上的作用；
3. 在配置空间 $\mathcal C$ 中，这等价于把缺陷位置从初态 $X_{\mathrm{i}}$ 通过若干短路径段拼接到终态 $X_{\mathrm{f}}$，其同伦类与原本的 $U_{\mathrm{top}}$ 对应的 braid / Dehn twist 相同。

在完全平坦的理想情况（例如精确可积点）下，任何两条同伦等价的路径给出相同的 holonomy；在近平坦的实际情况中，沿不同拼接方案得到的 holonomy 在“拓扑门”的 SU(2) 共轭类附近波动，其偏差可以用 Berry 曲率和 complexity 曲率来严格控制（见 kit-exp 6.4–6.5）。

**4.3 从时间演化到空间重构：等价性的论证方向**

要证明“空间诱导”的方案与传统的“时间中的 adiabatic 编织”在拓扑上等价，可以沿着以下路线：

1. 在 Hilbert 丛的语言中，把时间路径视为 $\mathcal C$ 中的一条曲线 $\gamma(t)$；
2. 构造一个仅用有限轮局域操作的“空间重构方案”，对应另一条曲线 $\gamma'(s)$，使得 $\gamma$ 与 $\gamma'$ 在 $\mathcal C$ 中同伦；
3. 在平坦/近平坦区域内，利用 $F\approx 0$ 的事实，证明（或数值上验证）
	 $$
	 U[\gamma'] \approx U[\gamma] \approx U_{\mathrm{top}}
	$$
	 且误差上线性/二次地受扰动强度、曲率上界和路径长度控制；
4. 通过 run\_dehn\_twist\_micro\_vs\_berry.py 等脚本的结果，给出“微观时间演化 vs Berry holonomy vs 空间重构 R‑word”三者之间的一致性数值证据。

严格的数学证明需要对 $F$ 的上界和路径同伦的细节做更多假设（如谱隙条件、adabatic 理论等），但我们已经有了一个明确的框架：

- 抽象类型论保证 F、R‑word 与拓扑门 $U_{\mathrm{top}}$ 的等价；
- R = e^{iH_P} 把这个 word 嵌入到具体的 lattice 哈密顿量中；
- Berry 曲率和平坦性保证“同伦等价路径”给出近似相同的 holonomy；
- 数值上，我们可以用 Dehn twist plateau、复杂度曲率 $K_{\mathrm{comp}}$ 等观测量来检验这一框架的有效性。

 

#### 5. 总结：从 YBE anyon 编织到空间诱导拓扑

通过上述三层结构，我们可以把原来的直觉归纳为一句话：

> 我们不再仅仅把 anyon 编织看作“随时间缓慢移动缺陷”的过程，而是通过 $R=e^{iH_P}$ 的指数表示，把抽象的半扭转 / Dehn twist 门嵌入到具体的晶格哈密顿量与缺陷配置空间中，从而在**空间重构**的层面上实现同一类拓扑幺正——其可靠性由 $H_P$ 空间中的 YBE / classical YBE 条件和 Berry 曲率的平坦性来量化。

这为后续工作提供了三个可扩展方向：

1. **曲率工程**：在 $H_P$ 空间中根据 classical YBE 和 $F$ 选择“好方向”与“坏方向”，为实验设计“拓扑更稳健”的路径；
2. **complexity 曲率**：用有限深局域电路拟合 Berry holonomy，定义 $K_{\mathrm{comp}}$ 来测量非可积性与拓扑稳健性的退化；
3. **空间重构方案**：系统构造一类常数深度（或对系统尺寸次线性）的空间重构电路，利用 R‑word 实现目标拓扑门，并在上述几何/复杂度框架下分析其鲁棒性与可实现性。

这些方向都建立在这里的“YBE anyon 编织 → R = e^{iH_P} → 配置空间几何扭转”这条思路之上，既整理和统一了我们之前的数值与分析，也为进一步的理论和实验探索提供了可操作的工具。

 

#### 6. 定理与引理框架（证明骨架草稿）

本节给出一个“准定理”的形式，把上面直觉性的等价关系写成精确的语句，并拆成若干可以分别攻击的引理。这里不追求完全严密的泛函分析细节，而是明确需要哪些假设、误差项如何分解，为后续严格化/写论文做准备。

**6.1 主命题（粗略版）**

设：

1. 给定一个有限维拓扑量子场论（例如 Ising TQFT），以及其中的某个“目标拓扑门” $U_{\mathrm{top}}$（如一个 braid word 或 Dehn twist）；

2. 选定一个满足常数 YBE 的 R‑矩阵 $R^{\sigma\sigma}$，并在两比特空间上选取一个生成元 $H_P^{(0)}$ 使得
$$
	R^{\sigma\sigma} = e^{iH_P^{(0)}}\big|_{\mathcal H_0}
$$

在零模/融合子空间 $\mathcal H_0$ 上与抽象的 R‑符号同构；

3. 在某个 2D 晶格模型（p+ip, 蜂窝等）中，引入缺陷配置空间 $\mathcal C$ 与 Hilbert 丛 $\mathcal E\to\mathcal C$，并假设存在一个带谱隙的参数区域 $\Omega\subset\mathcal C$，其 Berry 曲率满足 $\|F(X)\|\le \kappa$ 对所有 $X\in\Omega$ 成立；

4. 构造一条“理想时间路径” $\gamma:[0,1]\to\Omega$ 对应于在配置空间中的传统 adiabatic 编织方案，以及一条“空间重构路径” $\gamma':[0,1]\to\Omega$，对应于有限轮局域 $e^{iH_P}$ 门构成的 R‑word 电路，两条路径在 $\Omega$ 内同伦；

5. 假设整个演化过程中哈密顿量始终具有谱隙 $\Delta>0$，并满足足够的光锥/局域性条件（Lieb–Robinson 型）。

则在零模子空间上，存在一个整体相位 $e^{i\phi}$，使得
$$
\bigl\|U_{\mathrm{spatial}} - e^{i\phi}U_{\mathrm{top}}\bigr\|
\;\le\; C_1\,\varepsilon_{\mathrm{YBE}} + C_2\,\kappa\,\mathcal A_{\gamma,\gamma'} + C_3\,\varepsilon_{\mathrm{adiab}} + C_4\,\varepsilon_{\mathrm{Trotter}},
$$
其中：

- $U_{\mathrm{spatial}}$ 是由有限轮局域 $e^{iH_P}$ 门实现的“空间重构”幺正；
- $\varepsilon_{\mathrm{YBE}}$ 衡量所选 $H_P^{(0)}$ 偏离精确 YBE/经典 YBE 的程度；
- $\kappa$ 是 $\Omega$ 内 Berry 曲率的上界，$\mathcal A_{\gamma,\gamma'}$ 是 $\gamma$ 与 $\gamma'$ 夹住的面积（或更一般的“同伦 2‑链”体积）；
- $\varepsilon_{\mathrm{adiab}}$ 是 adiabatic 近似误差，依赖于演化速度与谱隙 $\Delta$；
- $\varepsilon_{\mathrm{Trotter}}$ 是用有限步局域门近似连续生成元的 Trotter/离散化误差；
- $C_i$ 是与具体模型/几何常数相关的数值系数。

换句话说，只要 YBE 偏差小、Berry 曲率有界且路径同伦“窄”、演化足够慢且离散近似足够精细，就可以把空间重构实现的门 $U_{\mathrm{spatial}}$ 看成与抽象拓扑门 $U_{\mathrm{top}}$ 在同一个 SU(2) 共轭类附近，且误差可以定量估计。

**6.2 引理 A：平坦/近平坦联络下的路径同伦稳定性**

设 Berry 联络 $A$ 在区域 $\Omega\subset\mathcal C$ 内的曲率满足 $\|F\|\le \kappa$。则对任意两条位于 $\Omega$ 内、具有相同起终点且同伦的闭合路径 $\gamma_1,\gamma_2$，有
$$
\bigl\|U[\gamma_1]-U[\gamma_2]\bigr\|\;\le\; C\,\kappa\,\mathcal A_{\gamma_1,\gamma_2},
$$
其中 $\mathcal A_{\gamma_1,\gamma_2}$ 是任意一个“填充”两条路径之间的 2‑链的面积上界，$C$ 为常数。

证明思路：利用非阿贝尔 Stokes 公式，将 holonomy 之差写成曲率在某个 2‑链上的有序指数，与 $\|F\|$ 和面积上界直接关联。平坦情形 $F\equiv 0$ 下则给出严格不变性 $U[\gamma_1]=U[\gamma_2]$。

**6.3 引理 B：adiabatic 时间演化与 Berry holonomy 的等价**

设 $H(X)$ 在路径 $\gamma(t)$ 上谱隙 $\Delta>0$ 下有界，且 $H$ 对参数变化足够光滑。令 $U_{\mathrm{phys}}(T)$ 为实际时间演化算符（以合适的慢速走完 $\gamma$），则在零模子空间上
$$
\bigl\|U_{\mathrm{phys}}(T) - e^{i\phi}\,U[\gamma]\bigr\|\;\le\; \varepsilon_{\mathrm{adiab}}(T,\Delta,\|\dot X\|),
$$
其中 $\varepsilon_{\mathrm{adiab}}\to 0$ 当演化时间 $T\to\infty$、速度 $\|\dot X\|\to 0$。

这是标准的 adiabatic 定理与 Berry 位相之间的关系，这里只需引用或略述证明思路：在带谱隙的情形，演化在瞬时本征子空间中的投影由 Berry 联络控制，非 adiabatic 跳跃概率在慢极限下指数/幂次地衰减。

**6.4 引理 C：抽象 F,R‑word 与 e^{iH\_P} 嵌入的一致性**

在理想可积点上，选择 $H_P^{(0)}$ 使得
$$
R^{\sigma\sigma}=e^{iH_P^{(0)}}\big|_{\mathcal H_0},
$$
且 R 满足常数 YBE。则由 YBE 的代数结构，可知将 R 放到多粒子 Hilbert 空间上定义的 braid group 表示，与 TQFT 中通过 F,R 给出的表示同构；因此，对任意拓扑门 $U_{\mathrm{top}}\simeq W(F,R^{\sigma\sigma})$，有
$$
U_{\mathrm{top}} \simeq \widetilde{W}\bigl(e^{iH_P^{(0)}}\bigr)\big|_{\mathcal H_0}.
$$

在小扰动 $H_P = H_P^{(0)}+V$ 下，只要谱隙不闭合，拓扑序的稳定性与 quasi‑adiabatic continuation 理论保证这一表示仅作连续变形，不跨越拓扑相：即存在一个局域幺正 $W$ 使得
$$
e^{iH_P}\big|_{\mathcal H_0} \approx W\,e^{iH_P^{(0)}}\big|_{\mathcal H_0}\,W^{-1},
$$
从而 F,R 表示仍与 TQFT 同属一个等价类，误差由 $\|V\|$ 和局域性常数控制（见 kit-exp 6.4 中的 Dyson 估计作为线性近似）。

**6.5 引理 D：空间重构电路对 quasi‑adiabatic 演化的逼近**

若存在一个准局域生成元 $K(X)$，使得 adiabatic 演化可写为
$$
U_{\mathrm{adiab}}[\gamma'] = \mathcal T\exp\Bigl(-i\int_0^1 K(X(s))\,ds\Bigr),
$$
且 $K$ 是由一簇局域 $H_P$ 函数构成，则可以通过 Trotter–Suzuki 分解和局域门近似构造一个有限深度电路 $U_{\mathrm{spatial}}$，满足
$$
\bigl\|U_{\mathrm{spatial}} - U_{\mathrm{adiab}}[\gamma']\bigr\|\;\le\; \varepsilon_{\mathrm{Trotter}}(N,\|K\|,\tau),
$$
其中 $N$ 是步数/电路深度，$\tau$ 是每步的有效时间片。$\varepsilon_{\mathrm{Trotter}}$ 在 $N\to\infty$、$\tau\to 0$ 时趋于 0，且可以用标准 Trotter 误差界来估计。

**6.6 主命题的拼接**

在上述引理成立的前提下，我们可以按如下顺序拼接误差：

1. 在理想可积点上，利用引理 C 得到
	$$
	U_{\mathrm{top}} \simeq \widetilde{W}\bigl(e^{iH_P^{(0)}}\bigr)\big|_{\mathcal H_0}.
	$$
2. 小扰动下，使用拓扑相稳定性与 Dyson 展开，得到
	$$
	e^{iH_P}\big|_{\mathcal H_0} \approx e^{iH_P^{(0)}}\big|_{\mathcal H_0} + O(\|V\|),
	$$
	从而 $\widetilde{W}(e^{iH_P})\big|_{\mathcal H_0}$ 与 $U_{\mathrm{top}}$ 的差异受 $\varepsilon_{\mathrm{YBE}}$ 控制。
3. 对应 R‑word 的配置空间路径 $\gamma'$ 与传统时间路径 $\gamma$ 在 $\Omega$ 内同伦，引理 A 给出
	$$
	U[\gamma'] \approx U[\gamma],
	$$
	误差 $\sim \kappa\,\mathcal A_{\gamma,\gamma'}$。
4. 引理 B 将 $U_{\mathrm{phys}}(T)$ 与 $U[\gamma]$ 联系起来，误差为 $\varepsilon_{\mathrm{adiab}}$。
5. 引理 D 说明用有限深空间重构电路可以逼近 $U_{\mathrm{adiab}}[\gamma']$，误差为 $\varepsilon_{\mathrm{Trotter}}$。

综合以上各项，得到 6.1 节的主不等式，这就把“空间重构实现的拓扑门”与“抽象 YBE/TQFT 中的拓扑门”之间的关系，系统地分解为四种可控误差源。进一步的严格化工作可以在每个引理内部补充更多技术细节与文献引用。

 

#### 7. 4‑Majorana toy 模型：从时间 half twist 到空间重构

为了把上述抽象框架具体化，这里构造一个最小的 4‑Majorana 模型，用来演示如何在不显式“移动缺陷”的情况下，通过局域的 $R=e^{iH_P}$ 门和配对图重写，在零模子空间上实现与 half twist 等价的幺正。

**7.1 模型与编码：4 个 Majorana 与逻辑子空间**

考虑 4 个 Majorana 算符
$$
\gamma_1,\gamma_2,\gamma_3,\gamma_4,\qquad \{\gamma_a,\gamma_b\}=2\delta_{ab}.
$$
用它们构造两组复费米子
$$
f_L = \frac{\gamma_1+i\gamma_2}{2},\qquad f_R = \frac{\gamma_3+i\gamma_4}{2}.
$$
总 Hilbert 空间是 4 维，约束总费米数为偶数后得到 2 维逻辑子空间（一个有效 qubit），可以取基
$$
|0_L\rangle = |0_L\rangle\otimes|0_R\rangle,\qquad |1_L\rangle = |1_L\rangle\otimes|1_R\rangle
$$
（下标 $L,R$ 表示左右复杂费米模的占据数，这里只示意逻辑编码）。

我们关心的是在这一逻辑 2 维空间上的幺正作用，特别是 half twist
$$
U_{\mathrm{ht}} = \exp\Bigl(\frac{\pi}{4}\,\gamma_2\gamma_3\Bigr),
$$
它在 Ising TQFT 里对应两个 $\sigma$ anyon 的一次交换。

**7.2 “时间中的 half twist”：标准连续演化**

传统的做法是：选取哈密顿量
$$
H_0 = \frac{i}{2}\,\gamma_2\gamma_3
$$
并在时间 $\tau=\pi/2$ 内 adiabatic 地演化：
$$
U_0(\tau) = e^{-iH_0\tau} = \exp\Bigl(\frac{\pi}{4}\,\gamma_2\gamma_3\Bigr) = U_{\mathrm{ht}}.
$$
这可以理解为在缺陷配置空间中，沿着一条使得 2 与 3 逐渐成对的路径缓慢移动涡旋 —— 典型的“时间中的编织”。

**7.3 配对图重写：三步空间重构方案**

我们改用“配对图重写”的方式，不显式移动缺陷，而是通过改变哪些 Majorana 成对耦合，来在空间上实现同一个 half twist。考虑三种配对模式：

1. 初态配对（A 型）：
$$
	H^{(A)} = \frac{iJ}{2}\,(\gamma_1\gamma_2 + \gamma_3\gamma_4),
$$
对应 $(1,2)$、$(3,4)$ 两条“键”被拉紧。

2. 中间配对（B 型）：
$$
	H^{(B)} = \frac{iJ}{2}\,(\gamma_1\gamma_3 + \gamma_2\gamma_4),
$$
对应 $(1,3)$、$(2,4)$ 成对（对角配对）。

3. 终态配对（C 型）：
$$
	H^{(C)} = \frac{iJ}{2}\,(\gamma_1\gamma_4 + \gamma_2\gamma_3),
$$
对应 $(1,4)$、$(2,3)$ 成对。

这三种配对可以视为“在平面上固定四个端点，通过重新连线重排配对”的三种拓扑等价类：A 是“横向”连线，C 是“纵向+中间”的连线，它们之间通过 B 作为中间形态平滑连接。

假设在每一种配对模式下，都存在一个足够大的谱隙，将我们限制在最低能的 2 维逻辑子空间内。在此子空间中，沿着路径
$$
H^{(A)} \to H^{(B)} \to H^{(C)}
$$
的 adiabatic 演化给出某个 holonomy $U_{ABC}$。直观上，这个序列等价于在平面上把内部的连线“翻转”一次，对应一次 half twist。

**7.4 用 R = exp(iH\_P) 表示局部步骤**

在我们的框架中，每一步“配对模式切换”可以分解成若干局域 $R$‑门。比如，从 A 切到 B 可以通过逐渐减小 $\gamma_1\gamma_2$、$\gamma_3\gamma_4$ 的耦合，同时打开 $\gamma_1\gamma_3$、$\gamma_2\gamma_4$ 的耦合；离散化后可以写成有限多步：
$$
U_{A\to B} \approx \prod_k \exp\bigl(i\,\delta t_k\,H_P^{(k)}\bigr),
$$
其中每个 $H_P^{(k)}$ 只在某两对 Majorana 上非零，对应两比特空间上的某个 $H_P$ 及其指数 $R^{(k)}=e^{iH_P^{(k)}}$。类似地，从 B 切到 C 也由一串 $R^{(k)}$ 组成：
$$
U_{B\to C} \approx \prod_l R^{(l)}.
$$

整个“空间重构”幺正就是
$$
U_{\mathrm{space}} = U_{B\to C}\,U_{A\to B},
$$
它在零模子空间上的作用 $U_{\mathrm{space}}\big|_{\mathcal H_0}$ 可以直接通过 4×4 或 2×2 矩阵算出，并与 $U_{\mathrm{ht}}$ 比较。在我们之前的数值脚本（如 verify/kit1/majorana\_braid\_check.py）中，已经验证了类似的配对图重写方案可以在数值上给出与 $\exp((\pi/4)\gamma_2\gamma_3)$ 共轭的 SU(2) 矩阵，这为“空间重构 = half twist”提供了具体例证。

**7.5 与 Dirac vortex 思路的对比**

在 Dirac vortex 的图像中，涡旋位置固定，而通过改变局域配对/耦合（或 branch cut 的走向），把零模“拖动”到不同的物理位置，最终实现一个拓扑门。上面的 4‑Majorana 模型正是这种思想在离散/有限维系统中的最简实现：

- 涡旋/端点（$\gamma_i$ 的物理位置）保持固定；
- 我们只是在“空间上重写”谁和谁成对（通过 $H^{(A)},H^{(B)},H^{(C)}$ 及其对应的 $R=e^{iH_P}$）；
- 在零模子空间上，这一重构与连续地打开 $H_0=(i/2)\gamma_2\gamma_3$ 的时间演化产生的 half twist 幺正等价（或在近平坦区域内仅差一个小扰动）。

因此，这个 4‑Majorana toy 模型为我们的主张提供了一个极简的“实验室”：

> 在带谱隙的拓扑系统中，通过 R = exp(iH\_P) 的局域门和配对图重写，可以在“空间上”实现与传统时间编织相同的拓扑门，而 Berry 曲率/YBE 结构为这一等价关系提供了几何和代数上的解释。



