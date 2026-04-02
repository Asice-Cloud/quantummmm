### Ising 相的试验：从 R 到任意子行为的最小例子

这一节专门用一个尽量小的 Ising 相例子，顺着

$$
R(a,b,c,d)\;\Rightarrow\;H_{\text{Kitaev}}\;\Rightarrow\;\text{2D+Z}_2\;\Rightarrow\;\text{Ising TQFT/MTC}\;\Rightarrow\;\text{Majorana/任意子行为}
$$

这条链走一遍，把前面抽象讨论的内容具体化。目标是：

- 选一组 $(a,b,c,d)$，其 2D 延拓在某参数区对应 Ising 类拓扑相；
- 写出对应的 1D Kitaev 链参数 $(t,\Delta,\mu)$ 与 2D BdG；
- 在一个非常小的几何上（比如两条 Kitaev 链形成的 T 形线网）找到 Majorana 零模；
- 在这个零模子空间上构造一个由 R 导出的“交换/半编织”算符，并与 Ising 任意子的 $R$ 矩阵进行对比，看它实现了怎样的粒子行为。

后续可以在本文件中逐步补全：先锁定具体的 $(a,b,c,d)$ 与 1D/2D 模型，再添加 Ising MTC 的 $F,R$ 数据和数值对比。

---

#### 1. 具体选择一组 R 与 Kitaev 链参数

在 1D 分析中，我们已经看过一个特别简单、物理直观的例子：
$$
a=0,\quad d=0,\quad b=1,\quad c=0.
$$
代入 1D 映射
$$
t\propto(b+c),\qquad \Delta\propto(b-c),\qquad \mu\text{ 由 }a,d\text{ 给出},
$$
可以直接取归一化为
$$
t=1,\qquad \Delta=1,\qquad \mu=0.
$$

这对应的 1D Kitaev 链哈密顿为
$$
H=\sum_j\Big[-(c_j^\dagger c_{j+1}+h.c.)+ (c_j c_{j+1}+h.c.)\Big],
$$
即标准形式中 $t=\Delta,\mu=0$ 的极限，处在拓扑相的中心（$|\mu|<2|t|$）。开边界时，两端各有一个 Majorana 零模，是典型的 1D Majorana 链。

在 Majorana 记号下，沿链的最近邻耦合可以重写为
$$
H= \frac i2\sum_j \big(\gamma_{2j}\gamma_{2j+1}\big),
$$
附近的局域双线性 $i\gamma_{2j}\gamma_{2j+1}$ 便是后续构造 braiding/交换算符的基本积木。

这一选择与前文从 R 导出的理想 braid 生成元
$$
U_{j,j+1}=\exp\Big(\frac\pi4\,\gamma_{2j}\gamma_{2j+1}\Big)
$$
是兼容的：在两零模子空间上，$U$ 实现 $\gamma_{2j}\mapsto\gamma_{2j+1},\ \gamma_{2j+1}\mapsto-\gamma_{2j}$，对应单对 Majorana 的一次“半交换”。

---

#### 2. 2D p+ip 扩展与 Ising 相

把上面的 1D 链扩展到 2D，可取一个标准的 spinless p+ip 超导格点模型（方格为例）：
$$
H=\sum_{\langle ij\rangle}\Big[-t(c_i^\dagger c_j+h.c.)+\Delta_{ij}(c_i c_j+h.c.)\Big]-\mu\sum_i(n_i-\tfrac12),
$$
其中
$$
\Delta_{i,i+\hat x}=+\Delta,\qquad \Delta_{i,i-\hat x}=-\Delta,\\
\Delta_{i,i+\hat y}=+i\Delta,\qquad \Delta_{i,i-\hat y}=-i\Delta,
$$
给出动量空间中典型的
$$
\Delta(\mathbf k)\propto\sin k_x + i\sin k_y,
$$
即 $p_x+ip_y$ 配对结构。取 $t=\Delta=1,\mu$ 落在 Chern=1 的参数区（例如 $-2<\mu<0$ 的一段），体系处于拓扑相，其低能极限由 Ising TQFT 描述：

- 粒子类型：$1,\psi,\sigma$；
- $\sigma$ 涡旋内束缚一个 Majorana 零模；
- $\psi$ 对应费米准粒子；
- 融合规则 $\sigma\times\sigma=1+\psi$ 等。

我们的 1D R→Kitaev 选择保证了沿 $x$ 方向的链在 $t=\Delta$、拓扑相中心，实现了每条“线”上的 Majorana 模式；在 2D 上连接成网络，并引入适当的涡旋/边界，即可得到标准的 Ising 任意子实现。

---

#### 3. 两个涡旋的 Majorana 与 braid 算符：微观 U 与 Ising 的 R

在上述 2D p+ip 模型中，考虑两个相距较远的 $\pi$‑涡旋，各自束缚一个零能 Majorana 模式，记为 $\gamma_1,\gamma_2$。在低能空间里，这两者生成一个 2 维 Hilbert 空间（总费米数守恒时等价于 Ising 中 $\sigma\times\sigma=1+\psi$ 的两条融合通道）。

在 Majorana 语言中，一个自然的“编织/交换”算符取为
$$
U_{12}=\exp\Big(\frac\pi4\,\gamma_1\gamma_2\Big).
$$
其对 Majorana 算符的共轭作用为
$$
U_{12}\,\gamma_1\,U_{12}^{-1}=\gamma_2,\qquad U_{12}\,\gamma_2\,U_{12}^{-1}=-\gamma_1.
$$

引入复费米算符
$$
f=\frac{\gamma_1+i\gamma_2}{2},\qquad f^\dagger=\frac{\gamma_1-i\gamma_2}{2},
$$
在占据数本征基 $\{|0\rangle,|1\rangle=f^\dagger|0\rangle\}$ 上，$U_{12}$ 的矩阵形式为
$$
U_{12}=e^{-i\pi/4}\begin{pmatrix}1&0\\0&i\end{pmatrix},
$$
即在偶/奇费米数子空间上分别赋予相位 $1$ 与 $i$（整体相位 $e^{-i\pi/4}$ 可忽略）。

另一方面，在 Ising 任意子范畴中，两 $\sigma$ 任意子在固定背景下的交换算符 $R^{\sigma\sigma}$ 在 fusion 基 $\{|1\rangle,|\psi\rangle\}$ 上通常写为
$$
R^{\sigma\sigma}=\begin{pmatrix}e^{-i\pi/8}&0\\0&e^{3i\pi/8}\end{pmatrix}.
$$
两者之间的关系是：

- 通过合适的基选择和整体相位重标度，$U_{12}$ 与 $R^{\sigma\sigma}$ 可以在 2 维子空间内视为“同一个编织”的两种实现；
- 更精确地，$U_{12}$ 实现的是 Ising 模型中 $\sigma$‑$\sigma$ 交换的某个规范选择下的表示，与标准 $R^{\sigma\sigma}$ 只差一个全局相位与可能的 $F$‑重排。

因此，在这一具体选择下：

- 从 $R(a,b,c,d)$ 映射到 1D/2D Majorana 模型得到的局域双线性 $\gamma_1\gamma_2$，给出的 $U_{12}=e^{\frac\pi4\gamma_1\gamma_2}$，
- 在两零模子空间上，等价于 Ising 任意子理论中 $\sigma$‑$\sigma$ 的一次编织操作，

这为 “R→Kitaev→Majorana→Ising 任意子” 的链条提供了一个具体、可计算的最小例子。

后续可以在此基础上：

- 引入第三个/第四个涡旋，研究多 $\sigma$ 任意子 fusion 空间上的 $F,R$ 矩阵和 R_to_Kitaev 给出的 $U$ 的对应关系；
- 在带边界/环面几何上考察 Dehn twist 与上文构造的 $U_\gamma$ 在零模子空间上的作用，进一步对接 Ising TQFT 的 mapping class 表示。

---

#### 4. 四个涡旋（四个 Majorana）上的 $F+R$ 结构

进一步考虑四个相距较远的 $\pi$‑涡旋，各自束缚一个零能 Majorana 模式 $\gamma_1,\gamma_2,\gamma_3,\gamma_4$。在拓扑极限下（忽略零模之间的小能级劈裂），这四个零模张成一个 4 维费米 Hilbert 空间；若固定总费米数奇偶（例如偶数），则有效拓扑子空间维数为 2，对应 Ising 范畴中四个 $\sigma$ 任意子的融合空间 $V_{\sigma\sigma\sigma\sigma}$ 的 2 维子空间：
$$
\sigma\times\sigma\times\sigma\times\sigma\to 1\oplus1\oplus\psi\oplus\psi,\qquad
\dim V^{\text{(fixed total)}}=2.
$$

在 Majorana 记号中，引入两个复费米算符
$$
f_1=\frac{\gamma_1+i\gamma_2}{2},\qquad f_2=\frac{\gamma_3+i\gamma_4}{2},
$$
令真空 $|0\rangle$ 满足 $f_1|0\rangle=f_2|0\rangle=0$。偶费米数子空间由
$$
|00\rangle=|0\rangle,\qquad |11\rangle=f_1^\dagger f_2^\dagger|0\rangle
$$
张成，是一个 2 维空间，可与 Ising 范畴中四 $\sigma$ 粒子在总拓扑荷为 $1$ 或 $\psi$ 的子空间同构（不同基选择对应不同 fusion 树展开）。

四个涡旋的局域 braid 生成元可以选为沿直线顺次交换相邻涡旋，对应的 Majorana 算符为
$$
\begin{aligned}
B_1 &= \exp\Big(\frac\pi4\,\gamma_1\gamma_2\Big),\\
B_2 &= \exp\Big(\frac\pi4\,\gamma_2\gamma_3\Big),\\
B_3 &= \exp\Big(\frac\pi4\,\gamma_3\gamma_4\Big),
\end{aligned}
$$
它们在算符层面的作用为
$$
\begin{aligned}
B_1:&\;\gamma_1\mapsto \gamma_2,\ \gamma_2\mapsto-\gamma_1,\quad \gamma_3,\gamma_4\text{ 不变},\\
B_2:&\;\gamma_2\mapsto \gamma_3,\ \gamma_3\mapsto-\gamma_2,\quad \gamma_1,\gamma_4\text{ 不变},\\
B_3:&\;\gamma_3\mapsto \gamma_4,\ \gamma_4\mapsto-\gamma_3,\quad \gamma_1,\gamma_2\text{ 不变},
\end{aligned}
$$
构成四 Majorana（四 $\sigma$ 任意子）的 braid 群 $B_4$ 的一个表示。

在偶费米数子空间基 $\{|00\rangle,|11\rangle\}$ 上，可以显式计算其矩阵（忽略整体相位）为：
$$
\begin{aligned}
B_1 &\sim \begin{pmatrix}1&0\\0&i\end{pmatrix},\\[4pt]
B_3 &\sim \begin{pmatrix}1&0\\0&i\end{pmatrix},\\[4pt]
B_2 &\sim \frac{1}{\sqrt2}\begin{pmatrix}1&-i\\-i&1\end{pmatrix},
\end{aligned}
$$
其中“$\sim$” 表示可差一个整体 U(1) 相位（不影响拓扑结构）。

另一方面，在 Ising 任意子范畴中，四 $\sigma$ 任意子的 braid 生成元在标准 fusion 基下可由 $F$ 与 $R$ 组合给出：

- 选定 fusion 树 $((\sigma_1\times\sigma_2)\times(\sigma_3\times\sigma_4))$，中间融合通道取基 $\{|1\rangle,|\psi\rangle\}$；
- 则 $B_1,B_3$ 在该基下是对角的：
	$$
	B_1^{\text{(Ising)}}\sim B_3^{\text{(Ising)}}\sim \begin{pmatrix}e^{-i\pi/8}&0\\0&e^{3i\pi/8}\end{pmatrix},
	$$
	与前面两 Majorana 的 $U_{12}$ 一致；
- 中间交换 $B_2$ 则通过 $F$‑变换实现：
	$$
	B_2^{\text{(Ising)}} = F^{-1}\,R^{\sigma\sigma}\,F,
	$$
	其中非平庸的 $F$ 矩阵为
	$$
	F^{\sigma\sigma\sigma}_\sigma = \frac{1}{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix},
	$$
	在此基下可计算出
	$$
	B_2^{\text{(Ising)}}\sim \frac{1}{\sqrt2}\begin{pmatrix}1&-i\\-i&1\end{pmatrix},
	$$
	与上面的 Majorana 结果完全一致（同样只差整体相位规范）。

因此：

- 由 R_to_Kitaev 映射得到的四 Majorana braid 生成元 $B_1,B_2,B_3$，在偶费米数子空间上给出了一组满足 braid 关系的 2×2 矩阵；
- 这组矩阵与 Ising 任意子范畴中四 $\sigma$ 任意子的 braid 生成元（经 $F,R$ 组合）一一对应；
- 从范畴角度看，R_to_Kitaev 提供的是 Ising MTC 在某个具体 physical fusion 基下的低能实现。

这就完成了从“局域 YBE 解 R → 1D/2D Majorana 模型 → Ising 任意子 MTC 的 F+R 数据”的一个具体四粒子例子，验证了我们之前的判断：在合适参数区，R 确实“提升”为一个完整的 TQFT/MTC（Ising），Majorana braid 算符正是该范畴编织结构的格点实现。

---

#### 5. Ising TQFT/MTC 的完整数据与字典

为后续讨论参数变化与上层空间结构变化，整理 Ising 拓扑序对应的模张量范畴（MTC）基本数据，并与前述 Majorana 例子建立明确字典。

1) **粒子类型与融合规则**

Ising 范畴的简单对象（粒子类型）为
$$
\{1,\psi,\sigma\},
$$
满足融合规则：
$$
\begin{aligned}
\psi\times\psi &= 1,\\
\psi\times\sigma &= \sigma,\\
\sigma\times\sigma &= 1+\psi.
\end{aligned}
$$
其中 $1$ 为真空，$\psi$ 为费米准粒子，$\sigma$ 为携带 Majorana 零模的涡旋。

2) **量子维数与总量子维数**

量子维数：
$$
d_1=1,\quad d_\psi=1,\quad d_\sigma=\sqrt2.
$$
总量子维数为
$$
\mathcal D=\sqrt{d_1^2+d_\psi^2+d_\sigma^2}=2.
$$

3) **F 矩阵（非平庸部分）**

除了平凡的单位融合外，唯一非平庸的 F 矩阵出现在三个 $\sigma$ 的融合中：
$$
F^{\sigma\sigma\sigma}_\sigma = \frac{1}{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix},
$$
作用在中间通道为 $\{1,\psi\}$ 的 2 维空间上。其它 F 要么为 $1\times1$ 的平凡相位，要么由此和单位、对偶性约束确定。

4) **R 矩阵与拓扑自旋**

非平庸的 R 矩阵为：
$$
R^{\sigma\sigma} = \begin{pmatrix}e^{-i\pi/8}&0\\0&e^{3i\pi/8}\end{pmatrix},
$$
在 fusion 基 $\{|1\rangle,|\psi\rangle\}$ 上，对应两 $\sigma$ 粒子在中间融合通道为 $1$ 或 $\psi$ 时的交换相位。其它 R 元为：
$$
R^{\psi\psi}=-1,\qquad R^{\psi\sigma}=R^{\sigma\psi}=-i.
$$

topological spin（自旋）$\theta_a$ 由
$$
	heta_a = e^{2\pi i h_a}
$$
给出，对 Ising 粒子为
$$
	heta_1=1,\qquad \theta_\psi=-1,\qquad \theta_\sigma=e^{i\pi/8},
$$
即共形权重 $h_1=0,\ h_\psi=1/2,\ h_\sigma=1/16$。这与 $R^{\sigma\sigma}$ 的本征相位一致（相差 F‑重排与整休系数）。

5) **S,T 矩阵与 torus 上的 mapping class 表示**

在环面 $T^2$ 上，Ising TQFT 的模变换由 S,T 矩阵给出（以粒子类型 $\{1,\psi,\sigma\}$ 为基底）：
$$
T=\mathrm{diag}(\theta_1,\theta_\psi,\theta_\sigma)=\mathrm{diag}(1,-1,e^{i\pi/8}),
$$
而 S 矩阵为
$$
S=\frac{1}{2}
\begin{pmatrix}
1 & 1 & \sqrt2\\
1 & 1 & -\sqrt2\\
\sqrt2 & -\sqrt2 & 0
\end{pmatrix}.
$$

其中：

- $T$ 对应 torus 上沿某一基本环（通常记为 a‑cycle）的 Dehn twist；
- $S$ 对应两个基本环之间的交换（$S: (a,b)\mapsto(b,-a)$），与拓扑熵、互易关系密切相关；
- 它们满足模群上的基本关系 $(ST)^3 = S^2$ 等。

6) **与 Majorana 模型的字典（总结）**

结合前述两、四涡旋的例子，可以给出一个工作中的字典：

- 粒子 ↔ 零模/缺陷：
	- 真空 $1$ ↔ 无涡旋或成对湮灭后的基态；
	- $\psi$ ↔ 费米准粒子（Majorana 配对成 Dirac 模式的占据态）；
	- $\sigma$ ↔ 单个携带 Majorana 零模的涡旋（或某类 Z2 通量缺陷）。
- 两涡旋（两 $\sigma$）子空间：
	- 由两个 Majorana $\gamma_1,\gamma_2$ 张成，偶/奇费米数对应 $\{1,\psi\}$ 通道；
	- 算符 $U_{12}=\exp(\frac\pi4\gamma_1\gamma_2)$ 的矩阵形式与 $R^{\sigma\sigma}$ 本征相位一致（到整体相位与基不变换）。
- 四涡旋（四 $\sigma$）子空间：
	- 偶费米数 2 维 Hilbert 空间与 $V_{\sigma\sigma\sigma\sigma}$ 的 $\{1,\psi\}$ 中间通道空间同构；
	- Majorana 上的 braid 生成元 $B_1,B_2,B_3$ 与 Ising 范畴中的 $R,F$ 组合给出的 braid 生成元矩阵一一对应。
- 上层空间与 mapping class：
	- 在简单几何（线网、平面上涡旋配置）中，由局域 R_to_Kitaev 导出的 Majorana braid 算符 $U,B_j$ 给出了 $B_n$ 的表示；
	- 在更复杂几何（如环面）上，加上 Z2 规范结构和全局操作（如 $U_\gamma$ 形式的 Wilson 线和平移），可构造 Dehn twist 及其在零模子空间上的作用，与 T 矩阵给出的拓扑自旋兼容；
	- 若考虑不同 Z2 背景（不同 Wilson 环），则对应 torus 上不同 spin 结构/通量扇区，Ising TQFT 在这些扇区上有相应的态空间与 mapping class 表示，这部分可在后续进一步展开。

这一节整理的 Ising TQFT/MTC 数据，为后续“参数变化 → 上层空间结构变化 → Majorana 行为/耦合变化”的分析提供了标准参照系：

- 在某个 R‑参数区被识别为 Ising 相时，Majorana 模型的所有拓扑行为应当在合适基底下与上述 F,R,S,T 数据同构；
- 参数变化若不跨越相界，则对应的是在同一个范畴内的规范重标度或几何变形；
- 若参数变化导致 bulk gap 关闭并进入另一个拓扑相，则上层范畴（MTC）本身发生改变，需要用新的粒子/融合/F,R,S,T 数据来描述。




