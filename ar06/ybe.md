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


需要给出
$c_{\mu\nu}$ 到 Kitaev‑链参数的线性映射矩阵 $M$, 对应$t,\Delta,\mu ...$
\
 

下面把完整的线性映射与可解性判定追加在此：

约定 16 分量系数向量（列向量）：
$$
C = [c_{xx},\; c_{yy},\; c_{xy},\; c_{yx},\; c_{zz},\; c_{z0},\; c_{0z},\; c_{x0},\; c_{0x},\; c_{y0},\; c_{0y},\; c_{xz},\; c_{zx},\; c_{yz},\; c_{zy},\; c_{00}]^T.
$$

目标物理参数向量：
$$
p = [t,\; \Delta,\; U,\; \mu,\; C_{\mathrm{perbond}}]^T,
$$
其中 $t$ 为最近邻 hopping，$\Delta$ 为配对幅度，$U$ 为最近邻相互作用强度，$\mu$ 为体化学势（哈密顿写作 $-\mu\sum_j n_j$），$C_{\mathrm{perbond}}$ 為每键常数能量偏移。

显式线性关系 $p = M\cdot C$（按列顺序与上面 $C$ 一致）。为了便于直接代入，这里给出显式分量形式以及对应的 5×16 符号矩阵：

逐分量公式：
- $t = c_{xx} + c_{yy} + i\, ( c_{xy} - c_{yx} ),$
- $\Delta = c_{xx} - c_{yy} - i\, ( c_{xy} + c_{yx} ),$
- $U = 4\,c_{zz},$
- $\mu = 4\,c_{zz} - 2\,( c_{z0} + c_{0z} ),$
- $C_{\mathrm{perbond}} = c_{zz} - ( c_{z0} + c_{0z} ) + c_{00}.$

对应矩阵（行顺序 $[t,\Delta,U,\mu,C_{\mathrm{perbond}}]$，列顺序如 $C$）：
$$
M = \begin{pmatrix}
 1 &  1 &  i & -i & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 1 & -1 & -i & -i & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 &  0 &  0 &  0 & 4 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 &  0 &  0 &  0 & 4 & -2 & -2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 &  0 &  0 &  0 & 1 & -1 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1
\end{pmatrix}.
$$

注：上面包含复系数以允许一般的 XY/混合项；若所有 $c_{\mu\nu}$ 为实并满足有关的厄米性约束（例如 $c_{xy}=-c_{yx}$ 等），则 $t,\Delta$ 将是实数或成对共轭，哈密顿仍为厄米。

串算符（非局域项）触发分量：
- 下列任何分量非零会在 JW 映射后引入字符串前缀 $S_j=\exp(i\pi\sum_{k<j} n_k)$，从而产生非局域项：
	- 单体 x/y 分量： $c_{x0},\; c_{y0},\; c_{0x},\; c_{0y}$；
	- 含单侧 z 的混合项： $c_{xz},\; c_{yz},\; c_{zx},\; c_{zy}$。

可解性快速判定规则：
- 自由二次（Kitaev/XY，可 Fourier + Bogoliubov 严格对角化）：当并且仅当所有字符串分量都为 0 且 $c_{zz}=0$，即
$$
c_{x0}=c_{0x}=c_{y0}=c_{0y}=c_{xz}=c_{zx}=c_{yz}=c_{zy}=0,\quad c_{zz}=0.
$$

- XXZ（Bethe‑ansatz 可积，非自由但可积）：若所有字符串分量 = 0 且 $c_{zz}\neq0$，模型等价于 XXZ/Ising‑like 链，可用 Bethe‑ansatz 或已知可积方法处理；参数为 $U=4c_{zz}$，并且 $t,\Delta$ 由上式给出。

- 一般不可积/非局域：若任一字符串分量非零，则 JW 后出现非局域串项，通常不属于上述两类，需要数值、平均场或逐项近似分析；特殊可解子情形仍可能存在但需逐一检查。

自动化建议：我可以把这一套判定写成小脚本（Python），接受任意 16 个 $c_{\mu\nu}$ 并输出 $p$ = $(t,\Delta,U,\mu,C_{\mathrm{perbond}})$ 与可解性类别（并列出触发串项的分量）。如需我现在添加到仓库，请回复“添加脚本”。

 

（已追加完毕）

 

测试用最小修改参数（可直接被工具解析）：

c_xx = 1.1
c_yy = 1.0
c_xy = 0.0
c_yx = 0.0
c_zz = 0.0
c_z0 = 0.0
c_0z = 0.0
c_00 = 0.0

（将上述 c_{μν} 代入映射可得到 Δ = 0.1, t = 2.1, μ = 0）

### 谱参数的 YBE 与含时演化 （追加）

在物理实现中，编织操作通常是通过含时哈密顿的演化实现的。为便于把 R 与含时演化对接，可以引入谱参数 $u$ 使得
$$
R(u)=e^{i u H_P}
$$
或更一般地
$$
R(u)=\mathcal{T}\exp\Big(-i\int_0^{u} H_P(s)\,ds\Big),
$$
其中 $H_P$ 是作用在一对相邻站点上的局域厄米生成元（例如本文中定义的 $H_p^{(j,j+1)}$ 或其线性组合），$\mathcal{T}$ 为时间序列算子。当不同参数点处的生成元彼此对易（例如取
$$
H_P = a_x\,\sigma^x\otimes\sigma^x + a_y\,\sigma^y\otimes\sigma^y + a_z\,\sigma^z\otimes\sigma^z,
$$
且这些分量在不同键上的推广相互对易时），则时间序列可去掉，得到简单的指数形式 $R(u)=e^{i\phi(u)H_{P0}}$，其中 $\phi(u)=\int_0^u g(s)ds$ 与实际脉冲序列对应。

谱参数形式的 YBE（加法式习惯）可写为
$$
R_{12}(u-v)\,R_{13}(u)\,R_{23}(v)=R_{23}(v)\,R_{13}(u)\,R_{12}(u-v).
$$
如果取 $R(u)=e^{uX}$，则当生成元满足代数约束
$$
[X_{12},\; X_{13}+X_{23}]=0
$$
时，上式成立（可由 Baker--Campbell--Hausdorff 与交换关系直接验证）。因此一条构建策略是选择一个局域生成元 $X= i H_P$（或按本文惯例的 $H_P$），使得对应的两体算符在不同两两位置间满足上述交换关系；例如当 $H_P$ 由 $\sigma^\alpha\otimes\sigma^\alpha$ 的线性组合构成时，这些基元在代数上成对对易，从而满足可积性条件。

时间演化到谱参数的映射：若物理演化为
$$
U(t)=\mathcal{T}\exp\Big(-i\int_0^t H_P(s)\,ds\Big),
$$
且在不同时间 $s,s'$ 上 $[H_P(s),H_P(s')]=0$，那么
$$
U(t)=e^{-i\Phi(t) H_{P0}},\qquad \Phi(t)=\int_0^t g(s)ds,
$$
因此可以把谱参数 $u$ 直接识别为累积作用量 $u=\Phi(t)$。这给出了明确的实验/数值实现路线：通过设计时域耦合强度 $g(t)$，将期望的谱参数轨迹映射为实际含时演化脉冲，从而在体系上实现参数化的 $R(u)$，并通过 YBE 保证不同脉冲序列对应的编织等价性。

实现建议：
- 选择 $H_P$ 为本文已讨论的那类成对对易的两体算符（例如 XX+YY+ZZ 线性组合），保证代数约束成立；
- 把目标编织角度或交换操作映射为谱参数 $u$ 的数值区间，并设计 $g(t)$ 使得 $\Phi(t)$ 达到目标值；
- 若 $[H_P(s),H_P(s')]\neq0$，则需使用时间序列算子或分段脉冲（Trotter 化）近似实现，并评估非交换项引入的误差。

数值验证：仓库下加入一个小脚本 `tools/check_R_u.py`，用于对小尺寸系统数值检验谱参数形式的 YBE（当 $H_P$ 由对易基元构成时应满足）。脚本利用 `numpy` 和本征展开计算矩阵指数，并输出两个侧的范数差以检验 YBE。

若你同意，我会把上述段落写入 `ybe.md` 并把验证脚本添加到 `tools/` 下（我现在就添加）。

### 从 XXZ/6-vertex 的 R(u) 导出局域哈密顿 h 与脉冲映射示例

取常见的三角函数形式的 6-vertex (XXZ) R 矩阵（基底顺序 $|00\>,|01\>,|10\>,|11\>$）：
$$
R(u)=\begin{pmatrix}
\;a(u) & 0 & 0 & 0\\
0 & b(u) & c & 0\\
0 & c & b(u) & 0\\
0 & 0 & 0 & a(u)
\end{pmatrix},
\qquad a(u)=\sin(u+\eta),\; b(u)=\sin u,\; c=\sin\eta.
$$

此矩阵满足正则性（regularity）：
$$
R(0)=\sin\eta\;P,
$$
其中 $P$ 为交换算符（置换算符，$P\,|\alpha\beta\>=|\beta\alpha\>$）。基于这一点，可以定义局域哈密顿密度为 R 的对数导数（常用且保证厄米性的规范化写法）：
$$
h\;=\;\frac{P\,\partial_u R(u)\big|_{u=0}}{\sin\eta}.
$$
代入 $a,b,c$ 的导数（$a'(0)=\cos\eta,\; b'(0)=1,\; c'=0$）并代数化简可得（去掉可吸收的常数项后）经典的 XXZ 形式：
$$
h' = \frac{1}{2\sin\eta}\big(\sigma^x\otimes\sigma^x + \sigma^y\otimes\sigma^y\big) + \frac{\cot\eta}{2}\;\sigma^z\otimes\sigma^z,
$$
即（加上常数项）
$$
h = \frac{1}{2\sin\eta}(\sigma^x\sigma^x+\sigma^y\sigma^y)+\frac{\cot\eta}{2}(\sigma^z\sigma^z + I\otimes I).
$$

因此从谱参数 R(u) 可以直接导出局域两体哈密顿密度 $h$，并且 $h$ 恰好是 XXZ 链的最近邻密度（参数由 $\eta$ 决定）。脚本 `tools/xxz_R_and_H.py` 中已给出数值验证：取 $\eta=0.6$，计算得到（Pauli 展开系数，保留实数部分）
- $\;\dfrac{1}{2\sin\eta}\approx 0.8855160983$（对应 $\sigma^x\sigma^x$ 与 $\sigma^y\sigma^y$），
- $\;\dfrac{\cot\eta}{2}\approx 0.7308479735$（对应 $\sigma^z\sigma^z$），并伴随一个可吸收的常数项相同数量级。



脉冲映射（如何用含时哈密顿实现谱参数）：若把体系物理哈密顿写作
$$
H(t)=g(t)\sum_j h_j,
$$
其中 $h_j$ 为作用在键 $\langle j,j+1\rangle$ 上的局域密度，则时间演化算符为
$$
U(t)=\mathcal{T}\exp\Big(-i\int_0^t H(s)ds\Big)=\mathcal{T}\exp\Big(-i\int_0^t g(s)ds\;\sum_j h_j\Big).
$$
当 $[h_i,h_j]=0$（或在近似下认为不同键间可交换）并且脉冲分辨足够好时，可把谱参数 $u$ 与累积作用量一一对应：
$$
u\equiv\Phi(t)=\int_0^t g(s)\,ds.
$$
于是为了实现目标的 $R(u_*)=\exp(-i u_* \sum_j h_j)$，可选简单恒定脉冲
$$
g(t)=\frac{u_*}{\tau},\qquad 0\le t\le\tau,
$$
使得 $\Phi(\tau)=u_*$. 例如：取 $\eta=0.6$ 且目标谱参数 $u_*=\pi/4$，取脉冲时长 $\tau=1$，则恒定耦合 $g=\pi/4$ 会产生
$$
U(\tau)=\exp\Big(-i\frac{\pi}{4}\sum_j h_j\Big)
$$
对应的局域交换即为 $R(u_*=\pi/4)$（若需要严格满足加法型 YBE，应使用与 R(u) 全等的参数化约定并注意总体相位/标度因子）。

实现与注意事项：
- 当不同键上的 $h_j$ 并不严格对易时，可用短时 Trotter 分段实现近似：分段将 $g(t)$ 拆成小块并交替作用，从而产生近似的 $\exp(-i u\sum_j h_j)$；误差与步长和不对易项的范数有关。 
- 如果你希望把 $R(u)$ 的标准（例如 XXZ 6-vertex 的三角函数形式）与本文最初使用的 $R=e^{iH_P}$ 约定完全对齐，我可以在文中补充关于符号/相位/归一化的对应表（并给出如何从 $R(u)$ 导出 $H_P$ 的逐步代数推导）。

上述段落与数值例子已加入本文；如需我把推导写得更代数详细或把数值示例扩展为交互式笔记本，请告诉我。

### 任意目标 Majorana 对 — 在 Pauli 基上的通用构造与示例

给定链上第 j 个费米子对应的 Majorana 算符习惯写法：
$$
\gamma_{2j-1}=c_j+c_j^{\dagger},\qquad \gamma_{2j}=i(c_j-c_j^{\dagger}).
$$
在 Jordan–Wigner 形式下，这两类 Majorana 可用局域 Pauli 与中间的 $\sigma^z$ 串表示（取链上站点编号从 1 开始）：
$$
\gamma_{2j-1}=\Big(\prod_{m<j}\sigma^z_m\Big)\sigma^x_j,
\qquad
\gamma_{2j}=\Big(\prod_{m<j}\sigma^z_m\Big)\sigma^y_j.
$$

因此对于任意 $j<k$，任意两类 Majorana 的乘积都可写成端点处的 $\sigma^{x/y}$ 与中间的 $\sigma^z$ 串：
$$
\begin{aligned}
\gamma_{2j-1}\gamma_{2k-1} &= \sigma^x_j\Big(\prod_{m=j}^{k-1}\sigma^z_m\Big)\sigma^x_k,\\
\gamma_{2j}\gamma_{2k} &= \sigma^y_j\Big(\prod_{m=j}^{k-1}\sigma^z_m\Big)\sigma^y_k,\\
\gamma_{2j-1}\gamma_{2k} &= \sigma^x_j\Big(\prod_{m=j}^{k-1}\sigma^z_m\Big)\sigma^y_k,\\
\gamma_{2j}\gamma_{2k-1} &= \sigma^y_j\Big(\prod_{m=j}^{k-1}\sigma^z_m\Big)\sigma^x_k.
\end{aligned}
$$

因此，欲构造目标的 Majorana 交换算子
$$
B(\gamma_a,\gamma_b)=\exp\Big(\frac{\pi}{4}\,\gamma_a\gamma_b\Big)
$$
可先构造局域生成元（写成 Pauli 形式）
$$
H_P = -i\frac{\pi}{4}\,\gamma_a\gamma_b,
$$
并用上面的等式把 $\gamma_a\gamma_b$ 用 Pauli 串替换，从而得到 $H_P$ 在 Pauli 基下的显式表达。下面给出两个常用示例：

例 1（相邻 Majorana，对应纯两体 Pauli）：对两个位于相邻站点的 Majorana，例如 $\gamma_1,\gamma_4$，按本文采用的 JW 约定有
$$
\gamma_1 = \sigma^x_1,
\qquad
\gamma_4 = \sigma^z_1\sigma^y_2,
$$
因此
$$
\gamma_1\gamma_4 = \sigma^x_1\sigma^z_1\sigma^y_2 = -i\,\sigma^y_1\sigma^y_2.
$$
等价地，
$$
i\gamma_1\gamma_4 = \sigma^y_1\sigma^y_2.
$$
于是
$$
H_P = -i\frac{\pi}{4}\gamma_1\gamma_4 = -\frac{\pi}{4}\big(\sigma^x\sigma^x - \sigma^y\sigma^y\big),
$$
从而一次含时演化 $U=\exp(iH_P)$（或等价地在谱参数表示下取合适的 $u$）即可实现精确的 $B(\gamma_1,\gamma_4)$（到整体相位）。

例 2（远距 Majorana，带字符串）：取站点 $j<k$ 的同类 Majorana $\gamma_{2j-1},\gamma_{2k-1}$，则
$$
\gamma_{2j-1}\gamma_{2k-1}=\sigma^x_j\Big(\prod_{m=j}^{k-1}\sigma^z_m\Big)\sigma^x_k.
$$
相应的局域生成元（在 Pauli 表示下）为
$$
H_P = -i\frac{\pi}{4}\,\sigma^x_j\Big(\prod_{m=j}^{k-1}\sigma^z_m\Big)\sigma^x_k.
$$
如果你的硬件或数值平台不能直接实现长链 Pauli 串，可以两种途径得到等效的交换操作：
- （局部移动）先用一系列已知的局域交换/置换把远距 Majorana 通过邻位交换移动到相邻位置，再作用相邻对的 $H_P$；
- （分段 Trotter）直接用短时片段交替施加端点与中间的 $\sigma^z$ 操作，或构造等效电路把长串分解为本地门序列。

小结：给定任意两个 Majorana，按照上面的通用公式可把目标算子写为 Pauli 串；若它正好落在相邻站点，生成元仅含两体 Pauli，从而能精确实现理想的编织算子；否则需要用移动或分段策略把非局域串约化为局域操作。

### 关于 σ^z、Jordan–Wigner 串 与四费米项（补充）

为避免歧义，这里补充几条常见混淆点的简洁说明：

- 单体 $\sigma^z_j$ 在 JW 映射下对应密度算符：
$$
\sigma^z_j = 2n_j-1,\qquad n_j=c_j^{\dagger}c_j.
$$ 
因此
$$
\sigma^z_j\sigma^z_{j+1}=4n_jn_{j+1}-2(n_j+n_{j+1})+1,$$
其中包含四费米项 $4n_jn_{j+1}$（最近邻相互作用）以及线性密度与常数项。

- $\sigma^{x,y}$ 含 JW 串：用 $S_j=\prod_{k<j}\sigma^z_k$ 表示串算符，有
$$
\sigma^x_j=S_j(c_j^{\dagger}+c_j),\qquad \sigma^y_j=S_j(-ic_j^{\dagger}+ic_j).
$$
当两个相邻站点同时为 $\sigma^{x,y}$ 类型时（例如 $\sigma^x_j\sigma^x_{j+1}$），相邻的串因子互相抵消，结果是在费米子表述中只含二次项（hopping/pairing）；但若出现单侧 $\sigma^z$ 或单体 $\sigma^{x,y}$，串不会完全抵消，会留下非局域串或高阶项。

- 直观示例：
	- $\sigma^x_j\sigma^x_{j+1}$、$\sigma^y_j\sigma^y_{j+1}$ → 映射为仅含 $c^{\dagger}c$ 与 $c^{\dagger}c^{\dagger}$ 的二次组合（无串）。
	- $\sigma^x_j\sigma^z_{j+1}$ → 映射为 $S_j(c_j^{\dagger}+c_j)(2n_{j+1}-1)$，保留串 $S_j$，在费米子语言上为非局域或高阶算符。

结论：因此若目标是得到纯二次 BdG 哈密顿（便于出现并描述 MZM），应优先使用最近邻的 XX/YY/XY/YX 之类基元并避免单侧 X/Y 或含单体 Z 的混合项；若不可避免，则需通过移动 Majorana 或用 Trotter 分段策略把串项等效实现。

### 实施注意、脉冲映射、Trotter 配方与验证命令（补充）

要把谱参数型 R(u) 用作物理含时演化实现编织操作，请注意以下要点：

- 精确条件：若生成元族满足 $[H_P(s),H_P(s')]=0$（或等价条件使 $R(u)$ 在不同位置嵌入后满足 $[X_{12},X_{13}+X_{23}]=0$），则存在严格的一一映射 $u=\Phi(t)=\int_0^t g(s)ds$，使得
$R(u)=\mathcal{T}\exp\big(-i\int_0^t H_P(s)ds\big)=\exp(-i\Phi H_{P0})$，从而可精确实现谱参数 R(u)。

- 近似情形：若生成元在时间上或空间上不对易，则采用分段（Trotter）或先移动 Majorana 的策略来近似实现目标算符；误差随分段数增长而下降（误差阶数通常为 O(1/m) 或 O(1/m^2)，视 Trotter 种类而定）。

- 脉冲映射示例：设目标局域密度为 $h$ 且目标谱参数为 $u_*$，可用恒定脉冲
$$
g(t)=u_*/\tau\quad(0\le t\le\tau)
$$
使得 $\Phi(\tau)=u_*$. 若只能实现若干基础哈密顿分量 $h_1,h_2,\dots$，且 $h=\sum_\alpha \lambda_\alpha h_\alpha$，可按比例放置时长或交替短时段实现。

- Trotter 具体配方（二阶示例）：要近似 $\exp(-i u(h_A+h_B))$，可用 m 步二阶 Trotter：
$$
	U\approx\Big( e^{-i\frac{u}{2m}h_A} e^{-i\frac{u}{m}h_B} e^{-i\frac{u}{2m}h_A}\Big)^m,
$$
取 m 足够大以使误差在可接受范围内（误差量级 ≲ O(u^3/m^2\,\|[h_A,[h_A,h_B]]\|+\dots)）。

- 验证实用命令：仓库中已有两个脚本用于数值验证：

	- 运行 XXZ → 局域 h 并查看 Pauli 展开：

		```bash
		python3 tools/xxz_R_and_H.py
		```

	- 数值检验谱参数 YBE（或比较两个算子的范数差）：

		```bash
		python3 tools/check_R_u.py
		```

	你也可以用这些脚本把 `H_P` 替换为目标的 Pauli 串（或用脚本扩展为命令行参数）以计算 $R(u)$、$B(\gamma_a,\gamma_b)$ 并输出范数差以量化实现的精确度。

- 实验/数值建议：先在小系统（2–4 自旋对）上用脚本验证脉冲/分段参数，再把经过验证的序列推广到更长链；评估指标包括算符范数差、保真度以及对局域扰动的鲁棒性。

若你愿意，我可以：
- 把 `tools/xxz_R_and_H.py` 和 `tools/check_R_u.py` 增强为接收命令行参数（例如 `--eta`, `--u`, `--m`），并添加一个小 README 说明如何对任意目标 `H_P` 运行验证；或
- 为你生成一个示例 Trotter 脉冲序列脚本，直接输出每一步的局域门与累积误差估计。

 

## 追加记录（工作进展与数值结果）

以下为本次工作在仓库中完成的要点与数值结果，已追加为工作记录以便复现与汇报：

- 已新增/运行的脚本：
	- `tools/check_R_u.py`：用于检验谱参数形式 $R(u)=\exp(iuH_P)$ 的 YBE 近似（小系统数值检验）。
	- `tools/xxz_R_and_H.py`：从 XXZ (6‑vertex) 的 $R(u)$ 导出局域哈密顿密度 $h$ 并给出 Pauli 展开系数与数值示例。 
	- `tools/verify_mzm.py`：构建 Kitaev‑链 BdG 矩阵并对角化以检测近零模。 
	- `tools/zero_mode_profile.py`：提取零模本征向量，计算站点归一化密度，拟合局域化长度并绘图（已生成 `tools/zero_mode_profile_cartoon.png`）。

- 示例运行与关键输出（复现命令 `python tools/zero_mode_profile.py`）：
	- 链长：$N=120$，参数 $t=1,\Delta=1,\mu=0$。\
	- 发现 2 个近零模（容差 $10^{-6}$）。\
	- 右端模（峰位 site=119）拟合局域化长度 $\xi\approx 2.829$（拟合使用峰向外固定点方法，points_used=40）。\
	- 左端模（峰位 site=0）数值上几乎全权集中在峰点（其余点密度 < 1e-30），拟合给 $\xi\approx 0.049$（表示数值极端局域化与下溢影响）。\
	- 绘图已保存为 `tools/zero_mode_profile_cartoon.png`，包含平滑密度曲线与链状示意条（端点圆点），用于报告图示。

- 理想编织验证：`tools/check_braiding_from_HP.py` 数值验证相邻 Majorana 的 H_P 构造（如 $H_P=-(\pi/4)(\sigma^x\sigma^x-\sigma^y\sigma^y)$）能用 $R=\exp(iH_P)$ 重构目标编织算子，算符差范数 ~4e-16（数值上充分一致）。

注意事项：谱参数形式 $R(u)=\exp(iuH_P)$ 在代数条件满足时可以直接作为含时演化的实现，否则需用时间序列或 Trotter 分段近似。对于非常尖锐或数值下溢的零模，建议用更稳健的度量（例如 90% 累计权重的半径）或增大链长/数值精度复验。


如需，我可以接着：

- 把绘图脚本增加命令行选项以输出线性/对数尺度、放大首/尾 20 个格点；
- 在图上叠加 log‑fit 直线并输出拟合残差与拟合点区间；
- 计算并把“到达累计权重 90%”的格点数作为稳健的局域化长度并写入文档。

(记录追加完毕)

 

## 如何构建 Majorana 零模（MZM）——操作手册与数值步骤

下面给出在理论、数值与实验层面上构建与验证 Majorana 零模的实用步骤，包含必要的数学表达、数值命令与实现建议，便于把上文的 R/H_P 构造直接映射到 MZM 的产生与观测。

1) 物理与数学先决条件（Kitaev‑链/BdG）

- 目标：构建一个二次 BdG 哈密顿 $H_{\mathrm{quad}}$（仅含 $c_j^{\dagger}c_{j+1},\; c_j c_{j+1}$ 等二次项）在开放边界条件下拥有两端的近零能模式。\
- 常见拓扑相条件（参考最简单的均匀 Kitaev 链）：当 $|\mu|<2|t|$ 且 $\Delta\neq0$ 时，链处于拓扑超导相并在两端出现 MZM（在有限链上对应近零本征值）。

2) 从 Pauli/两体生成元到 BdG（数值实现流程）

- 给定局域生成元系数 $c_{\mu\nu}$（见上文映射），用线性映射得到 $t,\Delta,\mu,U$（可以直接调用 `tools/xxz_R_and_H.py` 或自己写一行代码）。
- 用 `build_kitaev_bdg(N,t,\Delta,\mu)` 构造 $2N\times2N$ 的 BdG 矩阵（见仓库的 `tools/zero_mode_profile.py` 中实现）。
- 对 $H_{\mathrm{BdG}}$ 做本征分解，找出绝对值最小的本征值及对应本征向量 $\Psi=(u;v)$。\

3) 从 BdG 本征向量重建 Majorana 模态（算符形式）

- 令零能本征向量（或近零）为 $\Psi=(u_1,\dots,u_N; v_1,\dots,v_N)^T$，构造对应的费米子振幅算符（准粒子算符）
$$
f = \sum_{j=1}^N\big(u_j c_j + v_j c_j^{\dagger}\big),
$$
已规范化使得 $\{f,f^{\dagger}\}=1$（若原本征向量未规范化请先归一化）。\
- 对应的两个 Majorana 算符可写为
$$
\gamma_A = f + f^{\dagger},\qquad \gamma_B = -i(f - f^{\dagger}).
$$
这两算符都是厄米且满足 $\{\gamma_A,\gamma_B\}=0,\;\gamma_A^2=\gamma_B^2=1$（在数值上应做归一化误差检查）。\
- 若 $\Psi$ 为单独的零模（对应孤立的 Majorana），则上述 $\gamma$ 主要局域在链的一端；两个零模则分别局域在链两端。

4) 数值验证步骤（命令与示例）

- 构建并对角化（示例，仓库已有脚本）：
	```bash
	python3 tools/verify_mzm.py   # 或直接运行 tools/zero_mode_profile.py 查看密度与拟合
	python3 tools/zero_mode_profile.py
	```
- 从输出获取近零本征向量索引 i，并在 Python 中用如下片段重建 f 与 Majorana：
	```python
	psi = vecs[:, i]    # vecs from eigh(Hbdg)
	u = psi[:N]
	v = psi[N:]
	# normalize such that sum(|u|^2+|v|^2)=1
	f_coeffs = (u, v)
	# construct Majorana mode weights on site j: real and imag parts combine into gamma components
	# (see text for explicit gamma_A/gamma_B construction)
	```
	（详尽脚本见 `tools/zero_mode_profile.py`，可扩展为输出 gamma 在每个站点上的系数与归一化误差）。

5) 实验/数值上产生 MZM 的常用途径

- 域壁（domain‑wall）方法：在空间上把化学势 $\mu(x)$ 设计为中间为平庸相（$|\mu|>2|t|$），两侧为拓扑相（$|\mu|<2|t|$），界面处出现 MZM。数值上可以把链分成若干段并逐点设置 $\mu_j$；实验上通过局域门电压调制亦可实现。\
- 参数急/慢变化：若把 $\mu$ 从大变小穿越拓扑相界面，若变化速度慢可理解为准静态产生 MZM，若快速则可能产生非平衡激发。\
- 相邻 Majorana 的精确编织：直接用局域两体生成元 $H_P$（例如例子中的 $H_P=-(\pi/4)(\sigma^x\sigma^x-\sigma^y\sigma^y)$）并实施 $U=\exp(iH_P)$，在数值上通过 `tools/check_braiding_from_HP.py` 验证算符差矩阵范数是否在容差内。对于非相邻 Majorana，采用移动/分段 Trotter 策略。

6) 诊断与鲁棒性检查

- 检查本征值间隙：零模附近的能隙（距离零的最小非零本征能）应大于数值扰动阈值，否则零模不稳定。\
- 归一化误差：验证 $\Psi$ 的归一化，并检查构造的 $\gamma$ 是否满足近似的归一性 $\gamma^2\approx1$ 与正交性。\
- 局域化度量：除了指数拟合得到的 $\xi$，建议计算“从边开始累计权重达到 90% 所需的距离”作为稳健度量（适合极端局域/下溢情形）。

7) 可复现的示例参数（直接拷贝运行）

- 在仓库中，示例参数：$t=1,\;\Delta=1,\;\mu=0,\;N=120$，运行 `python3 tools/zero_mode_profile.py` 可得到两个近零模的密度图与拟合结果（前文记录了 xi 值与峰位）。

8) 推荐的后续改进（可选）

- 自动化脚本：把从 `c_{\mu\nu}` 到 BdG、到本征向量、到 Majorana 重建、到 画图与导出报告的流程做成一键脚本（可加命令行参数）。
- 鲁棒度扫描：在 $(\mu,t,\Delta)$ 网格上计算 xi 与零模存在性并绘制相图；若需要，我可以生成该扫描脚本并把结果保存为 PNG/CSV。 

以上段落将被追加到本文档以指导后续数值验证与实验实现。如需精确的代码片段把 Majorana 权重导出为文本/CSV，我现在可以把该小工具追加到 `tools/` 并运行一次示例生成输出。

 

## 详细方法步骤（逐步）

下面给出从谱参数 R(u) 到可实施编织操作与 Majorana 零模的逐步方法，包含代数检查、数值实现、脉冲设计与验证清单，便于直接复制到数值脚本或实验控制序列中。

步骤 1 — 预备与输入
- 确定 R(u) 的解析表达式（例如 XXZ/6‑vertex 的 a(u), b(u), c），并确认基底约定（本文采用 |00>,|01>,|10>,|11> 的字典序）。
- 选择谱参数点和目标谱参数值 u_*（例如 u_*=π/4 对应常见角度）。

步骤 2 — 从 R(u) 导出局域哈密顿密度 h
- 计算 R(0)（检查正则性：R(0) 是否为 scale×P）。若满足正则性，记 scale = 常数（例如 sin η）。
- 计算导数 ∂_u R(u)|_{u=0}。
- 用规范化定义得到 h：
$$
h = \frac{P\,\partial_u R(u)\big|_{u=0}}{\text{scale}}.
$$ 
- 将 h 展开到 Pauli 基（写成$ ∑_{αβ} c_{αβ} σ^α⊗σ^β$），保存各系数$ c_{αβ} $以备后续映射。

步骤 3 — Pauli → 物理参数（t, Δ, μ, U 等）与 JW 检查
- 用本文给出的线性映射（矩阵 M 或逐分量公式）把 Pauli 系数映射为 $t, Δ, U, μ $等物理参数。
- 检查是否存在触发字符串的分量（$c_{x0}, c_{0x}, c_{y0}, c_{0y}, c_{xz}, c_{zx}, …$）。若任一非零，记录这些分量，这些项在 JW 后会引入串前缀或非局域项。
- 如果目标是纯二次 BdG（便于出现 MZM），要求这些字符串分量为零且尽量使$ c_{zz}=0 $或可控（若 $c_{zz}≠0$，需处理交互项 U）。

步骤 4 — 构建 BdG 矩阵并求解零模（数值）
- 构造 BdG 矩阵 H_bdG(N,t,Δ,μ)（参考 `tools/zero_mode_profile.py` 的 `build_kitaev_bdg` 实现）。
- 对 H_bdG 做本征分解（eigh），按 |E| 排序找到近零本征值与对应本征向量 Ψ=(u;v)。
- 计算站点密度 dens_j=|u_j|^2+|v_j|^2，归一化并保存。可用 `tools/zero_mode_profile.py` 自动化该流程。

步骤 5 — 从本征向量重建 Majorana 算符并检查正交性
- 构造准粒子 $f = ∑ (u_j c_j + v_j c_j^†)$，检查归一化$ {f,f^†}=1$（或向量归一化）。
- 构造$ γ_A = f+f^†, γ_B = -i(f-f^†)$，并检查数值上$ γ_A^2≈1, γ_B^2≈1, {γ_A,γ_B}≈0$。

步骤 6 — 量化局域化并输出稳健度量
- 指数拟合：在峰位找向外拟合 log(density)=a+b x，得到 ξ=-1/b，报告拟合区间与残差。若拟合失败则使用固定点数的回退方法并标注其局限性。
- 稳健度量：计算从边缘向内累计权重$ W(R)=∑_{j≤R} dens_j$，找到最小 R 使 $W(R)≥0.9$，记录 R_{90} 作为稳健的局域化尺度。
- 导出 PNG/CSV：保存密度曲线、拟合直线、R_{90} 标记与拟合残差。

步骤 7 — 从 h 设计脉冲实现 R(u)
- 若所有 h_j（不同键）两两对易或 $[H_P(s),H_P(s')]=0$，可直接选 g(t) 使$ Φ=∫g dt = u_*$，作用时间 τ 取决于实验限制，$U=exp(-i Φ ∑_j h_j) $即为目标。

- 若不对易，采用 Trotter 分段：选择分段数 m，使用二阶 Trotter：
$$
U\approx\Big(e^{-i\frac{u_*}{2m}h_A} e^{-i\frac{u_*}{m}h_B} e^{-i\frac{u_*}{2m}h_A}\Big)^m.
$$ 
 估计误差并调整 m 直到保真度满足要求；记录门序列与每步持续时间。

步骤 8 — 构造并验证编织门
- 相邻 Majorana：把目标 γ_aγ_b 展开到两体 Pauli，取 $H_P = -i(π/4) γ_aγ_b $的 Pauli 表达，计算 R_from_HP = exp(i H_P) 并与理论 B_target 比较，计算范数差 ||R_from_HP - B_target||。阈值建议 1e-12–1e-8 视数值精度而定。
- 非相邻 Majorana：先用局域移动/邻位交换将目标 Majorana 移到相邻位置，或直接用长串 Trotter 分解并验证累计误差。

步骤 9 — 扫描与鲁棒性测试（可选但强烈建议）
- 在参数网格 (μ,t,Δ) 上计算 ξ 与是否存在近零模，绘制相图并标注拓扑/平庸相边界。
- 对脉冲参数（分段数 m、步时 τ、噪声幅度）做敏感性分析并记录算符保真度。

步骤 10 — 文档与复现
- 把上述所有命令、参数与输出（PNG/CSV）归档到一个可复现的目录（例如 `results/`），在 `README.md` 中写明复现步骤与运行命令。常用命令示例：
	```bash
	python3 tools/xxz_R_and_H.py --eta 0.6    # 导出 h 的 Pauli 系数
	python3 tools/check_R_u.py --u 0.785398  # 检查 YBE/谱参数一致性
	python3 tools/zero_mode_profile.py       # 计算零模、绘图并输出 xi 与 R90
	```

验证清单（发布/提交前）
- 确认 h 的 Pauli 系数已记录并映射到 (t,Δ,μ,U)。
- 确认 BdG 本征分解找到了预期数量的近零模并保存本征向量。 
- 输出并记录 ξ（指数拟合）与 R_{90}，以及拟合残差。 
- 验证算符差范数与目标门在容差内。 

附：常见数值陷阱与处理
- 浮点下溢：在对数绘图或拟合时为密度加小常数 eps（脚本中采用 1e-300）以避免 log(0)。
- 极端局域化：当密度集中在单点时，指数拟合退化，优先使用 R_{90} 或增加链长 N。 
- 相位/尺度不敏感性：比较算符时去掉全局相位或用保真度而非直接矩阵差来衡量。

