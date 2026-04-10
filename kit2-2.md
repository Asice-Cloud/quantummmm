### 从 1D R→Kitaev 链 到 2D Z2 规范推广（具体例子）

**本节目标**  
- 先回顾 1D 中从局域
    $$R_{i,i+1}=aI+b\,\sigma^x_i\sigma^x_{i+1}+c\,\sigma^y_i\sigma^y_{i+1}+d\,\sigma^z_i\sigma^z_{i+1}$$
    到 Kitaev 链参数 $(t,\Delta,\mu)$ 的映射；
- 再把同一个两体 $R$ 放到二维图的最近邻边 $(i,j)$ 上，用 Majorana+Z2 规范场写成链路变量 $u_{ij}$ 乘端点 Majorana 的二次项，构造 2D 拓扑超导模型；
- 最后说明在这个框架下，“远程算符”出现为：
    1) 由最近邻操作串成的复合幺正门 $R^{(\text{gate})}_{ij}$，以及  
    2) 在 2D+Z2 语言中带有路径串 $u_\gamma$ 的 Majorana 长程字符串算符 $B_{ij}(\gamma)=u_\gamma(i c_i c_j)$，与编织路径一一对应。

这一节把前面 1D 中的
$$
R_{i,i+1}=aI+b\,\sigma^x_i\sigma^x_{i+1}+c\,\sigma^y_i\sigma^y_{i+1}+d\,\sigma^z_i\sigma^z_{i+1}
$$
映射到 $t,\Delta,\mu$ 的结果，与本节的 Majorana+Z2 规范场框架结合，给出一个在二维方格上带 Z2 链变量的具体模型和动量空间形式。目的是说明：

- 每条键上的 R 决定这一键上的 $t_{ij},\Delta_{ij},\mu$；
- Z2 规范场通过链变量 $u_{ij}=\pm1$ 乘在这些耦合前面；
- 在简单的静态规范背景（例如 $u_{ij}=1$ 或均匀 $\pi$ 通量）下，可以写出 BdG 形式并分析拓扑性质。

#### 1. 1D 映射回顾（局域识别）

在 1D 链上，前文已得出：
$$
R_{i,i+1}=aI+b\,\sigma^x_i\sigma^x_{i+1}+c\,\sigma^y_i\sigma^y_{i+1}+d\,\sigma^z_i\sigma^z_{i+1}
$$
在 Jordan–Wigner 映射下对应于单键上的费米子算符
$$
\begin{aligned}
R_{i,i+1}
&=(b+c)\,(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i)
 +(b-c)\,(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i)\\
&\quad +d\,(4n_in_{i+1}-2(n_i+n_{i+1})+1)+aI.
\end{aligned}
$$

因此在“只看二次项”的层面上（忽略 $4n_in_{i+1}$ 或把其平均场重整化到 $\mu$ 中），可以局域识别
$$
t\propto (b+c),\qquad \Delta\propto (b-c),\qquad \mu\ \text{由 }a,d\text{ 的线性项给出}.
$$

#### 1.5 从相邻 $R_{i,i+1}$ 到图上一条边的 $R_{ij}$

**算符嵌入 vs. 模型的图结构**  
起点是 R_to_Kitaev 里定义在 $V\otimes V$ 上的
$$
R=aI + b\,\sigma^x\otimes\sigma^x + c\,\sigma^y\otimes\sigma^y + d\,\sigma^z\otimes\sigma^z.
$$
在 1D 自旋链上，底层“图”是一条线，只有边 $(i,i+1)$。把 $R$ 嵌入到链上时，
$$
R_{i,i+1}=aI + b\,\sigma^x_i\sigma^x_{i+1}+c\,\sigma^y_i\sigma^y_{i+1}+d\,\sigma^z_i\sigma^z_{i+1}
$$
表示“在这条最近邻边 $(i,i+1)$ 上作用 $R$”，这就是哈密顿量中真正出现的局域两体项。形式上当然可以在 Hilbert 空间里对任意 $i<j$ 定义
$$
R_{ij}=aI + b\,\sigma^x_i\sigma^x_j + c\,\sigma^y_i\sigma^y_j + d\,\sigma^z_i\sigma^z_j,
$$
作为“长程两体算符”，但在 1D Kitaev 链的物理模型中，我们的哈密顿量 $H=\sum_i R_{i,i+1}$ 只在图的最近邻边上取和，并不真正包含这些远程 $R_{ij}$。

**用最近邻操作实现等效远程作用（门序列的意义）**  
若只在“操作/量子门”的层面希望在远距 $(i,j)$ 上实现同样形式的两体变换，而底层硬件或模型只允许最近邻 $R_{k,k+1}$，可以利用交换算符序列。设 $S_{k,k+1}$ 是把第 $k$ 和 $k+1$ 个格点对调的某个局域算符（例如由一串局域 R 或 SWAP 门构成），定义
$$
S_{i\to j}=S_{j-1,j}\,S_{j-2,j-1}\cdots S_{i,i+1},
$$
它把第 $i$ 个格点搬运到第 $j$ 个位置。则门序列
$$
R^{\text{(gate)}}_{ij}=S_{i\to j}\,R_{j-1,j}\,S_{i\to j}^{\dagger}
$$
在算符代数层面等价于一个只作用在 $(i,j)$ 上的两体变换，但这是由若干最近邻操作拼出来的**复合幺正门**，而不是哈密顿量中的单个局域密度。

**费米/JW 视角与后文 Z2 规范的联系**  
若在 1D 中直接把上面的抽象 $R_{ij}$ 代入 Jordan–Wigner，则当 $i,j$ 不相邻时，中间会出现从 $i$ 到 $j$ 的长串 $e^{i\pi\sum_{k=i}^{j-1}n_k}$，等价于“沿路径的 Z2 链变量乘积”。这正是后文在二维中引入 Z2 链变量 $u_{ij}$、并把远程作用写成“路径上 $u$ 的乘积 × 端点的 Majorana 二次项”的动机：在 1D，可以把这类长程 $R_{ij}$ 看成沿链上一条固定路径的 JW 串；在 2D，则将这条路径推广为一般格上的路径 $\gamma:i\rightsquigarrow j$，并用 $\prod_{(kl)\in\gamma}u_{kl}$ 来封装路径依赖的符号与拓扑信息。

下面第 2 节起，我们就在二维方格/蜂窝格上显式构造带 Z2 链变量的版本：此时 $R_{ij}$ 只放在二维图的最近邻边 $(i,j)$ 上，用同一组 $(a,b,c,d)$ 或其导出的 $(t,\Delta,\mu)$ 在 2D 中分析拓扑 Majorana 模式。

#### 2. 在二维方格上加入 Z2 链变量

取一个二维方格晶格，格点用 $i,j$ 标记，最近邻键分为横向 $x$ 键和纵向 $y$ 键：

- 对每条有向键 $\langle ij\rangle$（例如从左到右、从下到上），放一个 Z2 链变量
    $$u_{ij}=\pm1,\qquad u_{ji}=u_{ij}.$$
- 定义 Z2 规范变换
    $$c_i\to s_i c_i,\qquad u_{ij}\to s_i u_{ij} s_j,\qquad s_i=\pm1,$$
    从而任何形如 $u_{ij}c_i^{\dagger}c_j$、$u_{ij}c_ic_j$ 的项都对规范变换不变。

现在把每条键上的 R 升级为带 Z2 链变量的版本。为了让规范场只耦合到“XY 方向”的跳跃/配对，我们选一个简单但常用的约定：
$$
R_{ij}(u_{ij}) = aI + u_{ij}\big(b\,\sigma^x_i\sigma^x_j + c\,\sigma^y_i\sigma^y_j\big) + d\,\sigma^z_i\sigma^z_j.
$$

在 Majorana+Z2 表示中（本文件前一节），各向异性键项可写为
$$
\sigma^a_i\sigma^a_j = u^{(a)}_{ij}\,(ic_i c_j),\qquad u^{(a)}_{ij}=-i b^a_i b^a_j,
$$
其中 $u^{(a)}_{ij}$ 本身就是一个 Z2 链变量。这里为了不引入过多符号，可以把“几何上引入的” $u_{ij}$ 和来自 Majorana 分解的 $u^{(a)}_{ij}$ 合并为一个总的 Z2 链变量
$$
ilde u_{ij}=u_{ij}\,u^{(a)}_{ij}=\pm1,
$$
并仍然记作 $u_{ij}$。这样在费米子语言中，横向/纵向最近邻耦合统一写成
$$
H_f = \sum_{\langle ij\rangle}\Big[ u_{ij}\,t_{ij}(c_i^{\dagger}c_j+h.c.)+u_{ij}\,\Delta_{ij}(c_ic_j+h.c.)\Big]
 -\mu\sum_i(n_i-\tfrac12)+\cdots.
$$
其中 $t_{ij},\Delta_{ij},\mu$ 由该键上选用的 $R_{ij}(a,b,c,d)$ 决定（例如对所有横向键用 $(a_x,b_x,c_x,d_x)$，纵向键用 $(a_y,b_y,c_y,d_y)$）。

规范场自身的哈密顿量写为标准的 Z2 规范项：
$$
H_g = -K\sum_p W_p,\qquad W_p=\prod_{\langle ij\rangle\in p} u_{ij},
$$
其中 $p$ 遍历所有 plaquette，$W_p=\pm1$ 是磁通（或 vison）变量。若要对电荷自由度也施加能量代价，可以再加星算符 $G_i$ 项（toric‑code 结构），本节暂不展开。

总哈密顿量即
$$
H = H_f[c,u]+H_g[u].
$$

当把 $u_{ij}$ 视作静态背景时（例如先固定某个 $\{u\}$ 扇区，再对 Majorana/费米子对角化），问题简化为“在某个给定 Z2 背景下的拓扑超导”。对每个给定的 $u$ 配置，可以像 1D 一样判定是否存在 Majorana 零模（边界、涡核处）并计算拓扑不变量；再在不同 $u$ 扇区之间比较总能量得到真实基态。

#### 3. 具体例子：方格上各向异性的 p‑wave 超导

为了给出一个完全具体、可计算的例子，考虑最简单的“各向异性 p‑wave” 模型：

- 横向键（$x$ 方向）上的 R 参数取
    $$a_x=0,\quad d_x=0,\quad b_x=1,\quad c_x=0,$$
    则
    $$t_x\propto b_x+c_x=1,\qquad \Delta_x\propto b_x-c_x=1.$$
- 纵向键（$y$ 方向）上的 R 参数取
    $$a_y=0,\quad d_y=0,\quad b_y=1,\quad c_y=0,$$
    同样得到
    $$t_y\propto1,\qquad \Delta_y\propto1.$$

在这一最简单情形中，我们先不追求精确的 $p_x+ip_y$ 结构，而是展示 Z2 规范的耦合方式和 BdG 形式。取所有链变量 $u_{ij}=+1$（无通量背景），则费米子哈密顿量为
$$
H_f = -t\sum_{\langle ij\rangle_x}(c_i^{\dagger}c_j+h.c.)
            -t\sum_{\langle ij\rangle_y}(c_i^{\dagger}c_j+h.c.)
            +\Delta\sum_{\langle ij\rangle_x}(c_ic_j+h.c.)
            +\Delta\sum_{\langle ij\rangle_y}(c_ic_j+h.c.)
            -\mu\sum_i(n_i-\tfrac12),
$$
其中为了简洁写成标量形式，取 $t_x=t_y=t,\;\Delta_x=\Delta_y=\Delta$，$\mu$ 由 $d_x,d_y,a_x,a_y$ 决定（在上述选择中 $d_x=d_y=a_x=a_y=0$，所以可以单独指定一个化学势）。

对该模型做 Fourier 变换，取自旋量子数省略的 Nambu 自旋量
$$
\Psi_{\mathbf{k}} = (c_{\mathbf{k}}, c^{\dagger}_{-\mathbf{k}})^T,
$$
得到标准的 BdG 形式
$$
H_f = \sum_{\mathbf{k}} \Psi^{\dagger}_{\mathbf{k}}\,\mathcal{H}(\mathbf{k})\,\Psi_{\mathbf{k}} + \text{const.},
$$
其中
$$
\mathcal{H}(\mathbf{k}) = \xi(\mathbf{k})\,\tau_z + \Re\Delta(\mathbf{k})\,\tau_x - \Im\Delta(\mathbf{k})\,\tau_y.
$$

对于上面的最近邻模型，在“所有键的配对相同相位”的简单约定下，得到
$$
\xi(\mathbf{k}) = -2t(\cos k_x+\cos k_y) - \mu,
$$
$$
\Delta(\mathbf{k}) = 2\Delta(\cos k_x+\cos k_y),
$$
对应的是一个各向异性但仍然具有 $s+ d$‑样结构的配对函数。若希望得到真正的 $p_x$ 或 $p_x\pm ip_y$ 结构，可以：

- 对配对项在反向键上引入负号，使得 $\Delta_x$ 在 $+\hat x$ 与 $-\hat x$ 键上相反，从而
    $$\Delta(\mathbf{k})\propto2i\Delta_x\sin k_x,$$
    类似地对 $y$ 方向做同样处理得到 $\Delta(\mathbf{k})\propto2i(\Delta_x\sin k_x+\Delta_y\sin k_y)$，对应 $p_x+p_y$；
- 若再允许在 $x,y$ 方向的配对之间加入相对相位（例如 $\Delta_y=\pm i\Delta_x$），即可得到 $p_x\pm ip_y$ 结构。这一步严格说需要 U(1) 级别的相位而非单纯 Z2 的 $\pm1$；在 Z2 框架中，这个相位可以通过更大空间中的规范场或有效低能理论产生（类似蜂窝格中磁场诱导的三自旋项）。

本节的重点是：

1. **局域 R → $(t,\Delta,\mu)$ 的识别完全可以在二维逐键复用**：只要指定每种键方向一组 $(a,b,c,d)$，就能在方格/蜂窝/任意图上得到对应的 $t_{ij},\Delta_{ij},\mu$；
2. **Z2 规范场自然以链变量 $u_{ij}=\pm1$ 的形式乘在这些耦合前面**，规范变换 $c_i\to s_i c_i, u_{ij}\to s_i u_{ij} s_j$ 保证了规范不变性；
3. 在固定的 $\{u\}$ 背景下，问题退化为标准的 BdG 拓扑超导，可以像 1D 那样用谱、能隙、Chern 数等判据来分析 Majorana 边界模和涡核零模；不同 $\{u\}$ 扇区之间的能量比较则给出有效的 Z2 规范动力学和拓扑序。

具体数值实现上，可以：

- 在一个有限的二维方格上固定若干 $\{u\}$ 构型（例如无通量、单漩涡、周期性 $\pi$ 通量等），按上式构造 BdG 矩阵 $\mathcal{H}(\mathbf{k})$ 或其实空间版本；
- 求谱并检查是否存在局域在边界/漩涡处的零能本征态；
- 改变 $(a,b,c,d)$（即改变每条键上的 R）观察零模和带结构的变化，从而直接在“R→Z2 规范→拓扑 Majorana” 这条链上进行扫描和对比。

#### 4. 2D 编织与路径算符 $U_\gamma$

在二维格上，“编织”可以理解为：把携带 Majorana 零模的拓扑缺陷（例如 vison/涡核）沿某条路径 $\gamma$ 移动，绕过或绕圈另一缺陷。时空里这对应它们的世界线有非平凡的绕缠数；在低能零模子空间中，这个过程实现一个幺正算符 $U_\gamma$。

- **局域一步：最近邻上的交换门**  
    若某一步是把零模从格点 $i$ 移到相邻格点 $j$，可以视为在一小段时间内只打开边 $(i,j)$ 上的耦合（由对应的 $R_{ij}(u_{ij})$ 决定），并选取演化时间使得在零模子空间中的有效演化近似为
    $$
    U_{ij}^{\text{(braid)}}\;\approx\;\exp\Big(\pm\frac{\pi}{4}\,\gamma_i\gamma_j\Big),
    $$
    这就是标准的 Majorana 交换算符形式。

- **整条路径：路径上局域步骤的有序乘积**  
    在二维格上选取一条最近邻路径
    $$
    \gamma: i_0\to i_1\to\cdots\to i_n,
    $$
    其中每对 $(i_k,i_{k+1})$ 是最近邻格点。把上述局域幺正按路径顺序相乘，得到整体路径算符
    $$
    U_\gamma\;=\;\prod_{k=0}^{n-1} U_{i_k i_{k+1}}^{\text{(braid)}}\;\approx\;\prod_{k=0}^{n-1}\exp\Big(\pm\frac{\pi}{4}\,\gamma_{i_k}\gamma_{i_{k+1}}\Big).
    $$
    这给出了“沿路径 $\gamma$ 把零模搬运一圈”的有效编织算符。

- **路径依赖与同伦类**  
    若两条路径 $\gamma,\gamma'$ 在不穿过其他缺陷的前提下可以连续变形到彼此（同伦等价），在理想拓扑极限下，对应的 $U_\gamma,U_{\gamma'}$ 在零模子空间中只相差整体相位；当一条路径把一个缺陷绕另一缺陷一整圈，而另一条路径不绕圈时，它们属于不同同伦类，对应的 $U_\gamma$ 在零模简并子空间中给出真正不同的线性变换，这就是 2D 中非阿贝尔任何子编织的本质。

在本节前面的 Z2+Majorana 框架中，沿路径 $\gamma$ 的“局域一步”总是通过某条边上的 $R_{ij}(u_{ij})$（或等价的 $u_{ij}(ic_ic_j)$）实现；因此 $U_\gamma$ 可以视作一串沿 $\gamma$ 的局域 R/链路耦合的时间演化所产生的整体幺正。不同的路径（尤其是绕圈与不绕圈）在零模子空间中给出不同的 $U_\gamma$，从而把“R→Kitaev→2D Z2 规范场” 与编织算符联系起来。


### Bulk 拓扑相图：动量空间 BdG + Chern 数

上面对 2D+Z2 模型给出了局域和实空间的构造。若只关心 bulk 拓扑相图（class D 的拓扑超导），可以在固定一个平移对称的 Z2 背景下，把问题完全化到动量空间并用 Chern 数来标记不同相位。本节给出一个与前文 R→Kitaev+Z2 一致的标准流程。

#### 1. 固定平移对称的 Z2 背景

为能做傅里叶变换，需要选择一个在晶格平移下周期的 Z2 链变量构型 $\{u_{ij}\}$。典型选择：

- **无通量背景**：对所有最近邻键取 $u_{ij}=+1$。这是最简单的情形，对应 $W_p=+1$ 的扇区；
- **均匀 $\pi$‑flux**：要求每个 plaquette 的 $W_p=\prod_{(ij)\in p}u_{ij}=-1$，并选用一个有限大小的周期单胞（例如 2×1 或 2×2）来重复。这种背景在某些参数下更容易实现非平庸拓扑相。

在给定的 $\{u\}$ 背景下，费米子/ Majorana 的哈密顿 $H_f[c,u]$ 是平移不变或具有有限周期的，可以对 $c_i$ 做 Bloch 展开。

#### 2. 写出动量空间 BdG 哈密顿

以二维方格上 spinless p‑wave 超导为例，在无通量背景下，取一格点/一轨道的布里渊区，自然的 Nambu 自旋量为
$$
\Psi_{\mathbf{k}}=(c_{\mathbf{k}},c^{\dagger}_{-\mathbf{k}})^T.
$$
那么 BdG 哈密顿可以写成
$$
H_f = \sum_{\mathbf{k}} \Psi^{\dagger}_{\mathbf{k}}\,\mathcal{H}(\mathbf{k})\,\Psi_{\mathbf{k}} + \text{const.},
$$
其中
$$
\mathcal{H}(\mathbf{k}) = \xi(\mathbf{k})\,\tau_z + \Re\Delta(\mathbf{k})\,\tau_x - \Im\Delta(\mathbf{k})\,\tau_y.
$$

从局域 R 参数得到 $t_x,t_y,\Delta_x,\Delta_y,\mu$ 后，有
$$
\xi(\mathbf{k}) = -2t_x\cos k_x - 2t_y\cos k_y - \mu,
$$
而配对函数取决于配对在不同方向上的相位约定：

- 若在 $\pm\hat x,\pm\hat y$ 上配对幅度相同并且实数，则
    $$\Delta(\mathbf{k}) = 2\Delta_x\cos k_x + 2\Delta_y\cos k_y;$$
- 若想得到真正的 p‑wave 结构，可令正负方向配对相反号，例如 $+\hat x$ 上为 $+\Delta_x$，$-\hat x$ 上为 $-\Delta_x$，则
    $$\Delta(\mathbf{k}) = 2i\Delta_x\sin k_x + 2i\Delta_y\sin k_y,$$
    这对应 $p_x+p_y$ 型；若再让 $\Delta_y=\pm i\Delta_x$，可得到 $p_x\pm ip_y$ 型复配对。

对蜂窝格或更一般双子格/多轨道情形，只需把 $\mathcal{H}(\mathbf{k})$ 扩展为 $2N_{\text{orb}}\times2N_{\text{orb}}$ 的 BdG 矩阵，并用前面给出的 $f(\mathbf{k})$ 或 $A(\mathbf{k})$ 构造出跳跃和配对块结构即可。

#### 3. 闭缝条件：相界上的解析方程

bulk 能谱为
$$
E_{\pm}(\mathbf{k}) = \pm \sqrt{\xi(\mathbf{k})^2 + |\Delta(\mathbf{k})|^2}.
$$
能隙闭合的必要条件是存在某个动量 $\mathbf{k}_c$ 使得
$$
\xi(\mathbf{k}_c) = 0,\qquad \Delta(\mathbf{k}_c) = 0.
$$
在方格上，上式往往在高对称点处有解（例如 $(0,0),(\pi,0),(0,\pi),(\pi,\pi)$）。具体做法：

- 对每个高对称点代入 $\xi(\mathbf{k}),\Delta(\mathbf{k})$，得到若干关于 $t_x,t_y,\Delta_x,\Delta_y,\mu$ 的解析方程；
- 这些方程给出的曲线/超曲面，就是不同拓扑相之间的“相界”（gap closing 线）。

例如在各向同性连续极限的 $p_x+ip_y$ 模型中，相界条件简化为 $|\mu|=4t$ 一类的关系；在晶格上会有若干不同高对称点贡献，得出一系列闭缝曲线。

#### 4. Chern 数与拓扑相的标记

在有能隙的区域（对所有 $\mathbf{k}$ 有 $E_-(\mathbf{k})<0,E_+(\mathbf{k})>0$ 且无零交叉），可以对占据的 BdG 带计算 Chern 数来标记拓扑相：
$$
C = \frac{1}{2\pi}\int_{\text{BZ}} d^2k\;\mathcal{F}_{xy}(\mathbf{k}),
$$
其中 $\mathcal{F}_{xy}$ 是 Berry 曲率。数值上常采用离散 Brillouin 区网格和 Berry 链接公式来计算（例如 Fukui–Hatsugai–Suzuki 算法）：

1. 在每个格点 $\mathbf{k}$ 上计算占据带的本征矢；
2. 沿小方块的边计算 Berry 链接（内积的相位）；
3. 把每个小方块的曲率求和得到离散近似的 Chern 数。

计算策略通常是：

- 先用闭缝条件找到相界的解析曲线；
- 在每个由相界围成的区域中任选若干代表点，计算一次 Chern 数；
- 把该 Chern 数赋给整个区域，从而得到 bulk 拓扑相图（不同整数对应不同拓扑相）。

#### 5. 从 R 参数空间到拓扑相图

将以上 BdG 结果与本工作中的 R→$(t,\Delta,\mu)$ 识别结合：

- 选定一个子空间的 R 参数（例如固定 $a,d$，只在 $(b,c)$ 平面上扫描；或者固定 $b,c$，扫描 $d$ 以改变 $\mu$）；
- 对每一点的 $(a,b,c,d)$，用 1D 局域公式得到对应的 $t_x,t_y,\Delta_x,\Delta_y,\mu$（按需要对不同方向选不同的 R）；
- 把这些 BdG 参数代入上面的 $\xi(\mathbf{k}),\Delta(\mathbf{k})$，检查是否闭缝并在 gap 区域内计算 Chern 数；
- 这样可在 R‑参数空间中画出“拓扑相图”：每一个 R‑区域对应一个确定的 Chern 数（例如 C=0,\pm1），从而直接把“YBE 解族 / R 矩阵族”与 2D Majorana 拓扑相联系起来。

在具体数值实现时，可以：

- 选定一个固定的 Z2 背景（如 $u_{ij}=1$），在二维 $(\mu/t,\Delta/t)$ 或由 R 归一化后的二维参数平面上取网格；
- 对每个网格点构造 $\mathcal{H}(\mathbf{k})$，用离散 BZ 上的最小能隙判定是否 gapless，并在 gapped 点上计算 Chern 数；
- 用不同颜色/整数在参数平面上标出 Chern 数，得到 bulk 拓扑相图；
- 如需考察 Z2 背景变化的影响，可在几种代表性的 $\{u\}$ 构型上重复以上步骤比较。