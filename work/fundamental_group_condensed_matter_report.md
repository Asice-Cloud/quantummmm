# 基本群与凝聚态物质的拓扑结构

---

## 摘要

基本群（fundamental group）$\pi_1$ 是代数拓扑中最基础的不变量，用于刻画一个空间的「一维孔洞」结构。在凝聚态物理中，基本群以三种层层递进的方式出场：实空间中拓扑缺陷的分类——涡旋、向错线的稳定性由序参量流形的基本群决定；动量空间中能带拓扑不变量——当 Bloch 哈密顿量分类空间的基本群为 $G$ 时，Brillouin 区点缺陷周围的量子态单值性由 Drinfeld 中心 $\mathcal{Z}(\text{Vec}_G)$ 描述 [13]；辫群作为构型空间的基本群与映射空间 $\text{Map}(\Sigma, S^{2k})$ 的基本群（Heisenberg 群）——二者共同为非阿贝尔任意子的编织统计与可观测代数提供了严格的数学根基 [14]。本文主要通过一些调研，结合了基本群理论和物理现象来探讨二者之间的关系。文章沿三条主线展开，辅以两项 2026 年的最新案例研究，阐述基本群如何为凝聚态物理的拓扑现象提供统一的数学语言，并重点讨论非阿贝尔基本群（如 $Q_8$）对应的物理实现。

---

## 1. 引言

凝聚态物理过去四十年的标志性进展——从整数量子霍尔效应（von Klitzing, 1980）到拓扑绝缘体（Kane & Mele, 2005），再到拓扑量子计算（Kitaev, 2003 [4]）——共享一条核心线索：拓扑不变量制约着量子多体系统的低能物理。陈数（Chern number）、$Z_2$ 不变量、辫群表示……这些拓扑量背后有一个共同的数学母体：同伦理论。

同伦理论最简单的对象是**基本群** $\pi_1(\mathcal{M})$，它描述从单位圆 $S^1$ 到一个空间 $\mathcal{M}$ 的连续映射在连续形变下的等价类。初看之下，$\pi_1$ 似乎过于朴素——如何与能带理论、分数量子霍尔效应这些量子现象产生关联？本报告在结合了近年来的论文基础上，通过调研，试图证明一个相反的观点：基本群恰恰是凝聚态拓扑现象中最广泛、最深刻的数学结构。只需将「空间 $\mathcal{M}$」换上不同的物理对象——序参量流形、Brillouin 区上的分类空间、多个任意子的构型空间或映射空间——基本群就在截然不同的物理语境中给出了拓扑分类和物理预言。

本文调研主要按以下三条主线展开。第一层（实空间）：$\mathcal{M}$ 为序参量流形，$\pi_1(\mathcal{M})$ 分类涡旋等线缺陷。当 $\pi_1(\mathcal{M})$ 是非阿贝尔群时，缺陷的编织呈现非交换结构，双轴向列液晶中的 $Q_8$ 涡旋是这一现象的典型代表。第二层（动量空间）：$\mathcal{M}$ 为含隙 Bloch 哈密顿量的分类空间，当该空间的基本群为有限群 $G$ 时，Brillouin 区中绕点缺陷的量子态单值性由 Drinfeld 中心 $\mathcal{Z}(\text{Vec}_G)$ 描述 [13]——基本群直接决定了分数拓扑绝缘体中的任意子融合规则，超越了传统周期表给出的整数不变量。第三层（多体空间）：一方面，$n$ 个全同粒子的构型空间具有基本群 $B_n$（辫群），其高维不可约表示对应非阿贝尔任意子；另一方面，映射空间 $\text{Map}(\Sigma_g, S^2)$ 的基本群为整数 Heisenberg 群，其群代数恰好是 FQH 阿贝尔任意子的量子可观测代数 [14]——这一结果将辫群统计推广到了高维。

---

## 2. 数学预备

设 $\mathcal{M}$ 为拓扑空间，取基点 $x_0 \in \mathcal{M}$。考虑所有以 $x_0$ 为起止点的闭合路径 $\gamma: [0,1] \to \mathcal{M}$，$\gamma(0) = \gamma(1) = x_0$。两条路径称为**同伦等价**，如果它们可以经连续形变相互转换。路径的等价类在拼接操作（先走 $\gamma_1$ 再走 $\gamma_2$）下构成群 $\pi_1(\mathcal{M}, x_0)$。对路径连通的 $\mathcal{M}$，基本群不依赖于基点选择，记为 $\pi_1(\mathcal{M})$——这是整个报告的核心数学对象 [12]。

在凝聚态物理中最常出现的流形及其基本群如下表所示。$S^1$ 的基本群 $\mathbb{Z}$ 决定了超流涡旋的整数量子化；$\mathbb{R}P^n$（$n \ge 2$）的 $\mathbb{Z}_2$ 分类了向错线的二值拓扑；$SO(3)$ 的 $\mathbb{Z}_2$ 解释了自旋 $1/2$ 波函数在 $2\pi$ 旋转下的符号反转是不可消除的拓扑性质；$S^3/Q_8$ 的 $Q_8$（四元数群）是一个非阿贝尔基本群的例子，对应双轴向列液晶中的非阿贝尔涡旋。$S^n$（$n \ge 2$）的基本群平凡，因此这些流形中没有稳定的线缺陷。

| 流形 $\mathcal{M}$ | $\pi_1(\mathcal{M})$ | 物理含义 |
|---|---|---|
| $S^1$ | $\mathbb{Z}$ | 涡旋量子数可加 |
| $\mathbb{R}P^n$ ($n \ge 2$) | $\mathbb{Z}_2$ | 向错线二值分类 |
| $SO(3)$ | $\mathbb{Z}_2$ | 旋量表示的双值性 |
| $S^3/Q_8$ | $Q_8$ (四元数群) | 非阿贝尔涡旋 |
| $U(1)$ | $\mathbb{Z}$ | 超流/超导涡旋 |
| $S^n$ ($n \ge 2$) | $0$ (平凡) | 无稳定线缺陷 |

基本群 $\pi_1$ 分类一维孔洞，更一般地，$\pi_n(\mathcal{M})$ 分类 $S^n \to \mathcal{M}$ 的映射类。在凝聚态中，$\pi_2$ 分类点缺陷（单极子、hedgehog），$\pi_3$ 分类瞬子/拓扑项（Wess-Zumino 项），而 $\pi_1$ 是缺陷分类的基础。在对称性破缺的朗道范式中，有序相由序参量流形 $\mathcal{M} = G/H$ 描述，其中 $G$ 是系统的对称群，$H$ 是未破缺子群。拓扑缺陷由 $\mathcal{M}$ 的同伦群分类——这正是基本群的第一层出场点 [1, 2, 12]。

---

## 3. 第一层：实空间拓扑缺陷的分类

### 3.1 缺陷分类的核心定理

在 $d$ 维空间中，$m$ 维拓扑缺陷由同伦群 $\pi_{d-m-1}(\mathcal{M})$ 分类 [1, 2]。取 $d=3$（三维实空间）：$m=0$（点缺陷）由 $\pi_2(\mathcal{M})$ 分类；$m=1$（线缺陷，如涡旋、向错线）由 $\pi_1(\mathcal{M})$ 分类；$m=2$（面缺陷，如畴壁）由 $\pi_0(\mathcal{M})$ 分类。也就是说，凝聚态中所有稳定线缺陷的分类就是基本群。

先看几个阿贝尔基本群的例子。超流 He-4 的序参量是复标量 $\psi = |\psi| e^{i\varphi}$，基态流形为 $S^1 = U(1)$。$\pi_1(S^1) = \mathbb{Z}$，每个涡旋携带一个缠绕数 $n \in \mathbb{Z}$，两个涡旋碰撞时缠绕数相加：$n_1 + n_2$（阿贝尔融合），环流量子化为 $\oint \mathbf{v}_s \cdot d\mathbf{l} = \frac{2\pi\hbar}{m}n$。向列型液晶的序参量是无方向的棍状分子，流形为 $\mathbb{R}P^2 = S^2/\mathbb{Z}_2$（球面上对径点等同）。Mermin [1] 在 1979 年的综述中系统阐述了这一分类：$\pi_1(\mathbb{R}P^2) = \mathbb{Z}_2$，只有两类不可约的向错线——平凡类 $[0]$ 和非平凡类 $[1]$，两个 $[1]$ 类向错线可以合并湮灭：$[1] + [1] = [0]$。三维旋转群 $SO(3)$ 的基本群也为 $\mathbb{Z}_2$，这意味着自旋 $1/2$ 波函数在 $2\pi$ 旋转下获得负号是拓扑不可消除的，这一「旋量双值性」直接源于 $\pi_1(SO(3)) = \mathbb{Z}_2$。

### 3.2 非阿贝尔基本群：双轴向列液晶中的 $Q_8$ 涡旋

双轴向列液晶（biaxial nematic）的分子有三个不等价的轴 $(\mathbf{n}, \mathbf{m}, \mathbf{l})$，构成正交标架。序参量流形为 $\mathcal{M} = SO(3)/D_2$，其中 $D_2 = \{1, R_x^\pi, R_y^\pi, R_z^\pi\}$ 是二面体群。这个商空间同胚于三维球面商去四元数群：$\mathcal{M} \cong S^3/Q_8$ [1, 3]。

由于 $S^3$ 是单连通的（$\pi_1(S^3) = 0$），由商空间的提升定理立即得到：

$$\pi_1(S^3/Q_8) = Q_8 = \{\pm 1, \pm i, \pm j, \pm k\}$$

$Q_8$ 的乘法规则为 $i^2 = j^2 = k^2 = -1$，$ijk = -1$，$ij = k$，$ji = -k$，其他群元同理。在双轴向列液晶中，这意味着存在 **8 种**基本涡旋类型，标记为 $\pm 1, \pm i, \pm j, \pm k$。涡旋的融合规则服从 $Q_8$ 乘法：$i$ 型涡旋与 $j$ 型涡旋融合产生 $k$ 型（$ij = k$）。由于 $ij \neq ji$，涡旋编织是**非交换的**：先绕 $i$ 涡旋再绕 $j$ 涡旋的结果不同于先绕 $j$ 再绕 $i$。这是实空间中纯经典的非阿贝尔辫群结构,完全不需要量子力学的表示。

Mermin [1] 最早系统建立了非阿贝尔缺陷的一般理论，Volovik 和 Mineev [3] 在超流 He-3 的语境中进行了推广。设 $\pi_1(\mathcal{M}) = \Gamma$ 是非阿贝尔群，则缺陷的「拓扑量子数」标记为 $\Gamma$ 中的共轭类，两个缺陷的融合结果为其共轭类在特定条件下的积，编织效应由 $\Gamma$ 在自身上的共轭作用描述。

超流 He-3 提供了另一个丰富的实例。A 相的序参量流形为 $(SO(3)/\mathbb{Z}_2) \times S^2$，$\pi_1 = \mathbb{Z}_4$（4 种涡旋，携带 $n/4$ 单位的环流）。B 相的序参量流形为 $SO(3)$，$\pi_1 = \mathbb{Z}_2$（2 种涡旋）。这些数据整理自 Volovik 和 Mineev [3]（另见 Mermin [1]）。

---

## 4. 第二层：动量空间能带拓扑

### 4.1 分类空间与周期表

能带理论将电子态按晶格动量 $\mathbf{k}$ 标记。含隙哈密顿量在对称类约束下的**分类空间**记为 $\mathcal{C}_q$（复对称类）或 $\mathcal{R}_q$（实对称类）。这一分类框架由 Kitaev [5] 提出，并由 Ryu、Schnyder、Furusaki 和 Ludwig [6] 系统化为「十倍道」（tenfold way）分类。拓扑不同的绝缘体由映射 $H: T^d \to \mathcal{C}_q$（或 $\mathcal{R}_q$）的同伦类分类。

Kitaev 和 Ryu 等人的周期表本质上计算了各种对称类下的同伦群。在一维系统中，基本群 $\pi_1$ 直接给出拓扑不变量：AIII 类（手征幺正）的 $\pi_1(\mathcal{C}_1) = \mathbb{Z}$ 给出 SSH 模型的缠绕数 [10]；D 类（粒子-空穴对称）的 $\pi_1(\mathcal{R}_2) = \mathbb{Z}_2$ 给出 Kitaev 链中 Majorana 零能模的存在性判据。二维陈绝缘体（A 类）由 $\pi_2(\mathcal{C}_0) = \mathbb{Z}$（第一陈数）分类；三维 $Z_2$ 拓扑绝缘体（AII 类）涉及稳定同伦群 $\pi_4$；一维 D 类的 $\pi_0(\mathcal{R}_2) = \mathbb{Z}_2$ 则分类涡旋束缚的 Majorana 零能模的存在与否。同一系统中不同低位的同伦群各有不同的物理含义：$\pi_0$ 给出二值判断，$\pi_1$ 给出缠绕数，$\pi_2$ 给出陈数。

### 4.2 案例研究：Drinfeld 中心与 Brillouin 区缺陷的单值性 [13]

上述周期表框架虽然系统，但基本群在其中扮演的角色只是「给出一个整数或 $Z_2$ 不变量」。Sati 与 Schreiber（2026）的最新工作 [13] 将这一范式推向了一个全新的层次：基本群直接生成融合范畴。

考虑一个含隙的 $d$ 维 Bloch 哈密顿量 $H(\mathbf{k})$，其分类空间记为 $\mathcal{B}$，并假设 $\mathcal{B}$ 的**基本群非平凡**：$\pi_1(\mathcal{B}) = G$（$G$ 为有限群）。在 Brillouin 区中考虑一个点缺陷（如 Weyl 点的投影），沿包围该缺陷的小环路 $\gamma \simeq S^1$，$H(\mathbf{k})|_{\mathbf{k} \in \gamma}$ 构成一个 $S^1 \to \mathcal{B}$ 的映射，其同伦类恰好是 $\pi_1(\mathcal{B}) = G$ 中的元素。将缺陷周围的参数空间视为穿孔圆盘 $D^2 \setminus \{0\}$（边界为 $\gamma$），绕缺陷一周的绝热输运给出基态空间上的一个单值性变换 $T_\gamma$。

Sati 与 Schreiber 证明了以下核心结果：当 Bloch 哈密顿量分类空间的基本群为 $G$ 时，穿孔圆盘上 gapped 基态的单值性——即缺陷附近的拓扑序——由 **Drinfeld 中心** $\mathcal{Z}(\text{Vec}_G)$ 描述，其任意子类型的融合规则与 $\mathcal{Z}(\text{Vec}_G)$ 的融合范畴完全一致。逻辑链如下：$\pi_1(\mathcal{B}) = G$ 意味着 $\mathcal{B}$ 的万有覆盖 $\widetilde{\mathcal{B}}$ 具有 Deck 变换群 $G$，绕缺陷的环路提升到万有覆盖后端点的差异由某个 $g \in G$ 给出；将绕缺陷的 $S^1$ 视为有效 (1+1)d 时空的「空间方向」，则该单值性 $g$ 定义了一个 $G$-规范通量；已知在 (2+1)d 中，纯 $G$ 规范理论的拓扑序由 Drinfeld 中心描述——其任意子为 $G$ 的共轭类标记的磁通和不可约表示标记的电荷的直积；这些融合规则直接给出了分数拓扑绝缘体材料中点缺陷附近涌现的准粒子统计。

以 $G = \mathbb{Z}_N$ 为例，$\mathcal{Z}(\text{Vec}_{\mathbb{Z}_N})$ 是 $\mathbb{Z}_N$ 量子双重模型（toric code 的推广），具有 $N^2$ 种任意子类型，标记为 $(a, \chi)$（$a \in \mathbb{Z}_N$ 为磁通，$\chi \in \widehat{\mathbb{Z}_N} \cong \mathbb{Z}_N$ 为电荷），融合规则为 $(a,\chi) \times (b,\xi) = (a+b, \chi \cdot \xi)$，编织相位为 $S_{(a,\chi),(b,\xi)} = \frac{1}{N} \chi(b) \xi(a)$。这给出了一套**直接从 $\pi_1(\mathcal{B})$ 计算分数拓扑绝缘体缺陷附近任意子谱**的完整方案。与经典的 SSH 缠绕数 $\pi_1(S^1) = \mathbb{Z}$ 相比，Sati-Schreiber 框架的深刻之处在于：它处理了有限非平凡基本群，而且基本群输出的不再是一个整数——而是一个完整的融合范畴（任意子类型、融合规则、编织统计）。

---

## 5. 第三层：辫群、Heisenberg 群与任意子统计

### 5.1 构型空间的基本群：辫群

$n$ 个全同不可区分粒子在二维平面 $\mathbb{R}^2$ 上的构型空间为 $\text{Conf}_n(\mathbb{R}^2) = (\mathbb{R}^{2n} - \Delta)/S_n$，其中 $\Delta$ 是对角线（排除粒子重合的情况），$S_n$ 是置换群（体现不可区分性）。**辫群恰好是这个构型空间的基本群** [11, 12]：$B_n = \pi_1(\text{Conf}_n(\mathbb{R}^2))$。

辫群 $B_n$ 的生成元为 $\sigma_1, \ldots, \sigma_{n-1}$（$\sigma_i$ 表示第 $i$ 与 $i+1$ 个粒子的逆时针交换），满足 Artin 辫关系 [11]：$\sigma_i\sigma_{i+1}\sigma_i = \sigma_{i+1}\sigma_i\sigma_{i+1}$，以及 $\sigma_i\sigma_j = \sigma_j\sigma_i$（$|i-j| \ge 2$）。辫群的非交换性来自第一个 Yang-Baxter 类型的辫关系。若加上 $\sigma_i^2 = 1$ 的约束，$B_n$ 退化为置换群 $S_n$（对应玻色子/费米子）；在非约束情况下，$\sigma_i$ 的特征值不在 $\{\pm 1\}$ 中——这就是**任意子**的数学根源 [7]。

任意子的统计由辫群 $B_n$ 的**酉表示** $R: B_n \to U(D)$ 决定 [7]：一维表示 $R(\sigma_i) = e^{i\theta}$ 给出阿贝尔任意子（$\theta = 0$ 为玻色子，$\theta = \pi$ 为费米子）；高维表示 $R(\sigma_i) = R^{ab}_c$（R-矩阵）给出非阿贝尔任意子。对于非阿贝尔任意子，粒子编织的最终量子态依赖于编织的顺序，因此携带集体量子信息——这正是拓扑量子计算的物理基础 [4, 7]。Ising 任意子（对应 $\nu=5/2$ Moore-Read 态 [8]）的融合规则 $\sigma \times \sigma = 1 + \psi$ 和 Fibonacci 任意子（对应 $\nu=12/5$ Read-Rezayi 态 [9]）的融合规则 $\tau \times \tau = 1 + \tau$ 都是辫群高维表示的经典实例。

### 5.2 案例研究：映射空间的基本群与 Heisenberg 群任意子代数 [14]

辫群 $B_n = \pi_1(\text{Conf}_n)$ 描述的是点状粒子的编织。Kallel、Sati 与 Schreiber（2026）[14] 揭示了基本群进入任意子理论的另一条路径——这一次，$\mathcal{M}$ 是**映射空间**，基本群是整数 Heisenberg 群，物理对象是通量插入算符的代数。

在封闭曲面 $\Sigma_g$（亏格 $g$）上的 FQH 系统中，阿贝尔任意子的量子可观测代数构成一个被长期忽视的结构：它们构成整数 Heisenberg 群的群代数。对于 Laughlin 态 $\nu = 1/m$，环绕 $a_i$ 和 $b_j$ 周期（$\Sigma_g$ 的 1-同调基）的磁通插入算符满足 $U_{a_i} U_{b_j} = e^{2\pi i \delta_{ij}/m} U_{b_j} U_{a_i}$——这是一组 Heisenberg 对易关系，对应的群为 $\text{Heis}(\mathbb{Z}^g, m) = \mathbb{Z}^g \times \mathbb{Z}^g \times \mathbb{Z}_m$。

关键洞察来自同伦论的经典结果（Thom, 1950s；被 Kallel-Sati-Schreiber [14] 重新挖掘并推广）：亏格 $g$ 的封闭曲面到 2-球面的映射空间 $\text{Map}(\Sigma_g, S^2)$ 的**基本群**的非挠部分恰好由整数 Heisenberg 群给出：

$$\pi_1\left(\text{Map}(\Sigma_g, S^2)\right)_{\text{non-torsion}} \cong \text{Heis}(\mathbb{Z}^{2g}, 2)$$

其中 level $2$ 来自 Hopf 不变量 $h(S^3 \to S^2) = 1$ 的倒数关系。完整的逻辑链为：FQH 系统中，复合费米子/通量附着将物理电子映射到「映射空间」中的点——$\Sigma_g$ 上的每个位形对应一个映射 $\Sigma_g \to S^2$（本质上是 2-Cohomotopy 理论中的通量量子化）；$\pi_1(\text{Map}(\Sigma_g, S^2))$ 中的环路对应绕 FQH 系统中各周期的通量插入操作，基本群的乘法规则（即 Heisenberg 对易关系）精确地给出了任意子量子可观测量之间的代数结构；Heisenberg 群的群代数 $\mathbb{C}[\text{Heis}]$ 恰好就是 Laughlin 态 $\nu = 1/m$ 的阿贝尔任意子的所有量子可观测量的代数。

Kallel、Sati 与 Schreiber 的核心贡献是将上述结果系统推广到高维。他们证明，对于 $k \in \{1, 2, 4\}$：

$$\pi_1\left(\text{Map}\big((S^{2k-1})^2, \; S^{2k}\big)\right)_{\text{non-torsion}} \cong \text{Heis}(\mathbb{Z}, 2)$$

其中 level $2$ 统一地由 $4k-1$ 维球面生成元的 Hopf 不变量决定：$\text{level} = 2 / h(\pi_{4k-1}(S^{2k})_{\text{gen}})$。$k=1$ 对应前述 FQH 在环面上的情况（$S^1 \times S^1 \to S^2$）；$k=2$ 对应 (4+1)d 的「高维 FQH 任意子」（$S^3 \times S^3 \to S^4$）；$k=4$ 对应 11 维超引力/M-理论中「Hypothesis H」预言的任意子型激发（$S^7 \times S^7 \to S^8$）。这三个 $k$ 值恰好对应 Hopf 不变量非平凡的唯一三个维度——实除法代数 $\mathbb{R}, \mathbb{C}, \mathbb{H}, \mathbb{O}$ 的 Hopf 纤维化——这是代数拓扑与凝聚态物理之间惊人的交叉。

传统辫群 $B_n = \pi_1(\text{Conf}_n)$ 和映射空间基本群是互补而非排斥的。前者描述点粒子的编织（物理粒子为点状任意子，基本群为无限非阿贝尔辫群，限于 (2+1)d），后者描述通量算符的代数结构（物理对象为通量插入算符，基本群为幂零 Heisenberg 群，可推广到任意 $4k-1 \to 2k$ 维度）。两者共同构成了完整的基本群-任意子对应。

---

## 6. 综合讨论

三个层次的统一性体现在：它们共享同一数学结构——基本群，区别仅在于 $\mathcal{M}$ 的物理含义。实空间中 $\mathcal{M}$ 为序参量流形 $G/H$，$\pi_1$ 分类线缺陷（涡旋、向错线）；动量空间中 $\mathcal{M}$ 为 Bloch 分类空间 $\mathcal{B}$，$\pi_1(\mathcal{B}) = G$ 决定了缺陷 monodromy 的 Drinfeld 中心 $\mathcal{Z}(\text{Vec}_G)$ [13]，或由分类空间 $\mathcal{C}_q/\mathcal{R}_q$ 给出缠绕数等不变量 [5,6]；多体空间中 $\mathcal{M}$ 既可以是构型空间 $\text{Conf}_n$（$\pi_1 = B_n$，辫群统计），也可以是映射空间 $\text{Map}(\Sigma, S^{2k})$（$\pi_1 = \text{Heis}$ 群，任意子可观测代数）[14]。

非阿贝尔性的起源同样有三种截然不同的物理来源。第一种是序参量流形的非阿贝尔 $\pi_1$（如 $S^3/Q_8$ 的 $Q_8$ 群）——涡旋融合的非交换性是纯实空间几何效应，属于经典拓扑缺陷。第二种是编织矩阵的非阿贝尔性——源于构型空间 $\text{Conf}_n$ 的基本群 $B_n$ 的高维不可约表示，编织顺序影响量子态，是拓扑量子计算的存储机制。第三种是 Berry 联络的非阿贝尔性——源于简并基态空间上 Wilson 环 $W = \mathcal{P}\exp(-\oint A)$ 的非交换holonomy。

综合以上，本报告提出以下统一的叙事框架。对称性破缺产生序参量流形 $\mathcal{M}$，其同伦群分类实空间的拓扑缺陷 [1,2,3]。能带结构产生 Bloch 分类空间 $\mathcal{B}$，其基本群 $G$ 通过 Drinfeld 中心 $\mathcal{Z}(\text{Vec}_G)$ 描述缺陷 monodromy 的融合范畴 [13]；同时，分类空间 $\mathcal{C}_q/\mathcal{R}_q$ 的同伦群给出传统拓扑不变量 [5,6]。多体波函数产生构型空间 $\text{Conf}_n$（$\pi_1 = B_n$，辫群统计 [11]）和映射空间 $\text{Map}(\Sigma, S^{2k})$（$\pi_1 = \text{Heis}$，可观测代数 [14]）。五条线在数学结构上同源于同伦理论，但在物理上对应截然不同但互补的现象。

---

## 7. 参考文献

[1] N. D. Mermin, *The topological theory of defects in ordered media*, Rev. Mod. Phys. **51**, 591 (1979).

[2] G. Toulouse and M. Kléman, *Principles of a classification of defects in ordered media*, J. Phys. Lett. **37**, L-149 (1976).

[3] G. E. Volovik and V. P. Mineev, *Investigation of singularities in superfluid He-3 in liquid crystals by the homotopic topology methods*, Sov. Phys. JETP **45**, 1186 (1977).

[4] A. Kitaev, *Anyons in an exactly solved model and beyond*, Ann. Phys. **321**, 2 (2006).

[5] A. Kitaev, *Periodic table for topological insulators and superconductors*, AIP Conf. Proc. **1134**, 22 (2009).

[6] S. Ryu, A. P. Schnyder, A. Furusaki, and A. W. W. Ludwig, *Topological insulators and superconductors: tenfold way and dimensional hierarchy*, New J. Phys. **12**, 065010 (2010).

[7] C. Nayak, S. H. Simon, A. Stern, M. Freedman, and S. Das Sarma, *Non-Abelian anyons and topological quantum computation*, Rev. Mod. Phys. **80**, 1083 (2008).

[8] G. Moore and N. Read, *Nonabelions in the fractional quantum Hall effect*, Nucl. Phys. B **360**, 362 (1991).

[9] N. Read and E. Rezayi, *Beyond paired quantum Hall states: Parafermions and incompressible states in the first excited Landau level*, Phys. Rev. B **59**, 8084 (1999).

[10] W. P. Su, J. R. Schrieffer, and A. J. Heeger, *Solitons in Polyacetylene*, Phys. Rev. Lett. **42**, 1698 (1979).

[11] J. S. Birman, *Braids, Links, and Mapping Class Groups*, Princeton University Press (1974).

[12] M. Nakahara, *Geometry, Topology and Physics*, 2nd ed., Taylor & Francis (2003).

[13] H. Sati and U. Schreiber, *Drinfeld Center as Quantum State Monodromy over Bloch Hamiltonians around Defects*, arXiv:2603.22029 (2026).

[14] S. Kallel, H. Sati, and U. Schreiber, *Higher-Dimensional Anyons via Higher Cohomotopy*, arXiv:2601.03150 (2026).
