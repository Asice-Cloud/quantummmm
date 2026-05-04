### kit2-exp：在 2D / Majorana + Z2 框架下把 R 改写成 exp 形式

本节是对 [kit2.md](kit2.md) 的补充，专门讨论：

1. 把局域两体算符从线性形式
   $$
   R_{ij}=a I + b\,\sigma^x_i\sigma^x_j + c\,\sigma^y_i\sigma^y_j + d\,\sigma^z_i\sigma^z_j
   $$
   改写为指数形式
   $$
   R_{ij}=e^{iK_{ij}},\qquad K_{ij}=J_x\,\sigma^x_i\sigma^x_j + J_y\,\sigma^y_i\sigma^y_j + J_z\,\sigma^z_i\sigma^z_j,
   $$
   以及与 2D Majorana+Z2 规范场表示之间的关系。
2. 在 2D 蜂窝/一般格子上，R=exp(iK) 如何自然地给出“局域平行移动元”和“Wilson 回路”的表述。
3. 对于之后在配置空间上讨论 holonomy（Dehn twist/half twist）时，2D 这一层的 exp 形式如何作为“空间方向上的联络”。


 

#### 1. 线性 R 与 exp 形式 R=exp(iK) 的对应

在 [R_to_kitaev2.md](R_to_kitaev2.md) 中，我们已经把 SU(2) 不变的常数 YBE 解改写为
$$
R=e^{iK},\qquad K=J_x\,\sigma^x\otimes\sigma^x + J_y\,\sigma^y\otimes\sigma^y + J_z\,\sigma^z\otimes\sigma^z.
$$
这里再强调几条对 2D 有用的事实：

- 运算子 $A=\sigma^x\otimes\sigma^x$、$B=\sigma^y\otimes\sigma^y$、$C=\sigma^z\otimes\sigma^z$ 两两对易，且 $A^2=B^2=C^2=I$，因此
  $$
  e^{iK}=e^{iJ_xA}e^{iJ_yB}e^{iJ_zC}
  $$
  可以逐项展开成 $I,A,B,C$ 的线性组合。
- 从而存在显式的三角表达式
  $$
  R = a(J) I + b(J)A + c(J)B + d(J)C
  $$
  （具体见 [R_to_kitaev2.md](R_to_kitaev2.md) 中的公式），给出指数参数 $(J_x,J_y,J_z)$ 和线性参数 $(a,b,c,d)$ 之间的对应。
- 对于 2D 来说，真正出现在哈密顿量中的仍然是 $K$ 本身：
  - 若取“哈密顿 = 生成元”的观点：$H=\sum_{\langle ij\rangle}K_{ij}$；
  - 若把 R 理解为“一步离散时间演化/门”，则 $U(\Delta t=1)=\prod_{\langle ij\rangle}R_{ij}$，在 Trotter 极限下等价于 $H$ 产生的时间演化。

因此：**在 2D Majorana+Z2 框架中，主要对象是 $K_{ij}$；改写成 exp(iK) 把边上的作用自然解读为“沿边的局域平行移动元”。**


 

#### 2. 在 2D 上的 K_{ij}：Majorana + Z2 映射

回顾 [kit2.md](kit2.md) 中的 Majorana+Z2 映射：在每个格点 $j$ 上引入 $b^x_j,b^y_j,b^z_j,c_j$，定义
$$
\sigma^a_j = i b^a_j c_j,\qquad u^{(a)}_{ij}=-i b^a_i b^a_j,\qquad \sigma^a_i\sigma^a_j = u^{(a)}_{ij} (i c_i c_j).
$$
对于单一类型的键 $(i,j)$（例如 $x$ 键），令
$$
K^{(x)}_{ij}=J_x\,\sigma^x_i\sigma^x_j.
$$
则在 Majorana 表示中
$$
K^{(x)}_{ij} = J_x\,u^{(x)}_{ij}\,(i c_i c_j).
$$
类似地，$y,z$ 键对应各自的 $u^{(y)}_{ij},u^{(z)}_{ij}$。于是一般的
$$
K_{ij}=J_x\,\sigma^x_i\sigma^x_j + J_y\,\sigma^y_i\sigma^y_j + J_z\,\sigma^z_i\sigma^z_j
$$
在 Majorana+Z2 语言下写成
$$
K_{ij}=J_x u^{(x)}_{ij}\,(i c_i c_j) + J_y u^{(y)}_{ij}\,(i c_i c_j) + J_z u^{(z)}_{ij}\,(i c_i c_j) + \text{(密度/四费米)}.
$$
若我们采取 Kitaev 蜂窝式的各向异性分配：每条边只带一种 `a` 型耦合，则对给定边 $(i,j)$ 实际上只有一个 $u_{ij}$ 出现，$K_{ij}$ 简化为
$$
K_{ij}=J_{a(ij)}\,u_{ij}\,(i c_i c_j) + \text{(可能的密度/四费米)}.
$$
在“自由 Majorana 主导 + 相互作用弱”的情形下，可以聚焦于二次部分：
$$
K_{ij}^{\text{(quad)}}\approx J_{a(ij)}\,u_{ij}\,(i c_i c_j).
$$

**关键点**：

- $K_{ij}^{\text{(quad)}}$ 在给定 Z2 背景 $u_{ij}=\pm1$ 下就是 Majorana 之间的最近邻双线性；
- 把 $R_{ij}=e^{iK_{ij}}$ 看成门/时间步时，其主要非平凡作用就是在对应边上的 Majorana 双线性旋转，幅度由 $J_{a(ij)}$ 和演化时间决定；
- 在 2D 上，沿某条空间路径 $\gamma$ 依次作用这些 $R_{ij}$，可以看作对经过的 Majorana 模式做一串“平行移动”。


 

#### 3. R=exp(iK) 作为 2D 上的“空间联络平行移动元”

在 Majorana+Z2 语言下，给定某个静态的 Z2 背景 $\{u_{ij}\}$，二次部分的哈密顿量可以统一写成
$$
H_c = \frac i4\sum_{p,q} A_{pq} c_p c_q
$$
其中 $A$ 是反对称实矩阵，$A_{pq}$ 由所有 $J_{a(ij)}u_{ij}$ 叠加而来。从 exp 形式的角度看，有两种自然的“联络”解释：

1. **空间方向上的离散联络**：
   - 对每条边 $(i,j)$，$K_{ij}^{\text{(quad)}}\sim J_{a(ij)}u_{ij}(i c_i c_j)$ 是一个 $so(2N)$ 代数元（Majorana 双线性）；
   - $R_{ij}=e^{iK_{ij}}\sim e^{\theta_{ij}\,c_i c_j}$ 则是 $SO(2N)$ 群元，是在 Majorana 模式空间中的一个“旋转”；
   - 把它理解成沿边 $(i,j)$ 的一次平行移动：
     $$
     \Gamma_{ij} \equiv R_{ij}=\exp\big(i K_{ij}\big).
     $$
2. **回路上的 Wilson 线/面元**：
   - 取一条空间闭合路径 $\gamma$（由边序列 $(i_1,i_2),\dots,(i_{n},i_1)$ 组成），定义
     $$
     U_\gamma = \prod_{(ij)\in\gamma}R_{ij} = \prod_{(ij)\in\gamma}\exp\big(iK_{ij}\big).
     $$
   - 在 $K_{ij}$ 之间“足够对易”或演化足够缓慢的局域极限下，可以近似为
     $$
     U_\gamma \approx \exp\Big(i\sum_{(ij)\in\gamma}K_{ij}\Big) = \exp\Big(\frac i4\sum_{p,q}\Omega^{(\gamma)}_{pq} c_p c_q\Big)
     $$
     其中 $\Omega^{(\gamma)}$ 是某个沿路径积分得到的 $so(2N)$ 角度矩阵。

与 [kit2.md](kit2.md) 中的 Z2 规范场 $u_{ij}$ 对比：

- $u_{ij}$ 与 plaquette flux $W_p=\prod_{(ij)\in p}u_{ij}$ 描述的是“真正的空间 Z2 规范场”；
- 而 $R_{ij}=e^{iK_{ij}}$ 中的 Majorana 双线性部分，则提供了一个值在 $SO(2N)$ 的“上层联络”；
- 二者叠加：实际的 Majorana 平行移动包含了 Z2 背景（通过 $u_{ij}$ 出现在 $K_{ij}$ 中）和 Majorana 自身的旋转（通过 $e^{iK}$）。

在之后 kit3-exp 中讨论配置空间/holonomy 时，2D 上的这一层可以被看成是“在给定缺陷/通量配置下，Majorana 模式在晶格上的平行移动联络”。


 

#### 4. 从边上的 exp(iK) 到路径上的 effective Majorana 双线性

为了后续与 Dehn twist/half twist 对接，这里把“边上的 exp(iK)”如何在低能 Majorana 子空间中等效为某对零模的双线性指数整理一下。

设某个 2D 构型中存在若干低能 Majorana 零模 $\{\gamma_a\}$（例如边界、涡旋/vison 端点处的零能态）。在给定参数 $J_x,J_y,J_z$ 与 Z2 背景 $u_{ij}$ 下，对应的有效低能哈密顿量可写为
$$
H_{\text{eff}} = \frac i2\sum_{a<b} A_{ab}\,\gamma_a\gamma_b.
$$
若我们沿某条空间路径 $\gamma$ 依次作用 $R_{ij}=e^{iK_{ij}}$，这些门在低能子空间上的累积作用可以压缩为
$$
U_\gamma^{\text{(low)}} = \exp\Big(\frac12\Theta_\gamma^{ab}\,\gamma_a\gamma_b\Big),
$$
其中 $\Theta_\gamma^{ab}$ 是某个 antisymmetric 角度矩阵，取决于：

- 路径 $\gamma$ 穿过了哪些 edge（以及对应的 $J_{a(ij)},u_{ij}$）；
- 演化的总“时间长度”（重复次数/脉冲面积）；
- 投影到零模子空间时的重叠结构。

在非常理想的极限下（只耦合到一对 Majorana，且其它耦合被抑制），可以得到
$$
U_\gamma^{\text{(low)}}\approx \exp\Big(\frac{\theta_\gamma}{2}\,\gamma_p\gamma_q\Big)
$$
这与 [R_to_Kitaev.md](R_to_Kitaev.md) 最后一节中典型 braid 生成子
$$
U=\exp\big((\pi/4)\,\gamma_2\gamma_3\big)
$$
的形式完全一致。

**要点**：exp 形式使得

- 从边上的门 $R_{ij}=e^{iK_{ij}}$ 到有效 braid gate $U^{\text{(low)}}=e^{(\theta/2)\gamma\gamma}$ 的连接变得非常直接；
- 与其事后对 R 取对数，不如从一开始就把 K 视为“生成元”，在 2D 上做路径积分、再投影到低能子空间，得到所需的 braid/holonomy；
- 这为 kit3-exp 中讨论“配置空间上的 Dehn twist / half twist = 某条路径上的 exp(iK) 的 holonomy”打下了基础。


 

#### 5. 小结：kit2 视角下 exp 形式带来的信息增量

相对于线性 R 形式，把 2D 上的边算符写成
$$
R_{ij}=e^{iK_{ij}},\qquad K_{ij}=J_{a(ij)}\,\sigma^{a(ij)}_i\sigma^{a(ij)}_j
$$
带来的主要“增量”是：

- 在 2D Majorana+Z2 语言下，$K_{ij}$ 自然被视为 Majorana 模式空间的 $so(2N)$ 联络元，$R_{ij}$ 则是对应的 $SO(2N)$ 平行移动元；
- 沿空间路径的门序列 $\prod_{(ij)\in\gamma}R_{ij}$ 可以统一看作某个有效 Majorana 双线性的指数 $e^{(1/2)\Theta_\gamma^{ab}\gamma_a\gamma_b}$，这为之后解释 braid/Dehn twist 提供了直接桥梁；
- 与 [kit2.md](kit2.md) 中的 Z2 链变量 $u_{ij}$ / plaquette flux $W_p$ 叠加后，2D 上实际的“联络”包含两个层次：
  1. Z2 规范场（$u_{ij}$ 决定 topological sector 与 vison 构型）；
  2. 在给定 Z2 背景下，由 K_{ij} 诱导的 Majorana 模式空间旋转（值在 $SO(2N)$）。

在这样的视角下，2D 模型既可以被看作“Majorana + Z2 规范场”的静态谱问题，也可以被看作“在空间与参数空间上携带平坦/近平坦联络的离散化 bundle”，这正好对接到 kit3 中关于配置空间与 holonomy 的讨论。


 

#### 6. 具体 2D Ising‑like 模型与“空间变换”的例子

这里选一个尽量贴近 Kitaev 蜂窝 Ising‑like 相的例子，说明：

1. 如何把 1D 中的 $K=J_x\sigma^x\sigma^x+J_y\sigma^y\sigma^y+J_z\sigma^z\sigma^z$ 推广到 2D 蜂窝格；
2. 在 exp 形式下，某些“空间上的操作”（如沿路径翻转 $u_{ij}$、施加局域 $R_{ij}=e^{iK_{ij}}$ 序列）怎样在 Majorana 描述里对应到拓扑缺陷/vison 的生成与移动，为 kit3 的配置空间和 Dehn twist/half twist 做准备。


**6.1 2D 蜂窝格上的 exp(iK) 模型**

- 取蜂窝格，按 [kit2.md](kit2.md) 的标准约定，将三类最近邻边标记为 $x,y,z$ 键，向量分别为 $\delta_x,\delta_y,\delta_z$；
- 在每条 $a$ 型键 $(i,j)$ 上放置同一个局域生成元
  $$
  K^{(a)}_{ij}=J_a\,\sigma^a_i\sigma^a_j,\qquad a\in\{x,y,z\};
  $$
  整体哈密顿量（生成元）为
  $$
  H_K=\sum_{\langle i,j\rangle}K_{ij}=\sum_{\langle i,j\rangle\in x}J_x\sigma^x_i\sigma^x_j+\sum_{\langle i,j\rangle\in y}J_y\sigma^y_i\sigma^y_j+\sum_{\langle i,j\rangle\in z}J_z\sigma^z_i\sigma^z_j.
  $$
- 若把一步离散演化写为
  $$
  U(\Delta t=1)=\prod_{\langle i,j\rangle}e^{iK_{ij}},
  $$
  在 Trotter 极限下 $U(\Delta t=1)\approx e^{iH_K}$，这就是一个典型的 Kitaev 蜂窝‑型二维模型的 exp(iK) 实现方式。

在 Majorana+Z2 语言下，取 $u_{ij}$ 为守恒 Z2 链路变量，在无磁通扇区（所有 $W_p=+1$）且忽略四费米项的极限中，
$$
H_K^{\text{(quad)}}=\frac i4\sum_{p,q}A_{pq}c_pc_q,
$$
其中 $A(k)$ 的动量空间形式与 [kit2.md](kit2.md) 中的蜂窝 Majorana 谱分析完全一致（$J_a$ 决定 $f(\mathbf k)$ 与 gapless/gapped 区域）。


**6.2 通过 exp(iK) 生成/移动 2D vison 的空间操作**

在 Z2 规范场视角下，plaquette flux
$$
W_p=\prod_{(ij)\in p}u_{ij}
$$
描述了每个小回路上的 Z2 磁通。vison（涡旋）对应 $W_p=-1$ 的 plaquette。标准的“生成一对 vison 并移动它们”的操作是：

- 选择一条路径 $\gamma$，为其上的边 $(i,j)\in\gamma$ 逐一翻转 $u_{ij}\to -u_{ij}$；
- 这样会使得路径两端的 $W_{p_\mathrm{start}},W_{p_\mathrm{end}}$ 翻转符号，其余 plaquette 保持不变。

在 exp 形式的微观实现中，可以用一串局域 $R_{ij}=e^{iK_{ij}}$ 来“实现”或“模拟”这样的操作。粗略图像是：

1. 某些特定的 $K_{ij}$ 作用可以有效翻转对应边上的 $u_{ij}$（例如在扩展模型中加入使 $u_{ij}$ 变成动力学自由度的项，再在适当时间尺度上施加脉冲）；
2. 沿路径 $\gamma$ 对各边作用合适的 $e^{iK_{ij}}$ 序列，相当于在 Majorana 模式空间中实施一次 Wilson 线操作 $U_\gamma$，并在起终端生成/移动 vison；
3. 在低能 Majorana 子空间中，这个 $U_\gamma^{\text{(low)}}$ 可以近似为某个 $\exp(\tfrac12\Theta_\gamma^{ab}\gamma_a\gamma_b)$，在包含端点附近零模的 Hilbert 空间里给出非平凡的拓扑相位。

因此，从 2D 空间的角度看：

- $R_{ij}=e^{iK_{ij}}$ 是“定义在边上的平行移动”；
- 通过合适的 $R_{ij}$ 序列，我们不仅可以沿 Majorana 链路做旋转，还可以改变 Z2 背景（$u_{ij}$）本身，从而在空间上移动拓扑缺陷（vison/涡旋）；
- 这些空间操作对应到 kit3 的“多体配置空间”里，就是在缺陷位置空间中的路径（braid、Dehn twist 等）。


**6.3 与 1D Ising 例子、kit3-exp 的衔接**

在 [kit3-exp.md](kit3-exp.md) 中，我们用一个 1D 例子说明了：

- 单键 $K\sim\sigma^x\sigma^x$ 在 JW+Majorana 之后给出 $i\gamma_2\gamma_3$；
- 选择演化时间 $\tau=\pi/4$，$U=\exp((\pi/4)\gamma_2\gamma_3)$ 的本征相位与 Ising 任意子 $R^{\sigma\sigma}$ 的本征值一致（up to overall phase）。

把这个局域构件嵌入 2D：

- 在某个 Ising‑like 相区内（例如 $J_x\approx J_y>0,J_z$ 适中且加小磁场使体系进入非阿贝尔相），边界或涡旋端点上会出现一串空间分布的 Majorana 零模；
- 在这组零模中，选出一对邻近的 $\gamma_a,\gamma_b$，对应于两个待交换的缺陷或边界端点；
- 通过沿一条空间路径 $\gamma$ 的若干 $R_{ij}=e^{iK_{ij}}$ 操作，把“1D 式”的 $U_{ab}=\exp((\pi/4)\gamma_a\gamma_b)$ 在 2D 上实现为一个具体的 braid 路径；
- 然后在低能 Hilbert 空间里检查这个 $U_{ab}$ 的本征值与 Ising 的 R 矩阵是否一致——这正是 kit3 系列数值实验的核心任务之一。

简而言之：kit2‑exp 给出了“如何在 2D 上用 exp(iK) + Majorana+Z2 描述空间上的平行移动与拓扑缺陷的生成/移动”，而 kit3‑exp 则在此基础上，把这些空间路径视为多体配置空间中的曲线，计算其 holonomy（Dehn twist、half twist）。

