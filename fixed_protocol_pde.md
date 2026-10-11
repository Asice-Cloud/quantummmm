# 固定协议下的几何状态 PDE：仅含不可控参数 E₁、t₁

本轮保留每段时长 \(\tau>0\)、门控编织幅度 \(t_c\) 和量子点能级幅度 \(E_0\) 的符号，把它们作为已知的可控设计量。门控的归一化形状和两轮六段协议保持既定设置。只有 \(E=E_1\) 与 \(\eta=|t_1|\) 的幅度不可控，并在一次演化中保持为定值。对每组给定的 \((\tau,t_c,E_0)\)，研究对象是参数—时间空间中的完整状态场

\[
\mathcal B:\Omega\times[0,6\tau]\longrightarrow
SO(5)/(SO(2)\times SO(3)),\qquad
(E,\eta,t)\mapsto\mathcal B(E,\eta;t\mid\tau,t_c,E_0).
\]

其中 \(\Omega\subset\mathbb R\times[0,\infty)\) 是待分析的参数区域。下文省略可控设计量的条件记号；所有 \(E,\eta\) 偏导都在 \(\tau,t_c,E_0\) 给定时计算。调整 \(\tau\) 会按比例拉伸门控时间轴，归一化门控形状不变。

本文记录时间演化、参数响应及其几何表示。带辅助连接场的相容系统保留完整状态，没有对 \(E\tau\)、\(\eta\tau\) 或状态角度作小量截断；但连接场仍是未知量，因此不能将它称为已经得到系数直接已知的独立参数平面 PDE。具体分量系数与混合偏导方程现已展开，见 [explicit_parameter_equations.md](fixed_protocol_pde/explicit_parameter_equations.md)。只含终点保真度 \(F_{\rm end}(E,\eta)\) 的闭合参数 PDE 尚未得到。

此前研究统一保存在 [previous_work](previous_work/README.md)。本轮从完整状态重新建立方程；旧的慢旋转近似及其数值误差是后续比较资料。

当前采用的具体路线见 [显式方程第 5 节](fixed_protocol_pde/explicit_parameter_equations.md#5-初值接口与保真度)：保留 \(\mathcal B(E,\eta,t)\) 的十个全局分量，把状态及两个参数方向的响应组成 30 分量初值问题，再通过精确的保真度变化率读取其影响。以下连接场表示作为几何补充，具体求解不要求引入连接场。

## 1. 把固定的模型系数写成已知函数

沿用 [model.md](model.md) 中 \((1,2,3,a,b)\) 的基底顺序，令

\[
L_{ij}=e_i e_j^T-e_j e_i^T,\qquad \hbar=1.
\]

每段局部时间 \(s\in[0,\tau]\)，既定包络为

\[
a(s)=\frac{1-\cos(\pi s/\tau)}2,\qquad b(s)=1-a(s).
\]

设 \(p=|t_2|\)、\(r=|t_3|\)、\(e_d=E_d\)，实际 \(t_1\) 耦合为 \(q=\eta g\)。已知控制系数为：

| 段类型 | 全局时间区间 | \(p\) | \(r\) | \(e_d\) | \(g\) |
|---|---|---|---|---|---|
| 1 | \([0,\tau]\)、\([3\tau,4\tau]\) | \(t_ca\) | \(0\) | \(E_0b\) | \(a\) |
| 2 | \([\tau,2\tau]\)、\([4\tau,5\tau]\) | \(t_cb\) | \(t_ca\) | \(0\) | \(b\) |
| 3 | \([2\tau,3\tau]\)、\([5\tau,6\tau]\) | \(0\) | \(t_cb\) | \(E_0a\) | \(0\) |

表中 \(s\) 始终从所在段的起点计时。采用与此前波函数相同的五维实向量表示，生成元为

\[
\boxed{
K(E,\eta,t)=K_c(t)+E K_E+\eta K_\eta(t),}
\]

\[
K_c=2(-pL_{2a}+rL_{3a}+e_dL_{ab}),\qquad
K_E=2L_{12},\qquad K_\eta(t)=2g(t)L_{1b}.
\]

因此 \(\partial_EK=K_E\)、\(\partial_\eta K=K_\eta\)，所有二阶不可控参数导数为零。这里 \(t\) 是物理演化坐标，\(\tau\) 是给定的可控段时长。

若采用每段归一化时间 \(v=s/\tau\)，则 \(\partial_v=\tau\partial_t\)，所以状态方程变为 \(\partial_v\mathcal B=[\tau K,\mathcal B]\)。物理时间生成元的因子仍为 \(2\)，归一化时间生成元的因子为 \(2\tau\)。

## 2. 完整状态及其时间方程

研究初态和目标态

\[
\psi_0=\frac{e_1-ie_2}{\sqrt2},\qquad
\psi_\star=\frac{e_1+ie_2}{\sqrt2}.
\]

对 \(\psi=(u-iv)/\sqrt2\)，其中 \(u,v\) 是实正交单位向量，定义

\[
\mathcal B=uv^T-vu^T,\qquad P=-\mathcal B^2,\qquad
\rho=\frac{P+i\mathcal B}{2}.
\]

\(\mathcal B\) 保留该波函数除整体相位之外的全部信息；本表示针对上述选定波函数，不把它解释为任意多体旋量态的完整描述。

由 \(\partial_t\psi=K\psi\) 或 \(u_t=Ku,v_t=Kv\)，直接得到

\[
\boxed{\partial_t\mathcal B=[K,\mathcal B],\qquad
\mathcal B(E,\eta;0)=\mathcal B_0=L_{12}.}
\tag{1}
\]

方程保持

\[
\mathcal B^T=-\mathcal B,\quad
\mathcal B^3=-\mathcal B,\quad
\operatorname{rank}\mathcal B=2,\quad
\|\mathcal B\|_F^2=2.
\]

因此它描述的是六维有向平面流形上的运动。十个反对称矩阵分量是一套带约束的全局表示。

仅式 (1) 本身是一个带二维参数的时间 ODE 族，也可视作参数—时间场上的退化 PDE；它没有说明参数偏导如何计算。下面补上这些方向的精确关系。

## 3. 先直接推导参数偏导方程

记

\[
S_E=\partial_E\mathcal B,\qquad
S_\eta=\partial_\eta\mathcal B.
\]

对式 (1) 求偏导，利用乘积法则：

\[
\boxed{
\begin{aligned}
\partial_t S_E&=[K,S_E]+[K_E,\mathcal B],\\
\partial_t S_\eta&=[K,S_\eta]+[K_\eta,\mathcal B],
\end{aligned}\qquad
S_E(E,\eta;0)=S_\eta(E,\eta;0)=0.}
\tag{2}
\]

这与式 (1) 构成闭合系统。等价地，它给出混合偏导 PDE

\[
\partial_{tE}\mathcal B=[K,\partial_E\mathcal B]+[K_E,\mathcal B],
\qquad
\partial_{t\eta}\mathcal B=[K,\partial_\eta\mathcal B]+[K_\eta,\mathcal B].
\]

这些导数在任意有限参数点定义，不是零参数展开。\(S_E,S_\eta\) 是流形切向量，例如对几何约束求导可得

\[
\mathcal B^2S_\lambda+\mathcal BS_\lambda\mathcal B
+S_\lambda\mathcal B^2=-S_\lambda.
\]

二阶导数 \(H_{\lambda\mu}=\partial_\lambda\partial_\mu\mathcal B\) 同样闭合：

\[
\boxed{
\partial_tH_{\lambda\mu}
=[K,H_{\lambda\mu}]+[K_\lambda,S_\mu]+[K_\mu,S_\lambda],
\qquad H_{\lambda\mu}(0)=0.}
\tag{3}
\]

当 \(\lambda=\mu\) 时，后两项合并为 \(2[K_\lambda,S_\lambda]\)。因为 \(K\) 对两个参数是仿射的，没有遗漏 \(K_{\lambda\mu}\) 项。

## 4. 写成三个方向统一的一阶几何 PDE

为了把参数方向也写成与时间方向相同的几何形式，引入两个反对称辅助场 \(\Xi_E,\Xi_\eta\)，由以下初值问题定义：

\[
\boxed{
\begin{aligned}
\partial_t\Xi_E&=[K,\Xi_E]+K_E,\\
\partial_t\Xi_\eta&=[K,\Xi_\eta]+K_\eta,
\end{aligned}\qquad
\Xi_E(E,\eta;0)=\Xi_\eta(E,\eta;0)=0.}
\tag{4}
\]

这些方程只含已知 \(K,K_E,K_\eta\)，并由初值唯一确定。辅助场描述参数改变在累计传播中形成的无穷小旋转；它们不是额外的物理控制。

现在证明

\[
\boxed{\partial_E\mathcal B=[\Xi_E,\mathcal B],\qquad
\partial_\eta\mathcal B=[\Xi_\eta,\mathcal B].}
\tag{5}
\]

对任一 \(\lambda\)，令

\[
D_\lambda=\partial_\lambda\mathcal B-[\Xi_\lambda,\mathcal B].
\]

利用 Jacobi 恒等式以及式 (1)、(4)，有

\[
\partial_t[\Xi_\lambda,\mathcal B]
=[K,[\Xi_\lambda,\mathcal B]]+[K_\lambda,\mathcal B].
\]

与式 (2) 相减，得到 \(\partial_tD_\lambda=[K,D_\lambda]\)、\(D_\lambda(0)=0\)。唯一性给出 \(D_\lambda\equiv0\)，从而证明式 (5)，无需把未知传播矩阵 \(R\) 留在方程右边。

因此主系统可以统一写为

\[
\boxed{
\begin{aligned}
\partial_t\mathcal B&=[K,\mathcal B],\\
\partial_E\mathcal B&=[\Xi_E,\mathcal B],\\
\partial_\eta\mathcal B&=[\Xi_\eta,\mathcal B],\\
\partial_t\Xi_E&=[K,\Xi_E]+2L_{12},\\
\partial_t\Xi_\eta&=[K,\Xi_\eta]+2g(t)L_{1b}.
\end{aligned}}
\tag{6}
\]

式 (6) 是带辅助场的精确相容系统，主未知量为 \(\mathcal B,\Xi_E,\Xi_\eta\)。它通过时间初值问题闭合，参数方向的系数 \(\Xi_E,\Xi_\eta\) 仍待求解，不能视为已消去时间历史的参数平面 PDE。若需要一、二阶响应，也可使用等价的式 (1)—(3)；不必同时把所有辅助量都作为独立主未知量。

### 4.1 参数方向的相容条件

参数偏导必须满足 \(\partial_E\partial_\eta\mathcal B=\partial_\eta\partial_E\mathcal B\)。对辅助场，定义

\[
\mathcal C=\partial_E\Xi_\eta-\partial_\eta\Xi_E-[\Xi_E,\Xi_\eta].
\]

对式 (4) 求导并使用 Jacobi 恒等式可得

\[
\partial_t\mathcal C=[K,\mathcal C],\qquad\mathcal C(0)=0.
\]

所以

\[
\boxed{\partial_E\Xi_\eta-\partial_\eta\Xi_E
-[\Xi_E,\Xi_\eta]=0.}
\tag{7}
\]

令 \(\Xi_t=K\)，则三对方向统一满足

\[
\partial_\mu\Xi_\nu-\partial_\nu\Xi_\mu
-[\Xi_\mu,\Xi_\nu]=0,\qquad \mu,\nu\in\{t,E,\eta\}.
\]

等价地，协变导数 \(\mathcal D_\mu=\partial_\mu-\operatorname{ad}_{\Xi_\mu}\) 满足 \(\mathcal D_\mu\mathcal B=0\) 且相互对易。

这叫辅助连接的零曲率条件。它表示参数与时间的微分描述相容，**不表示两个参数没有相互作用**。\([\Xi_E,\Xi_\eta]\) 一般非零，由连接场的参数偏导补偿；参数曲面在状态流形中的诱导度量也不因此成为平坦度量。

### 4.2 不可控参数保持固定意味着什么

实际时间轨迹为 \(t\mapsto(E,\eta,\mathcal B(E,\eta;t))\)，其中 \(\dot E=\dot\eta=0\)。式 (5) 比较不同固定参数的实验轨迹，不给参数引入时间漂移。

因此模型没有产生参数空间中的扩散项，也不需要假设噪声分布。式 (6) 的几何关系来自同一个确定性初值问题。它提供可分析的参数偏导和相容结构，但不能被理解为已经消除了时间求解的工作。

## 5. 初值、段间接口与参数边界

整个系统的初始数据在 \(t=0\) 给出：

\[
\mathcal B(E,\eta;0)=L_{12},\qquad
\Xi_E(E,\eta;0)=\Xi_\eta(E,\eta;0)=0.
\]

若使用响应方程，则 \(S_\lambda(0)=0,H_{\lambda\mu}(0)=0\)。这些零初值来自所有参数共用同一初态。

在 \(t_j=j\tau\)，\(j=1,\ldots,5\)，下一段接收上一段的完整输出：

\[
\boxed{
\mathcal B(t_j^+)=\mathcal B(t_j^-),\qquad
\Xi_\lambda(t_j^+)=\Xi_\lambda(t_j^-).}
\tag{8}
\]

同样有 \(S_\lambda(t_j^+)=S_\lambda(t_j^-)\)、\(H_{\lambda\mu}(t_j^+)=H_{\lambda\mu}(t_j^-)\)。每段只更换已知系数函数的表达式，不能重置状态或辅助场。

现有门控在接口处数值匹配、一阶时间导数为零，所以 \(K,K_E,K_\eta\) 连续且一阶时间导数连续；主状态及连接场具有相应的 \(C^2\) 时间正则性。不存在需要人为补入的跳跃项。

参数区域边缘只界定要查看的范围。对于由时间初值问题确定的场，不额外指定任意的 \(E\) 或 \(\eta\) 边界值；否则可能与式 (6) 不相容。\(\eta=0\) 的物理偏导可取右侧导数，系数本身在该处正常。

对给定的 \(\tau\)，终点为 \(T=6\tau\)。目标态只用于读出保真度，**不要求** \(\mathcal B(E,\eta;6\tau)\) 等于目标；否则就预先强制了希望研究的结果。

## 6. 从 PDE 读出保真度、梯度和参数几何

令

\[
\mathcal B_\star=-L_{12},\qquad P_\star=-\mathcal B_\star^2.
\]

保真度作为完整状态的函数为

\[
\boxed{
F(E,\eta;t)=\frac14\left[
\operatorname{tr}(P_\star P)-\operatorname{tr}(\mathcal B_\star\mathcal B)
\right],\qquad F_{\rm end}(E,\eta\mid\tau,t_c,E_0)=F(E,\eta;6\tau\mid\tau,t_c,E_0).}
\tag{9}
\]

本文 \(F\) 是平方重叠，论文的模长为 \(f=\sqrt F\)。等价地，

\[
F=\frac14\left[\sum_j(\mathcal B_{1j}^2+\mathcal B_{2j}^2)-2\mathcal B_{12}\right].
\]

对 \(\mu=t,E,\eta\)，已有

\[
\partial_\mu P=[\Xi_\mu,P],\qquad
\partial_\mu\mathcal B=[\Xi_\mu,\mathcal B].
\]

所以三个方向的保真度导数统一为

\[
\boxed{
\partial_\mu F=\frac14\left[
\operatorname{tr}(P_\star[\Xi_\mu,P])
-\operatorname{tr}(\mathcal B_\star[\Xi_\mu,\mathcal B])
\right].}
\tag{10}
\]

这里右边依赖完整状态及连接，不能未经消元证明就改写成仅依赖 \(F\) 的一阶标量方程。第三段的 \(K\) 与目标态对易，因而 \(\partial_tF=0\)，但 \(\mathcal B\) 仍然运动并影响下一轮；这说明标量保真度没有包含全部演化信息。

### 6.1 任意有限参数处的梯度和 Hessian

用 \(S_\lambda\) 表示响应时，

\[
P_\lambda=-(\mathcal BS_\lambda+S_\lambda\mathcal B),
\]

\[
P_{\lambda\mu}=-\left(
\mathcal BH_{\lambda\mu}+H_{\lambda\mu}\mathcal B
+S_\lambda S_\mu+S_\mu S_\lambda\right).
\]

于是

\[
F_\lambda=\frac14\left[\operatorname{tr}(P_\star P_\lambda)
-\operatorname{tr}(\mathcal B_\star S_\lambda)\right],
\]

\[
F_{\lambda\mu}=\frac14\left[\operatorname{tr}(P_\star P_{\lambda\mu})
-\operatorname{tr}(\mathcal B_\star H_{\lambda\mu})\right].
\]

这提供终点高保真等值线 \(F_{\rm end}=F_0\) 的法向 \((F_E,F_\eta)\)。当 \(F_\eta\ne0\)，其局部斜率为 \(d\eta/dE=-F_E/F_\eta\)。固定 \(E\) 的内部峰由 \(F_\eta=0,F_{\eta\eta}<0\) 判别。它们都是完整有限参数场的导数条件。

### 6.2 参数平面在状态流形中的诱导度量

在参数方向定义纯态的 Fubini–Study 度量

\[
\boxed{
G_{\lambda\mu}=\frac12\operatorname{tr}(\rho_\lambda\rho_\mu)
=\frac18\left[
\langle P_\lambda,P_\mu\rangle_F
+\langle S_\lambda,S_\mu\rangle_F\right].}
\tag{11}
\]

这里 \(\langle X,Y\rangle_F=\operatorname{tr}(X^TY)\)。\(G\) 半正定；若参数映射在某处降秩，它可以退化。其特征方向给出完整物理态对参数变化的敏感方向，仍需结合 \(\nabla F\) 才能判断保真度变化。

由式 (2) 及 \(P_t=[K,P]\)，可以进一步得到精确的度量演化式

\[
\begin{aligned}
\partial_tG_{\lambda\mu}=\frac18\big\{&
\langle[K_\lambda,P],P_\mu\rangle_F
+\langle P_\lambda,[K_\mu,P]\rangle_F\\
&+\langle[K_\lambda,\mathcal B],S_\mu\rangle_F
+\langle S_\lambda,[K_\mu,\mathcal B]\rangle_F\big\}.
\end{aligned}
\tag{12}
\]

纯粹由共同 \(K\) 引起的旋转项在内积中抵消；参数敏感性由不同参数产生的驱动差异建立。这使几何量和状态、响应一起闭合，不需要单独拟合一个参数空间度量。

## 7. 已完成的独立核对

此前数值核对取 \(\tau=100,t_c=E_0=0.3\)，作为上述符号推导的一个具体设计实例。这些取值只用于核对，不是推导成立的条件。代码在 \((E,\eta)=(0,0),(0.01,0.023),(0.027,0.001)\) 上传播六段，分别求解：

- 完整 \(\mathcal B\)、连接场、一二阶响应及连接的混合参数导数；
- 独立复波函数分量方程，以及它的一二阶参数导数。

两套实现的比较覆盖每段 101 个采样时刻。归一化残差采用固定单位下的 \(|A-B|/(1+|A|+|B|)\)；连接曲率的分母相应包含三个被相减项。

| 核对项 | 三个样本、全过程的最大误差 |
|---|---:|
| \(\mathcal B\) 与波函数重构 | \(2.99\times10^{-12}\) |
| 保真度与直接波函数重叠 | \(1.16\times10^{-12}\) |
| \(S_\lambda=[\Xi_\lambda,\mathcal B]\)，绝对残差 | \(1.09\times10^{-9}\) |
| 同一参数 PDE，归一化残差 | \(8.66\times10^{-10}\) |
| 连接零曲率，绝对残差 | \(1.54\times10^{-7}\) |
| 连接零曲率，归一化残差 | \(1.99\times10^{-9}\) |
| 保真度参数梯度，归一化残差 | \(8.06\times10^{-12}\) |
| 保真度 Hessian，归一化残差 | \(4.29\times10^{-9}\) |
| 参数度量，归一化残差 | \(1.95\times10^{-9}\) |

较大的连接绝对残差对应含很大参数导数的量，表中同时列出其归一化误差。完整矩阵二阶响应与波函数重构的最大归一化残差为 \(6.10\times10^{-8}\)。这些是样本上的数值一致性检验；精确性来自前述推导，不把三个样本当作整个空间的数值误差界。

例如 \(E=0.01,\eta=0.023\)，终点得到

\[
F_{\rm end}=0.8893735894738,\qquad
(F_E,F_\eta)=(-36.52301185,\ 0.01670052).
\]

同一点的两个参数方向具有很不相同的保真度梯度；这由全程响应给出，而不是通过修改门控或时长得到。

运行：

```sh
.venv/bin/python -B fixed_protocol_pde/verify_pde.py
```

代码：[verify_pde.py](fixed_protocol_pde/verify_pde.py)；结果：[verification.json](fixed_protocol_pde/verification.json)。该数值实例保留既有设置；符号系数核对见 [derive_explicit_equations.py](fixed_protocol_pde/derive_explicit_equations.py)。

## 8. 从完整状态到峰分支

以式 (6) 及式 (9)—(10) 为基础，对给定的可控设计 \((\tau,t_c,E_0)\)，在终点 \(6\tau\) 研究高保真区域 \(F_{\rm end}(E,\eta\mid\tau,t_c,E_0)\geq F_0\) 的边界和敏感方向。先解释条带的横向宽度如何由两个不可控参数的共同作用决定，再与归档中的完整参数图对照。推导保留可控量的符号；不可控参数的响应针对 \(E,\eta\)。

本轮进一步把状态方程与峰的驻点条件、分支延拓及几何高度联系起来，见 [peak_structure.md](fixed_protocol_pde/peak_structure.md)。为比较论文 Fig. 1(d)，另对给定的可控时长 \(\tau\) 求设计响应；这不增加不可控参数。新文档区分完整状态的精确预测、慢角度的近似解释及尚未对齐的论文作图约定。

图注参数现已统一为 \(E_1=0.01\,\mathrm{meV}\)，正文的 \(0.001\) 按用户确认的笔误处理。本次从原 PDF 提取色带并进行定量比较，见 [paper_alignment.md](fixed_protocol_pde/paper_alignment.md)。完整状态与独立波函数的对应已通过核对；当前 \(K=2A\) 在图注参数下仍未复现论文亮带，整体减半的候选也有累计偏移。原图网格的插值检验不足以解释低耦合区的差异，因此没有据此修改当前模型。

进一步从原始有效哈密顿量独立构造四元数旋量矩阵和费米子 BdG 矩阵，已严格验证五维生成元就是 \(K=2A\)，且当前十个分量方程的系数全部一致，见 [effective_hamiltonian_audit.md](fixed_protocol_pde/effective_hamiltonian_audit.md)。`model.md` 中三个旋量耦合标签已修正，PDE 实现本身采用的符号正确。旋量、BdG 与分段 PDE 的完整演化核对也通过；原图对齐的剩余问题不能通过任意更改这个已确定的因子来解决。

按该 PDE 绘制的完整状态与保真度图见 [state_evolution/README.md](fixed_protocol_pde/state_evolution/README.md)：包括 \(\tau=100\) 时的全程十分量状态、\(F\) 与论文重叠模长 \(\sqrt F\)，以及沿用 Fig. 1(d) 坐标、参数范围和调色板的终点参数图。该图的 \(E_1=0.01\)，没有使用此前的经验系数校准。
