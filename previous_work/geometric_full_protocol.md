# 用完整状态量描述六段演化及其边界

本文回到完整的有向平面状态量 \(\mathcal B\)，把三段协议及其第二次重复接成一个初值问题。每一段都接收上一段的完整状态。所有分段公式都保留有限 \(E_1,|t_1|\) 的作用，不要求 \(E_1\tau\) 或 \(|t_1|\tau\) 很小。

基本定义见 [geometric_state_analysis.md](geometric_state_analysis.md)；本文统一 [相干方程](geometric_coherence_dynamics.md) 和 [第二段的完整分块](geometric_second_stage.md)，补充重复协议所需的第一段状态分量、第三段的完整传递，以及段间匹配条件。

## 1. 主变量、初态和测量量

按 \((1,2,3,a,b)\) 排列基向量，令

\[
L_{ij}=e_i e_j^T-e_j e_i^T,\qquad
\kappa=\frac2\hbar,\qquad E=E_1,\qquad \eta=|t_1|_{\mathrm{amplitude}}.
\]

研究论文选定的 Majorana 波函数

\[
\psi(0)=\frac{e_1-ie_2}{\sqrt2},\qquad
\psi_\star=\frac{e_1+ie_2}{\sqrt2}.
\]

沿用与其对应的五维实向量表示，而非把以下状态量解释为任意多体旋量态。定义

\[
\mathcal B=uv^T-vu^T,\qquad
P=-\mathcal B^2,\qquad
\rho=\frac12(P+i\mathcal B).
\]

完整状态满足

\[
\boxed{\dot{\mathcal B}=[K,\mathcal B],\qquad
\mathcal B(0)=L_{12}.}
\]

状态始终位于六维有向 Grassmann 流形

\[
\mathcal M=SO(5)/(SO(2)\times SO(3)),
\]

并保持 \(\mathcal B^T=-\mathcal B\)、\(\mathcal B^3=-\mathcal B\)、\(\operatorname{rank}\mathcal B=2\)、\(\|\mathcal B\|_F^2=2\)。用十个反对称矩阵分量记录状态可以避免在分段边界更换局部坐标；这十个量受到上述约束，并非十个独立物理自由度。

目标是 \(\mathcal B_\star=-L_{12}\)。把保真度作为状态的函数：

\[
\boxed{
F=\Phi(\mathcal B)=\frac14\left[
\sum_j(\mathcal B_{1j}^2+\mathcal B_{2j}^2)-2\mathcal B_{12}\right].}
\]

本文 \(F\) 是平方重叠；论文的重叠模长是 \(f=\sqrt F\)。定义

\[
d=\mathcal B_{12},\qquad
\ell=1-\frac12(P_{11}+P_{22}),
\]

则 \(F=(1-\ell-d)/2\)。这三个量用于读出轨迹，下一段的输入仍是完整 \(\mathcal B\)。

## 2. 把控制协议写清楚

按照 [model.md](../model.md) 的时间约定，每段时长 \(\tau\)，一轮为 \(3\tau\)，两轮为 \(6\tau\)。每段局部时间为 \(s\in[0,\tau]\)，定义

\[
a(s)=\frac{1-\cos(\pi s/\tau)}2,\qquad b(s)=1-a(s).
\]

记 \(p=|t_2(t)|\)、\(q=|t_1(t)|\)、\(r=|t_3(t)|\)、\(e_d=E_d(t)\)。采用论文中 \(t_1,t_2\) 随同一门控变化的协议：

| 段 | 第一轮时间 | \(p(s)\) | \(q(s)\) | \(r(s)\) | \(e_d(s)\) |
|---|---|---|---|---|---|
| 1 | \(0\to\tau\) | \(t_c a\) | \(\eta a\) | \(0\) | \(E_0 b\) |
| 2 | \(\tau\to2\tau\) | \(t_c b\) | \(\eta b\) | \(t_c a\) | \(0\) |
| 3 | \(2\tau\to3\tau\) | \(0\) | \(0\) | \(t_c b\) | \(E_0 a\) |

第二轮重复同一张控制表，时间整体平移 \(3\tau\)。各段共用

\[
\boxed{K=\kappa\left(EL_{12}+qL_{1b}-pL_{2a}+rL_{3a}+e_dL_{ab}\right).}
\]

下面第三段的守恒结论使用了 \(q=0\)。若另选始终开启的 \(t_1\)，仍可使用完整对易子方程，但需要重新分析第三段。

## 3. 在同一套状态坐标中写全部演化

定义复相干和辅助分量

\[
Z_j=\mathcal B_{1j}+i\mathcal B_{2j},\quad j=3,a,b,
\qquad h_{3a}=\mathcal B_{3a},\quad
h_{3b}=\mathcal B_{3b},\quad h_{ab}=\mathcal B_{ab}.
\]

\(d\)、三个复 \(Z_j\) 和三个实 \(h\) 给出全部十个分量。对易子方程等价于

\[
\begin{aligned}
\dot d&=-\kappa(p\operatorname{Re}Z_a+q\operatorname{Im}Z_b),\\
\dot Z_3&=\kappa(-iEZ_3+rZ_a+ip h_{3a}-q h_{3b}),\\
\dot Z_a&=\kappa(-iEZ_a-rZ_3+e_dZ_b+pd-qh_{ab}),\\
\dot Z_b&=\kappa(-iEZ_b-e_dZ_a+iqd-ip h_{ab}),\\
\dot h_{3a}&=\kappa(e_dh_{3b}-p\operatorname{Im}Z_3),\\
\dot h_{3b}&=\kappa(rh_{ab}-e_dh_{3a}+q\operatorname{Re}Z_3),\\
\dot h_{ab}&=\kappa(-rh_{3b}+q\operatorname{Re}Z_a+p\operatorname{Im}Z_b).
\end{aligned}
\]

保真度及其变化率为

\[
F=\frac{d^2-d}{2}+\frac{|Z_3|^2+|Z_a|^2+|Z_b|^2}{4},
\]

\[
\boxed{\begin{aligned}
\dot F=\frac\kappa2\big\{&
p[(1-d)\operatorname{Re}Z_a
  +h_{3a}\operatorname{Im}Z_3-h_{ab}\operatorname{Im}Z_b]\\
&+q[(1-d)\operatorname{Im}Z_b
  -h_{ab}\operatorname{Re}Z_a-h_{3b}\operatorname{Re}Z_3]\big\}.
\end{aligned}}
\]

这些式子已经区分了参数的作用：\(\eta\) 通过 \(q\) 参与改变瞬时重叠；\(E\) 通过 \(-iEZ_j\) 旋转相干，其自身不直接改变瞬时 \(F\)，却会改变随后进入 \(\dot F\) 的分量。\(r,e_d\) 的单独作用也保持 \(F\)，但会改变后续演化所需的状态。

### 3.1 第一段：两个旋转通道，加模式 3 的四个分量

第一段 \(r=0\)。定义

\[
X_\pm=d\pm h_{ab},\qquad
C_\pm=Z_a\pm iZ_b,\qquad
V_\pm=p\mp q,\qquad \Delta_\pm=E\pm e_d.
\]

两个精确通道为

\[
\boxed{\dot X_\pm=-\kappa V_\pm\operatorname{Re}C_\pm,
\qquad
\dot C_\pm=\kappa V_\pm X_\pm-i\kappa\Delta_\pm C_\pm.}
\]

每个通道的 \(X_\pm^2+|C_\pm|^2\) 守恒。参数通过 \(V_\pm=(t_c\mp\eta)a(s)\) 和 \(\Delta_\pm=E\pm E_0b(s)\) 改变通道旋转。

还需要四个实量

\[
W=(\operatorname{Re}Z_3,\operatorname{Im}Z_3,h_{3a},h_{3b})^T,
\]

其精确演化是

\[
\boxed{\dot W=\kappa
\begin{pmatrix}
0&E&0&-q\\-E&0&p&0\\0&-p&0&e_d\\q&0&-e_d&0
\end{pmatrix}W.}
\]

该矩阵反对称，因此 \(\|W\|^2=P_{33}\) 在第一段保持不变。它对应模式 3 的波函数权重 \(\rho_{33}=P_{33}/2\)。由于 \(\mathcal B\) 是简单双向量，还有

\[
\boxed{X_+^2+|C_+|^2=X_-^2+|C_-|^2=1-P_{33}.}
\]

证明只需四指标 Plücker 关系
\(d h_{ab}-\mathcal B_{1a}\mathcal B_{2b}
+\mathcal B_{1b}\mathcal B_{2a}=0\)，以及
\(\sum_{i<j}\mathcal B_{ij}^2=1\)。两个通道平方范数之差等于四倍上述关系；平均值加 \(\|W\|^2\) 等于 1。

**第一轮初始时**，\(W=0\)，两个通道都从 \((X_\pm,C_\pm)=(1,0)\) 出发，位于单位球面。

**第二轮第一段开始时**，一般 \(W\ne0\)，通道初值继承第一轮第三段的输出，球面半径是 \(\sqrt{1-P_{33}}\)。两通道加 \(W\) 才能完整传递状态；它们的初值仍须共同满足几何约束。

### 3.2 第二段：全部十个分量参与传递

第二段 \(e_d=0\)。定义

\[
X=\begin{pmatrix}\mathcal B_{1a}\\\mathcal B_{23}\\
\mathcal B_{2b}\\\mathcal B_{3b}\end{pmatrix},\qquad
Y=\begin{pmatrix}d\\\mathcal B_{13}\\\mathcal B_{1b}\\
-\mathcal B_{2a}\\-\mathcal B_{3a}\\\mathcal B_{ab}\end{pmatrix}.
\]

精确方程为

\[
\boxed{\dot Y=\kappa LX,\qquad \dot X=-\kappa L^T Y,\qquad
L=\begin{pmatrix}
-p&0&-q&0\\r&E&0&-q\\0&0&E&0\\
E&r&0&0\\0&p&0&0\\q&0&p&-r
\end{pmatrix}.}
\]

上一段输出的每个分量都可以转换到这里的 \((X,Y)\)，并继续演化。没有对第二段使用慢球面近似，也没有把输入替换成理想暗态。第二轮的第二段使用 \(\mathcal B(4\tau)\)，它包含此前四段的全部状态影响。

### 3.3 第三段：沿保真度等值面运动

第三段 \(p=q=0\)。按 \((1,2)\mid(3,a,b)\) 分块：

\[
\mathcal B=\begin{pmatrix}dJ&C\\-C^T&H\end{pmatrix},\qquad
J=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
D(s)=\begin{pmatrix}0&r&0\\-r&0&e_d\\0&-e_d&0\end{pmatrix}.
\]

这里 \(C\) 的两行分别是 \(Z=(Z_3,Z_a,Z_b)\) 的实部和虚部，\(H\) 包含三个 \(h\)。定义

\[
O(s)=e^{\kappa EJs},\qquad
\dot Q(s)=\kappa D(s)Q(s),\quad Q(0)=I_3.
\]

完整传递是

\[
\boxed{d(s)=d_{\mathrm{in}},\qquad
C(s)=O(s)C_{\mathrm{in}}Q(s)^T,\qquad
H(s)=Q(s)H_{\mathrm{in}}Q(s)^T.}
\]

等价地，

\[
Z(s)=e^{-i\kappa Es}Z_{\mathrm{in}}Q(s)^T.
\]

由于 \(\|C\|_F\) 不变，\(d,\ell,F\) 都保持不变。但 \(C,H\) 一般会改变：第三段沿目标重叠的等值面运动。该运动属于保持目标态射线的 \(SO(2)\times SO(3)\) 作用。

第一轮第三段的输出要用于下一轮，因此这里的相干旋转不可丢弃。第二轮第三段之后没有更多控制时，终点保真度已在 \(5\tau\) 确定，即 \(F(6\tau)=F(5\tau)\)；若还要继续第三轮，则仍需传递完整状态。

## 4. 每一段之间怎样匹配

### 4.1 物理状态直接连续传递

令第 \(j\) 段的局部状态为

\[
\mathcal B_j(s)=\mathcal B((j-1)\tau+s),\quad j=1,\ldots,6.
\]

唯一所需的状态接口条件是

\[
\boxed{\mathcal B_{j+1}(0)=\mathcal B_j(\tau).}
\]

切换控制不会重新准备波函数。固定基底中，\(\mathcal B\)、\(P\)、\(\rho\) 及所有相干分量都连续，不存在额外的边界旋转或投影。若计算时更换局部正交标架 \(S_j\)，坐标才需要转换：

\[
\widetilde{\mathcal B}_{j+1}(0)=
(S_{j+1}^TS_j)\widetilde{\mathcal B}_j(\tau)(S_{j+1}^TS_j)^T.
\]

这表示同一个物理状态，不是一次额外演化。

### 4.2 余弦门控使接口处足够光滑

三个控制接口的生成元为

| 边界 | 左右共同的 \(K/\kappa\) | 对状态的含义 |
|---|---|---|
| \(\tau,4\tau\) | \(EL_{12}+\eta L_{1b}-t_cL_{2a}\) | 第一段的完整输出成为第二段输入 |
| \(2\tau,5\tau\) | \(EL_{12}+t_cL_{3a}\) | \(p,q\) 已降为零，开始保真度守恒段 |
| \(3\tau\) | \(EL_{12}+E_0L_{ab}\) | 控制恢复起点值，状态继续传给下一轮 |

\(0\) 和 \(6\tau\) 的控制也为第三行所列值。控制恢复初始值并不要求状态恢复初始值。

因为 \(a',b'\) 在两端都为零，各接口满足

\[
K(t_j^-)=K(t_j^+),\qquad
\dot K(t_j^-)=\dot K(t_j^+)=0.
\]

利用

\[
\dot{\mathcal B}=[K,\mathcal B],\qquad
\ddot{\mathcal B}=[\dot K,\mathcal B]+[K,[K,\mathcal B]],
\]

得到 \(\mathcal B\) 及其一、二阶时间导数连续，即全轨迹至少为 \(C^2\)。控制的二阶导数一般会跳变，所以不能一概要求三阶状态导数连续。

在进入第三段的 \(2\tau,5\tau\)，\(\dot d=\dot F=0\)，其二阶时间导数也为零；这并不意味着 \(\dot{\mathcal B}=0\)。在 \(3\tau\) 开始下一轮时，\(p,q\) 由零以二次阶开启，\(d,F\) 的初始变化一般从局部时间的三次阶开始。状态相干可以同时以非零速度旋转。

### 4.3 控制接口不是几何瓶颈

时间接口是 \(t=j\tau\)，由实验控制预先指定。几何瓶颈是状态空间中的集合

\[
\Sigma=\{\mathcal B\in\mathcal M:\mathcal B_{12}=0\}.
\]

此前已证明：在 \(\Sigma\) 上，\(\ell\ge1/2\)、\(F\le1/4\)。有限参数时轨迹可以在段内穿过 \(\Sigma\)，也可以在某段多次穿过；段间匹配并不要求 \(d=0\)。

穿过瓶颈后朝哪一侧继续，由当时的完整状态决定。其法向速度是

\[
\left.\dot d\right|_\Sigma=-\kappa\left(p\operatorname{Re}Z_a+q\operatorname{Im}Z_b\right).
\]

右侧为负，轨迹局部进入 \(d<0\)；为正，则进入 \(d>0\)。等于零时可能是切触，需检查更高阶导数。第三段本身保持 \(d\)，不能完成从正侧到负侧的穿越。

这给出了可供几何分析的两个不同对象：时间接口上的状态传递，以及状态流形中瓶颈集合的穿越。

## 5. 全过程是分段流的复合

令 \(\mathscr T_j(s;E,\eta)\) 是第 \(j=1,2,3\) 种控制从局部时间 0 到 \(s\) 对 \(\mathcal B\) 的精确作用，并简记 \(\mathscr T_j=\mathscr T_j(\tau)\)。它们都是对反对称矩阵空间的线性作用，限制到 \(\mathcal M\) 后保持合法状态。

一轮作用为

\[
\mathscr C=\mathscr T_3\circ\mathscr T_2\circ\mathscr T_1.
\]

两轮的边界状态满足

\[
\begin{aligned}
\mathcal B_0&=L_{12},&
\mathcal B_1&=\mathscr T_1\mathcal B_0,&
\mathcal B_2&=\mathscr T_2\mathcal B_1,\\
\mathcal B_3&=\mathscr T_3\mathcal B_2,&
\mathcal B_4&=\mathscr T_1\mathcal B_3,&
\mathcal B_5&=\mathscr T_2\mathcal B_4,\\
\mathcal B_6&=\mathscr T_3\mathcal B_5=\mathscr C^2\mathcal B_0.
\end{aligned}
\]

这里 \(\mathcal B_j=\mathcal B(j\tau)\)，下标表示边界状态。最终参数依赖为

\[
\boxed{F(E,\eta,6\tau)=\Phi\!\left(\mathscr C(E,\eta)^2L_{12}\right).}
\]

这给出了任意有限参数下的完整结构；一般仍需计算分段流，尚未给出所有参数下的初等闭式解。

本问题是给定初态的演化问题。\(\mathcal B_\star=-L_{12}\) 是用于评价结果的目标，不能另外强行设为终点边界条件。哈密顿演化也没有使 \(d\approx\pm1\) 成为耗散吸引子的机制；我们研究的是有限时间内的到达及返回。

### 5.1 参数影响也在各接口完整传递

令 \(\lambda=E,\eta\)，\(S_\lambda=\partial_\lambda\mathcal B\)。对精确方程求导：

\[
\boxed{\dot S_\lambda=[K,S_\lambda]+[\partial_\lambda K,\mathcal B],
\qquad S_\lambda(0)=0,}
\]

\[
\partial_E K=\kappa L_{12},\qquad
\partial_\eta K=\kappa g(t)L_{1b},
\]

其中 \(g=a,b,0\) 分别对应三段。接口处

\[
\boxed{S_{\lambda,j+1}(0)=S_{\lambda,j}(\tau).}
\]

因此每一段的参数影响都有两个来源：输入态已经随参数改变，以及当前段的生成元随参数改变。用分段流表示即

\[
\partial_\lambda\mathcal B_j=
(\partial_\lambda\mathscr T_j)\mathcal B_{j-1}
+\mathscr T_j\partial_\lambda\mathcal B_{j-1},
\]

对第二轮使用对应重复的 \(\mathscr T_j\)。这些导数是在任意参数点的精确灵敏度方程，求它们不需要先假设 \(E\tau,\eta\tau\ll1\)。

最终保真度导数可以直接读出：

\[
\partial_\lambda F=\frac14\left[
-\operatorname{tr}\{P_0(\mathcal B S_\lambda+S_\lambda\mathcal B)\}
+\operatorname{tr}(L_{12}S_\lambda)\right],
\qquad P_0=-L_{12}^2.
\]

第三段虽然保持 \(F\)，仍保留输入态的参数依赖，并产生与 \(E\) 有关的相干相位。第一轮的这些影响会经下一轮传递到最终保真度。

## 6. 分段拼接的独立检验

验证脚本：[verify_full_protocol.py](geometric_validation/verify_full_protocol.py)。结果：[full_protocol.json](geometric_validation/full_protocol.json)。

验证只演化状态量：第一段用两个通道和四个 \(W\) 分量，第二段用 \(4+6\) 分块，第三段用上面的精确旋转公式。六段完整传递输入，与此前独立积分并保存的复波函数轨迹比较。没有重新拟合参数或使用慢球面近似。

五个已有案例覆盖 \(\tau=100,500\)、\(E=0,0.01,0.02\) 和不同 \(\eta\)。全部通过；所有记录误差的最大值小于 \(8\times10^{-12}\)。检验包括完整 \(\mathcal B\)、保真度、简单平面约束、第一段不变量和接口的一、二阶状态导数。符号计算另核对了四分量子系统、通道半径恒等式、\(\dot F\) 及第三段分块。

以 \(\hbar=1\)、\(t_c=E_0=0.3\,\mathrm{meV}\)、\(\tau=100/\mathrm{meV}\)、\(E=0.01\,\mathrm{meV}\)、\(\eta=0.023\,\mathrm{meV}\) 为例：

| 边界 \(t/\tau\) | \(d\) | \(\ell\) | \(F\) |
|---|---:|---:|---:|
| 0 | 1.000000 | 0.000000 | 0.000000 |
| 1 | 0.001066 | 0.500520 | 0.249207 |
| 2 | -0.237678 | 0.471755 | 0.382962 |
| 3 | -0.237678 | 0.471755 | 0.382962 |
| 4 | -0.000054 | 0.761421 | 0.119317 |
| 5 | -0.886132 | 0.107385 | 0.889374 |
| 6 | -0.886132 | 0.107385 | 0.889374 |

这个例子揭示了三个接口机制：

1. 第一段终点接近 \(d=0\)，但不是严格落在瓶颈上；第二轮的第一段输出也不是理想输入。
2. 第一轮第三段保持 \(F=0.382962\)，但 \(\|\mathcal B(3\tau)-\mathcal B(2\tau)\|_F\approx1.94170\)，完整状态显著改变。
3. 第二轮第一段开始时，\(P_{33}\approx0.943509\)，两通道球面半径仅约 \(0.23768\)。保留 \(W\) 对整条轨迹的描述是必要的。

第三段的状态变化与保真度守恒并不矛盾：它沿等值面运动，为下一轮准备了不同的相干输入。

## 7. 在这个完整框架上继续研究什么

这条链的后续分析已整理到 [geometric_second_cycle_chain.md](geometric_second_cycle_chain.md)：第一段的精确几何上界、反向传回的所需输入平面、最终换号边界，以及参数响应在前期输入和随后两段之间的分解，并完成独立波函数检验。

这项研究集中在 **\(\mathcal B(3\tau)\to\mathcal B(4\tau)\to\mathcal B(5\tau)\)** 这一条链：第一轮结束时的相干，怎样在第二轮第一段被旋转，再怎样决定第二段的定向穿越和保真度增长。

\(X_\pm,C_\pm,W\) 在 \(3\tau\to4\tau\) 的旋转给出第一段不变量和能力上界。进一步把最终目标反向传回接口，就能用完整输入与所需输入的对齐关系精确评价终点，并用 \(S_E,S_\eta\) 和反向梯度区分参数改变输入和改变段内演化的贡献。

下一项研究可在有限参数区域内追踪实际输入与所需输入的相对定向，分析最终换号曲线、高保真度等值线及接口球面解释的适用范围。
