# PRB 111 模型的微扰、几何结构与对称性抵消

本文保留两条已完成的推导：由辅助快模消去得到几何输运与有效动力学；由该动力学得到原文读出的二阶、四阶公式，以及保留编织门的群对称性抵消条件。一般累计相位下的收敛公式、解析延拓、图像验证及误差预算统一见 [解析延拓与保真度](analytic_continuation.md)。

**核心结论：**可以得到具体的微扰保真度公式；存在使指定读出的二阶损失消失的比例，但它不等于整个门无误差。利用编织共轭和严格反向路径，可以在一阶投影模型内同时抵消门误差的一阶、二阶项。原模型一般含低能劈裂，不能把其演化直接称为纯几何门。

## 1. 模型、协议与读出约定

来源为 [原始模型记录](../model.md)、[PRB 111 论文](../arts/prb111.pdf)。几何分析思路来自 [2008 年论文](../arts/2008.pdf)，其组会解读保留在 [文献笔记](literature/2008.md)。

取 Majorana 归一化 $\{\gamma_i,\gamma_j\}=2\delta_{ij}$，固定原协议中的相位与耦合符号：

$$
H=i\gamma_a(d\gamma_b+u_2\gamma_2-u_3\gamma_3)
+i\gamma_1(u_1\gamma_b+e\gamma_2),
\qquad e=E_1>0,\quad \eta=|t_1|_{\max}.
$$

辅助控制等幅 $E_0=t_c=C$。每段持续 $\tau$，局部时间 $s\in[0,1]$，定义

$$
a=\frac{1-\cos\pi s}{2},\quad b=1-a,\quad
q=\sqrt{a^2+b^2},\quad g=\frac a q.
$$

| 段 | $d$ | $u_2$ | $u_3$ | $u_1$ |
|---|---:|---:|---:|---:|
| I | $Cb$ | $Ca$ | $0$ | $\eta a$ |
| II | $0$ | $Cb$ | $Ca$ | $\eta b$ |
| III | $Ca$ | $0$ | $Cb$ | $0$ |

这里 $u_1(t)$ 随门控开关；不是全程恒定的 $t_1$。两圈连续执行，总时长 $6\tau$，中间不重置状态。

小参数有两种，不能混用：

$$
\epsilon=\frac{\sqrt{e^2+\eta^2}}C,\qquad
L=\frac{C\tau}{\hbar};\qquad
r=\frac\eta e,\quad \lambda=\frac{e\tau}{\hbar},\quad
\nu=\frac{\eta\tau}{\hbar}=r\lambda.
$$

辅助模式投影要求 $\epsilon\ll1,L\gg1$；本篇的累计角展开还要求 $|\lambda|\sqrt{1+r^2}\ll1$。例如可以在 $1\ll L\ll1/\epsilon$ 的窗口同时成立。另一篇的完整标量函数不要求累计角小。

统一记原文波函数重叠为

$$
f=|\langle\psi_1^+(0)|\psi_1^-(6\tau)\rangle|,
\qquad F=f^2.
$$

它不是平均门保真度。完整单粒子传播采用 $K=2A/\hbar$，$\dot{\mathcal R}=K\mathcal R$。在 $(1,2,3,a,b)$ 顺序下，令 $\zeta=(e_1-ie_2)/\sqrt2$，则

$$
f=|\zeta^T\mathcal R\zeta|,
\quad
F=\frac{(\mathcal R_{11}-\mathcal R_{22})^2+
(\mathcal R_{12}+\mathcal R_{21})^2}{4}.
$$

这里是转置而非厄米共轭，因为目标态是初态的粒子—空穴伙伴。

## 2. 精确能谱与“低能简并”条件

令 $\boldsymbol\gamma=(\gamma_b,\gamma_2,\gamma_3)$，

$$
\boldsymbol v=(d,u_2,-u_3),\quad
\boldsymbol w=(u_1,e,0),\quad
R=|\boldsymbol v|,\ W=|\boldsymbol w|,\ p=\boldsymbol v\cdot\boldsymbol w.
$$

二次型的非零奇异值由两行矩阵 $M=(\boldsymbol v;\boldsymbol w)$ 决定：

$$
s_\pm^2=\frac{R^2+W^2\pm\sqrt{(R^2-W^2)^2+4p^2}}2.
$$

加入不参与耦合的旁观 Majorana 并固定宇称后，四个能级可写为

$$
-s_+-s_-,\quad -s_++s_-,\quad s_+-s_-,\quad s_++s_-.
$$

低能内部劈裂为 $2s_-$，与高能子空间的最小间隔为 $2(s_+-s_-)$。因此“辅助能隙大”不意味着“逻辑子空间保持简并”。

$$
s_-=0\iff\boldsymbol v\times\boldsymbol w=0
\iff eu_3=u_1u_3=0,\quad de=u_2u_1.
$$

在 $e>0$、$u_3\ne0$ 的原路径上，简并条件不成立。这正是与理想暗态 holonomy 构造的关键差别：存在几何联络，也存在非平凡的内部动力学。

## 3. 理想几何门与快模投影

### 3.1 快模与慢模

**快模是能量较高、内禀振荡频率较大的辅助费米子模式。** 这里的“快”指冻结瞬时控制参数后，模式自身的动力学时间尺度；控制沿路径变化的速度则由段时长 $\tau$ 决定。

先取 $e=u_1=0$ 的理想参考模型。令 $\boldsymbol n=\boldsymbol v/R$，并定义归一化的 Majorana 组合

$$
\gamma_c=\boldsymbol n\cdot\boldsymbol\gamma
=\frac{d\gamma_b+u_2\gamma_2-u_3\gamma_3}{R},
\qquad R=\sqrt{d^2+u_2^2+u_3^2}.
$$

强耦合部分为 $H_0=iR\gamma_a\gamma_c$。将这两个 Majorana 组成普通费米子

$$
c_f=\frac{\gamma_a+i\gamma_c}{2},\qquad
\boxed{H_0=R(2c_f^\dagger c_f-1).}
$$

其空占据、满占据的能量分别为 $-R,+R$，激发能为 $2R$。冻结参数后，Majorana 算符满足

$$
\dot\gamma_a=\frac{2R}{\hbar}\gamma_c,\qquad
\dot\gamma_c=-\frac{2R}{\hbar}\gamma_a,
\qquad \omega_{\rm fast}=\frac{2R}{\hbar}.
$$

当前三段协议均有 $R=Cq$，因此

$$
\frac{C}{\sqrt2}\le R\le C.
$$

快模始终有量级为 $C$ 的大能隙，但它的组成随路径转动：初始时 $\gamma_c=\gamma_b$，第一段末为 $\gamma_2$，第二段末为 $-\gamma_3$，第三段末回到 $\gamma_b$。

理想参考模型中，$\gamma_1$ 以及在 $(\gamma_b,\gamma_2,\gamma_3)$ 空间内与 $\boldsymbol n$ 正交的两个 Majorana 组合都是零模。加入较小的 $e,u_1$ 后，它们形成低能内部自由度，其动力学能量尺度由弱耦合决定，故称慢模。严格说，非零扰动也会改变快、慢模式的组成；第 2 节的精确奇异值 $s_+,s_-$ 给出频率 $2s_+/\hbar,2s_-/\hbar$，此外还有一个零频模式。弱耦合下 $s_+\simeq R$，$s_-\ll R$。

“消去快模”是将演化约化到快费米子保持低能占据的子空间。理想投影为

$$
P_0=1-c_f^\dagger c_f=\frac{1-i\gamma_a\gamma_c}{2}.
$$

所需条件是 $e,\eta\ll C$、$C\tau/\hbar\gg1$，并匹配初末低能制备。投影子空间沿路径变化，仍会产生几何联络；辅助空间的虚跃迁则产生下一节的有效能量修正。因此约化后既保留几何输运，也保留慢模内部动力学。

**辅助快模绝热不要求慢模内部也完全绝热。** 标量保真度函数保留低能内部演化的全部累计相位；完整五维矩阵函数还保留辅助模式及快慢转移，可用于计算有限速度泄漏和检验投影误差。

### 3.2 移动标架与理想 holonomy

选逻辑 Pauli 轴

$$
X=-i\gamma_2\gamma_3,\quad
Y=-i\gamma_3\gamma_1,\quad
Z=-i\gamma_1\gamma_2.
$$

用 $\boldsymbol n=(\cos\theta,\sin\theta\cos\phi,\sin\theta\sin\phi)$ 参数化方向，可取移动标架

$$
Q=e^{-\phi\gamma_2\gamma_3/2}
  e^{-\theta\gamma_b\gamma_2/2}
  e^{+\phi\gamma_2\gamma_3/2}.
$$

将其导数投影到 $P_b=(1-i\gamma_a\gamma_b)/2$，得

$$
\mathcal A=iP_bQ^\dagger dQP_b
=-\frac{1-\cos\theta}{2}X\,d\phi,
\quad
\mathcal F=-\frac12X\sin\theta\,d\theta\wedge d\phi.
$$

原三段路径包围有向立体角 $\Omega=-\pi/2$，去除共同动态相位后

$$
B=e^{i\pi X/4},\qquad
B^\dagger YB=Z,\quad B^\dagger ZB=-Y.
$$

这给出了编织门与球面几何的对应，不依赖于每段速度；辅助模式必须能跟随该子空间。

## 4. Schrieffer–Wolff 消去及几何修正

分解弱耦合：

$$
\omega=\boldsymbol n\cdot\boldsymbol w,\qquad
\boldsymbol w_\perp=\boldsymbol w-\omega\boldsymbol n,
\qquad q_\perp=|\boldsymbol w_\perp|.
$$

扰动的块内、块外部分分别为

$$
D=i\gamma_1\boldsymbol w_\perp\cdot\boldsymbol\gamma,
\quad V_{\rm off}=i\omega\gamma_1\gamma_c.
$$

$D$ 与 $i\gamma_a\gamma_c$ 对易，而 $V_{\rm off}$ 反对易。取反厄米生成元

$$
S=\frac{\omega}{2R}\gamma_a\gamma_1,\qquad [S,H_0]=-V_{\rm off}.
$$

展开 $e^SHe^{-S}$，投影后得到

$$
H_{\rm eff}=-\left(R+\frac{\omega^2}{2R}\right)I
+P_0DP_0+O(W^3/R^2).
$$

二阶静态项是共同相位；影响逻辑读出的下一阶静态能量通常是三阶。又因为 $P_0\gamma_a\gamma_1P_0=0$，几何联络的一阶修正消失。

还可直接由精确投影看出几何修正。对 $R>W$ 的连续分支，令

$$
\vartheta=\tfrac12\operatorname{atan2}(2p,R^2-W^2),
\quad \gamma_+=\cos\vartheta\gamma_a+\sin\vartheta\gamma_1,
$$

$$
\boldsymbol q_+=\cos\vartheta\boldsymbol v+\sin\vartheta\boldsymbol w,
\quad \boldsymbol n_+=\boldsymbol q_+/s_+,
\quad P=\frac{1-i\gamma_+\boldsymbol n_+\cdot\boldsymbol\gamma}{2}.
$$

中心 Majorana 转动的投影导数为零；联络仍是上一节的球面联络，但方向改为 $\boldsymbol n_+$。展开得

$$
\boldsymbol n_+=\boldsymbol n+\frac{\omega\boldsymbol w_\perp}{R^2}
+O(\epsilon^4),\qquad
s_-=q_\perp\left(1-\frac{\omega^2}{2R^2}\right)+O(\epsilon^5R).
$$

立体角变分为

$$
\delta\Omega=\oint\frac\omega{R^2}
\boldsymbol n\cdot(\boldsymbol w_\perp\times d\boldsymbol n)+O(\epsilon^4).
$$

当前门控路径只有第二段贡献。用 $\beta$ 表示该段方向角，积分化为

$$
\delta\Omega=\frac{e\eta}{C^2}
\int_0^{\pi/2}\cos^2\beta(\sin\beta+\cos\beta)d\beta
+O(\epsilon^4)
=\boxed{\frac{e\eta}{C^2}+O(\epsilon^4)}.
$$

因此原路径的几何部分也有二阶偏差。后面的组合门只消去指定有效模型中的动力学项，并未自动消去这个几何修正。

## 5. 原协议的有效动力学

取每段平行输运的暗平面基矢 $\boldsymbol f_1,\boldsymbol f_2$：

| 段 | $\boldsymbol n$ | $\boldsymbol f_1$ | $\boldsymbol f_2$ |
|---|---|---|---|
| I | $(\cos\theta,\sin\theta,0)$ | $(-\sin\theta,\cos\theta,0)$ | $(0,0,1)$ |
| II | $(0,\cos\beta,-\sin\beta)$ | $(-1,0,0)$ | $(0,\sin\beta,\cos\beta)$ |
| III | $(\sin\zeta,0,-\cos\zeta)$ | $(-\cos\zeta,0,-\sin\zeta)$ | $(0,1,0)$ |

一阶块内扰动为

$$
h_K=-\boldsymbol w\cdot\boldsymbol f_1 Z+
\boldsymbol w\cdot\boldsymbol f_2Y,
$$

故三段分别为

$$
h_I=\frac{\eta a^2-eb}{q}Z,\quad
h_{II}=\eta bZ+egY,\quad h_{III}=eY.
$$

定义椭圆积分常数，参数约定为 $m=1/2$：

$$
k=\frac{\mathsf K(1/2)}\pi,\quad j=\frac{\mathsf E(1/2)}\pi,
\qquad \mathcal A=\frac12+j,\quad \mathcal B=1+2k.
$$

因为 $\int g=\int b/q=k$、$\int a^2/q=j$，单圈传播为

$$
U_{\rm cyc}=BD,\quad D=D_3D_2D_1,
\quad D_1=e^{-i\lambda(rj-k)Z},\quad D_3=e^{-i\lambda Y},
$$

$$
D_2=\mathcal T\exp\left[-i\lambda\int_0^1(rbZ+gY)ds\right].
$$

这里的时间排序将在下一篇中用收敛积分、微分方程和系数递推实际求出，不作为最终“未知传播子”。两圈的误差矩阵必须包括几何共轭：

$$
U_{\rm cyc}^2=B^2E,\qquad E=(B^\dagger DB)D.
$$

它不是 $B^2D^2$。写

$$
E=c_0I-i(c_xX+c_yY+c_zZ),
$$

无泄漏的 $SO(3)$ 与 $SU(2)$ 对应给出原波函数读出

$$
\boxed{f=1-c_x^2-c_y^2.}
$$

这是为什么本文的 $f$ 虽然定义为单粒子重叠模，却在逻辑自旋表示中是一个概率；不能再任意取一次平方根。

## 6. 二阶保真度及抵消比例

令

$$
y=1+k,\quad z=r\mathcal A-k,\quad
v=y-z=\mathcal B-\mathcal Ar,\quad w=y+z=1+\mathcal Ar.
$$

单圈一阶生成元为 $zZ+yY$。经过两圈共轭，

$$
\log E=-i\lambda(vY+wZ)+O_r(\lambda^2).
$$

于是

$$
\boxed{f=1-(\mathcal B-\mathcal Ar)^2\lambda^2+O_r(\lambda^4)},
\qquad F=1-2(\mathcal B-\mathcal Ar)^2\lambda^2+O_r(\lambda^4).
$$

其中 $\mathcal A\simeq0.929923$、$\mathcal B\simeq2.18034$。指定读出的抵消比例是

$$
\boxed{r_* =\mathcal B/\mathcal A\simeq2.3446.}
$$

它使 $Y$ 方向误差消失，但 $Z$ 方向仍在。若用去掉目标门后的过程保真度

$$
F_{\rm proc}=|\operatorname{tr}E|^2/4,
$$

则

$$
1-F_{\rm proc}=\lambda^2(v^2+w^2)+O_r(\lambda^4).
$$

因此 $r_*$ 不是普适量子门的二阶无误差点。

## 7. 任意比例的四阶公式：系数与推导

### 7.1 只由脉冲决定的常数

设

$$
P(s)=\int_0^s b(u)du=\frac s2+\frac{\sin\pi s}{2\pi},
\quad G(s)=\int_0^s g(u)du,
$$

$$
J_p=\int_0^1g(1-2P)ds>0,\quad
T_p=\frac12\int_0^1g(2P-\tfrac12)^2ds,\quad
U_p=\frac12\int_0^1b(2G-k)^2ds.
$$

下标用于避免与传播矩阵、椭圆积分常数混淆。再定义

$$
\ell_p=\frac k2-J_p,\qquad
 a_4=U_p-\frac{k^2}{4}-(1+k)J_p,\quad
 b_4=\frac k8+jJ_p-T_p.
$$

这些是一维已知积分，无未知态或传播子。

### 7.2 第二段展开到三阶

写 $D_2=d_{20}I-i(d_{2x}X+d_{2y}Y+d_{2z}Z)$，Dyson 展开得

$$
\begin{aligned}
d_{20}&=1-\frac{\lambda^2}{2}(k^2+r^2/4)+O_r(\lambda^4),\\
d_{2x}&=r\ell_p\lambda^2+O_r(\lambda^4),\\
d_{2y}&=k\lambda-(k^3/6+r^2T_p)\lambda^3+O_r(\lambda^5),\\
d_{2z}&=r\lambda/2-(r^3/48+rU_p)\lambda^3+O_r(\lambda^5).
\end{aligned}
$$

例如两次 $Z$、一次 $Y$ 的排列满足 $ZZY=YZZ=Y$、$ZYZ=-Y$。固定 $Y$ 所在时间为 $s$，其余积分变成

$$
\int_0^1g\left[\frac{P^2}2+\frac{(1/2-P)^2}2-P(1/2-P)\right]ds
=T_p.
$$

$U_p$ 同理。三重有序积分由此降成一维常数。

### 7.3 乘回第一、三段与第二圈

令 $t=rj-k$、$N=y^2+z^2$、$m=r\ell_p+kt+z$。单圈 $D=d_0I-i\boldsymbol d\cdot\boldsymbol\sigma$ 的系数为

$$
\begin{aligned}
d_0&=1-N\lambda^2/2+O_r(\lambda^4),&d_x&=m\lambda^2+O_r(\lambda^4),\\
d_y&=y\lambda-Y_3\lambda^3+O_r(\lambda^5),&
d_z&=z\lambda-Z_3\lambda^3+O_r(\lambda^5),
\end{aligned}
$$

$$
Y_3=y^3/6+z^2/2+kt^2/2+r\ell_pt+r^2T_p,
\quad
Z_3=z^3/6+z/2+t(k^2/2+k)+r(U_p+\ell_p).
$$

两圈乘法精确给出

$$
c_x=2d_0d_x-d_y^2-d_z^2,\qquad c_y=(d_0+d_x)(d_y-d_z).
$$

记

$$
C_2=2m-N=-(v^2+2rJ_p),\quad
\Delta_3=Z_3-Y_3=-v^3/6+ra_4+r^2b_4.
$$

则 $c_x=C_2\lambda^2+O_r(\lambda^4)$，
$c_y=v\lambda+(\Delta_3+C_2v/2)\lambda^3+O_r(\lambda^5)$。代入读出，得到全部四阶系数：

$$
\boxed{f=1-v^2\lambda^2+\mathcal P_4(r)\lambda^4+O_r(\lambda^6)},
$$

$$
\boxed{\mathcal P_4(r)=\frac{v^4}{3}-2rJ_pv^2-4r^2J_p^2
-2rv(a_4+b_4r).}
$$

若用平方读出：

$$
F=1-2v^2\lambda^2+[v^4+2\mathcal P_4(r)]\lambda^4+O_r(\lambda^6).
$$

固定 $r$ 下仅有偶次项，源于 Pauli 分级与 $\lambda\to-\lambda$ 对称。$O_r$ 不对 $r\to\infty$ 一致；它也不包含投影和有限速度误差。

### 7.4 抵消点与单参数极限

在 $r_*$：

$$
\boxed{1-f=4r_*^2J_p^2\lambda^4+O(\lambda^6)},\quad
1-F=8r_*^2J_p^2\lambda^4+O(\lambda^6).
$$

$J_p>0$ 说明四阶剩余严格非零。两个单参数极限在投影模型内则有全阶三角函数形式：

$$
\eta=0:\quad f=\cos^2(\mathcal B\lambda),\qquad
e=0:\quad f=\cos^2(\mathcal A\nu).
$$

二者的小角展开分别还原四阶式。有限 $r$ 时，不对易项不能用一条经验余弦替代。

与更早的 [PRB 105, 054507](../arts/prb105.pdf) 对照，其单动态相位结果为 $f=\cos^2(\theta_2/2)$；取 $\theta_2=2\mathcal A\nu$，即与本文 $e=0$ 结果及二阶、四阶展开一致。但原文式 (3) 保留 $i\epsilon_d(t)\gamma_1\gamma_2$，其耦合标号和控制路径也与本文不同，不能把它直接等同于 $E_1=0$ 的完整模型。逐项模型映射、全部读出模长与验证边界见 [解析延拓文档第 6.5 节](analytic_continuation.md#65-与-prb-105-054507-2022-的单动态相位结果对照)。

更高阶可递推：对 $H_s=rbZ+gY$，令 $V_0=I$，

$$
V_{n+1}(s)=-i\int_0^sH_uV_n(u)du,\quad D_2=\sum_{n\ge0}\lambda^nV_n(1).
$$

以 $M=|\lambda|(r/2+k)$，矩阵尾项满足

$$
\left\|D_2-\sum_{n=0}^N\lambda^nV_n\right\|
\le e^M\frac{M^{N+1}}{(N+1)!}.
$$

这说明高阶系数可以构造，但大累计相位更适合另一篇的统一有序积分展开。

## 8. 从群结构寻找门误差抵消

### 8.1 对合分级：为什么一阶在平面内、二阶沿法向？

全模型属于 $\mathfrak{so}(5)$。取

$$
\mathscr C=-\gamma_a\gamma_b\gamma_2\gamma_3,
\quad \mathscr C H_0\mathscr C=H_0,\quad
\mathscr C V\mathscr C=-V.
$$

在低能空间它等价于 $X$ 共轭：$Y,Z$ 为奇，$X$ 为偶。由 $[Y,Z]=2iX$，奇数阶 Magnus 项位于 $YZ$ 平面，偶数阶沿 $X$。这解释了仅调 $r_*$ 无法同时消掉二阶误差生成元。

单圈写成

$$
D=\exp[-i\lambda K_1-i\lambda^2K_2+O_r(\lambda^3)],
\quad K_1=zZ+yY,\quad K_2=c(r)X,
$$

$$
c(r)=(1+k)(rj-k)+r/2+r\ell_p.
$$

### 8.2 编织共轭形成 $C_4$ 平均

定义 $\mathscr R(K)=B^\dagger KB$。在 $YZ$ 平面它是四分之一转动，

$$
\sum_{j=0}^3\mathscr R^j(K_1)=0.
$$

保持脉冲形状、仅缩放每圈时长 $s_j\tau$：

| 顺向圈数 | 权重 $s_j$ | 净几何门（忽略整体相位） | 总时长 |
|---:|---|---|---:|
| 5 | $(1,1,2,1,1)$ | $B$ | $18\tau$ |
| 6 | $(1,1,2,2,1,1)$ | $B^2$ | $24\tau$ |

两组的一阶向量和为零，对任意固定 $e,\eta$ 成立。二阶剩余生成元的 $X$ 系数分别为

$$
-4\lambda^2(v^2+2rJ_p),\qquad -6\lambda^2(v^2+2rJ_p).
$$

因此过程不保真度降到四阶，但不是全阶抵消。更一般地，同向正时长序列的二阶系数由

$$
\left(\sum_js_j^2\right)c(r)+Q_{\boldsymbol s}(y^2+z^2)
$$

组成，$Q_{\boldsymbol s}$ 仅由时长与方向决定。前者关于 $r$ 为一次式，后者含二次项；要对所有 $r$ 抵消，先要 $Q_{\boldsymbol s}=0$，随后还需正数 $\sum s_j^2$ 乘 $c(r)$ 恒为零，不能成立。这是该类同向设计的限制。

### 8.3 严格反向路径提供二阶符号翻转

反向必须同时反转三段次序和所有门控函数，包括 $u_1(t)$。对总长 $3\tau$ 的单圈，

$$
h_-(t)=B h_+(3\tau-t)B^\dagger,
\quad K_{1,-}=BK_{1,+}B^\dagger,
\quad K_{2,-}=-K_{2,+}.
$$

令第 $j$ 圈方向为 $\sigma_j=\pm1$，前面净方向数为 $q_{j-1}=\sum_{k<j}\sigma_k$，

$$
\ell_j=q_{j-1}-\frac{1-\sigma_j}{2}.
$$

一阶、二阶对任意固定弱耦合抵消，可由三个有限代数约束保证：

$$
\boxed{\sum_js_je^{-i\ell_j\pi/2}=0},
$$

$$
\boxed{\sum_j\sigma_js_j^2=0},\qquad
\boxed{\sum_{j>k}s_js_k\sin[(\ell_j-\ell_k)\pi/2]=0}.
$$

第一条消掉平面内向量和，第二条消掉各圈内部的二阶项，第三条消掉跨圈 BCH 交换子。还需指定净方向数，以保留目标几何门。

### 8.4 九圈二阶校正门：完整构造

取

$$
\boldsymbol\sigma=(+,+,+,+,-,-,-,-,+),\qquad
\boldsymbol\ell=(0,1,2,3,3,2,1,0,0).
$$

令 $a_0,t_0,i_0>0$，$h_0=\sqrt{a_0^2+i_0^2}$，$c_0=(a_0+h_0+i_0)/2$，权重为

$$
\boldsymbol s=(a_0,t_0,c_0,t_0,t_0,c_0,t_0,h_0,i_0).
$$

一阶实部为 $a_0-2c_0+h_0+i_0=0$，虚部为 $-t_0+t_0+t_0-t_0=0$；内部二阶权重为 $a_0^2-h_0^2+i_0^2=0$。跨圈交换子若暂将两处 $c_0$ 分别记为 $c_0,f_0$，求和为 $2t_0(f_0-c_0)|\boldsymbol k|^2$，故同样为零。

具体可取

$$
\boxed{\boldsymbol s=(3/4,1,3/2,1,1,3/2,1,5/4,1)}.
$$

总时长 $30\tau$、净几何门 $B$。在上述一阶投影动力学模型中，设 $\Lambda=\sqrt{e^2+\eta^2}$，则

$$
U_9\sim B\exp[O((\Lambda\tau/\hbar)^3)],\qquad
1-F_{\rm proc,9}=O((\Lambda\tau/\hbar)^6).
$$

重复两次得到目标 $B^2$。这是明确的存在性构造，不声称圈数最优；九圈之间幅度必须固定，累计弱角要小。它没有自动补偿几何偏移、三阶静态修正、有限速度泄漏或参数漂移。

## 9. 何时能称为几何门？以及对 Fig. 1(d) 的解释范围

精确低能投影的 Kato 输运可写为

$$
U_{\rm low}=U_g\,\mathcal T\exp\left[-\frac i\hbar\int h_K(t)dt\right].
$$

几何联络存在不等于动力学因子为标量。全程简并使 $h_K$ 无迹部分为零，是纯几何演化的充分条件；也可以通过额外控制让整体动力学因子为 $\pm I$。上面的组合序列给的是有阶数保证的近似抵消。

若在同一误差标架中能实现严格反对称 $h(T-t)=-h(t)$，后半段是前半段的逆，可全阶消去动力学。但物理路径反向本身也会反转 holonomy，不能据此直接宣称留下一个非平凡几何门；必须另外验证净几何路径。

对 Fig. 1(d)，微扰式解释左侧亮区在 $r_*$ 附近变平、两侧单参数主导时出现余弦振荡的原因。图中 $\lambda$ 延伸到 12，不能用四阶多项式解释整幅条纹。有限比例下的峰位漂移、峰高下降、长时间相干抵消，须用 [解析延拓文档](analytic_continuation.md) 的全累计相位公式。

本篇 $O_r(\lambda^6)$ 只描述微扰级数截断。物理投影误差、有限速度误差，以及理想峰顶误差为何能再高一阶，均在该文档的“假设与误差预算”中单独列出。
