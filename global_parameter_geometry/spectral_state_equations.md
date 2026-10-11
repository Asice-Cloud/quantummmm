# 全域状态方程与误差控制

本文从已有的完整状态方程出发，推导整个有限参数区域上的统一演化表示。主要结构是：参数生成元仿射依赖、Legendre 相邻阶耦合，以及反对称生成元给出的稳定性和误差界。

数学推导保留 \(\tau,t_c,E_0\)；只有 \(E=E_1\)、\(\eta=|t_1|\) 作为不可控参数。本文没有对 \(E\tau\)、\(\eta\tau\) 作小量截断。

## 1. 完整状态与固定协议

按 \((1,2,3,a,b)\) 排列模式，定义

\[
L_{ij}=e_ie_j^T-e_je_i^T,\qquad
K=2\bigl(EL_{12}+\eta gL_{1b}-pL_{2a}+rL_{3a}+e_dL_{ab}\bigr).
\]

每段局部时间为 \(v=(t-j\tau)/\tau\in[0,1]\)，\(j=0,\ldots,5\)。令

\[
a(v)=\frac{1-\cos\pi v}{2},\qquad c(v)=1-a(v).
\]

| 段 | 时间区间 | \(p\) | \(r\) | \(e_d\) | \(g\) |
|---|---|---|---|---|---|
| 1、4 | \([0,\tau]\)、\([3\tau,4\tau]\) | \(t_ca\) | \(0\) | \(E_0c\) | \(a\) |
| 2、5 | \([\tau,2\tau]\)、\([4\tau,5\tau]\) | \(t_cc\) | \(t_ca\) | \(0\) | \(c\) |
| 3、6 | \([2\tau,3\tau]\)、\([5\tau,6\tau]\) | \(0\) | \(t_cc\) | \(E_0a\) | \(0\) |

完整状态满足

\[
\partial_t\mathcal B=[K,\mathcal B],\qquad \mathcal B(E,\eta,0)=L_{12}.
\tag{1}
\]

它在六维有向平面流形

\[
\mathcal M=SO(5)/(SO(2)\times SO(3))
\]

上演化，满足

\[
\mathcal B^T=-\mathcal B,\quad \mathcal B^3=-\mathcal B,\quad
\operatorname{rank}\mathcal B=2,\quad \|\mathcal B\|_F^2=2.
\]

这套状态描述对应此前选定的五维波函数轨道，不扩展为任意多体量子态的完整描述。

使用十个全局分量

\[
\mathbf b=(d,x_1,x_2,x_3,y_1,y_2,y_3,h_1,h_2,h_3)^T,
\]

\[
d=\mathcal B_{12},\quad
\mathbf x=(\mathcal B_{13},\mathcal B_{1a},\mathcal B_{1b}),\quad
\mathbf y=(\mathcal B_{23},\mathcal B_{2a},\mathcal B_{2b}),\quad
\mathbf h=(\mathcal B_{3a},\mathcal B_{3b},\mathcal B_{ab}).
\]

真实解满足 \(\|\mathbf b\|_2=1\)，初值为 \(\mathbf b_0=(1,0,\ldots,0)^T\)。

## 2. 对两个参数的仿射分解

式 (1) 等价于

\[
\boxed{\partial_t\mathbf b=M(E,\eta,t)\mathbf b,
\qquad M=M_c(t)+E M_E+\eta M_\eta(t).}
\tag{2}
\]

三个矩阵的作用可以完全写出：

\[
M_c\mathbf b=2\begin{pmatrix}
-p x_2\\
r x_2\\
p d-rx_1+e_dx_3\\
-e_dx_2\\
r y_2+p h_1\\
-r y_1+e_dy_3\\
-e_dy_2-p h_3\\
-p y_1+e_dh_2\\
-e_dh_1+r h_3\\
p y_3-r h_2
\end{pmatrix},\qquad
M_E\mathbf b=2\begin{pmatrix}
0\\y_1\\y_2\\y_3\\-x_1\\-x_2\\-x_3\\0\\0\\0
\end{pmatrix},\qquad
M_\eta\mathbf b=2g\begin{pmatrix}
-y_3\\-h_2\\-h_3\\0\\0\\0\\d\\0\\x_1\\x_2
\end{pmatrix}.
\tag{3}
\]

它们满足

\[
M_c^T=-M_c,\quad M_E^T=-M_E,\quad M_\eta^T=-M_\eta,
\qquad \|M_E\|_2=2,\quad \|M_\eta\|_2=2g(t).
\tag{4}
\]

最后一个等式使用本协议的 \(g\ge0\)。这里的矩阵范数是诱导 Euclidean 范数。式 (3) 已包含物理时间生成元的因子 2；改用局部 \(v\) 时，整个右端再乘 \(\tau\)。

## 3. 任意阶参数响应及全域上界

记 \(D_{mn}=\partial_E^m\partial_\eta^n\mathbf b\)。因为 \(M\) 对参数仿射，直接求导得到

\[
\boxed{
\partial_tD_{mn}=M D_{mn}
+m M_E D_{m-1,n}+n M_\eta D_{m,n-1}.
}
\tag{5}
\]

下标为负的项省略，\(D_{00}(0)=\mathbf b_0\)，其余初值为零。已有的一阶和二阶响应是这套精确层级的低阶部分。

对固定实参数，\(M\) 的传播算子是正交矩阵。利用常数变易公式与式 (4)，令

\[
\alpha(t)=2t,\qquad \beta(t)=2\int_0^t g(s)\,ds,
\]

可归纳得到与参数点无关的上界

\[
\boxed{\|D_{mn}(E,\eta,t)\|_2\le\alpha(t)^m\beta(t)^n.}
\tag{6}
\]

归纳步骤的积分上界恰好是 \(\alpha^m\beta^n\) 的时间导数；\(m=n=0\) 时使用 \(\|\mathbf b\|_2=1\)。由于六段中 \(\int_0^{6\tau}g\,dt=2\tau\)，终点有

\[
\|\partial_E^m\partial_\eta^n\mathbf b(E,\eta,6\tau)\|_2
\le (12\tau)^m(4\tau)^n.
\tag{7}
\]

有限时间内，有限维线性方程对仿射参数的解具有解析依赖；也可由收敛的时间有序积分级数看出这一点。因此在任意给定紧矩形上，可以逐步提高多项式阶数逼近完整参数场。

式 (6) 是一个全域的响应上界，通常比较保守。它不直接给出实际条带间距，也不单独决定峰数。它说明参数域宽度和累计演化时间都会影响所需的展开阶数。

## 4. 整个参数区域上的 Legendre 展开

取

\[
\Omega=[E_c-\Delta E,E_c+\Delta E]
\times[\eta_c-\Delta\eta,\eta_c+\Delta\eta],
\]

其中 \(\Delta E,\Delta\eta>0\)，\(\eta_c\ge\Delta\eta\)。定义

\[
x=\frac{E-E_c}{\Delta E},\qquad
y=\frac{\eta-\eta_c}{\Delta\eta}.
\]

令 \(P_m\) 是 Legendre 多项式，取

\[
\phi_m(x)=\sqrt{2m+1}\,P_m(x),\qquad
\int_{-1}^1\phi_m\phi_n\,\frac{dx}{2}=\delta_{mn}.
\]

此归一化只是函数空间内积的约定，不假定不可控参数服从某个概率分布。乘以坐标的作用为

\[
x\phi_m=a_{m+1}\phi_{m+1}+a_m\phi_{m-1},\qquad
a_m=\frac{m}{\sqrt{(2m-1)(2m+1)}},\quad a_0=0.
\tag{8}
\]

展开完整状态场：

\[
\mathbf b(E,\eta,t)=\sum_{m,n\ge0}\mathbf c_{mn}(t)\phi_m(x)\phi_n(y),
\qquad
\mathbf c_{mn}=\int_{-1}^1\int_{-1}^1
\mathbf b\phi_m\phi_n\,\frac{dx\,dy}{4}.
\tag{9}
\]

这是有限区域上的全局表示；\(\mathbf c_{mn}\) 并非某个参数点的 Taylor 导数。

## 5. 系数演化、初值和接口

令 \(\bar M=M_c+E_cM_E+\eta_cM_\eta\)。将式 (9) 代入式 (2)，利用正交性和递推式 (8)，得到

\[
\boxed{
\begin{aligned}
\dot{\mathbf c}_{mn}={}&\bar M\mathbf c_{mn}\\
&+\Delta E\,M_E
\bigl(a_m\mathbf c_{m-1,n}+a_{m+1}\mathbf c_{m+1,n}\bigr)\\
&+\Delta\eta\,M_\eta
\bigl(a_n\mathbf c_{m,n-1}+a_{n+1}\mathbf c_{m,n+1}\bigr).
\end{aligned}}
\tag{10}
\]

负下标系数取零。无限阶系统精确地表示原始参数场；有限阶系统需要截断误差分析。

由于初态与参数无关且 \(\phi_0=1\)，

\[
\mathbf c_{00}(0)=\mathbf b_0,\qquad
\mathbf c_{mn}(0)=0\quad(m+n>0).
\tag{11}
\]

各段接口 \(t_j=j\tau\) 满足

\[
\mathbf c_{mn}(t_j^+)=\mathbf c_{mn}(t_j^-),\qquad j=1,\ldots,5.
\tag{12}
\]

全部状态系数传递到下一段；在 \(3\tau\) 不重置。参数域的四条边没有额外的 Dirichlet 或 Neumann 边界条件，因为不同参数标签之间没有物理输运。终点目标态用于读出保真度，也不是时间终边界条件。

式 (10) 保留了时间演化。它提供全参数域的系数描述，没有增加独立的参数方向动力学。

## 6. 有限阶截断保留什么结构

保留 \(0\le m\le P\)、\(0\le n\le Q\)，在式 (10) 中把区域外系数置零，这对应正交 Galerkin 截断。未知实数个数为 \(10(P+1)(Q+1)\)。

定义对称三对角矩阵 \((J_P)_{m,m+1}=(J_P)_{m+1,m}=a_{m+1}\)。按 \((m,n,\text{状态分量})\) 展平系数，截断系统的生成元是

\[
\mathscr M_{PQ}=
I_{P+1}\otimes I_{Q+1}\otimes\bar M
+\Delta E J_P\otimes I_{Q+1}\otimes M_E
+\Delta\eta I_{P+1}\otimes J_Q\otimes M_\eta.
\tag{13}
\]

由 \(J^T=J\) 和式 (4)，

\[
\mathscr M_{PQ}^T=-\mathscr M_{PQ},\qquad
\frac{d}{dt}\sum_{m=0}^P\sum_{n=0}^Q\|\mathbf c_{mn}\|_2^2=0.
\tag{14}
\]

这里说的是截断系统在精确时间积分下保持系数范数；实际积分仍有数值误差。

式 (14) 等价于保持重建状态的参数域平均平方范数。它不保证每一个参数点的 \(\|\widetilde{\mathbf b}\|_2=1\)，也不保证 \(\widetilde{\mathcal B}^3=-\widetilde{\mathcal B}\)。重建矩阵的反对称性可以由分量表示自动保证；其余几何约束必须检查。若对近似状态作额外投影，投影后的场需要重新计算残差。

## 7. 保真度与参数导数的读取

对真实状态使用

\[
F=\frac{d^2-d}{2}+\frac{|\mathbf x|^2+|\mathbf y|^2}{4}.
\tag{15}
\]

对近似状态也采用这一多项式来定义 \(\widetilde F\)，便于直接分析误差。真实状态上的等价公式为

\[
F=\frac{(1-d)^2-|\mathbf h|^2}{4}.
\]

在有限阶近似中，这两个读数的差恰好是

\[
\frac{\|\widetilde{\mathbf b}\|_2^2-1}{4}.
\tag{16}
\]

所以，不能一边违反几何约束，一边把不同的精确读出公式当作数值上必然相同。

参数导数直接由多项式求出，例如

\[
\partial_E\widetilde{\mathbf b}
=\frac1{\Delta E}\sum_{m,n}\mathbf c_{mn}\phi_m'(x)\phi_n(y),
\quad
\partial_\eta\widetilde{\mathbf b}
=\frac1{\Delta\eta}\sum_{m,n}\mathbf c_{mn}\phi_m(x)\phi_n'(y).
\tag{17}
\]

二阶导数同理；每次求物理参数导数都要保留相应的尺度因子。由式 (15)，对 \(\lambda,\mu\in\{E,\eta\}\)，

\[
F_\lambda=(d-\tfrac12)d_\lambda
+\tfrac12(\mathbf x\cdot\mathbf x_\lambda+\mathbf y\cdot\mathbf y_\lambda),
\tag{18}
\]

\[
\begin{aligned}
F_{\lambda\mu}={}&(d-\tfrac12)d_{\lambda\mu}+d_\lambda d_\mu\\
&+\tfrac12\bigl(
\mathbf x_\lambda\cdot\mathbf x_\mu+\mathbf x\cdot\mathbf x_{\lambda\mu}
+\mathbf y_\lambda\cdot\mathbf y_\mu+\mathbf y\cdot\mathbf y_{\lambda\mu}\bigr).
\end{aligned}
\tag{19}
\]

研究峰时优先使用 \(F\)。在 \(F>0\) 的区域，\(\sqrt F\) 有相同的局部峰位置；\(F=0\) 附近不应通过对平方根求导判断光滑性。

## 8. 状态残差与全域误差界

对任意在各段内可微、接口连续的近似场，定义

\[
r=\partial_t\widetilde{\mathbf b}-M\widetilde{\mathbf b},\qquad
e=\mathbf b-\widetilde{\mathbf b}.
\]

则 \(\partial_t e=Me-r\)。反对称项对 \(\|e\|_2\) 不作贡献，故

\[
\boxed{\|e(E,\eta,t)\|_2
\le\|e(E,\eta,0)\|_2+\int_0^t\|r(E,\eta,s)\|_2\,ds.}
\tag{20}
\]

若初值准确，且确实有 \(R_0(t)\ge\sup_\Omega\|r\|_2\)，则

\[
\sup_\Omega\|e(t)\|_2\le\varepsilon_0(t),\qquad
\varepsilon_0(t)=\int_0^tR_0(s)\,ds.
\tag{21}
\]

若数值近似在接口有跳跃，还要加上各次状态跳跃范数；精确传递接口值可以避免这类附加项。

对精确时间积分的 Galerkin 系统，残差只出现在截断边缘向外的一层：

\[
r_{P+1,n}=-\Delta E\,a_{P+1}M_E\mathbf c_{P,n},\quad 0\le n\le Q,
\]

\[
r_{m,Q+1}=-\Delta\eta\,a_{Q+1}M_\eta\mathbf c_{m,Q},\quad 0\le m\le P.
\tag{22}
\]

没有同时增加两个下标的角点项，因为原生成元没有 \(E\eta\) 项。实际积分还要计入保留模式内部的时间求解残差。

若 \(r=\sum r_{mn}\phi_m\phi_n\)，由 \(|P_m|\le1\) 得到一种保守的全域包络：

\[
\sup_\Omega\|r\|_2
\le\sum_{m,n}\|r_{mn}\|_2\sqrt{(2m+1)(2n+1)}.
\tag{23}
\]

这使得空间上的上确界有可计算的上界。时间上的残差积分仍需要误差控制；少量时间采样或求解器容差本身不等于严格认证。

将式 (15) 视为十维变量上的二次多项式，真实状态满足

\[
\|\nabla_{\mathbf b}F\|_2\le\tfrac32,\qquad
\|\nabla_{\mathbf b}^2F\|_2=1.
\]

因此有

\[
\boxed{|\widetilde F-F|\le\tfrac32\varepsilon_0+\tfrac12\varepsilon_0^2.}
\tag{24}
\]

这个界不要求近似状态精确位于物理流形上。

## 9. 梯度和曲率也必须控制误差

对误差方程求参数导数：

\[
\partial_t e_\lambda=M e_\lambda+M_\lambda e-r_\lambda,
\]

\[
\partial_t e_{\lambda\mu}
=M e_{\lambda\mu}+M_\lambda e_\mu+M_\mu e_\lambda-r_{\lambda\mu}.
\tag{25}
\]

在所有初始响应准确为零时，可取

\[
\varepsilon_\lambda(t)=\int_0^t
\bigl(\|M_\lambda(s)\|_2\varepsilon_0(s)+R_\lambda(s)\bigr)\,ds,
\]

\[
\varepsilon_{\lambda\mu}(t)=\int_0^t
\bigl(\|M_\lambda\|_2\varepsilon_\mu
+\|M_\mu\|_2\varepsilon_\lambda+R_{\lambda\mu}\bigr)(s)\,ds,
\tag{26}
\]

其中 \(R_\lambda,R_{\lambda\mu}\) 分别界定残差导数在参数域上的范数。它们可以通过对残差多项式求导再作区间包络获得。

例如，由式 (18) 可得保真度梯度误差界

\[
|F_\lambda-\widetilde F_\lambda|
\le\tfrac32\varepsilon_\lambda
+\varepsilon_0\|\widetilde{\mathbf b}_\lambda\|_2.
\tag{27}
\]

二阶误差可取

\[
\begin{aligned}
|F_{\lambda\mu}-\widetilde F_{\lambda\mu}|
\le{}&\tfrac32\varepsilon_{\lambda\mu}
+\varepsilon_0\|\widetilde{\mathbf b}_{\lambda\mu}\|_2\\
&+\varepsilon_\lambda\|\widetilde{\mathbf b}_\mu\|_2
+\varepsilon_\mu\|\widetilde{\mathbf b}_\lambda\|_2
+\varepsilon_\lambda\varepsilon_\mu.
\end{aligned}
\tag{28}
\]

需要全域界时，右端的近似导数范数也取域上上界。式 (24)、(27)、(28) 分别连接保真度、峰位置和峰曲率的可靠性。

## 10. 本文已经得到与尚待验证的内容

已得到的是原模型的一套全域表示、精确无限层级、有限截断的结构性质，以及在残差已被界定时成立的误差估计。Legendre 递推、矩阵分解和截断生成元的反对称性已作代数核对。

尚待实施的是六段谱演化的数值求解、实际收敛阶数的测量、峰线提取和全域误差包络。上述不等式并不意味着已有一次完成的数值认证，也没有直接给出所有峰的闭式坐标。

下一项任务见 [峰线几何与第一项验证](peak_geometry_and_validation.md)。

进一步的推导已写入 [从全域系数到峰线](global_peak_bounds.md)：包括有限谱系统与 Gauss 节点演化的等价关系、可在求解前计算的全域误差界、从边缘系数产生六个标量误差包络的方法，以及峰线与状态曲面几何的关系。
