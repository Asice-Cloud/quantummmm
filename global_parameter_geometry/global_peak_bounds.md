# 从全域系数到峰线：先验阶数、可计算误差与几何判据

本篇继续前两篇的推导，解决一个具体缺口：已有的峰位置误差公式含有未知的梯度、曲率误差界，怎样把这些界变成由参数区域、时间或实际谱系数决定的量？

本次得到三项结果。第一，仿射参数模型的有限阶 Legendre 系统，等价于 Gauss 节点上的原始演化及其插值。第二，在任何给定有限参数矩形上，可以先验界定状态、梯度和 Hessian 的截断误差。第三，峰的倾斜与高度变化可以直接写成状态曲面的度量和弯曲项。

这些推导保留 \(\tau,t_c,E_0\)，不要求 \(E\tau\) 或 \(\eta\tau\) 小。下文的数值表仅计算误差公式，没有运行完整六段演化，也没有计算出一条新的峰线。

## 1. 有限阶谱系统可以精确对角化到参数节点

沿用 [状态方程](spectral_state_equations.md) 的记号：

\[
\partial_t\mathbf b=(\bar M+\Delta E\,xM_E+\Delta\eta\,yM_\eta)\mathbf b,
\qquad (x,y)\in[-1,1]^2.
\tag{1}
\]

参数方向分别保留到 \(P,Q\) 阶。取 \(x_i\) 为 \(P_{P+1}(x)\) 的零点，\(w_i\) 为归一化测度 \(dx/2\) 下的 Gauss 权重，所以 \(\sum_i w_i=1\)。令

\[
U^x_{mi}=\sqrt{w_i}\,\phi_m(x_i),\qquad 0\le m,i\le P.
\]

由 Gauss 求积对次数不超过 \(2P+1\) 的多项式精确，以及 \(\phi_m\) 的正交性，

\[
U^x(U^x)^T=I.
\]

由三项递推及 \(\phi_{P+1}(x_i)=0\)，

\[
J_PU^x=U^x\operatorname{diag}(x_i).
\tag{2}
\]

\(y\) 方向同理，记节点、权重为 \(y_j,v_j\)，正交矩阵为 \(U^y\)。定义变换后的系数和节点状态：

\[
\widehat{\mathbf c}_{ij}=
\sum_{m,n}U^x_{mi}U^y_{nj}\mathbf c_{mn},\qquad
\mathbf z_{ij}=\frac{\widehat{\mathbf c}_{ij}}{\sqrt{w_iv_j}}.
\]

有限谱系统变成

\[
\boxed{
\partial_t\mathbf z_{ij}
=M(E_c+\Delta E x_i,\eta_c+\Delta\eta y_j,t)\mathbf z_{ij},
\qquad \mathbf z_{ij}(0)=\mathbf b_0.
}
\tag{3}
\]

变换矩阵与时间无关，因此该等价关系覆盖全部六段及连续接口，而不只是在某个时刻成立。

令 \(\ell_i(x),\ell_j(y)\) 为两组节点的 Lagrange 基函数。若时间演化精确求解，则重建的有限阶谱解恰好是

\[
\boxed{
\mathbf b_{PQ}(x,y,t)=
\sum_{i=0}^P\sum_{j=0}^Q
\ell_i(x)\ell_j(y)\mathbf b(x_i,y_j,t)
=I_P^xI_Q^y\mathbf b.
}
\tag{4}
\]

这有两个含义：节点上的状态精确位于物理流形上，节点之间的多项式一般不精确满足流形约束；谱表示也仍然需要求解时间演化。它提供的是完整参数场的可控表示，没有消除原来的动力学求解。

## 2. 参数区域与时间给出无量纲响应尺度

由已有的全阶参数响应上界，定义

\[
\chi_E(t)=2\Delta E\,t,\qquad
\chi_\eta(t)=2\Delta\eta\int_0^t g(s)\,ds.
\]

则在整个参数区域内

\[
\boxed{
\|\partial_x^r\partial_y^s\mathbf b(x,y,t)\|_2
\le\chi_E(t)^r\chi_\eta(t)^s.
}
\tag{5}
\]

本协议终点满足

\[
\chi_E(6\tau)=12\tau\Delta E,\qquad
\chi_\eta(6\tau)=4\tau\Delta\eta.
\tag{6}
\]

上界不依赖参数矩形的中心，也不显式含 \(t_c,E_0\)。原因是中心参数及控制生成的真实传播保持范数。实际解和实际峰线仍依赖这些量；此处只是一个保守的参数响应上界。

## 3. 一维 Gauss 插值的状态与导数误差

先考虑一个向量函数 \(f(x)\)，假设所有需要的阶数满足 \(\|f^{(k)}\|_2\le\chi^k\)。使用 \(n\) 个 Gauss 节点，插值次数为 \(n-1\)。节点多项式是

\[
\omega_n(x)=\prod_{i=1}^{n}(x-x_i)=\frac{P_n(x)}{k_n},
\qquad k_n=2^{-n}\binom{2n}{n}.
\]

插值余项可写成

\[
f-I_{n-1}f=\omega_n(x)D_n(x),\qquad
D_n(x)=f[x_1,\ldots,x_n,x].
\]

括号表示差商，节点重合时按连续延拓解释。差商的单纯形积分表示给出

\[
\|D_n^{(j)}(x)\|_2\le\frac{j!}{(n+j)!}\chi^{n+j}.
\tag{7}
\]

其中系数 \(j!/(n+j)!\) 来自积分中最后一个重心坐标的 \(j\) 次矩。向量情形可对任意单位向量投影后使用相同公式，再取上确界。

Legendre 多项式满足

\[
d_{n,r}:=\max_{[-1,1]}|P_n^{(r)}|
=\frac{(n+r)!}{2^r r!(n-r)!},\qquad 0\le r\le n;
\quad d_{n,r}=0\ (r>n).
\]

对余项求导，得到可显式计算的界

\[
\boxed{
\mathcal E_{n,r}(\chi)=
\frac1{k_n}\sum_{j=0}^{r}\binom rj d_{n,r-j}
\frac{j!\,\chi^{n+j}}{(n+j)!},\qquad
\|(f-I_{n-1}f)^{(r)}\|_\infty\le\mathcal E_{n,r}(\chi).
}
\tag{8}
\]

其中最常用的三项为

\[
q_n(\chi)=\mathcal E_{n,0}(\chi)
=\frac{2^n n!}{(2n)!}\chi^n,
\]

\[
\mathcal E_{n,1}=q_n\left[\frac{n(n+1)}2+\frac{\chi}{n+1}\right],
\]

\[
\mathcal E_{n,2}=q_n\left[
\frac{(n-1)n(n+1)(n+2)}8
+\frac{n(n+1)\chi}{n+1}
+\frac{2\chi^2}{(n+1)(n+2)}\right].
\tag{9}
\]

这些是插值余项的全区间界，没有将物理演化作小参数展开。对固定有限 \(\chi\)，随阶数提高，阶乘使误差上界最终下降。

## 4. 二维参数场的六个误差界

为控制插值算子的导数，定义

\[
C_{P,r}=\left[\sum_{m=0}^P(2m+1)d_{m,r}^2\right]^{1/2}.
\]

由正交矩阵 \(U\) 与 Cauchy–Schwarz 不等式，

\[
\|\partial_x^r I_P f\|_\infty
\le C_{P,r}\max_i\|f(x_i)\|_2,
\qquad C_{P,0}=P+1.
\tag{10}
\]

将插值误差分解成

\[
I-I_P^xI_Q^y=(I-I_P^x)+I_P^x(I-I_Q^y),
\]

或交换 \(x,y\) 的顺序，结合式 (5)、(8)、(10)，可得

\[
\boxed{
\begin{aligned}
\varepsilon_{rs}^{\rm par}(t)=\min\bigl\{&
\chi_\eta^s\mathcal E_{P+1,r}(\chi_E)
+C_{P,r}\mathcal E_{Q+1,s}(\chi_\eta),\\
&\chi_E^r\mathcal E_{Q+1,s}(\chi_\eta)
+C_{Q,s}\mathcal E_{P+1,r}(\chi_E)\bigr\},
\end{aligned}}
\tag{11}
\]

\[
\sup_\Omega\|\partial_x^r\partial_y^s(\mathbf b-\mathbf b_{PQ})\|_2
\le\varepsilon_{rs}^{\rm par}(t).
\]

取 \((r,s)=(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)\)，就得到状态、两个一阶响应和三个二阶响应的界。全部系数由 \(P,Q,\Delta E,\Delta\eta,\tau\) 和已知门控决定，不需要先知道峰的位置。

## 5. 转换成保真度、梯度与曲率的界

把保真度写成固定的二次多项式

\[
F(\mathbf b)=\frac12\mathbf b^TH\mathbf b+\ell^T\mathbf b,
\quad H=\operatorname{diag}(1,\tfrac12,\tfrac12,\tfrac12,
\tfrac12,\tfrac12,\tfrac12,0,0,0),\quad
\ell=(-\tfrac12,0,\ldots,0)^T.
\]

真实状态上 \(\|\nabla F\|_2\le3/2\)，\(\|H\|_2=1\)。令 \(e_0=\varepsilon_{00}\)、\(e_x=\varepsilon_{10}\)、\(e_y=\varepsilon_{01}\)、\(e_{xx}=\varepsilon_{20}\)、\(e_{xy}=\varepsilon_{11}\)、\(e_{yy}=\varepsilon_{02}\)。利用式 (5)，得到

\[
\delta_0=\tfrac32e_0+\tfrac12e_0^2,
\qquad
\delta_i=\tfrac32e_i+e_0(\chi_i+e_i),\quad i\in\{x,y\},
\tag{12}
\]

\[
\boxed{
\delta_{ij}=\tfrac32e_{ij}
+e_0(\chi_i\chi_j+e_{ij})
+\chi_j e_i+\chi_i e_j+e_i e_j.
}
\tag{13}
\]

这里 \(\chi_x=\chi_E,\chi_y=\chi_\eta\)，分别有

\[
|F-\widetilde F|\le\delta_0,\qquad
|F_i-\widetilde F_i|\le\delta_i,\qquad
|F_{ij}-\widetilde F_{ij}|\le\delta_{ij}.
\]

它们在整个参数矩形上成立。换回物理导数时，例如

\[
\delta_{F_\eta}=\delta_y/\Delta\eta,\quad
\delta_{F_{\eta\eta}}=\delta_{yy}/\Delta\eta^2,\quad
\delta_{F_{E\eta}}=\delta_{xy}/(\Delta E\Delta\eta).
\tag{14}
\]

## 6. 对拟定参数区域，阶数需要多高

前一篇拟定

\[
E\in[0.005,0.015],\qquad\eta\in[0.010,0.040]\quad(\mathrm{meV}),
\]

所以 \(\Delta E=0.005\)、\(\Delta\eta=0.015\,\mathrm{meV}\)。仅在本次数值示例取 \(\tau=100/\mathrm{meV}\)，得到终点的

\[
\chi_E=\chi_\eta=6.
\]

令 \(P=Q=N\)。以下为式 (11)–(13) 的先验上界，导数均相对于无量纲 \(x,y\)，并假设时间求解精确。表中数值向上保留三位有效数字。

| \(N\) | 状态误差 | \(F\) 误差 | \(F_x,F_y\) 误差 | \(F_{xx},F_{yy}\) 误差 | \(F_{xy}\) 误差 |
|---:|---:|---:|---:|---:|---:|
| 8 | \(2.93\) | \(8.67\) | \(3.46\times10^2\) | \(1.19\times10^4\) | \(1.92\times10^4\) |
| 12 | \(2.32\times10^{-2}\) | \(3.50\times10^{-2}\) | \(1.96\) | \(7.78\times10^1\) | \(1.63\times10^2\) |
| 16 | \(4.82\times10^{-5}\) | \(7.22\times10^{-5}\) | \(6.58\times10^{-3}\) | \(4.12\times10^{-1}\) | \(9.25\times10^{-1}\) |
| 20 | \(3.69\times10^{-8}\) | \(5.53\times10^{-8}\) | \(7.50\times10^{-6}\) | \(6.93\times10^{-4}\) | \(1.62\times10^{-3}\) |
| 24 | \(1.27\times10^{-11}\) | \(1.90\times10^{-11}\) | \(3.60\times10^{-9}\) | \(4.63\times10^{-7}\) | \(1.10\times10^{-6}\) |

低阶上界可能很松。表中某阶未达到目标，不证明该阶的实际误差一定很大；24 阶达到目标，也不是说它是最小必要阶数。

20 阶的保真度上界已经低于首轮 \(10^{-7}\) 目标，但二阶导数的上界尚未达到 \(10^{-5}\)。24 阶的参数截断界同时达到这两个目标。这为首次求解选择 24 阶、并与较低阶比较提供了依据。

计算脚本为 [derive_global_error_bounds.py](derive_global_error_bounds.py)，原始数值为 [bounds_at_tau100.json](bounds_at_tau100.json)。从项目根目录运行：

```bash
.venv/bin/python -B global_parameter_geometry/derive_global_error_bounds.py
```

脚本使用 70 位十进制精度计算公式；JSON 中存储的是舍入后的数值字符串，不是区间算术结果。它还核对了 Gauss 正交变换、完整十分量生成元的相似变换和初值变换，并用已知导数的旋转函数检查余项公式的实现。抽样检查不替代上面的解析证明。

## 7. 时间积分误差必须单独进入总误差

设实际节点计算 \(\widehat{\mathbf z}_{ij}\) 的误差已被界定为

\[
\max_{i,j}\|\widehat{\mathbf z}_{ij}(t)-\mathbf b(x_i,y_j,t)\|_2\le\nu(t).
\]

由张量插值的算子界，总误差可以取

\[
\boxed{
\varepsilon_{rs}^{\rm total}
\le\varepsilon_{rs}^{\rm par}+C_{P,r}C_{Q,s}\nu(t).
}
\tag{15}
\]

例如 \(P=Q=24\) 时，\(C_{24,0}=25\)，\(C_{24,1}\approx4503.33\)。最坏情况下，重建混合导数的时间误差包络会显著放大。这是保守的算子界，实际误差可能更小且具有参数相关性。

因此不能直接把既有求解器的 \(\mathrm{rtol}\) 代入 \(\nu\)，也不能把第 6 节的参数截断表称作最终数值解的完整认证。首次计算仍沿用已定时间设置；需要测量并界定其影响。

## 8. 求解之后，用边缘系数形成六个标量误差方程

先验上界用于选阶；实际系数可提供另一条误差控制路线。令

\[
w_{m,r}=\sqrt{2m+1}\,d_{m,r}.
\]

由前一篇的边缘残差公式，对精确时间求解的谱系统，可取无量纲残差导数包络

\[
\begin{aligned}
R_{rs}^{\rm edge}={}&
\Delta E\,a_{P+1}w_{P+1,r}
\sum_{n=0}^Qw_{n,s}\|M_E\mathbf c_{P,n}\|_2\\
&+\Delta\eta\,a_{Q+1}w_{Q+1,s}
\sum_{m=0}^Pw_{m,r}\|M_\eta\mathbf c_{m,Q}\|_2.
\end{aligned}
\tag{16}
\]

实际时间重建若有模式内部残差 \(r_{mn}^{\rm time}\)，再加上

\[
\sum_{m=0}^P\sum_{n=0}^Q
w_{m,r}w_{n,s}\|r_{mn}^{\rm time}\|_2
\]

作为 \(R_{rs}^{\rm time}\)。记 \(R_{rs}=R_{rs}^{\rm edge}+R_{rs}^{\rm time}\)。取

\[
\gamma_E=2\Delta E,\qquad\gamma_\eta(t)=2\Delta\eta g(t).
\]

准确初值和连续接口下，六个误差包络可由以下非负标量系统产生：

\[
\boxed{
\begin{aligned}
\dot\varepsilon_{00}&=R_{00},\\
\dot\varepsilon_{10}&=\gamma_E\varepsilon_{00}+R_{10},\\
\dot\varepsilon_{01}&=\gamma_\eta\varepsilon_{00}+R_{01},\\
\dot\varepsilon_{20}&=2\gamma_E\varepsilon_{10}+R_{20},\\
\dot\varepsilon_{11}&=\gamma_E\varepsilon_{01}+\gamma_\eta\varepsilon_{10}+R_{11},\\
\dot\varepsilon_{02}&=2\gamma_\eta\varepsilon_{01}+R_{02},
\end{aligned}\qquad\varepsilon_{rs}(0)=0.}
\tag{17}
\]

式 (16)–(17) 将“需要一个误差界”推进为“从边缘谱系数、时间残差和已知门控计算误差界”。空间包络已经显式给出；连续时间上的残差和积分仍需受到控制。若只是时间采样，就应报告为误差估计，而不是严格保证。

## 9. 峰线的几何含义：靠近目标、相切与弯曲

将真实状态嵌入 Euclidean 矩阵空间：

\[
Z(E,\eta,t)=\frac1{\sqrt8}(P,\mathcal B),\qquad
Z_\star=\frac1{\sqrt8}(P_\star,\mathcal B_\star),\qquad
\Delta Z=Z-Z_\star.
\]

由原保真度的精确距离公式，

\[
F=1-\|\Delta Z\|^2.
\]

定义

\[
G_{ij}=\langle Z_i,Z_j\rangle,\qquad
W_{ij}=\langle\Delta Z,Z_{ij}\rangle,
\quad i,j\in\{E,\eta\}.
\]

直接求导得到

\[
\boxed{F_i=-2\langle\Delta Z,Z_i\rangle,
\qquad F_{ij}=-2(G_{ij}+W_{ij}).}
\tag{18}
\]

沿 \(\eta\) 扫描的峰满足 \(\langle\Delta Z,Z_\eta\rangle=0\)：状态轨迹在该点与以目标为中心的等距离球面相切。曲率条件为 \(G_{\eta\eta}+W_{\eta\eta}>0\)。这里 \(W\) 是状态曲面在所选坐标下的二阶变化与目标偏差的配对，不能在一般有限保真度处任意丢掉。

峰线的倾斜因此满足

\[
\boxed{\frac{d\eta_*}{dE}
=-\frac{G_{E\eta}+W_{E\eta}}
{G_{\eta\eta}+W_{\eta\eta}}.}
\tag{19}
\]

它把峰线走向与两个参数如何改变完整状态直接联系起来。沿峰线的高度二阶变化为

\[
H''(E)=-2\left[
G_{EE}+W_{EE}
-\frac{(G_{E\eta}+W_{E\eta})^2}{G_{\eta\eta}+W_{\eta\eta}}
\right].
\tag{20}
\]

特别地，在一个真正的 \(F=1\) 点，\(\Delta Z=0\)，因而 \(W=0\)、\(\operatorname{Hess}F=-2G\)。若 \(G_{\eta\eta}>0\)，

\[
\boxed{
\frac{d\eta_*}{dE}=-\frac{G_{E\eta}}{G_{\eta\eta}},\qquad
H''=-2\frac{\det G}{G_{\eta\eta}}.
}
\tag{21}
\]

当 \(\det G>0\)，两个参数引起的状态响应线性独立。调节 \(\eta\) 可以沿条件峰线作最佳的一阶补偿，但不能完全抵消 \(E\) 引起的状态变化，峰高随离开这个完美点而下降。若 \(\det G=0\)，只说明存在一阶状态补偿方向，尚不足以证明存在一整条 \(F=1\) 曲线。

式 (18)–(20) 在任意正则有限参数点成立。式 (21) 是在实际 \(F=1\) 点的结论，不是在 \(E=\eta=0\) 附近的小参数假设，也不宣称本模型已找到这样的点。

## 10. 怎样保证整段峰线而不是几个离散峰

使用无量纲坐标，在矩形带

\[
S=[x_a,x_b]\times[y_-,y_+]
\]

上，若已有覆盖此区域的 \(\delta_y,\delta_{yy}\)，并能验证

\[
\widetilde F_y(x,y_-) > \delta_y,\qquad
\widetilde F_y(x,y_+) < -\delta_y
\quad\text{对所有 }x\in[x_a,x_b],
\]

\[
\widetilde F_{yy}(x,y)+\delta_{yy}\le-\kappa<0
\quad\text{对所有 }(x,y)\in S,
\tag{22}
\]

则每个 \(x\) 都恰有一个真实峰，形成一条光滑分支 \(y_*(x)\)。对同一带内的近似导数根 \(\widetilde y_*\)，

\[
|\eta_*-\widetilde\eta_*|
\le\Delta\eta\frac{\delta_y}{\kappa}.
\tag{23}
\]

数值求根残差需要另加到分子中。峰高界可使用前一篇的二次位置误差公式。

有限阶重建下，\(\widetilde F\) 是次数至多为 \((2P,2Q)\) 的双变量多项式。把它的导数在每个矩形带上转换为 Bernstein 基底，基函数非负且和为 1，因此多项式值被所有 Bernstein 系数的最小值和最大值夹住。这为式 (22) 的连续区域符号检查提供了一个具体办法；界太松时，可细分矩形。

浮点计算 Bernstein 系数时仍需计入舍入误差。式 (22) 是待验证的充分条件，目前没有宣称本模型在整个拟定矩形上都满足它。遇到退化点时，应转而检查 \(F_y=F_{yy}=0\)，而不是继续使用带有小分母的正则延拓公式。

## 11. 推导后的下一项实施任务

现在已有一条从已知模型到可检查峰线的完整链条：节点演化给出谱系数，先验界辅助选阶，边缘残差给出实际误差包络，多项式导数和曲率符号确定峰线，状态度量与目标偏差解释倾斜和高度。

首次实施可优先比较 20 阶和 24 阶的六段求解，沿用既有时间设置。在同一参数域内，分别检查保真度与二阶导数的实际收敛，再寻找一段能够满足式 (22) 的正则峰线。时间误差、峰线存在性和位置仍需通过这次实施来确定；第 6 节已经得到的是求解之前的参数截断估计。
