# 时间与参数的偏微分：完整 SO(5) 保真度的精确响应

本文件讨论参数响应和局部峰线。仅含 F 的高阶消元方程、全局群轨道输运 PDE，以及三变量标量闭合的边界，见 [so5_global_scalar_pde.md](so5_global_scalar_pde.md)。

可以把现有群演化写成时间和参数上的偏微分系统。对本模型最有用的是精确的参数变分方程：在任意实际 E₁、|t₁| 上，传播完整保真度，同时求梯度、Hessian 和高保真度峰线的移动。

这不要求 E₁τ 或 |t₁|τ 小，也不要求暗模投影。求到二阶偏导仅限制了所报告的灵敏度阶数；作为基础的物理传播没有作二阶截断。

需要区分：静态参数之间没有物理扩散，单纯把所有参数的传播称为 PDE 不会自动简化计算。真正提供分析能力的是下面的参数微分关系、群几何约束和峰线方程。

## 1. 完整参数场

记 E=E₁、η=|t₁|，固定协议长度 τ 和控制包络 g(t)。在参数平面内部 η>0 上研究光滑场

\[
\mathcal B(t,E,\eta)=R(t,E,\eta)L_{12}R(t,E,\eta)^T.
\]

令

\[
K(t,E,\eta)=K_c(t)+E K_E+\eta K_\eta(t),
\quad K_E=\frac2\hbar L_{12},\quad K_\eta(t)=\frac2\hbar g(t)L_{1b}.
\]

其中 Kc 包含 t₂、t₃、Ed，采用前文和论文的门控。完整演化为

\[
\boxed{\partial_t\mathcal B=[K,\mathcal B],\qquad
\mathcal B(0,E,\eta)=L_{12}.}
\]

保真度直接由这个场读取：

\[
\boxed{F(t,E,\eta)=\frac14\left[
\sum_{i=1}^2\sum_{j=1}^5\mathcal B_{ij}^2-2\mathcal B_{12}
\right].}
\]

取终点 t=6τ 就是论文双交换的平方重叠，模长为 f=√F。该参数场保留完整的非对易、有限角度和有限控制速度。

η=0 是幅度域的边界，可取右导数，或者先将 η 作为有符号耦合系数延拓再求导；不能把复隧穿振幅的所有实分量与这个幅度参数混同。

## 2. 精确的一阶与二阶响应系统

定义

\[
\mathcal B_E=\partial_E\mathcal B,\qquad
\mathcal B_\eta=\partial_\eta\mathcal B,
\]

以及 \(\mathcal B_{EE},\mathcal B_{E\eta},\mathcal B_{\eta\eta}\)。对完整主方程求偏导，得到

\[
\boxed{
\begin{aligned}
\partial_t\mathcal B_E
&=[K,\mathcal B_E]+[K_E,\mathcal B],\\
\partial_t\mathcal B_\eta
&=[K,\mathcal B_\eta]+[K_\eta,\mathcal B].
\end{aligned}}
\]

由于 K 对两个幅度是仿射的，二阶参数导数的源项尤其简单：

\[
\boxed{
\begin{aligned}
\partial_t\mathcal B_{EE}
&=[K,\mathcal B_{EE}]+2[K_E,\mathcal B_E],\\
\partial_t\mathcal B_{E\eta}
&=[K,\mathcal B_{E\eta}]+[K_E,\mathcal B_\eta]+[K_\eta,\mathcal B_E],\\
\partial_t\mathcal B_{\eta\eta}
&=[K,\mathcal B_{\eta\eta}]+2[K_\eta,\mathcal B_\eta].
\end{aligned}}
\]

所有参数导数字段的初值都是零。每个字段是一个实反对称矩阵，因此只需十个实坐标；主字段与五个响应字段共六十个独立坐标。

这些是混合偏导的演化关系，也可在固定参数点作为一组耦合 ODE 求解。在参数空间中，每一点的 K 都取该点的实际参数，不用零扰动点的 K₀ 替代。

若初始态也随参数变化，则相应的导数字段初值不再是零；本分析沿用论文指定的固定初态和目标态。

## 3. 从字段偏导精确得到保真度偏导

定义双线性配对

\[
\langle X,Y\rangle_{12}=\sum_{i=1}^2\sum_{j=1}^5X_{ij}Y_{ij}.
\]

对 p,q∈{E,η}，有

\[
\boxed{F_p=\frac12\left[
\langle\mathcal B,\mathcal B_p\rangle_{12}-(\mathcal B_p)_{12}
\right],}
\]

\[
\boxed{F_{pq}=\frac12\left[
\langle\mathcal B_p,\mathcal B_q\rangle_{12}
+\langle\mathcal B,\mathcal B_{pq}\rangle_{12}-(\mathcal B_{pq})_{12}
\right].}
\]

它们可以回答：

- FE、Fη：在当前参数处增加哪个耦合使 F 上升或下降。
- FEE、Fηη：各自的斜率变化有多快。
- FEη：改变 t₁ 后，E₁ 的作用如何变化，以及相反方向的变化。
- Hessian 的主方向：在同一能量单位下，哪些联合参数变化对应较大的局部曲率。

混合导数是联合作用的量度，但不能单凭它非零就断言全部来源是生成元不对易；观测量的几何和相位干涉也会产生混合响应。

Hessian 在非驻点不能单独判定极大或极小。即使参数本身不小，导数也只是该点的局部信息；若用该点的有限阶 Taylor 多项式代替整个平面，仍会回到先前的局部近似问题。

## 4. 如何从局部导数研究任意有限参数变化

对参数平面中的任意光滑曲线 p(s)=(E(s),η(s))，有精确恒等式

\[
\boxed{F(p(1))-F(p(0))
=\int_0^1\left[F_E(p(s))E'(s)+F_\eta(p(s))\eta'(s)\right]ds.}
\]

这不要求两个端点相近，但导数必须沿曲线取其真实值。因此参数响应方程可以用来延拓整个完整保真度，而不是只在原点展开一次。

固定保真度的规则等值线满足

\[
\frac{d\eta}{dE}=-\frac{F_E}{F_\eta}
\qquad(F_\eta\neq0).
\]

对固定 E、固定 τ 的局部峰分支 η*(E)，条件是

\[
F_\eta(E,\eta_*(E))=0,\qquad F_{\eta\eta}<0.
\]

隐函数定理给出峰线方程

\[
\boxed{\frac{d\eta_*}{dE}=-\frac{F_{E\eta}}{F_{\eta\eta}}.}
\]

它能够描述高保真度峰的位置随 E 移动。若 Fηη 接近零，峰线可能变平、转折或发生分支变化，应更换延拓方式，不能继续除以这个接近零的量。

这条峰线仅保证沿 η 方向极大；要得到整个二维参数平面的局部极大，还需要 FE=0 并检查完整 Hessian。

## 5. 群参数联络：真正的参数偏微分形式

定义右对数导数

\[
J_E=(\partial_E R)R^T,\qquad J_\eta=(\partial_\eta R)R^T.
\]

有完整的时间—参数联立系统

\[
\boxed{
\begin{aligned}
\partial_t\mathcal B&=[K,\mathcal B],\\
\partial_E\mathcal B&=[J_E,\mathcal B],\\
\partial_\eta\mathcal B&=[J_\eta,\mathcal B],\\
\partial_tJ_p&=K_p+[K,J_p],\qquad J_p(0)=0.
\end{aligned}}
\]

参数混合偏导的相容性是 Maurer–Cartan 关系

\[
\boxed{\partial_EJ_\eta-\partial_\eta J_E-[J_E,J_\eta]=0.}
\]

它保留非对易作用，同时保证同一个光滑 R(E,η) 的两种偏导次序一致。整个未投影群联络是纯规范的平坦联络；不能把 [JE,Jη] 单独非零称为非零参数曲率。投影到某个瞬时子空间后的 Berry 曲率是另一种对象。

这套形式更适合研究参数变化怎样旋转有向平面，也避免直接跟踪旋转角在简并点的标签交换和分支跳变。

## 6. 为什么不能仅给 F 加一个通用标量 PDE

主方程使 F 依赖于有向平面的位置及其相干分量。两个不同的平面可以给出相同的 F，却具有不同的后续 F 变化。因此，仅知道当前标量 F，一般不能构成对所有态都成立的一阶闭合时间方程。这个论证并不排除高阶标量方程、带记忆的方程或含参数偏导的特殊闭合；不能据此断言所有 F 的标量 PDE 都不可能。

本模型也没有要求相邻 E、η 之间存在扩散；未经推导加入 FEE+Fηη 等项会改变问题。正确做法是保留 \(\mathcal B\) 或等价的波函数投影，建立上述精确偏导系统，再提取 F。只有额外指定了参数噪声、随机过程或优化问题，才会出现具有其他物理来源的 PDE。

## 7. 一个有限参数点的核对

使用论文门控、tc=E₀=0.3 meV、ℏ=1，完整六段传播给出：τ=100 meV⁻¹、E=η=0.01 meV 时，

\[
F=0.009487011275,
\quad \nabla F=(-25.20181844,\ 10.67079341)\ \mathrm{meV}^{-1},
\]

\[
H_F=\begin{pmatrix}
48414.25118&-20482.39333\\
-20482.39333&8660.58066
\end{pmatrix}\ \mathrm{meV}^{-2}.
\]

在这个实际参数点，增加 E 降低 F，增加 η 提高 F；负的混合导数表示增加 η 会使 E 的斜率进一步下降。这不是关于零扰动点的展开。

对同一个 E、τ，在原来二阶抵消线附近求 Fη=0，得到

\[
\boxed{\eta_*=0.02300051168\ \mathrm{meV},\quad
\eta_*/E=2.300051168,\quad F_*=0.889373594.}
\]

此处 Fηη≈−32638.80390 meV⁻²，FEη≈73060.03551 meV⁻²，因此

\[
\boxed{d\eta_*/dE\simeq2.238440959.}
\]

峰线位置的比值 2.30005 与峰线斜率 2.23844 不同，直接表明峰线在这个区域不是简单的固定比例直线。这里 FE≈−36.48563 meV⁻¹，因此它也不是整个二维平面的驻点。

上述方程由 `so5_parameter_response.py` 实现，原始数据写入 `so5_parameter_response.json`。核对方式包括：

- 有向平面的一阶响应与独立五维传播矩阵的一阶响应一致，最大残差约 1.4×10⁻¹² meV⁻¹。
- 二阶混合响应与独立一阶响应的中心差分一致，Hessian 最大绝对残差约 3.3×10⁻⁶ meV⁻²。
- 峰线微分方程的斜率与在 E±10⁻⁷ meV 处重新求局部峰得到的斜率一致。
- 验证了反对称性、双向量范数与 \(\mathcal B^3=-\mathcal B\) 的群轨道约束。

这些是有限参数响应及局部峰线的核对，不代表已经延拓出整个参数平面的所有峰分支，也不代表得到了一般余弦协议的初等闭式解。

## 8. 在 Fig. 1(d) 坐标中的导数

论文以 y=lg(η/E) 为纵轴。在 E>0、η>0 时，η=E10ʸ，坐标变换给出

\[
\left.\partial_yF\right|_E=\ln10\,\eta F_\eta,
\qquad
\left.\partial_y^2F\right|_E=(\ln10)^2(\eta F_\eta+\eta^2F_{\eta\eta}),
\]

\[
\left.\partial_EF\right|_y=F_E+(\eta/E)F_\eta.
\]

因此固定 t₁ 和固定 t₁/E₁ 研究 E₁ 的影响，是两种不同的变化路径。偏导分析必须先说明固定的是哪一个坐标。
