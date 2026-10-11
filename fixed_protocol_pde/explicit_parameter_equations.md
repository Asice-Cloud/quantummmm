# E₁、t₁ 的显式偏导方程与闭合范围

上一版的 \(\partial_E\mathcal B=[\Xi_E,\mathcal B]\) 需要先求未知连接场；它没有给出仅在 \((E,\eta)\) 平面上、系数直接已知的独立 PDE。下面消去这层辅助表示，直接列出已知系数的分量方程。得到的是**混合时间—参数偏导系统**；单独终点保真度的参数平面 PDE 仍未得到。

当前采用的路线是保留真实演化时间，以完整 \(\mathcal B\) 为主变量。第 5 节将分量方程接成统一初值问题，并给出保真度及其参数响应与状态的精确耦合；这是接下来实际分析的系统。

从原始有效哈密顿量、四元数旋量矩阵和费米子 BdG 矩阵到本页系数的独立核对已完成，见 [effective_hamiltonian_audit.md](effective_hamiltonian_audit.md)。五项耦合的符号及物理时间因子 2 均与本页一致。

## 1. 固定设置与完整状态

取 \(\hbar=1\)，保留每段时长 \(\tau>0\) 及控制幅度 \(t_c,E_0\) 的符号。它们是给定的可控设计量；所有 \(E,\eta\) 偏导均在这些设计量给定时计算。每段使用 \(v=s/\tau\in[0,1]\)，定义

\[
A(v)=\frac{1-\cos\pi v}{2},\qquad C(v)=1-A(v).
\]

令 \(E=E_1,\eta=|t_1|\)，实际门控耦合为 \(q=\eta g(v)\)。全部已知系数为

| 段类型 | \(p\) | \(r\) | \(e_d\) | \(g\) |
|---|---|---|---|---|
| 第一段 | \(t_cA\) | \(0\) | \(E_0C\) | \(A\) |
| 第二段 | \(t_cC\) | \(t_cA\) | \(0\) | \(C\) |
| 第三段 | \(0\) | \(t_cC\) | \(E_0A\) | \(0\) |

三段重复两次。下面定义完整状态列向量

\[
\mathbf b=(d,x_1,x_2,x_3,y_1,y_2,y_3,h_1,h_2,h_3)^T,
\]

\[
\begin{aligned}
d&=\mathcal B_{12},&
(x_1,x_2,x_3)&=(\mathcal B_{13},\mathcal B_{1a},\mathcal B_{1b}),\\
(y_1,y_2,y_3)&=(\mathcal B_{23},\mathcal B_{2a},\mathcal B_{2b}),&
(h_1,h_2,h_3)&=(\mathcal B_{3a},\mathcal B_{3b},\mathcal B_{ab}).
\end{aligned}
\]

它就是完整 \(\mathcal B\) 的十个反对称分量，没有新增未知传播矩阵。

## 2. 系数完全展开的时间方程

\[
\boxed{\partial_v\mathbf b=2\tau N(E,\eta,v)\mathbf b,}
\]

其中

\[
\boxed{
N=\begin{pmatrix}
0&0&-p&0&0&0&-\eta g&0&0&0\\
0&0&r&0&E&0&0&0&-\eta g&0\\
p&-r&0&e_d&0&E&0&0&0&-\eta g\\
0&0&-e_d&0&0&0&E&0&0&0\\
0&-E&0&0&0&r&0&p&0&0\\
0&0&-E&0&-r&0&e_d&0&0&0\\
\eta g&0&0&-E&0&-e_d&0&0&0&-p\\
0&0&0&0&-p&0&0&0&e_d&0\\
0&\eta g&0&0&0&0&0&-e_d&0&r\\
0&0&\eta g&0&0&0&p&0&-r&0
\end{pmatrix}.}
\tag{A}
\]

所有 \(p,r,e_d,g\) 按表代入。例如第二段只有 \(p=t_cC,r=t_cA,e_d=0,g=C\)。式 (A) 中的 \(N\) 由不可控参数、已知设计量及 \(v\) 明确给出，不需要先算真实演化。物理时间满足 \(\partial_t\mathbf b=2N\mathbf b\)，换成 \(v\) 才出现因子 \(\tau\)。

## 3. 对两个参数的具体 PDE

直接求偏导，得到

\[
\boxed{
\left(\partial_v-2\tau N\right)\partial_E\mathbf b
=2\tau\begin{pmatrix}
0\\y_1\\y_2\\y_3\\-x_1\\-x_2\\-x_3\\0\\0\\0
\end{pmatrix},}
\tag{B}
\]

\[
\boxed{
\left(\partial_v-2\tau N\right)\partial_\eta\mathbf b
=2\tau g\begin{pmatrix}
-y_3\\-h_2\\-h_3\\0\\0\\0\\d\\0\\x_1\\x_2
\end{pmatrix}.}
\tag{C}
\]

式 (A)—(C) 就是系数明确的闭合混合偏导系统，未知量为 \(\mathbf b,\partial_E\mathbf b,\partial_\eta\mathbf b\)。同样可以仅把 \(\mathbf b(E,\eta,v)\) 看作未知场，把 (B)、(C) 中的导数直接理解为它的偏导。

这里没有 \(\Xi\)，也没有把参数作用留在未展开的对易子里。它对有限 \(E,\eta\) 精确成立。

### 3.1 把与定向最相关的分量单独写出

例如，第一行是

\[
\partial_vd=-2\tau(p x_2+\eta g y_3).
\]

两个参数方向因此满足

\[
\boxed{
\begin{aligned}
\partial_{vE}d
&=-2\tau\left(p\,\partial_Ex_2+\eta g\,\partial_Ey_3\right),\\
\partial_{v\eta}d
&=-2\tau\left(p\,\partial_\eta x_2+\eta g\,\partial_\eta y_3+g y_3\right).
\end{aligned}}
\tag{D}
\]

相干分量 \(x_2\) 的方程为

\[
\partial_vx_2=2\tau(p d-rx_1+e_dx_3+E y_2-\eta g h_3),
\]

所以

\[
\boxed{
\begin{aligned}
\partial_{vE}x_2=2\tau\big(&p d_E-r(x_1)_E+e_d(x_3)_E
+E(y_2)_E-\eta g(h_3)_E+y_2\big),\\
\partial_{v\eta}x_2=2\tau\big(&p d_\eta-r(x_1)_\eta+e_d(x_3)_\eta
+E(y_2)_\eta-\eta g(h_3)_\eta-g h_3\big).
\end{aligned}}
\tag{E}
\]

例如第二段把 \(p=t_cC,r=t_cA,e_d=0,g=C\) 代入 (D)：

\[
\begin{aligned}
\partial_{vE}d&=-2\tau C\left[t_c(x_2)_E+\eta(y_3)_E\right],\\
\partial_{v\eta}d&=-2\tau C\left[t_c(x_2)_\eta+\eta(y_3)_\eta+y_3\right].
\end{aligned}
\]

这明确显示两个参数的区别：\(\eta\) 在定向响应中有直接源项，\(E\) 对定向的影响通过相干响应传来。其余分量由 (A)—(C) 同时演化，不能只保留 (D) 就声称定向已经闭合。

## 4. 二阶参数偏导也具有明确系数

令 \(N_E=\partial_EN\)、\(N_\eta=\partial_\eta N\)。它们是 (B)、(C) 右端所列线性映射；\(N_E\) 为常矩阵，\(N_\eta\) 仅通过已知 \(g(v)\) 变化。

因为所有二阶参数导数 \(N_{EE},N_{E\eta},N_{\eta\eta}\) 为零，精确得到

\[
\boxed{
\begin{aligned}
(\partial_v-2\tau N)\mathbf b_{EE}&=4\tau N_E\mathbf b_E,\\
(\partial_v-2\tau N)\mathbf b_{E\eta}
&=2\tau(N_E\mathbf b_\eta+N_\eta\mathbf b_E),\\
(\partial_v-2\tau N)\mathbf b_{\eta\eta}
&=4\tau N_\eta\mathbf b_\eta.
\end{aligned}}
\tag{F}
\]

这些偏导方程可以计算实际参数点的梯度、曲率与峰线变化。使用导数本身不要求 \(E\tau,\eta\tau\) 小；若之后再用有限阶 Taylor 多项式替代整个参数区域，那才引入另外的局部近似。

## 5. 初值、接口与保真度

第一段初始数据为

\[
\mathbf b(E,\eta,0)=(1,0,0,0,0,0,0,0,0,0)^T,
\quad \mathbf b_E=\mathbf b_\eta=\mathbf b_{EE}
=\mathbf b_{E\eta}=\mathbf b_{\eta\eta}=0.
\]

每段 \(v=1\) 的全部状态和参数导数，连续传给下一段的 \(v=0\)。第二轮不重新置初态。六段的控制均按第 1 节表格取值。

从状态读出

\[
\boxed{F=\frac{d^2-d}{2}
+\frac14\sum_{j=1}^{3}(x_j^2+y_j^2),}
\]

\[
F_\lambda=(d-\tfrac12)d_\lambda
+\frac12\sum_{j=1}^{3}\left(x_j(x_j)_\lambda+y_j(y_j)_\lambda\right).
\]

终点为第六段 \(v=1\)，即全局时间 \(6\tau\)。完整记号为 \(F_{\rm end}(E,\eta\mid\tau,t_c,E_0)=F(E,\eta,6\tau\mid\tau,t_c,E_0)\)；下文省略可控量的条件记号。

### 5.1 统一的 30 分量初值问题

为避免在接口重新计时，现在改用全局物理时间 \(t\in[0,6\tau]\)。各段局部 \(v=(t-j\tau)/\tau\)，其中 \(j=0,\ldots,5\)。接口位于 \(t_j=j\tau\)，\(j=1,\ldots,5\)。定义

\[
M(E,\eta,t)=2N(E,\eta,v),\quad
M_E=2N_E,\quad M_\eta=2N_\eta.
\]

系数按段选取；这里的换元不改变每段时长或门控。记 \(u=\partial_E\mathbf b\)、\(w=\partial_\eta\mathbf b\)，则

\[
\boxed{
\partial_t
\begin{pmatrix}\mathbf b\\u\\w\end{pmatrix}
=\begin{pmatrix}
M&0&0\\
M_E&M&0\\
M_\eta&0&M
\end{pmatrix}
\begin{pmatrix}\mathbf b\\u\\w\end{pmatrix},\qquad
\begin{pmatrix}\mathbf b\\u\\w\end{pmatrix}_{t=0}
=\begin{pmatrix}\mathbf b_0\\0\\0\end{pmatrix}.}
\tag{G}
\]

三个 \(10\times10\) 系数块均由 (A) 显式给定。式 (G) 是状态演化及参数响应的一种统一表示；等价的混合偏导写法是

\[
\partial_t\mathbf b=M\mathbf b,\qquad
\partial_{tE}\mathbf b=M\mathbf b_E+M_E\mathbf b,\qquad
\partial_{t\eta}\mathbf b=M\mathbf b_\eta+M_\eta\mathbf b.
\]

它作为时间初值问题已经闭合。\(u,w\) 不是任意附加的状态：由对主方程求导、零初值及解的唯一性，它们正是同一状态场的两个参数偏导。不能在每个段的起点重新置零。

严格地说，主方程是参数化的 ODE 族，联合响应是混合偏导系统。这里不假定参数平面中存在另外的物理输运或扩散。

若使用全局归一化时间 \(\zeta=t/\tau\in[0,6]\)，则各段局部 \(v=\zeta-j\)，状态主方程的系数通过 \(\tau E,\tau\eta,\tau t_c,\tau E_0\) 进入；对原幅度 \(E,\eta\) 求导时，响应源项相应带有因子 \(\tau\)。因而保留 \(\tau\) 可以明确区分可控的时间尺度与两个不可控幅度，同时保持归一化门控形状。

反对称性 \(M^T=-M\) 给出 \(\partial_t\|\mathbf b\|^2=0\)。结合原始对易子方程，它保持完整有向平面约束。十个全局分量对应六个独立状态自由度；在定向穿过零值的过程中不必强制使用某一个局部角度坐标图。

### 5.2 保真度的时间与参数响应方程

定义两个由当前状态直接计算的量

\[
\mathcal A(\mathbf b)=(1-d)x_2+h_1y_1-h_3y_3,
\qquad
\mathcal C(\mathbf b)=(1-d)y_3-h_3x_2-h_2x_1.
\]

将 (A) 代入保真度导数，全部 \(E,r,e_d\) 的显式项抵消，精确得到

\[
\boxed{\partial_tF=p(t)\mathcal A(\mathbf b)
+\eta g(t)\mathcal C(\mathbf b),\qquad F(E,\eta;0)=0.}
\tag{H}
\]

\(F\) 也可始终用第 5 节的代数式读出，无需增加独立的保真度未知量。消失的只是瞬时变化率中的显式项：\(E,r,e_d\) 仍通过 \(\mathbf b\) 改变随后演化的保真度。

对任意状态方向 \(s\)，两个量的微分为

\[
\begin{aligned}
D\mathcal A[s]={}&-s_d x_2+(1-d)s_{x_2}
+s_{h_1}y_1+h_1s_{y_1}-s_{h_3}y_3-h_3s_{y_3},\\
D\mathcal C[s]={}&-s_d y_3+(1-d)s_{y_3}
-s_{h_3}x_2-h_3s_{x_2}-s_{h_2}x_1-h_2s_{x_1}.
\end{aligned}
\]

因此还有完全明确的混合偏导方程

\[
\boxed{
\begin{aligned}
\partial_{tE}F&=p\,D\mathcal A[u]+\eta g\,D\mathcal C[u],\\
\partial_{t\eta}F&=p\,D\mathcal A[w]+\eta g\,D\mathcal C[w]+g\mathcal C.
\end{aligned}}
\tag{I}
\]

式 (G) 确定 \(\mathbf b,u,w\)，式 (H)—(I) 给出保真度及其响应。这是一套具体的、保留完整状态的联合描述，没有未指定的函数或待拟合系数。第三段 \(p=g=0\)，所以 \(F,F_E,F_\eta\) 均保持不变，而完整状态继续演化并传给下一周期。

分析时可以比较源项 \(M_E\mathbf b\)、\(M_\eta\mathbf b\) 怎样建立不同相干响应，再用 (I) 判断这些响应怎样影响保真度。最后在给定设计的终点 \(t=6\tau\) 读出 \(F_{\rm end}\) 及其梯度。若要研究条带的局部曲率，再接入 (F)；建立主系统不需要先做角度约化。

## 6. 它与“只在参数平面上的 PDE”相差什么

式 (B)、(C) 含有 \(\partial_{vE}\)、\(\partial_{v\eta}\)，而非一个独立的终点参数空间算子。它们由时间方程求导得到，刻画参数响应；并未获得一个只含

\[
F_{\rm end},\quad\partial_EF_{\rm end},\quad\partial_\eta F_{\rm end},
\quad\partial_{EE}F_{\rm end},\quad\partial_{E\eta}F_{\rm end},\ldots
\]

且所有系数都由固定设计直接给出的闭合关系。

因此需要修正此前的表述：**已经得到精确的状态演化和显式参数响应 PDE；尚未得到独立闭合的 \(F_{\rm end}(E,\eta)\) 参数平面 PDE。** 不能把二者当作同一成果，也不能通过增加未知连接场就声称后一个问题已解决。

若继续要求后一种方程，需要从原系统消去时间历史和隐藏相干分量，证明所选参数偏导之间存在有限的闭合关系。当前没有完成这个消元。保真度丢失状态信息会阻碍直接的一阶闭合，但这不足以证明任何高阶参数 PDE 都不可能。

## 7. 系数的代数核对

[derive_explicit_equations.py](derive_explicit_equations.py) 用五维完整对易子独立生成十个分量，保留正符号 \(\tau\)，核对了 (A) 的全部矩阵元素、反对称性、(B)/(C) 的全部参数源项及二阶系数导数为零。另外直接对保真度代数式求导，核对了 (H)、归一化时间关系 \(\partial_vF=\tau\partial_tF\) 及 (I) 的显式参数源项。结果见 [explicit_equations.json](explicit_equations.json)。

```sh
.venv/bin/python -B fixed_protocol_pde/derive_explicit_equations.py
```

此前 [verify_pde.py](verify_pde.py) 对完整状态及一、二阶参数响应所做的独立波函数核对，取 \(\tau=100,t_c=E_0=0.3\)，是这里符号系统的一个具体设计实例。符号推导不依赖这些数值；此次只恢复可控量的符号，不重新扫描参数或更改归一化门控。
