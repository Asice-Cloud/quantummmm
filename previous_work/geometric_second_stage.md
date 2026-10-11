# 第二段：瞬时旋转几何、相干记忆与定向返回

承接 [geometric_coherence_dynamics.md](geometric_coherence_dynamics.md)。这里单独分析第二段，重点是第二周期的 \([4\tau,5\tau]\)：为什么进入这一段时都接近 \(d=0\)，不同 \(E_1,|t_1|\) 却可以使终点定向相反？

第二段可以得到两类精确结果：从哈密顿量直接求出瞬时零轴和两种旋转频率；将十个相干分量分为两组，得到核完全已知的六分量记忆方程。后者包含实际状态和保真度，不要求先计算完整 \(R\)。以下没有小参数展开，也没有假设状态绝热跟随零轴。

## 1. 第二段的控制与输入态

把这一段的局部物理时间记为 \(0\leq t\leq\tau\)，定义

\[
E=E_1,\qquad \eta=|t_1|\text{ 的控制幅度},\qquad
\kappa=2/\hbar,
\]

\[
g(t)=\frac{1+\cos(\pi t/\tau)}2,\qquad
p(t)=t_cg(t),\quad q(t)=\eta g(t),\quad
r(t)=t_c[1-g(t)],\quad E_d=0.
\]

生成元是

\[
K=\kappa A,\qquad A=EL_{12}+qL_{1b}-pL_{2a}+rL_{3a}.
\]

状态仍用 \(\mathcal B\) 表示，满足 \(\dot{\mathcal B}=[K,\mathcal B]\)。进入第二段的 \(\mathcal B_{\rm in}\) 由此前的演化决定，不能只保留 \(d_{\rm in}\)，也不能把其他相干重新置零。

直接定向驱动为

\[
\boxed{\dot d=-\kappa g(t)
[t_c\mathcal B_{1a}(t)+\eta\mathcal B_{2b}(t)].}
\]

这一段的控制在两次周期中相同，但输入态通常不同。本文中的输入交换对照都在第二周期进行。

## 2. 可以直接求出的瞬时旋转几何

将基底重新排列为 \((1,a)\mid(2,3,b)\)，则

\[
A=\begin{pmatrix}0&M\\-M^T&0\end{pmatrix},\qquad
M=\begin{pmatrix}E&0&q\\p&-r&0\end{pmatrix}.
\]

所有耦合均跨越这两个基底子空间。它们并非两个独立演化子系统。

### 2.1 瞬时零轴

在 \(\operatorname{rank}M=2\) 时，\(A\) 有一条零轴。在原基底 \((1,2,3,a,b)\) 中可以直接写成

\[
\boxed{
n(t)=\frac{(0,qr,pq,0,-Er)^T}
{\sqrt{q^2(p^2+r^2)+E^2r^2}},\qquad A(t)n(t)=0.
}
\]

零轴是无定向的直线，可以用投影 \(N=nn^T\) 描述。这是生成元的几何结构，并不等于实际波函数或有向状态平面 \(\mathcal B\)。

对于 \(q,r\neq0\)，其分量比为

\[
\frac{n_3}{n_2}=\frac p r,\qquad
\boxed{\frac{n_b}{n_2}=-\frac E q=-\frac{E}{\eta g(t)}.}
\]

第一项由既定门控决定；第二项显式揭示有限参数改变零轴的方式。若 \(\eta>0\)，段初零轴为模式 3；当 \(E\neq0\)，段末零轴为模式 \(b\)。若 \(E=0\)，段内零轴沿 \((2,3)\) 平面旋转，趋向模式 2；恰好在段末，零空间扩为三维。

因此零轴的两个极限不交换：

\[
\lim_{E\to0}\lim_{t\to\tau^-}N(E,t)=e_be_b^T,
\qquad
\lim_{t\to\tau^-}\lim_{E\to0}N(E,t)=e_2e_2^T
\quad(\eta>0).
\]

这说明不能把 \(E=0\) 的瞬时零轴图像直接用于有限 \(E\) 的整个第二段。它不意味着有限时长的实际保真度在 \(E=0\) 不连续；原动力学系数光滑，有限时间的传播及保真度仍光滑。

### 2.2 两种旋转频率

\[
MM^T=\begin{pmatrix}E^2+q^2&Ep\\Ep&p^2+r^2\end{pmatrix}.
\]

其两个本征值为

\[
\boxed{
\mu_\pm^2=\frac12\left[S\pm
\sqrt{(E^2+q^2-p^2-r^2)^2+4E^2p^2}\right],
\qquad S=E^2+q^2+p^2+r^2.
}
\]

\(K\) 的谱是 \(0,\pm i\kappa\mu_+,\pm i\kappa\mu_-\)。有限 \(E\) 通过 \(Ep\) 混合两种瞬时旋转；\(\eta\) 通过 \(q=\eta g\) 改变频率与零轴方向。\(E=0\) 时，\((1,b)\) 与 \((2,3,a)\) 两个固定子空间的旋转确实解耦。

一般参数下瞬时旋转平面随时间变化，\(A(t)\) 与 \(A(s)\) 通常不对易。因此只积分 \(\mu_\pm\)，或只画 \(n(t)\)，还不能得到终点保真度。实际状态的输运需要下面的相干方程。

## 3. 实际状态的四分量与六分量分解

定义列向量

\[
X=\begin{pmatrix}
\mathcal B_{1a}\\\mathcal B_{23}\\\mathcal B_{2b}\\\mathcal B_{3b}
\end{pmatrix},\qquad
Y=\begin{pmatrix}
d\\\mathcal B_{13}\\\mathcal B_{1b}\\-\mathcal B_{2a}\\-\mathcal B_{3a}\\\mathcal B_{ab}
\end{pmatrix}.
\]

\(X\) 包含前述两个基底子空间内部的平面分量，\(Y\) 包含跨越两个子空间的分量。注意这里的 \(X,Y\) 与第一段的 \(X_\pm\) 不是同一个量。

由对易子逐项得到

\[
\boxed{\dot Y=\kappa L(t)X,\qquad
\dot X=-\kappa L(t)^TY,}
\]

\[
\boxed{
L(t)=\begin{pmatrix}
-p&0&-q&0\\
r&E&0&-q\\
0&0&E&0\\
E&r&0&0\\
0&p&0&0\\
q&0&p&-r
\end{pmatrix}.}
\]

这里仍包含十个矩阵分量；原有的纯态平面约束保持，独立状态自由度仍为六。\(X,Y\) 不是两个独立物理系统，也不是任意十个自由变量。

几何范数保持：

\[
\frac d{dt}(\|X\|^2+\|Y\|^2)=0,
\qquad \|X\|^2+\|Y\|^2=1.
\]

它们互相交换分量，而非产生耗散。保真度也能直接恢复：

\[
\boxed{F=\frac{d^2-d}{2}
+\frac{X_1^2+X_2^2+X_3^2+Y_2^2+Y_3^2+Y_4^2}{4}.}
\]

下标从 1 开始，\(Y_1=d\)。仍取 \(F=f^2\)，其中 \(f\) 是此前对应论文的重叠模。

## 4. 消去四个分量：核全部已知的精确记忆方程

设输入为 \(X_0,Y_0\)。积分第二个方程：

\[
X(t)=X_0-\kappa\int_0^tL(s)^TY(s)\,ds.
\]

代回第一个方程，得到

\[
\boxed{
\dot Y(t)=\kappa L(t)X_0
-\int_0^t\mathscr K(t,s)Y(s)\,ds,
\qquad
\mathscr K(t,s)=\kappa^2L(t)L(s)^T,
\qquad Y(0)=Y_0.
}
\]

这是第二段完整的六分量 Volterra 方程。核直接由 \(E,\eta,g,p,r\) 给定，不包含未知 \(R\) 或待求传播矩阵。解出 \(Y\) 后，用上一式恢复 \(X\)，从而恢复完整 \(\mathcal B\) 与 \(F\)。

它是含历史的闭合方程，不是六个任意变量的自治常微分方程，也不是仅对 \(F\) 的闭合标量方程。已知核不等于已经得到一般闭式解，但现在可以直接分析已知积分算子如何随参数改变。

\(\mathscr K(t,t)\) 半正定，不过这不表示单调衰减：这里的结构是保范数的相干反馈，而且 \(\mathscr K(t,s)\) 在不同时刻一般不对称。

## 5. 定向分量的精确记忆来源

取上一式的第一行，并使用 \(p=t_cg,q=\eta g\)，得到

\[
\boxed{
\begin{aligned}
\dot d(t)={}&-\kappa g(t)
[t_c\mathcal B_{1a}(0)+\eta\mathcal B_{2b}(0)]\\
&-\kappa^2g(t)\int_0^t\Big\{
g(s)[(t_c^2+\eta^2)d(s)-2t_c\eta\mathcal B_{ab}(s)]\\
&\hspace{35mm}-t_cr(s)\mathcal B_{13}(s)
+E[t_c\mathcal B_{2a}(s)-\eta\mathcal B_{1b}(s)]
\Big\}\,ds.
\end{aligned}}
\]

这条式子把实际机制分成四类：

1. 输入相干 \(\mathcal B_{1a}(0),\mathcal B_{2b}(0)\) 给出初始驱动。
2. \(d\) 与 \(\mathcal B_{ab}\) 通过 \(t_c^2+\eta^2\) 和 \(2t_c\eta\) 持续反馈；不是只有一个 \(\eta^2\) 的“误差强度”。
3. 打开的 \(r\) 通过 \(\mathcal B_{13}\) 接入模式 3 的相干。
4. \(E\) 把另一相干象限 \(\mathcal B_{2a},\mathcal B_{1b}\) 转入定向驱动，包含 \(E\eta\) 的交叉项。

右侧各个历史分量必须共同满足第 4 节的六分量方程；这一行单独不是 \(d\) 的闭合标量方程。

段初 \(g(0)=1,g'(0)=r(0)=0\)，因此还有精确局部关系

\[
\dot d(0)=-\kappa[t_c\mathcal B_{1a}(0)+\eta\mathcal B_{2b}(0)],
\]

\[
\ddot d(0)=-\kappa^2\left[
(t_c^2+\eta^2)d_0-2t_c\eta\mathcal B_{ab}(0)
+E(t_c\mathcal B_{2a}(0)-\eta\mathcal B_{1b}(0))\right].
\]

所以即使 \(d_0\approx0\)，输入相干仍决定初始速度与曲率。不过它们只描述段初，不能替代整个历史积分。

## 6. 终点定向的精确积分账本

定义剩余门控权重

\[
G(s)=\int_s^\tau g(t)dt
=\frac{\tau-s}{2}-\frac{\tau}{2\pi}\sin(\pi s/\tau)\geq0.
\]

交换积分次序可得

\[
\boxed{d_{\rm out}-d_{\rm in}
=D_{\rm in}+D_d+D_{ab}+D_3+D_E,}
\]

\[
\begin{aligned}
D_{\rm in}&=-\frac{\kappa\tau}{2}
(t_c\mathcal B_{1a}(0)+\eta\mathcal B_{2b}(0)),\\
D_d&=-\kappa^2(t_c^2+\eta^2)\int_0^\tau Gg\,d\,ds,\\
D_{ab}&=2\kappa^2t_c\eta\int_0^\tau Gg\,\mathcal B_{ab}\,ds,\\
D_3&=\kappa^2t_c\int_0^\tau Gr\,\mathcal B_{13}\,ds,\\
D_E&=-\kappa^2E\int_0^\tau G
(t_c\mathcal B_{2a}-\eta\mathcal B_{1b})\,ds.
\end{aligned}
\]

这是精确分解，允许检查哪几条相干反馈共同形成净驱动。但各项不是独立因果贡献：全部历史分量都受 \(E,\eta\) 影响，且各项可能大幅抵消。例如已有高保真度样本 \((E,\eta,\tau)=(0.01,0.023,100)\) 中，

\[
D_d\approx37.429421,
\qquad D_3\approx-37.078634,
\]

\[
D_{\rm in}\approx-0.004415,
\quad D_{ab}\approx-1.681043,
\quad D_E\approx+0.448593,
\]

总和才是 \(d_{\rm out}-d_{\rm in}\approx-0.886078\)。因此不能把某个正、负系数或某个单项的绝对值当成最终成功判据。

## 7. 有限参数对照：输入态与段内旋转分别做了什么

全协议改变一个参数，同时改变进入第二段的状态和第二段的生成元。为了区分这两种作用，可以固定一个样本的实际输入态，再单独更换第二段参数。以下均取第二周期、\(\tau=100,t_c=0.3,\hbar=1\)；能量单位 meV。

### 7.1 改变 \(E\)：这组翻转主要随输入相干改变

固定 \(\eta=0.023\)。表内是第二段终点 \(d\)：

| 输入态来自哪套完整协议 | 段内 \(E=0.01\) | 段内 \(E=0.02\) |
|---|---:|---:|
| 全程 \(E=0.01\) 的 \(\mathcal B_{\rm in}\) | −0.886132 | −0.766303 |
| 全程 \(E=0.02\) 的 \(\mathcal B_{\rm in}\) | +0.403502 | +0.368447 |

对角线对应原完整协议。非对角线是固定输入的段内对照，不是另一组全协议保真度样本。

这组参数中，更换段内 \(E\) 没有改变定向符号；更换输入态则在两个段内 \(E\) 下都改变了符号。此前阶段积累的相干对这次翻转至关重要。两套实际输入的 \(d\) 都很小，分别为 \(-0.0000544,-0.002625\)，但不能因此认为输入态几乎相同。

还有一个明确反例：原 \(E=0.02\) 样本段初 \(\dot d\approx-0.0004913\)，起初向负定向运动，段末却为 \(d\approx+0.368447\)。只看段初驱动会判断错终点。

### 7.2 改变 \(\eta\)：输入与段内作用都能改变定向

固定 \(E=0.01\)：

| 输入态来自哪套完整协议 | 段内 \(\eta=0.01\) | 段内 \(\eta=0.023\) |
|---|---:|---:|
| 全程 \(\eta=0.01\) 的 \(\mathcal B_{\rm in}\) | +0.805180 | −0.304304 |
| 全程 \(\eta=0.023\) 的 \(\mathcal B_{\rm in}\) | −0.207929 | −0.886132 |

在这组对照中，改变输入或改变段内 \(\eta\) 都可以使正定向变为负定向；两者共同改变时负定向更接近目标。具体平方保真度为对角线 \(0.009487,0.889374\)，非对角线 \(0.425298,0.364772\)。这说明 \(d<0\) 仍不等于高保真度，还需要平面位置正确。

两张表是对具体样本的有限参数诊断，不证明一般参数下 \(E\) 只改变输入，也不证明 \(\eta\) 增大总会提高保真度。

## 8. 下一步如何继续做参数分析

第二段的核可显式分解为

\[
L(t)=L_0(t)+E L_E+\eta g(t)L_\eta,
\]

\[
L_E=\begin{pmatrix}0&0&0&0\\0&1&0&0\\0&0&1&0\\1&0&0&0\\0&0&0&0\\0&0&0&0\end{pmatrix},
\quad
L_\eta=\begin{pmatrix}0&0&-1&0\\0&0&0&-1\\0&0&0&0\\0&0&0&0\\0&0&0&0\\1&0&0&0\end{pmatrix}.
\]

因此核的全部参数依赖是已知的 \(E^2,\eta^2,E\eta\) 等项。**这不是保真度的二阶展开**：解方程时积分算子反复作用，解包含参数的所有阶，并可发生完整的振荡和返回。

例如对任意有限参数点，令 \(\lambda=E\) 或 \(\eta\)，精确参数导数满足

\[
\dot Y_\lambda=\kappa(L_\lambda X_0+LX_{0,\lambda})
-\int_0^t(\mathscr K_\lambda Y+\mathscr K Y_\lambda)ds,
\qquad Y_\lambda(0)=Y_{0,\lambda},
\]

\[
\mathscr K_\lambda(t,s)=\kappa^2[
L_\lambda(t)L(s)^T+L(t)L_\lambda(s)^T].
\]

固定输入的分析取 \(X_{0,\lambda}=Y_{0,\lambda}=0\)；分析全协议必须保留它们。这不是在 \(E=\eta=0\) 附近截断，而是在任意参数点描述完整解的敏感性。

当前最具体的推进方向是：把前面阶段得到的输入相干，接入这一已知核，研究净驱动何时换号；同时用完整 \(X,Y\) 恢复保真度。瞬时零轴和频率提供几何结构，记忆方程提供实际状态输运。尚未得到整个参数空间中成功区域的一般解析边界。

为进一步降低分析难度，已选定并初步验证 [移动标架与慢球面动力学](geometric_moving_frame.md)：区分累计相位和几何混合，保留慢部分的完整有限参数动力学。[固定输入的换号位置预测](geometric_boundary_validation.md) 已与完整波函数独立核对，并从快慢返回反馈得到有效修正。该简化的适用性需要单独验证，完整方程仍作为基准。

## 9. 独立核对

代码：[verify_second_stage.py](geometric_validation/verify_second_stage.py)；数据：[second_stage.json](geometric_validation/second_stage.json)。

符号计算验证了分块生成元、零轴、频率特征多项式、四加六相干方程、范数守恒与记忆核第一行。将新分块方程独立传播，与已有五组完整复波函数轨迹在两个周期的第二段比较，最大状态、保真度及约束误差低于 \(4.7\times10^{-12}\)。

用独立稠密波函数解，在 \(\tau=100,500\) 的有限参数样本上核对六分量记忆方程及定向记忆行，归一化时间导数误差低于 \(9.6\times10^{-11}\)。五组样本的终点积分账本误差低于 \(1.4\times10^{-11}\)。输入交换表也直接使用完整波函数传播。

复现命令（使用已保存的 `geometric_validation/bottleneck.npz`）：

```sh
.venv/bin/python -B previous_work/geometric_validation/verify_second_stage.py
```
