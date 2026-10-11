# 相干通道怎样决定定向返回：完整方程与参数分析入口

本文承接 [geometric_state_analysis.md](geometric_state_analysis.md) 的有向平面状态与几何瓶颈。当前问题是：\(E_1,|t_1|\) 怎样改变相干，使经过瓶颈的轨迹最终走向 \(d\approx-1\) 或 \(d\approx+1\)？

先得到三项可以直接从哈密顿量推导的结构：完整相干方程；第一段两个精确的球面旋转通道及其标量记忆方程；第三段保持保真度、但传递相干相位的精确映射。第二段把这些通道与模式 3 耦合，是继续分析定向返回的重点。

这些结果不需要先求出完整 \(R\)，不要求 \(E_1\tau,|t_1|\tau\) 小，也不使用绝热投影。它们仍是动力学方程；尚未给出整个协议的终点保真度一般闭式解。

第二段的进一步推导见 [geometric_second_stage.md](geometric_second_stage.md)：瞬时零轴与旋转频率、核全部已知的六分量记忆方程，以及区分输入相干和段内参数作用的有限参数对照。

## 1. 直接推动定向的是哪几个量

记

\[
E=E_1,\quad q(t)=|t_1(t)|=\eta g(t),\quad
p(t)=|t_2(t)|,\quad r(t)=|t_3(t)|,\quad \kappa=2/\hbar.
\]

定义三个复相干分量和三个实分量：

\[
Z_j=\mathcal B_{1j}+i\mathcal B_{2j}\quad(j=3,a,b),
\qquad h_{3a}=\mathcal B_{3a},\quad
h_{3b}=\mathcal B_{3b},\quad h_{ab}=\mathcal B_{ab}.
\]

连同 \(d=\mathcal B_{12}\)，它们恰好包含十个实矩阵分量。纯态轨道约束仍使独立状态自由度为六。

定向方程为

\[
\boxed{\dot d=-\kappa J(t),\qquad
J(t)=p(t)\operatorname{Re}Z_a(t)+q(t)\operatorname{Im}Z_b(t).}
\]

\(J>0\) 时 \(d\) 减小，\(J<0\) 时 \(d\) 增大。任取一个瓶颈穿越时刻 \(d(t_\times)=0\)，都有

\[
\boxed{d(T)=-\kappa\int_{t_\times}^{T}J(s)\,ds.}
\]

所以终点定向由穿越之后的净相干驱动决定。驱动可以反复变号，单看穿越瞬间的方向并不足以决定终点。高保真度还需要满足平面位置条件。

## 2. 完整的相干动力学

从 \(\dot{\mathcal B}=[K,\mathcal B]\) 展开得到

\[
\begin{aligned}
\dot Z_3&=\kappa[-iEZ_3+rZ_a+ip h_{3a}-q h_{3b}],\\
\dot Z_a&=\kappa[-iEZ_a-rZ_3+E_dZ_b+pd-qh_{ab}],\\
\dot Z_b&=\kappa[-iEZ_b-E_dZ_a+iqd-ip h_{ab}],\\
\dot h_{3a}&=\kappa[E_dh_{3b}-p\operatorname{Im}Z_3],\\
\dot h_{3b}&=\kappa[rh_{ab}-E_dh_{3a}+q\operatorname{Re}Z_3],\\
\dot h_{ab}&=\kappa[-rh_{3b}+q\operatorname{Re}Z_a+p\operatorname{Im}Z_b].
\end{aligned}
\]

加上第 1 节的 \(\dot d\)，这是一套完整闭合的十实分量线性方程。初始条件为 \(d=1\)，其余分量为零。

这里已经明确了参数作用：\(E\) 以 \(-i\kappa E Z_j\) 旋转相干相位；\(q\) 参与相干的生成与反馈；\(r\) 把与模式 3 有关的分量接入动力学。

可以直接从这些变量恢复保真度：

\[
\boxed{
\ell=1-d^2-\frac12(|Z_3|^2+|Z_a|^2+|Z_b|^2),\qquad
F=\frac{d^2-d}{2}+\frac{|Z_3|^2+|Z_a|^2+|Z_b|^2}{4}.
}
\]

于是 \(h\) 分量虽然不直接出现在保真度中，仍会改变下一时刻的相干及其驱动。

## 3. 第一段：两个精确的球面旋转通道

定义

\[
X_\pm=d\pm h_{ab},\qquad C_\pm=Z_a\pm iZ_b,
\qquad V_\pm=p\mp q,\qquad \Delta_\pm=E\pm E_d.
\]

在完整协议中也有精确关系

\[
\boxed{
\begin{aligned}
\dot X_\pm&=-\kappa V_\pm\operatorname{Re}C_\pm
\mp\kappa r h_{3b},\\
\dot C_\pm&=-i\kappa\Delta_\pm C_\pm
+\kappa V_\pm X_\pm-\kappa r Z_3.
\end{aligned}}
\]

第一段 \(r=0\)，两套方程完全解耦。令 \(C_\pm=Y_\pm+iW_\pm\)，得到

\[
\frac d{dt}\begin{pmatrix}X_\pm\\Y_\pm\\W_\pm\end{pmatrix}
=\kappa\begin{pmatrix}
0&-V_\pm&0\\
V_\pm&0&\Delta_\pm\\
0&-\Delta_\pm&0
\end{pmatrix}
\begin{pmatrix}X_\pm\\Y_\pm\\W_\pm\end{pmatrix}.
\]

这是围绕轴 \(\kappa(-\Delta_\pm,0,V_\pm)\) 的三维旋转，并保持

\[
X_\pm^2+|C_\pm|^2=\text{常数}.
\]

第一周期第一段的初值是 \((X_\pm,C_\pm)=(1,0)\)，两条轨迹均在单位球面上。这个形式与 Bloch 向量旋转相同，但不把它们解释为两个独立物理量子比特。

论文门控满足 \(p=t_cg,q=\eta g\)。因此两个旋转轴由

\[
\boxed{\Delta_\pm=E\pm E_d(t),\qquad
V_\pm=(t_c\mp\eta)g(t)}
\]

直接确定。\(E\) 改变相位旋转；\(\eta\) 对两条通道的驱动以相反符号出现。是否使终点保真度提高，仍取决于累计相位及后续混合，不能从单个系数推出单调性。

回到原变量，\(d=(X_++X_-)/2\)、\(h_{ab}=(X_+-X_-)/2\)。这使两条通道的联合旋转与定向变化直接联系起来。

## 4. 第一段的标量记忆方程：系数全部已知

定义完整相位

\[
\Theta_\pm(t,s)=\kappa\int_s^t\Delta_\pm(u)\,du.
\]

在任意 \(r=0\) 的区间 \([t_0,t]\)，复相干可以精确消去：

\[
C_\pm(t)=e^{-i\Theta_\pm(t,t_0)}C_\pm(t_0)
+\kappa\int_{t_0}^t e^{-i\Theta_\pm(t,s)}
V_\pm(s)X_\pm(s)\,ds.
\]

因此

\[
\boxed{
\begin{aligned}
\dot X_\pm(t)={}&-\kappa V_\pm(t)
\operatorname{Re}[e^{-i\Theta_\pm(t,t_0)}C_\pm(t_0)]\\
&-\int_{t_0}^t\mathscr K_\pm(t,s)X_\pm(s)\,ds,\\
\mathscr K_\pm(t,s)={}&\kappa^2V_\pm(t)V_\pm(s)
\cos\Theta_\pm(t,s).
\end{aligned}}
\]

第一周期第一段 \(C_\pm(0)=0\)，第一行的初始相干项消失。对于余弦门控 \(g=f_-\)、\(E_d=E_0f_+\)，在 \(0\leq s\leq t\leq\tau\) 上

\[
\Theta_\pm(t,s)=\kappa\left[
E(t-s)\pm\frac{E_0}{2}\left(t-s+
\frac\tau\pi[\sin(\pi t/\tau)-\sin(\pi s/\tau)]\right)\right],
\]

\[
\boxed{\mathscr K_\pm(t,s)=\kappa^2(t_c\mp\eta)^2
g(t)g(s)\cos\Theta_\pm(t,s).}
\]

这是两个标量 Volterra 方程，其核不含未知 \(R\) 或待求传播矩阵。\(E\) 的全有限参数依赖表现为余弦内的相位，\(\eta\) 的依赖表现为两条通道不同的权重。没有截断累计相位。

这个表示允许研究已知积分算子的相位抵消、参数导数和受控近似。它仍需要解积分方程；不能把余弦核的符号直接当作终点保真度的符号，也不能将第一段结果无条件推广到整个协议。

第二周期第一段必须使用上一周期传来的 \(X_\pm(t_0),C_\pm(t_0)\)。将相干初值重新置零会丢失第三段传递的相位。

## 5. 一个不依赖数值传播的参数判据

第一周期第一段若 \(\eta=t_c\)，则 \(V_+=0\)。由初值 \((X_+,C_+)=(1,0)\) 得到

\[
X_+(t)=1,\qquad C_+(t)=0
\]

对所有 \(E\)、所有 \(\tau\) 精确成立。另一通道在单位球上，\(-1\leq X_-\leq1\)，因此

\[
\boxed{0\leq d(t)=\frac{1+X_-(t)}2\leq1,
\qquad F(t)\leq\frac14\quad\text{（第一周期第一段）}.}
\]

这是直接从参数得到的定向限制，不依赖 \(E\tau\) 小量。它只约束这一段；第二段的 \(r\) 耦合会破坏 \(X_+\) 的这个不变量，不能据此判断整个双周期一定失败。

## 6. 第二段混合通道，第三段传递相位

第二段有 \(E_d=0\)、\(r\neq0\)，两个通道的相位项均为 \(E\)，而 \(rZ_3,rh_{3b}\) 将它们与模式 3 耦合。需要保留第 2 节中 \(Z_3,h_{3a},h_{3b}\) 的反馈；第一段的两个独立球面方程不再单独闭合。

第三段 \(p=q=0\)，因此 \(d\) 恒定。令复行向量 \(Z=(Z_3,Z_a,Z_b)\)，并定义只作用于 \((3,a,b)\) 的矩阵

\[
D(t)=\begin{pmatrix}0&r&0\\-r&0&E_d\\0&-E_d&0\end{pmatrix},
\qquad
\dot Q=\kappa DQ,\quad Q(t_0)=I.
\]

此时有精确映射

\[
\boxed{Z(t)=e^{-i\kappa E(t-t_0)}Z(t_0)Q(t)^T.}
\]

\(Q\) 仅依赖既定控制与协议时长，和 \(E,\eta\) 无关。对于实际控制它本身也可能需要积分，但可在固定协议下单独求一次。相干范数保持，故 \(\ell,F\) 也恒定；相干相位继续变化。

一段长度 \(\tau\) 给出相位因子 \(e^{-2iE\tau/\hbar}\)。这会改变下一周期输入相干的实部和虚部，进而改变 \(J\)。改变全协议 \(E\) 同时改变进入第三段的态，不能把全部参数影响归结为这一因子。

## 7. 独立核对与下一步的具体对象

核对代码为 [verify_coherence_dynamics.py](geometric_validation/verify_coherence_dynamics.py)，记录为 [coherence_dynamics.json](geometric_validation/coherence_dynamics.json)。符号计算逐项验证完整相干方程、通道变换和球面范数；新相干方程的完整传播与原五组复波函数数据比较，十个实分量及保真度一致，所有报告误差均低于 \(1.5\times10^{-11}\)。

用独立的稠密波函数解，核对 \(\tau=100,500\) 两组有限参数的第一段标量记忆方程，局部归一化时间导数误差低于 \(9\times10^{-11}\)。另取 \(E=0,0.01,0.08\)、\(\eta=t_c=0.3\)、\(\tau=100\)，核对 \(X_+=1\) 不变量，误差低于 \(3.2\times10^{-12}\)。

下表显示原有限参数样本在第二周期第一、第二段结束时的定向。能量单位 meV，\(\tau\) 单位 \(\mathrm{meV}^{-1}\)。

| \(E\) | \(\eta\) | \(\tau\) | \(d(4\tau)\) | \(d(5\tau)=d(6\tau)\) |
|---|---|---|---|---|
| 0.01 | 0.01 | 100 | +0.001191 | +0.805180 |
| 0.01 | 0.023 | 100 | −0.0000544 | −0.886132 |
| 0.02 | 0.023 | 100 | −0.002625 | +0.368447 |
| 0.01 | 0.02344645627 | 500 | −0.002564 | +0.997814 |

这些样本在 \(4\tau\) 附近都接近瓶颈投影，但随后第二段走向不同的定向；第三段保持该结果。这个样本观察不证明一般参数下 \(d(4\tau)\) 总接近零。

因此下一步集中分析 **第二周期第二段的净驱动 \(\int_{4\tau}^{5\tau}J(t)dt\)**：它怎样依赖进入这一段的 \(C_\pm,Z_3,h_{3a},h_{3b}\)，以及已知的 \(E,(t_c\mp\eta)g,r\)。第一段的积分核与第三段的相位映射用于组织输入态，第二段的完整耦合用于判断正、负定向返回。

当前已经明确了分析变量、参数进入的方式、分段结构和一个参数依赖的不变量。一般参数下这段净驱动的符号及终点高保真度区域，仍是待解决的动力学问题。

复现命令（需要上一轮生成的 `geometric_validation/bottleneck.npz`）：

```sh
.venv/bin/python -B previous_work/geometric_validation/verify_coherence_dynamics.py
```
