### $E_1, t_1$ and $SO(5)$

思路：

将 $E_1,t_1$ 的影响视为两部分：一是实际演化与理想演化流形之间的误差；二是误差在非交换演化中的累积方式，以及它是否引起峰漂移等现象。主要就是尝试能否直接从哈密顿量的结构预测演化结果，尽可能避开数值积分，或者找到可以进行展开估计的量。



#### 1.

研究 Majorana 模式在耦合 $E_1,t_1$ 存在时是否仍实现 $\gamma_2$-$\gamma_3$ 编织。思路依次为：

1. 用五个 Majorana 双线性生成元构造 $SO(5)$ 演化；
2. 在群空间中定义理想编织和无泄漏编织；
3. 将误差拆成子空间泄漏与子空间内部旋转；
4. 用相对演化、非交换展开和流形观点解释 $E_1,t_1$ 的影响

---

#### 2. 物理模型与 $SO(5)$ 演化

系统模式为

\[
\boldsymbol\gamma=(\gamma_1,\gamma_2,\gamma_3,\gamma_a,\gamma_b)^T.
\]

有效哈密顿量为

\[
H(t)=iE_d\gamma_a\gamma_b+iE_1\gamma_1\gamma_2
 +i|t_2|\gamma_a\gamma_2-i|t_1|\gamma_b\gamma_1
 -i|t_3|\gamma_a\gamma_3.
\]

定义 $X_{ij}=i\gamma_i\gamma_j$，所有双线性生成元闭合为 $\mathfrak{so}(5)$。因此

\[
\dot R(t)=\Omega(t)R(t),
\qquad R(t)\in SO(5).
\]

当前归一化必须使用

\[
\Omega_{ij}=2h_{ij}.
\]

三段门控函数为 $f_\pm(s)=[1\pm\cos(\pi s/\tau)]/2$，论文图对应两个三段周期，总时间 $T=6\tau$。

---

#### 3. 理想编织和无泄漏条件

##### 3.1 理想编织

一次交换的算符变换为

\[
\gamma_1\mapsto\gamma_1,\qquad
\gamma_2\mapsto\gamma_3,\qquad
\gamma_3\mapsto-\gamma_2.
\]

对应 $SO(3)$ 矩阵

\[
B_{23}=\begin{pmatrix}1&0&0\\0&0&-1\\0&1&0\end{pmatrix},
\qquad B_{23}^2=\operatorname{diag}(1,-1,-1).
\]

由于量子点内部最终旋转无关紧要，理想结果不是单个 $SO(5)$ 矩阵，而是目标集合

\[
\mathcal C_B=\{\operatorname{diag}(B,A_Q):A_Q\in SO(2)\}.
\]

##### 3.2 无泄漏编织

取

\[
P=\operatorname{diag}(1,1,1,0,0),\qquad Q=I_5-P,
\]

并分块 $R=\begin{pmatrix}M&L\\N&A_Q\end{pmatrix}$。

无泄漏的充要条件为

\[
\boxed{QRP=0}\quad(\text{即 }N=0).
\]

正确编织还要求

\[
\boxed{PRP=B}.
\]

因此必须区分“子空间闭合”和“闭合后的内部旋转正确”。三维子空间属于

\[
\mathrm{Gr}(3,5)\simeq SO(5)/(SO(3)\times SO(2)).
\]

~~让AI整理了各种条件，稳定子群、对易性、时间回文、群回波、参数交换、暗态平行输运和新增生成元等可能产生约束曲线的机制。这些是理论设计条件，当前并未证明固定协议存在全部精确曲线。感觉没什么可行性~~

---

#### 4. 群空间保真度

##### 4.1 子空间存活率

由 $R^TR=I$，有 $M^TM+N^TN=I_3$。定义

\[
\boxed{F_{\rm surv}=\frac13\operatorname{Tr}(M^TM)=1-\frac13\|N\|_F^2}.
\]

#### 4.2 内部旋转保真度

若存在泄漏，$M$ 不再严格是 $SO(3)$，不再是纯旋转，而是“旋转 + 收缩/形变”。这是因为部分 MZM 信息已经流入量子点子空间，此时一般有

\[
M^TM<I_3.
\]

对 $M$ 做极分解

\[
M=HU,
\qquad H=(MM^T)^{1/2},
\qquad U\in SO(3).
\]

定义

\[
\boxed{F_{\rm rot}^{\rm cond}=\frac{1+\operatorname{Tr}(B^TU)}4}.
\]

若相对旋转角为 $\vartheta$，则 $F_{\rm rot}^{\rm cond}=\cos^2(\vartheta/2)$。

例如若 $M=aB$ 且 $0<a<1$，则 $U=B$，说明旋转方向完全正确；但 $a<1$ 表示存在泄漏，因而 $F_{\rm surv}<1$，最终 $F_{\rm group}<1$。所以必须同时报告 $F_{\rm surv}$ 和 $F_{\rm rot}^{\rm cond}$，不能只看其中一个指标。

##### 4.3 联合指标

\[
\boxed{F_{\rm group}=F_{\rm surv}F_{\rm rot}^{\rm cond}}.
\]

这个指标同时惩罚泄漏和内部编织错误，直接作用于 Majorana 群变换，不依赖某一个波函数。

代码 [so5_group_fidelity_map.py](/home/asice-cloud/projects/pyyy/quantumsss/ss3/so5_group_fidelity_map.py) 实现了 $F_{\rm group}$ 的二维扫描，并生成 `.png`/`.npz` 数据。它使用 $E_1=0.01$ meV、归一化因子 2 和两个三段周期。



![so5_group_fidelity_map](/home/asice-cloud/projects/pyyy/quantumsss/ss3/so5_group_fidelity_map.png)



---

#### 6. 综合结论与验证边界

##### 6.1 (E_1,t_1) 独立性的直接检验

为直接检验能否把两个参数分离，脚本 [verify_nonseparable.py](/home/asice-cloud/projects/pyyy/quantumsss/ss3/verify_nonseparable.py)，计算四点差分

\[
\Delta_{\oplus}F
=F(E,t)-F(E,0)-F(0,t)+F(0,0),
\]

以及归一化乘积残差

\[
\Delta_{\otimes}F
=F(E,t)-\frac{F(E,0)F(0,t)}{F(0,0)}.
\]

若参数效应可加，则 $Delta_{\oplus}F=0$；若可按基准值归一化后相乘，则 $Delta_{\otimes}F=0$。



在 $\tau=1\,\mathrm{meV}^{-1}$、$T=6\tau$、$Omega_{ij}=2h_{ij}$下，取

\[
E=2\times10^{-4}\,\mathrm{meV},
\qquad t=3\times10^{-4}\,\mathrm{meV},
\]

得到：

\[
\begin{array}{c|c|c}
 &\Delta_{\oplus}F&\Delta_{\otimes}F\\ \hline
F_{\rm pq}&3.1878933\times10^{-9}&3.1877423\times10^{-9}\\
F_{\rm group}&1.6263352\times10^{-8}&1.6263336\times10^{-8}
\end{array}
\]

两种残差均明显非零，因此当前模型中 $E_1,t_1$ 既不能用简单加和表示，也不能用归一化乘积表示。



进一步将 $(E,t)$ 同时缩放为 $(sE,st)$，得到：

\[
\begin{array}{c|c|c}
s&\Delta_{\oplus}F_{\rm pq}/s^2&\Delta_{\oplus}F_{\rm group}/s^2\\ \hline
1&3.1878933\times10^{-9}&1.6263352\times10^{-8}\\
0.5&3.1673197\times10^{-9}&1.6264857\times10^{-8}\\
0.25&3.0846632\times10^{-9}&1.6258217\times10^{-8}
\end{array}
\]

除去数值积分误差后，$\Delta F\propto s^2$，说明该残差由二阶 $E_1t_1$ 混合项主导。
