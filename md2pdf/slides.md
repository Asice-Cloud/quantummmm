---
theme: default
title: E1,t1 与 SO(5) 编织误差
info: |
  $E_1,t_1$ 耦合下 Majorana 编织的群空间保真度分析与数值验证
class: text-center
transition: slide-left
mdc: true
highlighter: shiki
lineNumbers: false
drawings:
  persist: false
---

# $E_1,\;t_1$ 与 $SO(5)$ 编织误差

### 群空间保真度与参数独立性分析

<div class="pt-12 opacity-60 text-sm">阶段性汇报</div>

---
layout: default
---

# 研究思路

将 $E_1,\,t_1$ 的影响视为两部分：

- **误差来源**：实际演化与理想演化流形之间的误差
- **累积方式**：误差在非交换演化中如何累积，是否引起峰漂移等现象

核心目标：尝试**直接从哈密顿量的结构预测演化结果**，尽可能避开数值积分，
或找到可以进行展开估计的量。

---
layout: default
---

# 研究目标

研究 Majorana 模式在耦合 $E_1,\,t_1$ 存在时是否仍实现 $\gamma_2$-$\gamma_3$ 编织：

1. 用五个 Majorana 双线性生成元构造 $SO(5)$ 演化；
2. 在群空间中定义**理想编织**和**无泄漏编织**；
3. 将误差拆成**子空间泄漏**与**子空间内部旋转**；
4. 用相对演化、非交换展开和流形观点解释 $E_1,t_1$ 的影响及峰漂移。

---
layout: default
---

# 物理模型

系统模式为

$$\boldsymbol\gamma=(\gamma_1,\gamma_2,\gamma_3,\gamma_a,\gamma_b)^T$$

有效哈密顿量

$$
H(t)=\begin{aligned}
&iE_d\gamma_a\gamma_b+iE_1\gamma_1\gamma_2\\[2pt]
&+i|t_2|\gamma_a\gamma_2-i|t_1|\gamma_b\gamma_1-i|t_3|\gamma_a\gamma_3
\end{aligned}
$$

其中 $E_d$ 为量子点能级，$E_1$ 为内部耦合，$t_{1,2,3}$ 为隧穿幅度。

---
layout: default
---

# $SO(5)$ 演化

定义双线性生成元

$$X_{ij}=i\gamma_i\gamma_j$$

所有双线性生成元闭合为 $\mathfrak{so}(5)$，因此

$$\dot R(t)=\Omega(t)R(t),\qquad R(t)\in SO(5)$$

当前归一化必须使用

$$\Omega_{ij}=2h_{ij}$$

三段门控函数 $f_\pm(s)=\dfrac{1\pm\cos(\pi s/\tau)}{2}$，论文图对应两个三段周期，
总时间 $T=6\tau$。

---
layout: default
---

# 理想编织

一次交换的算符变换为

$$
\gamma_1\mapsto\gamma_1,\qquad
\gamma_2\mapsto\gamma_3,\qquad
\gamma_3\mapsto-\gamma_2
$$

对应 $SO(3)$ 矩阵

$$
B_{23}=\begin{pmatrix}1&0&0\\0&0&-1\\0&1&0\end{pmatrix},
\qquad B_{23}^2=\operatorname{diag}(1,-1,-1)
$$

由于量子点内部最终旋转无关紧要，理想结果不是单个 $SO(5)$ 矩阵，而是目标集合

$$\mathcal C_B=\{\operatorname{diag}(B,A_Q):A_Q\in SO(2)\}$$

---
layout: default
class: text-sm
---

# 无泄漏编织

取投影算符

$$
P=\operatorname{diag}(1,1,1,0,0),\qquad Q=I_5-P,\qquad
R=\begin{pmatrix}M&L\\N&A_Q\end{pmatrix}
$$

**无泄漏的充要条件**

$$\boxed{QRP=0}\quad(\text{即 }N=0)$$

**正确编织还要求**

$$\boxed{PRP=B}$$

因此必须区分「子空间闭合」与「闭合后的内部旋转正确」。
三维子空间属于

$$\mathrm{Gr}(3,5)\simeq SO(5)/(SO(3)\times SO(2))$$

---
layout: default
class: text-sm
---

# 关于约束曲线机制的备注

> ~~让 AI 整理了各种条件，稳定子群、对易性、时间回文、群回波、参数交换、
> 暗态平行输运和新增生成元等可能产生约束曲线的机制。这些是理论设计条件，
> 当前并未证明固定协议存在全部精确曲线。感觉没什么可行性~~

<br>

<div class="opacity-70">

结论：该方向暂不作为主线，转而用**群空间保真度**定量刻画误差。

</div>

---
layout: default
---

# 群空间保真度 · 子空间存活率

由 $R^TR=I$，有 $M^TM+N^TN=I_3$。定义**子空间存活率**

$$\boxed{F_{\rm surv}=\frac13\operatorname{Tr}(M^TM)=1-\frac13\|N\|_F^2}$$

它只衡量泄漏，不区分旋转方向是否正确。

---
layout: default
class: text-sm
---

# 群空间保真度 · 内部旋转保真度

若存在泄漏，$M$ 不再是纯旋转，而是「旋转 + 收缩/形变」，一般有 $M^TM<I_3$。

对 $M$ 做极分解

$$M=HU,\qquad H=(MM^T)^{1/2},\qquad U\in SO(3)$$

定义

$$\boxed{F_{\rm rot}^{\rm cond}=\frac{1+\operatorname{Tr}(B^TU)}4}$$

若相对旋转角为 $\vartheta$，则 $F_{\rm rot}^{\rm cond}=\cos^2(\vartheta/2)$。

例如 $M=aB$ 且 $0<a<1$ 时 $U=B$，旋转方向完全正确；但 $a<1$ 表示存在泄漏，
$F_{\rm surv}<1$。**所以必须同时报告两个指标。**

---
layout: default
---

# 联合指标

$$\boxed{F_{\rm group}=F_{\rm surv}\cdot F_{\rm rot}^{\rm cond}}$$

该指标同时惩罚**泄漏**和**内部编织错误**，直接作用于 Majorana 群变换，
不依赖某一个波函数。

<br>

<div class="text-sm opacity-80">

代码 `so5_group_fidelity_map.py` 实现了 $F_{\rm group}$ 的二维扫描，
并生成 `.png` / `.npz` 数据。参数：$E_1=0.01$ meV、归一化因子 2、两个三段周期。

</div>

---
layout: default
---

# 保真度扫描结果

<img src="/so5_group_fidelity_map.png" class="block h-80 mx-auto rounded shadow" />

<div class="text-xs opacity-60 text-center pt-2">
$F_{\rm group}$ 二维扫描（$E_1=0.01$ meV，归一化因子 2，两个三段周期）
</div>

---
layout: default
class: text-sm
---

# $(E_1,t_1)$ 独立性检验

为直接检验能否把两个参数分离，差分检验脚本计算四点差分

$$\Delta_{\oplus}F=F(E,t)-F(E,0)-F(0,t)+F(0,0)$$

以及归一化乘积残差

$$\Delta_{\otimes}F=F(E,t)-\frac{F(E,0)\,F(0,t)}{F(0,0)}$$

- 若参数效应**可加**，则 $\Delta_{\oplus}F=0$
- 若可按基准值**归一化后相乘**，则 $\Delta_{\otimes}F=0$

---
layout: default
class: text-sm
---

# 数值结果

在 $\tau=1\,\mathrm{meV}^{-1}$、$T=6\tau$、$\Omega_{ij}=2h_{ij}$ 下，取

$$E=2\times10^{-4}\,\mathrm{meV},\qquad t=3\times10^{-4}\,\mathrm{meV}$$

$$
\begin{array}{c|c|c}
 &\Delta_{\oplus}F&\Delta_{\otimes}F\\ \hline
F_{\rm pq}&3.1878933\times10^{-9}&3.1877423\times10^{-9}\\
F_{\rm group}&1.6263352\times10^{-8}&1.6263336\times10^{-8}
\end{array}
$$

两种残差均**明显非零**，因此当前模型中 $E_1,t_1$ 既不能用简单加和表示，
也不能用归一化乘积表示。

---
layout: default
class: text-sm
---

# 缩放分析与结论

进一步将 $(E,t)$ 同时缩放为 $(sE,st)$，得到

$$
\begin{array}{c|c|c}
s&\Delta_{\oplus}F_{\rm pq}/s^2&\Delta_{\oplus}F_{\rm group}/s^2\\ \hline
1&3.1878933\times10^{-9}&1.6263352\times10^{-8}\\
0.5&3.1673197\times10^{-9}&1.6264857\times10^{-8}\\
0.25&3.0846632\times10^{-9}&1.6258217\times10^{-8}
\end{array}
$$

除去数值积分误差后 $\Delta F\propto s^2$，说明该残差由**二阶 $E_1t_1$ 混合项**主导。

<div class="pt-4 opacity-80 text-xs">

结论：$E_1,t_1$ 的影响不可分离；用 $F_{\rm group}$ 作为统一量化指标，
可同时捕捉泄漏与编织错误。

</div>
