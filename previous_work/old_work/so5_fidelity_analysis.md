# 用 SO(5) 研究论文保真度及 E₁、t₁ 的影响

本文件保留先前的局部扰动分析。完整有限角度、保留五模和非对易作用的分析见 [so5_full_fidelity_analysis.md](so5_full_fidelity_analysis.md)。本文件的二阶抵消线不能推广为全时间的高保真度线。

本文直接分析论文图中的波函数重叠。主要结果是：

\[
\boxed{
1-F\simeq\frac{2\tau^2}{\hbar^2}
\left[2.180340599\,E_1-0.929923300\,|t_1|\right]^2.
}
\]

这里 \(F=f_{\rm paper}^2\)，\(f_{\rm paper}\) 是论文画出的重叠模长。该显式公式采用论文的 t₁ 门控规律、两次三段交换，并要求理想控制等幅、绝热且累积扰动小。下面同时给出不要求绝热的定义与参数展开，以便区分精确结论和近似结论。

## 1. 对准论文的观测量

论文第二节（205411-3 页）取

\[
c_-=(e_1-ie_2)/\sqrt2,\qquad
c_+=(e_1+ie_2)/\sqrt2,
\]

对应 \(\psi_1^-(0)=(u_1-iu_2)/\sqrt2\) 与目标 \(\psi_1^+(0)=(u_1+iu_2)/\sqrt2\)。双交换总时间 \(T=6\tau\)。

在 \((1,2,3,a,b)\) 的初始 Majorana 波函数基中，以 \(R\in SO(5)\) 表示有效模型的传播，则

\[
f_{\rm paper}=|c_+^\dagger R(T)c_-|,
\qquad
\boxed{F=\frac{(R_{11}-R_{22})^2+(R_{12}+R_{21})^2}{4}.}
\]

这在有效模型内是精确恒等式，不要求绝热、弱扰动或二能级投影。未耦合的 \(\gamma_4\) 是旁观者，不影响这里的 \(\psi_1\) 重叠。完整纳米线模型若有向其他能级的泄漏，则应扩大传播空间。

这是 Majorana/BdG 波函数的指标；固定宇称多体量子态的传播需使用 Spin(5) 旋量表示，不能把不同表示的矩阵直接互换。

## 2. 保真度的精确几何含义

### 2.1 两种失效：回到原分支与离开原平面

定义

\[
F_{\rm return}=|c_-^\dagger R c_-|^2
=\frac{(R_{11}+R_{22})^2+(R_{12}-R_{21})^2}{4},
\]

\[
\mathcal L=\frac12\sum_{j\in\{3,a,b\}}
\left(R_{j1}^2+R_{j2}^2\right).
\]

\(\mathcal L\) 是波函数离开初始 \((1,2)\) 平面的总权重。由于 \(c_+,c_-,e_3,e_a,e_b\) 构成复正交基，

\[
\boxed{F+F_{\rm return}+\mathcal L=1.}
\]

令 \(T_{12}\) 为 R 的左上 \(2\times2\) 块，定义有向平面的重叠

\[
\mathcal O=\det T_{12}
=\langle e_1\wedge e_2,\,R e_1\wedge R e_2\rangle.
\]

则

\[
\boxed{F=\tfrac12(1-\mathcal L-\mathcal O).}
\]

理想双交换有 \(\mathcal L=0\)、\(\mathcal O=-1\)：原来的二维平面保留下来，但其定向被翻转，因而 \(c_-\) 转换成 \(c_+\)，允许一个整体相位。

### 2.2 相关的商空间是有向二维平面的群轨道

选择参考旋转

\[
R_B=\operatorname{diag}(1,-1,-1,1,1),\qquad
W=R_B^T R,
\]

有 \(F=|c_-^\dagger W c_-|^2\)。使 F 保持为 1 的误差满足

\[
Wc_-=e^{i\varphi}c_-.
\]

其稳定子为 \(SO(2)\times SO(3)\)：在 \((1,2)\) 平面内旋转只增加相位，正交的三维补空间也可以任意旋转。相应的几何空间是

\[
\boxed{SO(5)/(SO(2)\times SO(3)),}
\]

即有向二维平面的 Grassmann 流形，维数为 \(10-1-3=6\)。因此论文的指标关注二维平面的定向和混合，而不等同于整个群元的距离，也不只取决于 \(\gamma_1\) 的方向。

## 3. 哪些李代数方向降低保真度？

若 \(W=\exp\Omega\) 接近单位元，\(\Omega^T=-\Omega\)，则

\[
1-F=\|(I-c_-c_-^\dagger)Wc_-\|^2
\simeq\|(I-c_-c_-^\dagger)\Omega c_-\|^2.
\]

定义 \(L_{ij}=e_i e_j^T-e_j e_i^T\)。分解

\[
\Omega=\alpha L_{12}
+\sum_{j\in\{3,a,b\}}(p_jL_{1j}+q_jL_{2j})
+\Omega_{\{3,a,b\}}.
\]

二阶的保真度损失为

\[
\boxed{1-F\simeq\frac12\sum_{j\in\{3,a,b\}}(p_j^2+q_j^2).}
\]

\(L_{12}\) 和补空间内部的旋转在这个指标中不产生直接损失；六个跨平面方向产生损失。这是保真度诱导的局部度量。

## 4. E₁、t₁ 的精确扰动参考系

采用 \(\{\gamma_i,\gamma_j\}=2\delta_{ij}\)、\(H=(i/2)\gamma^T A\gamma\)，则 \(\dot R=(2A/\hbar)R\)。令 \(E=E_1\)、\(\eta=|t_1|\) 为幅度，实际时变耦合为 \(\eta g(t)\)，有

\[
K=K_0+\frac2\hbar(E L_{12}+\eta g(t)L_{1b}).
\]

无扰动演化写成 \(R_0=\operatorname{diag}(1,Q)\)，\(Q\in SO(4)\)，四维子空间按 \((2,3,a,b)\) 排序。令 \(S=R_0^T R\)，则

\[
\dot S=M(u)S,\qquad
M(u)=\begin{pmatrix}0&u^T\\-u&0\end{pmatrix},
\quad
u=\frac2\hbar Q^T(Ee_2+\eta g e_b).
\]

这个误差方程对任意控制速度都是精确的。其一阶对数为 \(\Omega_1=M(z)\)，其中

\[
z=Ea+\eta b,\qquad
a=\frac2\hbar\int_0^{6\tau}Q^T e_2\,dt,
\quad b=\frac2\hbar\int_0^{6\tau}g Q^T e_b\,dt.
\]

当无扰动终点实现目标波函数时，\(R_0(T)\) 与 \(R_B\) 的区别属于态的稳定子。因此第 3 节给出

\[
\boxed{
1-F=\frac12\left(z_{13}^2+z_{1a}^2+z_{1b}^2\right)
+O(\epsilon^4),
}
\]

这里 z 的第一个分量是 \(z_{12}\)，它在二阶损失中不出现；\(\epsilon\) 表示小的累积扰动尺度。于是敏感性由三维响应向量

\[
a_\perp=(a_3,a_a,a_b),\qquad b_\perp=(b_3,b_a,b_b)
\]

决定：

\[
1-F\simeq\frac12\left[E^2\|a_\perp\|^2
+2E\eta\,a_\perp\cdot b_\perp+\eta^2\|b_\perp\|^2\right].
\]

如果 \(b_\perp\neq0\)，该二次式的无约束最小点是

\[
\eta_*=-E\frac{a_\perp\cdot b_\perp}{\|b_\perp\|^2}.
\]

实际还要求 \(\eta_*\ge0\)。只有两个响应向量反向共线时，一个参数才足以消除全部二阶损失。不同协议下的交叉项符号不能仅由裸耦合强度判断。

四阶余项来自本模型的对称性：以 \(P=\operatorname{diag}(-1,1,1,1,1)\) 共轭可将 \((E,\eta)\) 同时变号，而 F 不变。因此固定初态和控制下，F 的参数展开没有总次数为奇数的项。

## 5. 非绝热终点：保留真实的无扰动基线

若 \(F_0(\tau)=|c_+^\dagger R_0(T)c_-|^2<1\)，不能直接把上一节的正定二次式加到 \(1-F_0\) 上，因为扰动也可能改变原有的非绝热误差。

一个不要求绝热的计算工具是传播矩阵的参数展开。记

\[
R=R_0+E R_E+\eta R_\eta
+E^2 C_{EE}+E\eta C_{E\eta}+\eta^2 C_{\eta\eta}+\cdots,
\]

\(K_E=2L_{12}/\hbar\)、\(K_\eta(t)=2g(t)L_{1b}/\hbar\)，则

\[
\begin{aligned}
\dot R_E&=K_0R_E+K_E R_0,\\
\dot R_\eta&=K_0R_\eta+K_\eta R_0,\\
\dot C_{EE}&=K_0C_{EE}+K_E R_E,\\
\dot C_{E\eta}&=K_0C_{E\eta}+K_E R_\eta+K_\eta R_E,\\
\dot C_{\eta\eta}&=K_0C_{\eta\eta}+K_\eta R_\eta.
\end{aligned}
\]

除 \(R_0(0)=I\) 外，所有初值为零。对各矩阵取 \(c_+^\dagger(\cdot)c_-\)，得到振幅系数 \(a_0,a_E,a_\eta,a_{EE},a_{E\eta},a_{\eta\eta}\)。于是

\[
F=F_0+E^2 h_{EE}+2E\eta h_{E\eta}+\eta^2 h_{\eta\eta}+O(\epsilon^4),
\]

\[
\begin{aligned}
h_{EE}&=|a_E|^2+2\operatorname{Re}(a_0^*a_{EE}),\\
h_{E\eta}&=\operatorname{Re}(a_E^*a_\eta+a_0^*a_{E\eta}),\\
h_{\eta\eta}&=|a_\eta|^2+2\operatorname{Re}(a_0^*a_{\eta\eta}).
\end{aligned}
\]

这只需对每个 \(\tau\) 求一次响应系数，就得到两个参数在零扰动附近的完整二次依赖，不需要对 E、t₁ 做二维扫描。非绝热基线附近 h 不必负半定。

## 6. 代入论文的实际门控

取 \(t_c=E_0=\lambda_0\) 的等幅情况，局部时间 \(v\in[0,1]\)，\(f_\pm=(1\pm\cos\pi v)/2\)。论文中 t₁ 与 t₂ 共用 G₁，因而

| 段 | \(|t_2|\) | \(|t_3|\) | \(E_d\) | \(g\) |
|---|---|---|---|---|
| 1 | \(\lambda_0 f_-\) | 0 | \(\lambda_0 f_+\) | \(f_-\) |
| 2 | \(\lambda_0 f_+\) | \(\lambda_0 f_-\) | 0 | \(f_+\) |
| 3 | 0 | \(\lambda_0 f_+\) | \(\lambda_0 f_-\) | 0 |

完整计算重复这个三段协议两次。E₁ 恒定，而 t₁ 不是全程恒定。

理想耦合方向为 \(n=(|t_2|,-|t_3|,E_d)/\lambda\)，其垂直平面给出两个暗模。沿球面的平行输运基在各段的扰动投影为

| 段 | \(h\cdot d_1\) | \(h\cdot d_2\) |
|---|---|---|
| 1 | \(E\cos\alpha-\eta f_-\sin\alpha\) | 0 |
| 2 | \(-\eta f_+\) | \(E\sin\alpha\) |
| 3 | 0 | E |

其中 \(\alpha(v)=\operatorname{atan2}(f_-,f_+)\)。定义常数

\[
C=\int_0^1\sin\alpha\,dv
=\frac{\mathcal K(1/2)}\pi=0.5901702995,
\]

\[
A=\int_0^1f_-\sin\alpha\,dv
=\frac{\mathcal E(1/2)}\pi=0.4299233001,
\qquad B=A+\frac12=0.9299233001.
\]

\(\mathcal K(m)\)、\(\mathcal E(m)\) 分别是参数约定的第一、第二类完全椭圆积分。

一次三段协议的一阶暗模响应为

\[
Z=\frac{2\tau}{\hbar}(CE-B\eta,(1+C)E)^T.
\]

一次理想编织的暗模 holonomy 是

\[
J_d=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]

第二次协议在初始参考系中的响应是 \(J_d^T Z\)。因此双编织的总响应为

\[
\boxed{
z_{\rm dark}=Z+J_d^T Z
=\frac{2\tau}{\hbar}
\begin{pmatrix}-E-B\eta\\(1+2C)E-B\eta\end{pmatrix}.
}
\]

第一个分量沿 \(L_{12}\)，只增加初态相位；第二个分量沿 \(L_{13}\)，进入保真度损失。忽略绝热极限中被抑制的快模分量，得到

\[
\boxed{
1-F\simeq\frac{2\tau^2}{\hbar^2}
\bigl[(1+2C)E-B\eta\bigr]^2,
\qquad
1-f_{\rm paper}\simeq\frac{\tau^2}{\hbar^2}
\bigl[(1+2C)E-B\eta\bigr]^2.
}
\]

对 \(E>0\)，二阶损失的局部抵消线是

\[
\boxed{\eta/E=(1+2C)/B\simeq2.344645627.}
\]

对负 E 和非负 \(\eta\)，这个近似下两者增强损失，不存在正幅度的抵消线。等幅条件或门控路径改变后，以上常数和比值也会改变。

这个结果替代了 `213.md` 中恒定 t₁、单次三段协议、总混合量的最小点；两者分析的协议和观测量不同，不能混用。

## 7. 单参数的显式振荡：不要求累计角度小

在绝热暗模投影模型内，慢子空间的演化属于 SO(3)，按 \((1,d_1,d_2)\) 排序。定义

\[
x=2\tau E/\hbar,\qquad y=2\tau\eta/\hbar.
\]

一次协议的误差传播可以因子化成

\[
S_1=e^{(Cx-Ay)L_{12}},\qquad S_3=e^{xL_{13}},
\]

\[
S_2=\mathcal T\exp\int_0^1
[-y f_+(v)L_{12}+x\sin\alpha(v)L_{13}]dv.
\]

令 \(J=\operatorname{diag}(1,J_d)\)，则双编织慢传播为

\[
R_{\rm slow}=(J S_3 S_2 S_1)^2.
\]

对于任意这样的 SO(3) 矩阵，重叠公式进一步简化为

\[
\boxed{f_{\rm paper}=(1-R_{{\rm slow},33})/2.}
\]

当只有一个扰动参数时，每段所需生成元彼此共线，可直接积分并相乘：

\[
\boxed{
\begin{aligned}
\eta=0:\quad &f_{\rm paper}=\cos^2\left[\frac{\tau}{\hbar}(1+2C)E\right],\\
E=0:\quad &f_{\rm paper}=\cos^2\left[\frac{\tau}{\hbar}B\eta\right].
\end{aligned}
}
\]

相应的 F 是上式的平方，即 cos⁴。它们在投影模型内对任意累计角度成立，给出单参数导致的相干振荡及随耦合强度缩短的周期。完整五模模型仍会有有限速度、快模混合等修正。

当两个参数同时非零时，S₂ 中的 \(L_{12}\)、\(L_{13}\) 不对易，因此不能把小扰动公式中的线性组合直接放入 cos²，并当成一般精确解。

## 8. 抵消线为何仍有高阶损失

对慢模型做二阶 Magnus 展开，得到 \(L_{23}\) 方向的有向面积旋转。定义

\[
I=2\int_0^1\left[\frac v2+\frac{\sin\pi v}{2\pi}\right]\sin\alpha(v)dv-\frac C2
=0.2346331992,
\qquad H=A(1+C)+\frac12+I=1.4182844620.
\]

一次协议的面积系数为 \([C(1+C)x^2-Hxy]/2\)。两次协议的面积相加，并计入前后两次一阶位移之间的 BCH 对易项，得到

\[
\Theta_2=\frac12\left[(1+2C)^2x^2-2(H+CB)xy+B^2y^2\right].
\]

在抵消线 \(y=((1+2C)/B)x\) 上，

\[
\Theta_2=0.1417384015\,x^2,
\qquad
\boxed{1-F=0.0100448872\,x^4+O(x^6)}
\]

成立于慢投影模型。因一阶只剩不影响指标的相位，二阶的 \(L_{23}\) 旋转产生四阶损失。完整五模的快模响应不一定同时相消，可能保留较小的二阶损失。因此局部抵消不等于精确恢复编织。

## 9. 适用条件与验证结果

精确定义、几何恒等式以及第 5 节有限速度参数展开不依赖绝热投影。第 6 节简式要求

\[
\hbar\|\dot n\|/(2\lambda)\ll1,
\qquad
\frac2\hbar\int_0^{6\tau}\sqrt{E^2+\eta^2g(t)^2}\,dt\ll1.
\]

单参数的非小角度公式只去掉了第二个小角度要求，仍以绝热、弱于快模尺度的扰动和可忽略的投影修正为条件。它们不能用于推断无限长时间或非绝热极限。对固定有界控制，\(\tau\to0\) 时精确传播趋近单位元，论文目标重叠趋近零。

使用 `verify_so5_fidelity.py` 进行了小规模独立核对，数据写入 `so5_fidelity_validation.json`。取 \(t_c=E_0=0.3\) meV、\(\hbar=1\)，并以 \(x=2\tau E\)、\(y=2\tau\eta\) 表示参数。有限时间展开为 \(F=F_0+(x,y)H_{xy}(x,y)^T+O((|x|+|y|)^4)\)：

| \(\tau\)（meV⁻¹） | \(F_0\) | \(H_{xx}\) | \(H_{xy}\) | \(H_{yy}\) |
|---|---|---|---|---|
| 100 | 0.9999844193 | −2.370087 | 1.011951 | −0.432069 |
| 300 | 0.9999999879 | −2.376191 | 1.013576 | −0.432346 |
| 1000 | 0.9999999999 | −2.376875 | 1.013757 | −0.432376 |
| 绝热解析极限 | 1 | −2.376943 | 1.013775 | −0.432379 |

还核对了：

- 复波函数重叠与四个实矩阵元素的表达式完全一致，三个权重之和为 1。
- 有限速度参数展开的残差在扰动减半时缩小约 15.97、15.99 倍，符合四阶余项。
- 单参数 cos² 公式与独立的慢模型传播一致；在 \(\tau=1000\)、\(x=0.8\) 的 E-only 例子中，与完整五模重叠的差约为 \(1.24\times10^{-5}\)。
- 抵消线上慢模型的四阶系数趋近 0.0100448872。

这些核对支持上述局部响应和绝热解析结果；没有据此宣称它们能覆盖论文整张参数图。若要研究大累计角度下两个参数同时存在的行为，可以在第 7 节的 SO(3) 因子化基础上继续分析 S₂ 的非对易旋转。

源文件：`arts/prb111.pdf` 第二节、图 1，`arts/pq.txt` 的对应段落，以及 `model.md`。
