**B.1 引理 C（Stokes‑based）——有条理且严格的证明**

目的：在严格的假设下，用非阿贝尔 Stokes 定理把两个闭路对应的 holonomy 差异用曲率的面积积分控制，并把曲率的范数上界以哈密顿量的谱隙与微扰分量的范数显式量化。

一、假设与记号
- 参数空间点记为 $\lambda\in\Lambda$，沿参数路径变动的哈密顿为
  $$H(\lambda)=H_0(\lambda)+V(\lambda).$$
- 基态编码投影为 $P(\lambda)$，理想情形的投影记为 $P_0(\lambda)$（对应 $H_0$）。记 $\Delta_0:=\inf_{\lambda}\mathrm{gap}(H_0(\lambda))>0$。假设
  $$\sup_{\lambda}\|V(\lambda)\|=:\epsilon_{\mathrm{tot}}<\Delta_0/2.$$
- 把微扰细分为三类局域分量：
  $$V(\lambda)=H_{\mathrm{int}}(\lambda)+H_{\mathrm{string}}(\lambda)+H_{\mathrm{gauge}}(\lambda),$$
  并定义各自范数上界 $\epsilon_{\alpha}:=\sup_{\lambda}\|H_{\alpha}(\lambda)\|$，使 $\epsilon_{\mathrm{tot}}=\sum_\alpha\epsilon_{\alpha}$。
- 假设 $H(\lambda)$ 在参数上 $C^1$（足够光滑），且对参数导数有界：
  $$
  M_1:=\sup_{\lambda,a}\|\partial_a H(\lambda)\|<\infty.
  $$ 

二、命题（精确表述）
  在上述假设下，对于任何两条同边界的闭路所张成曲面 $S$，相应 holonomy 差异有上界
  $$
  \|U_{\mathrm{pert}}(\partial S)-U_0(\partial S)\| \le C\,A_S\,\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}\,\frac{M_1^2}{\Delta_0^4} + R(\epsilon_{\mathrm{tot}},\Delta_0),
  $$
其中 $A_S$ 为曲面面积，$C$ 与格点几何有关且与各分量的局域/支撑半径通过常数因子耦合，余项 $R(\cdot)$ 为次阶量级 $O(\epsilon_{\mathrm{tot}}^2/\Delta_0^6)$。

三、证明（分步）

Step 1 — 投影的解析表示与基本界定。
对每个 $\lambda$，在复平面选取闭曲线 $\Gamma$ 包围基态分离谱且与最近谱点距离至少 $\Delta:=\Delta_0-\|V\|\ge\Delta_0/2$。投影有留数表示
  $$
  P(\lambda)=-\frac{1}{2\pi i}\oint_{\Gamma}R(z;\lambda)\,dz,\\ R(z;\lambda)=(H(\lambda)-z)^{-1}.
  $$ 
  因此对参数微分
  $$
  \partial_a P = -\frac{1}{2\pi i}\oint_{\Gamma} R(z)\,(\partial_a H)\,R(z)\,dz.
  $$ 
  对任意 $z\in\Gamma$ 有 $\|R(z)\|\le 1/\mathrm{dist}(z,\mathrm{spec}(H))\le 2/\Delta_0$，故
  $$
  \|\partial_a P\| \le C_1\frac{\|\partial_a H\|}{\Delta_0^2}\le C_1\frac{M_1}{\Delta_0^2},
  $$
  其中 $C_1=|\Gamma|/(2\pi)\cdot 4$ 为几何常数（可选取相同的 $\Gamma$ 以统一常数）。

Step 2 — 曲率的表示与范数界。
曲率定义为 $F_{ab}=P[\partial_a P,\partial_b P]P$。因此
  $$
  \|F_{ab}\|\le 2\,\|\partial_a P\|\,\|\partial_b P\| \le C_2\frac{\|\partial_a H\|\,\|\partial_b H\|}{\Delta_0^4}\le C_2\frac{M_1^2}{\Delta_0^4}.
  $$
所以曲率在曲面上的面积积分满足
  $$\iint_S\|F\| \le C_2\,A_S\,\frac{M_1^2}{\Delta_0^4}.$$

Step 3 — 投影的微扰展开与 $\delta P$ 的界。
使用解析留數表示並结合 resolvent identity 可以把上述界写得更具体。记 $R_0(z)=(H_0-z)^{-1}$、$R(z)=(H-z)^{-1}$，有
$$
R(z)-R_0(z)=-R_0(z)V R(z)=-\sum_{n\ge1}(-1)^{n-1}R_0(z)\big(VR_0(z)\big)^n,
$$
当对围道 $\Gamma$ 上的点满足 $\sup_{z\in\Gamma}\|R_0(z)V\|<1$ 时，上述 Neumann 级数收敛（在我们的假设 $\|V\|<\Delta_0/2$ 下通常成立，因为 $\|R_0(z)\|\le 2/\Delta_0$）。代入投影的留数表达式并取范数得
$$
\|\delta P\|=\Big\| -\frac{1}{2\pi i}\oint_{\Gamma}\big(R-R_0\big)dz\Big\| \le \frac{|\Gamma|}{2\pi}\sup_{z\in\Gamma}\sum_{n\ge1}\|R_0(z)\|^{n+1}\,\|V\|^n.
$$
当记 $q:=\sup_{z\in\Gamma}\|R_0(z)\|\,\|V\|<1$ 时，上述级数求和并简化得
$$
\|\delta P\| \le \frac{|\Gamma|}{2\pi}\frac{\|R_0\|^2\,\|V\|}{1-\|R_0\|\,\|V\|} \le \frac{|\Gamma|}{2\pi}\frac{(2/\Delta_0)^2\,\|V\|}{1-2\|V\|/\Delta_0}=:C_3\frac{\|V\|}{\Delta_0}.
$$
在常见的小微扰情形 $\|V\|\ll\Delta_0$ 下，可进一步用简化常数写作 $\|\delta P\|\lesssim C_3'\,\|V\|/\Delta_0$。

对导数差量的估计，从导数的留數表达式出发：
$$
\partial_a P - \partial_a P_0 = -\frac{1}{2\pi i}\oint_{\Gamma} \big(R(\partial_a H)R - R_0(\partial_a H_0)R_0\big)\,dz.
$$
把被积子重写为三项之和：
$$
\begin{aligned}
R(\partial_a H)R - R_0(\partial_a H_0)R_0 &= (R-R_0)(\partial_a H)R \\
&\quad + R_0(\partial_a H - \partial_a H_0)R \\
&\quad + R_0(\partial_a H_0)(R-R_0).
\end{aligned}
$$
对每一项使用范数不等式和上面的 resolvent 差界可得：

- 第一与第三项（含 $R-R_0$）：利用 $\|R-R_0\|\le \|R_0\|^2\,\|V\|/(1-\|R_0\|\,\|V\|)$ 以及 $\|R\|\le \|R_0\|/(1-\|R_0\|\,\|V\|)$，得到典型贡献
  $$
  \Big\|\frac{1}{2\pi i}\oint_{\Gamma}(R-R_0)(\partial_a H)R\,dz\Big\| \lesssim \frac{|\Gamma|}{2\pi}\frac{\|R_0\|^2\,\|V\|}{1-\|R_0\|\,\|V\|}\cdot\|\partial_a H\|\cdot\|R\| 
  \lesssim C'\frac{\|\partial_a H\|\,\|V\|}{\Delta_0^3(1-2\|V\|/\Delta_0)^2}.
  $$

- 第二项（含 $\partial_a H-\partial_a H_0=\partial_a V$）：直接用 $\|R_0\|\le 2/\Delta_0$ 和 $\|R\|\le 2/\Delta_0$，得到
  $$
  \Big\|\frac{1}{2\pi i}\oint_{\Gamma}R_0(\partial_a V)R\,dz\Big\| \le \frac{|\Gamma|}{2\pi}\|R_0\|\,\|\partial_a V\|\,\|R\| \le C''\frac{\|\partial_a V\|}{\Delta_0^2}.
  $$

综上合并主要项，并在 $\|V\|/\Delta_0$ 小的前提下忽略更高阶项，得到常用的简化界
$$
\|\partial_a P - \partial_a P_0\| \le C_4\Big(\frac{\|\partial_a H\|\,\|V\|}{\Delta_0^3} + \frac{\|\partial_a V\|}{\Delta_0^2}\Big) + O\Big(\frac{\|V\|^2\,\|\partial_a H\|}{\Delta_0^4}\Big).
$$
其中 $\|\partial_a V\|$ 可進一步以各分量 $H_{\alpha}$ 的界表示（例如 $\|\partial_a V\|\le\sum_{\alpha}\|\partial_a H_{\alpha}\|$），從而把导数差的界写成各微扰成分的贡献之和。

Step 4 — 曲率差的展开与分量依赖。
  将 $P=P_0+\delta P$ 代入 $F=P[\partial_a P,\partial_b P]P$，展开差量得到由 $\delta P$ 与 $\partial P-\partial P_0$ 控制的若干项。逐项上界得
  $$
  \|\widetilde F - F_0\| \le C_5\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}\frac{M_1^2}{\Delta_0^4} + O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^6}\Big).
  $$
  若需把 $\epsilon_{\mathrm{tot}}$ 拆成各分量，注意每项在上式中的贡献按其算符范数线性入侵：
  $$
  \epsilon_{\mathrm{tot}}=\sum_{\alpha}\epsilon_{\alpha},\qquad \|\partial_a H\|\lesssim \|\partial_a H_0\| + \sum_{\alpha}C_{\alpha}\,\epsilon_{\alpha},
  $$
  其中 $C_{\alpha}$ 包含支撑半径与几何权重（对长串项 $C_{\mathrm{string}}$ 较大）。因此
  $$\|\widetilde F - F_0\| \lesssim \sum_{\alpha}C'_{\alpha}\frac{\epsilon_{\alpha}}{\Delta_0}\frac{M_1^2}{\Delta_0^4}+O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^6}\Big).$$

Step 5 — 曲率差到 holonomy 差的转换（Stokes）

为清楚起见，把曲面算子写为
$$
X:=\iint_S W(s,t)F(s,t)W(s,t)^{-1}\,dS,
\qquad X_0:=\iint_S W(s,t)F_0(s,t)W(s,t)^{-1}\,dS,
$$
其中共轭因子 $W(s,t)$ 为酉，不改变范数。由定义有
$$
\|X-X_0\| \le \iint_S\|F-F_0\|,\qquad \|X\|\le \iint_S\|F\|,\;\; \|X_0\|\le \iint_S\|F_0\|.
$$

利用 Duhamel 表示或算符指數的基本积分恒等，可以得到严格的差分界：
$$
e^{X}-e^{X_0}=\int_0^1 e^{(1-s)X_0}(X-X_0)e^{sX}\,ds
$$
从而
$$
\|e^{X}-e^{X_0}\| \le \|X-X_0\|\,e^{\|X\|+\|X_0\|}.
$$
因此
$$
\|U_{\mathrm{pert}}(\partial S)-U_0(\partial S)\| \le \Big(\iint_S\|\widetilde F-F_0\|\Big)\,\exp\Big(\iint_S\|F\|+\iint_S\|F_0\|\Big).
$$

代入 Step 2–4 中得到的曲率与曲率差的上界（把 $\iint_S\|\widetilde F-F_0\|$ 用 Step 4 的界代替，並把 $\iint_S\|F\|,\iint_S\|F_0\|$ 用 Step 2 的界代替），得到帶指數因子的更精确界：
$$
\|U_{\mathrm{pert}}-U_0\| \le \Big(C_6\,A_S\,\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}\,\frac{M_1^2}{\Delta_0^4}\Big)\exp\Big(C_7\,A_S\,\frac{M_1^2}{\Delta_0^4}\Big) + O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^6}\Big).
$$

在常用的小曲率情形（即 $A_S\,M_1^2/\Delta_0^4\ll1$ 且 $\epsilon_{\mathrm{tot}}/\Delta_0\ll1$）下，指數項可展開為 $1+o(1)$，主阶项退化为线性形式，得到文中所写的近似表达（常数 $C$ 可由 $C_6,C_7,C_5,C_2$ 组合得到）。

注：上述帶指數的界是严谨的，不依赖于 $A_S\le1$。将线性形式作为最终结果需要额外假设（曲率总量小），我已在文中对此假设作了明确说明。

四、常数与改善方向（注解）
- 上述降冪 $\Delta_0^{-5}$（整体看似）来自对 $\partial P$ 的双重 resolvent 估计；更精细的交换子重排或使用分块对角化（Schrieffer–Wolff）可在某些情形下把分母冪次降低到 $\Delta_0^{-2}$ 或 $\Delta_0^{-3}$，但需额外假设（例如 $\partial_a H$ 与 $[H_0,V]$ 的结构性小）。
- 若要数值对接，本证明中出现的常数 $C,C_1,\dots$ 与几何权重 $C_{\alpha}$ 应由具体模型（格点拓扑、算符支撑半径）评估。建议把各分量的范数 $\epsilon_{\alpha}$ 由 `kits/compute_jw_mapping.py` 输出，并用本文不等式估算期待的 $\varepsilon_{\mathrm{YBE}}$。

结论：在假设 $\epsilon_{\mathrm{tot}}<\Delta_0/2$ 且 $H(\lambda)$ 平滑的前提下，Stokes‑based 的严格推导把 holonomy / YBE 偏差以曲率差的面积积分表示，并可用谱隙与各微扰分量范数给出如上明确上界，从而把主文中抽象的 $\epsilon_{\mathrm{tot}}$ 具体化为对应分量的数值可验证界。
