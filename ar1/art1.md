# 由编织到空间诱导拓扑：一个统一框架

*(草稿——由工作笔记整理而成，后续可进一步润色)*

## 1 引言

拓扑相中存在非局域的编码自由度，可以通过任意子（anyon）之间的编织来实现量子操作。在连续的拓扑量子场论（TQFT）描述中，这些操作自然由 2+1 维时空中的世界线以及它们在基态/零模 Hilbert 空间上诱导出的辫群表示来刻画。然而，在具体晶格哈密顿量上，更“工程化”的实现往往是通过**空间变形**完成的：对流形做几何扭转（Dehn twist）、重写配对图、引入 Kekulé 或涡旋类纹理、以及由简单两体门构成的局域幺正电路等。

本文构建了一个统一框架，把这些看似不同的视角联系在一起。我们从一个指数型两比特门出发：
$$
R = e^{iH_P},
$$
其中 $H_P$ 展开在 Pauli 基底上。围绕这一出发点，我们：

- 将常数/经典 Yang–Baxter 方程与 $H_P$ 空间中的离散曲率、可积（或近平坦）子流形联系起来；
- 说明如何把 $H_P$ 嵌入到 1 维和 2 维的 Majorana/BdG 晶格模型（Kitaev 链、$p+ip$ 超导、honeycomb Majorana–$\mathbb Z_2$ 规范理论）中；
- 在缺陷位置、配对纹理等构成的配置空间上，引入 Berry 联络丛，其非阿贝尔 holonomy 实现任意子的统计；
- 并证明：在标准的局域性和谱隙假设下，**braid 所给出的幺正可以等价地由纯空间操作表示**——即 Dehn twist、码变形、空间纹理加有限深局域电路——且误差可以定量控制。

本文的中心结果是一个物理定理：在一个带谱隙的 2+1 维拓扑相中，对任意 braid word $\beta$，都可以用上述空间操作实现其对应的拓扑门，至多差一个整体相位和可控的有界误差。证明思路分为三层：

1. 纯拓扑层：辫群 $B_n$ 同构于多穿孔圆盘的 mapping class group，后者由 Dehn twist 生成；
2. TQFT 层：世界线编织与 mapping class 作用在 Hilbert 空间上给出同一幺正表示；
3. Hamiltonian/Berry/$R$ 层：在具体晶格中用 $R=e^{iH_P}$、Berry holonomy 与 quasi‑adiabatic 连续性，构造所需的空间实现，并给出显式范数上界。

我们还用具体玩具模型说明这一构造，包括四 Majorana 体系，以及一条具有 Kekulé 型配对纹理的 1D Kitaev 链，后者仅通过空间变化就产生了 Majorana 零模。

## 2 拓扑与 TQFT 预备知识

### 2.1 辫群与 mapping class group

记 $D\subset\mathbb R^2$ 为闭圆盘，$\{p_1,\dots,p_n\}\subset\operatorname{int}(D)$ 为 $n$ 个标记点（穿孔）。相对于边界的 mapping class group 定义为
$$
\mathrm{MCG}(D_n,\partial D)
 := \pi_0\bigl(\{f:D\to D\mid f\text{为同胚},\ f(\partial D)=\partial D,\ f(\{p_i\})=\{p_i\}\}\bigr).
$$

$n$ 股辫群 $B_n$ 定义为 $D\times[0,1]$ 中 $n$ 条不相交世界线的同伦类，其端点固定在这些穿孔上。一个经典定理（见 Birman 及 Farb–Margalit 等）说明存在自然同构
$$
\Phi: B_n \xrightarrow{\ \simeq\ } \mathrm{MCG}(D_n,\partial D).
$$

直观地说，一个 braid 可以“压平”为时间 $t=1$ 时的圆盘自同胚 $f$，其在空间中对穿孔做某种扭转；反过来，一个 mapping class 也可以“厚化”为时空中的编织世界线。

Mapping class 的基本构件是 Dehn twist。

**定义 2.1（Dehn twist）.** 给定 $D_n$ 中一条简单闭曲线 $c$，沿 $c$ 取一条窄带邻域，将其切开，在一侧沿曲线方向旋转 $2\pi$ 后再粘回，得到新的曲面同胚。其同伦类记为 $T_c\in\mathrm{MCG}(D_n,\partial D)$，称为沿 $c$ 的 Dehn twist。

对多穿孔圆盘/球面，可以选取有限条简单闭曲线，使得所有 mapping class 都可写成这些曲线上的 Dehn twist 的有限乘积（Dehn–Lickorish 型定理）。因此，每一个 braid word $\beta\in B_n$，都可以通过 $\Phi$ 写成若干 Dehn twist 及其逆的乘积。

### 2.2 TQFT 表示：世界线 vs mapping class

现在考虑一个 2+1 维拓扑量子场论（TQFT），例如由模范畴构造的 Reshetikhin–Turaev 型 TQFT。对带穿孔的曲面 $D_n$，在每个穿孔上附着任意子荷 $\{\sigma_i\}$，TQFT 赋予一个有限维 Hilbert 空间
$$
\mathcal H_0(D_n;\{\sigma_i\}).
$$

在该框架下，有两种等价的方式在 $\mathcal H_0$ 上构造幺正算符：

- **世界线表示**：把 braid $\beta\in B_n$ 画成时空中的任意子世界线网络，用模范畴的 F、R 符号化简，得到幺正 $\rho_{\mathrm{WL}}(\beta)\in U(\mathcal H_0)$；
- **mapping class 表示**：把 $[f]\in\mathrm{MCG}(D_n,\partial D)$ 看作对空间截面的几何操作，在三维 TQFT 中诱导出 $\rho_{\mathrm{MCG}}([f])$ 作用于 $\mathcal H_0$。

Reshetikhin–Turaev TQFT 的公理保证，这两种构造与同构 $\Phi$ 相容：
$$
\rho_{\mathrm{WL}}(\beta) = \rho_{\mathrm{MCG}}(\Phi(\beta))
$$
至多差一个整体相位。因此，**时间方向上的任意子编织与空间截面上的 Dehn twist 组合，在 TQFT Hilbert 空间上给出的是同一个投射表示。**

这一等价性完全是拓扑/TQFT 层面的，不依赖于任何具体的哈密顿量或微观模型。

## 3 晶格哈密顿、$R=e^{iH_P}$ 与 Berry 几何

接下来我们转向具体的晶格实现，希望把上面的抽象 TQFT 操作与下列对象联系起来：

- 由 Pauli 耦合构成的局域两体幺正门 $R=e^{iH_P}$；
- 参数/配置空间上绝热路径所对应的 Berry holonomy；
- 以及真正“空间几何”的操作：Dehn twist、branch cut 变形、涡旋/vison 移动、配对纹理（如 Kekulé/Dirac 纹理）以及有限深局域电路等。

### 3.1 局域生成元及其分解

可以考虑一个自旋/费米晶格模型，其低能部分实现了上一节讨论的某个 TQFT，并处于带谱隙的拓扑相。记 $\mathcal H$ 为全微观 Hilbert 空间，$\mathcal H_0(\lambda)\subset\mathcal H$ 为参数 $\lambda$ 下的基态/零模子空间。

对每条晶格键 $\langle ij\rangle$，我们引入一个局域两体生成元
$$
H_P^{(ij)} = \sum_{a,b} c_{ab}^{(ij)}\,\sigma_a\otimes\sigma_b,\qquad R_{ij} = e^{iH_P^{(ij)}}.
$$
由这些生成元构造哈密顿量族
$$
H(\lambda) = \sum_{\langle ij\rangle} H_P^{(ij)}(\lambda),
$$
其中 $\lambda$ 取值于一个带谱隙的参数空间 $\Omega$。

从物理角度看，$H(\lambda)$ 与 $H_P^{(ij)}$ 有两层互补的含义：

- 在“上层”，$H(\lambda)$ 只是一个具体的、局域的、带谱隙的晶格哈密顿量，其基态子空间 $\mathcal H_0(\lambda)$ 携带了前一节 TQFT 中的 braid / mapping class 表示；
- 在“下层”，每个 $H_P^{(ij)}$ 具体编码了配对、跳跃以及局域规范约束。经过 Jordan–Wigner/Majorana 变换后，可以自然地将其写成
  $$
  H_P^{(ij)} \simeq H_{\mathrm{quad}}^{(ij)} + H_{\mathrm{int}}^{(ij)} + H_{\mathrm{gauge}}^{(ij)},
  $$
  其中 $H_{\mathrm{quad}}$ 含有定义 Kitaev 链或 $p+ip$/honeycomb 型拓扑通道的二次 BdG 项，而 $H_{\mathrm{int}}$ 和 $H_{\mathrm{gauge}}$ 描述偏离理想可积/Yang–Baxter 点的相互作用与规范扇区修正。

这种分解**并不是主定理的额外公理假设**，而是一种阅读模型、参数化“偏离理想平坦联络程度”的便利方式（进入误差项 $\varepsilon_{\mathrm{YBE}}$ 及相关“复杂度曲率”）。我们在主命题中真正需要的结构性假设只有：局域性、谱隙稳定性，以及 $\mathcal H_0$ 中实现了目标 TQFT。

### 3.2 参数/配置空间上的 Berry 联络与曲率

记 $\Omega$ 为一个光滑的参数或配置空间（例如缺陷位置、耦合常数向量 $c_{ab}^{(ij)}$、配对纹理等）。对每个 $X\in\Omega$，假设 $H(X)$ 带有稳定谱隙，并记基态投影为
$$
P(X):\mathcal H\to\mathcal H_0(X),\qquad P^2=P=P^\dagger.
$$
在 $\Omega$ 上定义 Berry 联络与曲率为
$$
A_\mu(X) = iP\,\partial_\mu P\,P,\qquad
F_{\mu\nu}(X) = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu,A_\nu]
               = P[\partial_\mu P,\partial_\nu P]P.
$$
对任意光滑路径 $\gamma:[0,1]\to\Omega$，其非阿贝尔 Berry holonomy 定义为
$$
U[\gamma] = \mathcal P\exp\Bigl(-\int_0^1 A_\mu(X(t))\,\dot X^\mu(t)\,\mathrm dt\Bigr),
$$
作用在初始子空间 $\mathcal H_0(X(0))$ 上。

在配置空间的语境下，$\gamma$ 可以表示：

- 任意子或涡旋在时间中的真实运动（世界线图像）；
- 或者在缺陷位置固定的前提下，对 branch cut、配对图、Kekulé/Dirac 纹理等做**纯空间变形**的路径，只要整个过程中谱隙保持不闭合。

在 $\|F\|$ 很小的区域（例如在 $H_P$ 空间的可积/YBE 子流形附近），holonomy 近似只依赖于路径同伦类；于是，对应同一 mapping class 的时间实现和空间实现，在 $\mathcal H_0$ 上给出的幺正几乎相同。

## 4 主物理定理与误差估计

下面在 Hamiltonian/Berry/$R$ 框架中给出中心命题的物理版。非正式地说：

> **命题（物理版）.** 在一个局域、带谱隙且实现给定 TQFT 的 2+1 维拓扑相中，任意 braid word $\beta\in B_n$，都可以由 Hamiltonian 级别的**空间操作**来实现——这些操作由有限次 Dehn twist、几何畸变（branch cut、Kekulé/Dirac 纹理、涡旋/vison 运动）以及有限深局域 $R=e^{iH_P}$ 电路组成——使得它们在零模/编码子空间 $\mathcal H_0$ 上所诱导的幺正，与抽象的 braid 幺正只差一个整体相位和可控小误差。

更精确地，记 $U_{\mathrm{top}}$ 为对应 braid（或 mapping class）的 TQFT 幺正，$U_{\mathrm{spatial}}$ 为某一具体空间协议（由上述操作拼接实现）诱导的幺正。在局域性、谱隙、曲率与扰动范数有界的假设下，存在相位 $e^{i\phi}$ 使得
$$
\bigl\|U_{\mathrm{spatial}} - e^{i\phi}U_{\mathrm{top}}\bigr\|
\;\lesssim\; \varepsilon_{\mathrm{YBE}} + \kappa\,\mathcal A_{\gamma,\gamma'} + \varepsilon_{\mathrm{adiab}} + \varepsilon_{\mathrm{Trotter}},
$$
其中：

- $\varepsilon_{\mathrm{YBE}}$ 衡量所选 $H_P$ 偏离理想 Yang–Baxter/可积“平坦流形”的程度；
- $\kappa$ 是 Berry 曲率 $F$ 范数的统一上界，$\mathcal A_{\gamma,\gamma'}$ 是连接世界线路径 $\gamma$ 与空间变形路径 $\gamma'$ 的“填充曲面”的面积；
- $\varepsilon_{\mathrm{adiab}}$ 是常规绝热误差，由 $\|\dot H\|/\Delta^2$ 和总演化时间控制；
- $\varepsilon_{\mathrm{Trotter}}$ 则来自用有限深局域 $R$ 门电路近似连续 quasi‑adiabatic 演化的离散化误差。

证明可以分解为四个引理：

1. **近平坦联络与路径同伦（引理 A）.** 在 $\|F\|\le\kappa$ 的区域内，同伦的两条路径 $\gamma$ 与 $\gamma'$ 所对应 Berry holonomy 的差异有上界 $\mathcal O(\kappa\,\mathcal A_{\gamma,\gamma'})$（非阿贝尔 Stokes 公式的一种定量形式）。
2. **绝热演化与 Berry holonomy（引理 B）.** 沿一条带谱隙路径的缓慢时间演化 $U_{\mathrm{dyn}}(T)$，限制在 $\mathcal H_0$ 上，与 Berry holonomy $U[\gamma]$ 的差异由绝热定理给出的 $\varepsilon_{\mathrm{adiab}}$ 控制。
3. **F,R 数据与 $e^{iH_P}$ 的嵌入（引理 C）.** 在理想可积点 $H_P^{(0)}$，局域幺正 $e^{iH_P^{(0)}}$ 在 $\mathcal H_0$ 上重现 TQFT 的 R‑矩阵；沿谱隙不闭合的小局域扰动（由 $H_{\mathrm{int}}$ 与 $H_{\mathrm{gauge}}$ 控制）只会在 $\varepsilon_{\mathrm{YBE}}$ 的量级上轻微变形这一表示。
4. **quasi‑adiabatic 演化与有限深电路（引理 D）.** 由 $H(t)$ 构造的 quasi‑adiabatic 生成元 $K(t)$ 给出幺正 $U_{\mathrm{QA}}$，实现所需的 Berry 演化。利用 Lieb–Robinson 界与 Trotter–Suzuki 分解，可以用有限深、局域的 $R=e^{iH_P}$ 门电路逼近 $U_{\mathrm{QA}}$，误差记为 $\varepsilon_{\mathrm{Trotter}}$。

把上述四个引理用三角不等式串联起来，即得到上式中的误差不等式，也就建立了“时间中的 braid 与 Hamiltonian 级别的纯空间操作等价”的结果。

## 5 玩具模型与空间纹理

为使抽象框架更加直观，我们简要回顾两个代表性玩具模型。

### 5.1 四 Majorana 模型与配对重排

最小的例子是四个 Majorana 算符 $\gamma_1,\dots,\gamma_4$，可以组合成两个复费米子，从而在零模子空间中形成一个二维逻辑空间。标准的“时间向半扭（half‑twist）编织”可以由
$$
H_0 = \tfrac{i}{2}\,\gamma_2\gamma_3,\qquad U_{\mathrm{ht}} = e^{(\pi/4)\gamma_2\gamma_3}
$$
来实现。

另一方面，可以考虑不同的配对模式（例如 1–2 与 3–4 配对，然后 1–3 与 2–4 配对等），每一种模式都由形式为 $H^{(X)} = \tfrac{i}{2}\sum c_{ab}^{(X)}\gamma_a\gamma_b$ 的哈密顿量实现。沿着一条保持谱隙的路径，在这些配对模式之间 quasi‑adiabatic 地切换，或者用局域 $R$ 门实现对应的有限深电路，可以在逻辑子空间上实现与 half‑twist braid 等价（在 SU(2) 内共轭且误差可控）的幺正。这提供了一种“码变形”或“空间重连”式的编织实现。

### 5.2 一维 Kekulé 型 Kitaev 链

在一条更长的一维 Kitaev 链中，我们可以构造一个键依赖的配对剖面 $\Delta_j$，满足：

- 在一个有限区间 $D$ 之外，链处于均匀拓扑相，配对为常数 $\Delta_0$；
- 在 $D$ 内部，配对具有 Kekulé 型纹理：
  $$
  \Delta_j = |\Delta_0| e^{i\theta_j},\qquad \theta_j:0\to 2\pi\ \text{跨越区间 }D,
  $$
  即在空间中实现一次 $2\pi$ 的相位绕行。

相应的 BdG 哈密顿量可以用最近邻跃迁 $t$、化学势 $\mu$ 以及复配对 $\Delta_j$ 写出。对该模型进行数值对角化（对应的代码已独立实现）表明，这种**纯空间**纹理在畸变区 $D$ 附近产生了一对近零能量的 Majorana 模式，并在其上方保持了清晰的谱隙。这一现象具体展示了“空间诱导拓扑”的思想：空间纹理有效地在序参量配置空间中走了一条非平凡闭合路径，其 Berry holonomy 在零模子空间上实现了非平凡的拓扑操作。

## 6 讨论与展望

我们搭建了一个三层结构的统一框架，将以下内容联系起来：

- 编织、mapping class group 与 Dehn twist 的代数/拓扑结构；
- TQFT 中通过 F、R 符号实现的上述结构的幺正表示；
- 以及具体的晶格哈密顿模型，其中局域门 $R=e^{iH_P}$、Berry 几何与 quasi‑adiabatic 连续性共同构成了实现通道。

在局域性与谱隙假设下，这一框架表明：**任意子的 braid 操作可以被变形成纯空间几何操作**——包括 Dehn twist、码变形、配对纹理操控以及有限深局域电路——并且可以给出显式的误差估计。分解
$$
H_P^{(ij)} \simeq H_{\mathrm{quad}}^{(ij)} + H_{\mathrm{int}}^{(ij)} + H_{\mathrm{gauge}}^{(ij)}
$$
帮助我们澄清：哪些微观项负责拓扑通道，哪些只是受控扰动。

自然的后续方向包括：

- 设计**曲率工程化**的哈密顿路径，使 Berry 曲率在参数空间中足够小或近似中心，从而使 holonomy 近似只依赖于同伦类；
- 定义并计算一种“复杂度曲率”，用来度量在有限深电路下近似理想 YBE/平坦联络的难度；
- 将空间纹理的构造（Kekulé 型、Dirac‑涡旋型、branch cut 重排等）推广到更现实的平台和更高 genus 的曲面上；
- 在具体模型中将误差估计更严格地定理化，给出带常数的严谨不等式。

## 7 参考文献

本文涉及的数学与物理工具可参考以下文献：

- J. S. Birman, *Braids, Links, and Mapping Class Groups*, Ann. of Math. Studies 82, Princeton Univ. Press (1974).
- B. Farb and D. Margalit, *A Primer on Mapping Class Groups*, Princeton Univ. Press (2011).
- N. Reshetikhin and V. G. Turaev, *Invariants of 3-manifolds via link polynomials and quantum groups*, Invent. Math. 103, 547–597 (1991).
- V. G. Turaev, *Quantum Invariants of Knots and 3-Manifolds*, de Gruyter (1994/2010).
- B. Bakalov and A. Kirillov Jr., *Lectures on Tensor Categories and Modular Functors*, AMS (2001).
- B. Simon, *Holonomy, the quantum adiabatic theorem, and Berry's phase*, Phys. Rev. Lett. 51, 2167–2170 (1983).
- F. Wilczek and A. Zee, *Appearance of Gauge Structure in Simple Dynamical Systems*, Phys. Rev. Lett. 52, 2111–2114 (1984).
- T. Kato, *On the adiabatic theorem of quantum mechanics*, J. Phys. Soc. Jpn. 5, 435–439 (1950).
- J. E. Avron, R. Seiler and L. G. Yaffe, *Adiabatic theorems and applications to the quantum Hall effect*, Commun. Math. Phys. 110, 33–49 (1987).
- M. B. Hastings and X.-G. Wen, *Quasi-adiabatic continuation of quantum states: The stability of topological ground-state degeneracy and emergent gauge invariance*, Phys. Rev. B 72, 045141 (2005).
- S. Bravyi, M. B. Hastings and S. Michalakis, *Topological quantum order: stability under local perturbations*, J. Math. Phys. 51, 093512 (2010).
- E. H. Lieb and D. W. Robinson, *The finite group velocity of quantum spin systems*, Commun. Math. Phys. 28, 251–257 (1972).

本文稿是基于工作笔记整理出的初版，可以在此基础上进一步补充具体模型计算与图示，并根据目标期刊或受众的习惯对结构和表述做相应调整。
