### 命题：任意 braid 在拓扑相中可等价为空间几何操作（Dehn twist 组合）

本节在纯拓扑 → TQFT → 本文 R = exp(iH\_P) / Berry 框架三个层次上，给出以下命题的严谨表述与推导骨架。

```
2+1 维：指二维空间加一维时间（即在平面上的二维晶格系统，随时间演化的量子体系）。
带谱隙（gapped）：系统的能谱在基态（或基态简并子空间）与最低能激发之间存在严格的能量差 Δ>0（在考虑的参数区间或热力学极限下不趋于 0）。数学上可表为：对哈密顿量 H，基态能量集合与其余谱之间有下界 Δ>0。
拓扑相（topological phase）：一种零温下的物质相，其低能性质不是由局域对称破缺描述，而由全局/拓扑不变量刻画（例如基态简并依赖基底拓扑、存在任意子/非阿贝尔统计、边缘态与拓扑量子数等）。基态之间无法由局域算符区分（长期纠缠结构）。
```


```
指存在一组（family）哈密顿量 H(λ)（或一系列局域哈密顿量的连续/离散变化），参数 λ 沿某条路径变化，使得对应的低能态投影 P(λ) 在零模/编码子空间上产生所需的幺正变换（即实现 Dehn‑twist/space‑reconstruction）。
“Hamiltonian‑level” 强调用修改或局域调控哈密顿量本身来实现操作（例如改变耦合/配对图），而不是仅靠外加测量或把任意子在物理空间中移动的显式时间演化。
“一族”隐含两种常见实现方式：
连续路径 H(λ), 用绝热/ quasi‑adiabatic 演化实现（产生 Berry holonomy）；
离散序列 H_k 或 局域门 R=e^{iH_P} 的有限深电路，作为对连续生成元的 Trotter 近似。
为什么需要“家族”：单个哈密顿量对应静态模型，无法产生非平凡 holonomy；通过参数化家族并沿路径变形，低能子空间会累积几何相或投影算符的路径有序指数，从而得到所需的拓扑幺正。
物理要求（与误差项关联）：该族需保持体态谱隙（Δ>0）与局域性，以保证绝热近似、quasi‑adiabatic 构造和局域电路逼近的误差可控（误差由 ε_YBE, ε_adiab, ε_Trotter, κ·面积 等项给出）。
举例（直观）：给定一条闭合配置空间路径 γ，构造 H(λ)（λ∈γ）＝在不同 λ 下局域改变配对/Kekulé 模式；缓慢沿 λ 变化或用序列 R 门逼近，就能在编码子空间上实现与 braid/Dehn‑twist 等价的幺正。
```


> **命题（物理版）.** 在一个带谱隙的 2+1 维拓扑相中，
> 给定任意 braid word $\beta\in B_n$，存在一族 Hamiltonian‑level 的纯空间操作（由 Dehn twist / 几何畸变 / 局域 $R=e^{iH_P}$ 门重构给出），
> 使得在零模/编码子空间上，其非阿贝尔 Berry holonomy 与沿时间编织任意子的 holonomy 等价（在 SU(2) 内共轭，误差由曲率/YBE 偏差、非绝对慢走与离散化精度控制）。

下面依次在三个层次上论证：

1. 纯拓扑：$B_n \cong \mathrm{MCG}(D_n,\partial D)$，且 $\mathrm{MCG}$ 由 Dehn twist 生成；
2. TQFT：世界线 braid 表示与 mapping class 表示一致；
3. R/Berry：在 R = exp(iH\_P) 框架中构造空间操作并用 Berry 几何控制误差。

---

#### 1. 纯拓扑层：辫群与 Dehn twist 生成的 mapping class group 等价

记 $D\subset\mathbb R^2$ 为闭圆盘，$\{p_1,\dots,p_n\}\subset\mathrm{int}(D)$ 为 n 个标记点。我们考虑保持边界点集 $\partial D$ 不动的同胚的同伦类：
$$
\mathrm{MCG}(D_n,\partial D)
 := \pi_0\bigl(\{f:D\to D\mid f\text{为同胚},\ f(\partial D)=\partial D,\ f(\{p_i\})=\{p_i\}\}\bigr).
$$

另一方面，n 股辫群 $B_n$ 定义为 n 个点在区间 $[0,1]$ 上的编织同伦类。

**定理 1.1（经典同构，见 Birman 等）.** 有自然同构
$$
\Phi: B_n \simeq  \mathrm{MCG}(D_n,\partial D).
$$


接下来考虑 mapping class group 的生成元。

**定义 1.2（Dehn twist）.** 对一条简单闭曲线 $c\subset D_n$，沿着 $c$ 取一个窄带邻域，将其切开，在一侧沿曲线方向旋转 $2\pi$ 后再粘回，得到新的曲面同胚。其同伦类记为 $T_c\in\mathrm{MCG}(D_n,\partial D)$，称为沿 $c$ 的 Dehn twist。

**定理 1.3（Dehn–Lickorish 生成定理的球面多穿孔版本）.** 对于多穿孔圆盘（或等价的球面带标记点），存在有限族简单闭曲线 $\{c_j\}$，使得
$$
\mathrm{MCG}(D_n,\partial D)
 = \langle T_{c_j}^{\pm1}\mid j\in J\rangle
$$
即任意 mapping class 都可写成有限个 Dehn twist 的乘积。

同构 $\Phi$ 与生成定理合并即得：

> **推论 1.4（纯拓扑版主命题）.** 对任意 braid word $\beta\in B_n$，存在简单闭曲线 $c_1,\dots,c_m$ 及符号 $\varepsilon_k=\pm1$，使得
> $$
> \Phi(\beta) = [f] = T_{c_1}^{\varepsilon_1}\cdots T_{c_m}^{\varepsilon_m}\in\mathrm{MCG}(D_n,\partial D).
> $$
> 换言之，每一个 braid 都等价于有限个 Dehn twist 的组合，而每个 Dehn twist 是一个**纯空间几何操作**。

这一步纯拓扑论证不涉及任何量子物理，仅说明“时间中的编织”和“空间中的扭转”在 mapping class group 意义下是同一个对象的不同代表。

---

#### 2. TQFT 层：世界线 braid 表示与 mapping class 表示一致

现在假设给定一个 2+1 维拓扑量子场论（如 Reshetikhin–Turaev 型 TQFT），在一个带 n 个穿孔的空间截面上，其 Hilbert 空间记为
$$
\mathcal H_0(D_n; \{\sigma_i\}),
$$
其中 $\sigma_i$ 是附着在各穿孔上的 anyon 类型（例如 Ising TQFT 中的 $\sigma$）。

在这样的 TQFT 中，存在两种等价的幺正表示构造：

1. **世界线表示**：
    - 把 braid $\beta\in B_n$ 画成 2+1 维时空中的 worldlines；
     - 按照 F,R‑符号的 Reshetikhin–Turaev 规则计算振幅，得到

$$
		 \rho_{\mathrm{WL}}(\beta)\in U(\mathcal H_0).
$$

2. **mapping class 表示**：

    - 把 mapping class $[f]\in\mathrm{MCG}(D_n,\partial D)$ 看作“空间手术”操作：切‑扭‑粘；

     - 由 3D TQFT 的公理，该操作诱导
        $$
         \rho_{\mathrm{MCG}}([f])\in U(\mathcal H_0).
        $$

标准 TQFT 结构（见 Reshetikhin–Turaev、Walker 等）告诉我们：

**定理 2.1（表示一致性）.** 在上述构造下，有
$$
\rho_{\mathrm{WL}}(\beta) = \rho_{\mathrm{MCG}}(\Phi(\beta))
$$
至多差一个整体 U(1) 相位。这意味着：

- “时间中沿 worldlines 编织 anyon” 与
- “在空间截面上做对应的 mapping class（Dehn twist word）”

在 Hilbert 空间上给出的是**同一个幺正表示**。

结合推论 1.4 即得：

> **推论 2.2（TQFT 版主命题）.** 在任意 2+1 维 TQFT 中，对任意 braid word $\beta\in B_n$，存在简单闭曲线 $c_k$ 及符号 $\varepsilon_k$，使得在 $\mathcal H_0$ 上有
> $$
> \rho_{\mathrm{WL}}(\beta) = \rho_{\mathrm{MCG}}\Bigl(\prod_k T_{c_k}^{\varepsilon_k}\Bigr)
> $$
> （至多差一整体相位）。
>
> 也就是说，在 TQFT 的层面 **“时间中的 braid” 与若干“空间 Dehn twist” 的组合本来就是同一类幺正操作。**

---

#### 3. R = exp(iH\_P) / Berry 层：从抽象表示到具体 Hamiltonian 空间操作

我们现在转到具体的晶格拓扑相，实现上述 TQFT 的 Hilbert 空间与幺正表示。假设：

1. 有一个局域自旋/费米体系，其低能有效理论实现了某个拓扑序，对应的零模/编码子空间为 $\mathcal H_0$；

2. 在每条键 $\langle ij\rangle$ 上有一个两体生成元

$$
	 H_P^{(ij)} = \sum_{a,b} c_{ab}^{(ij)}\,\sigma_a\otimes\sigma_b,\qquad R_{ij}=e^{iH_P^{(ij)}};
$$

3. 总哈密顿量 $H(\lambda)$ 由这些 $H_P^{(ij)}(\lambda)$ 组成，参数 $\lambda$ 取值在某个带谱隙的参数空间 $\Omega$ 内。

这里的哈密顿量有两层物理意义：

- “上层”：$H(\lambda)$ 只是一个具体晶格模型，用来实现给定的 2+1 维 TQFT。它的零模/编码子空间 
    $\mathcal H_0(\lambda)$ 携带了前两节讨论的 braid/MCG 幺正表示；
- “下层”：每个局域块 $H_P^{(ij)}$ 描述的是微观上的配对、跳跃、规范约束等局域物理。我们在别处把它写成

$$
	H_P^{(ij)} \simeq H_{\mathrm{quad}}^{(ij)} + H_{\mathrm{int}}^{(ij)} + H_{\mathrm{gauge}}^{(ij)},
$$

并不是为命题额外增加一个假设，而是一种“解读”模型的方式：$H_{\mathrm{quad}}$ 选出 Kitaev/p+ip/honeycomb 那个拓扑通道，$H_{\mathrm{int}}, H_{\mathrm{gauge}}$ 则是偏离理想可积/YBE 点的局域扰动，它们用来定义 $\varepsilon_{\mathrm{YBE}}$、复杂度曲率等。**braid 可以变成空间几何操作** 这一结论只要求 $H(\lambda)$ 局域、带谱隙并实现目标 TQFT，与这三部分的具体形式无关；分解本身只是为了在具体模型里看清“谁在做拓扑事、谁只是扰动”。





前面的 kit-new/kit-new2/kit-new3 已经说明：

- 每个 $H_P^{(ij)}$ 在 JW/Majorana 映射后自然分解为:
    $H_{\mathrm{quad}} + H_{\mathrm{int}} + H_{\mathrm{gauge}}$，前者给出 BdG 拓扑结构，后者描述偏离可积/平坦的扰动；

- YBE / classical YBE 在 $H_P$ 空间中刻画出可积/近平坦子流形，使 Berry 曲率 $F$ 小，保证 holonomy 对路径同伦“刚性”。

对参数空间 $\Omega$ 中的路径 $\gamma$，其 Berry holonomy 定义为
$$
U[\gamma] = \mathcal P\exp\Bigl(-\int_\gamma A\Bigr),\qquad
 A_\mu(X) = iP\,\partial_\mu P\,P,\qquad
$$
（$P$ 为投影到 $\mathcal H_0(\lambda)$）。在配置空间的语境下，$\gamma$ 可以理解为：

- A. 随时间移动缺陷/序参量纹理的路径（worldline 图像）；
- B. 在空间上重写配对图 / Kekulé 畸变 / branch cut 的路径（空间操作）。

第 6 节中已给出如下形式的主不等式（略去常数）：
$$
\bigl\|U_{\mathrm{spatial}} - e^{i\phi}U_{\mathrm{top}}\bigr\|
\;\lesssim\; \varepsilon_{\mathrm{YBE}} + \kappa\,\mathcal A_{\gamma,\gamma'} + \varepsilon_{\mathrm{adiab}} + \varepsilon_{\mathrm{Trotter}},
$$
其中：

- $U_{\mathrm{top}}$ 对应抽象 TQFT 中的 $\rho_{\mathrm{WL}}(\beta)$ 或 $\rho_{\mathrm{MCG}}([f])$；
- $U_{\mathrm{spatial}}$ 是由有限轮局域 $R=e^{iH_P}$ 门、配对图重写或 Kekulé/Dirac 畸变构成的空间操作；
- $\varepsilon_{\mathrm{YBE}}$ 测量所选 $H_P$ 偏离可积/平坦子流形的程度；
- $\kappa$ 为 Berry 曲率上界，$\mathcal A_{\gamma,\gamma'}$ 为两条同伦路径间“填充面积”；
- $\varepsilon_{\mathrm{adiab}}$ 来自时间演化不完全 adiabatic；
- $\varepsilon_{\mathrm{Trotter}}$ 来自用离散 R 门近似连续生成元的误差。

这一不等式的证明骨架依赖于四个引理：

1. **引理 A（平坦/近平坦联络下的路径同伦稳定性）**：在 $\|F\|\le\kappa$ 区域内，同伦路径给出 holonomy 之差 $\sim\kappa\times$ 填充面积；
2. **引理 B（adiabatic 时间演化与 Berry holonomy 等价）**：在带谱隙的情形，真实时间演化与 Berry holonomy 只差 $\varepsilon_{\mathrm{adiab}}$；
3. **引理 C（F,R‑word 与 e^{iH_P} 嵌入一致性）**：在理想可积点上，$e^{iH_P^{(0)}}$ 在 $\mathcal H_0$ 上等同于 TQFT 的 R‑符号，小扰动下拓扑表示稳定；
4. **引理 D（空间重构电路逼近 quasi‑adiabatic 演化）**：quasi‑adiabatic 生成元可由有限轮局域 R 门逼近，误差 $\varepsilon_{\mathrm{Trotter}}$ 可控。

下面将上述四个引理的内容和主不等式的拼接过程展开说明。

##### 3.1 设定与记号

令参数/配置空间为光滑流形 $\Omega$，对每个 $X\in\Omega$ 有带谱隙的哈密顿量 $H(X)$，其零模/编码子空间投影为
$$
P(X): \mathcal H \to \mathcal H_0(X),\qquad P^2=P=P^\dagger.
$$
在 $\Omega$ 上定义 Berry 联络与曲率
$$
 A_\mu(X) = iP\,\partial_\mu P\,P,\qquad
F_{\mu\nu}(X) = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu]
						 = P[\partial_\mu P,\partial_\nu P]P.
$$
对任意光滑路径 $\gamma:[0,1]\to\Omega$，Berry holonomy 定义为
$$
U[\gamma] = \mathcal P\exp\Bigl(-\int_0^1 A_\mu(X(t))\,\dot X^\mu(t)\,\mathrm dt\Bigr).
$$
我们假设在考虑的区域内存在统一界
$$
\|F_{\mu\nu}(X)\| \le \kappa,\qquad \forall X\in \Omega.
$$

设 $\gamma$ 是“理想世界线”路径（anyon 真实编织或 TQFT 抽象路径），$\gamma'$ 是通过空间重构/Kekulé/Dirac 畸变等实现的“空间操作路径”。二者端点相同，且在保持谱隙的区域内同伦。记 $S$ 为两条路径围成的有向曲面，其面积为 $\mathcal A_{\gamma,\gamma'}$。

##### 3.2 引理 A：近平坦联络下的路径同伦稳定性

**引理 A.** 若 $\|F\|\le\kappa$，则存在常数 $C_A$，使得对任意同伦的两条路径 $\gamma,\gamma'$，有
$$
\bigl\|U[\gamma']-U[\gamma]\bigr\| \le C_A\,\kappa\,\mathcal A_{\gamma,\gamma'}.
$$

*证明思路.* 这是非阿贝尔 Stokes 公式的一个定量化版本。

- 由 Berry 联络的 Wilson 线定义和 path‑ordering，可将 $U[\gamma']U[\gamma]^{-1}$ 写成沿着由 $\gamma$、$\gamma'$ 和若干“连接线”组成的闭合回路 $\partial S$ 的 Wilson 环；

- 非阿贝尔 Stokes 定理在此可写作更具体的表面序形式：取基点 $x_0\in S$ 并定义把 $x_0$ 平行传輸到 $x$ 的 Wilson 線

$$
U(x,x_0)=\mathcal P\exp\Bigl(-\int_{x_0}^x A\Bigr).
$$

那麼曲率 $F=dA+A\wedge A$ 的表面序表示為
$$
U[\partial S]=\mathcal P_S\exp\Bigl(-\iint_S U(x,x_0)^{-1}F(x)U(x,x_0)\,dS_x\Bigr),
$$
其中 $\mathcal P_S$ 表示對面元按某種固定順序進行排列（離散化為小面元後按該順序相乘並取極限）。

離散化的直觀近似為：把 $S$ 劃分為小面元 $\{\Delta S_k\}$（中心點 $x_k$），選參考路徑 $p_k$ 把基點運到 $x_k$，則
$$
U[\partial S]\approx\prod_k^{\mathcal P_S}\exp\big(-U(p_k)^{-1}F(x_k)U(p_k)\,\Delta S_k\big).
$$

對表面有序指數做 Dyson 級數展開（令 $\tilde F(x)=U(x,x_0)^{-1}F(x)U(x,x_0)$），前兩項為
$$
U[\partial S]=I -\iint_S \tilde F(x)\,dS_x - \iint_{S_>}\tilde F(x_1)\tilde F(x_2)\,dS_{x_1}dS_{x_2} + O(\mathrm{Area}^3),
$$
其中 $S_>$ 表示按表面序的有序二重積分；高階項包含非對易的嵌套積分與對易子。

取算符範數並假設 $\sup_{x\in S}\|F(x)\|\le\kappa$，可得主階估計
$$
\|U[\partial S]-I\| \le C\,\kappa\,\mathrm{Area}(S) + O(\kappa^2\mathrm{Area}(S)^2),
$$
因此在小面積或曲率受界情況下有
$$
\|U[\partial S]-I\| \lesssim \kappa\,\mathrm{Area}(S).
$$

由於 $U[\partial S]=U[\gamma']U[\gamma]^{-1}$，所以
$$
\|U[\gamma']-U[\gamma]\| = \|U[\partial S]-I\| \lesssim \kappa\,\mathcal A_{\gamma,\gamma'},
$$
從而導出引理 A 中所需的不等式（常數 $C_A$ 可選取為上述 $C$ 的合適放大）。

物理上，这说明在曲率有界的近平坦联络下，同伦的两条路径给出的 holonomy 之差被“曲率上界 × 填充面积”控制，从而确保在 YBE/近平坦子流形附近，空间变形路径与世界线路径的 Berry 相位/矩阵几乎相同。

##### 3.3 引理 B：adiabatic 时间演化与 Berry holonomy

**引理 B.** 设 $H(X(t))$ 是一条光滑、谱隙 $\Delta>0$ 统一有界的时间参数路径，对应的时间演化算符为 $U_{\mathrm{dyn}}(t)$。若演化足够慢（总演化时间 $T$ 足够大），则存在常数 $C_B$ 使得
$$
\bigl\|U_{\mathrm{dyn}}(T)P(X(0)) - U[\gamma]P(X(0))\bigr\| \le C_B\,\varepsilon_{\mathrm{adiab}},
$$
其中 $\varepsilon_{\mathrm{adiab}} \sim \mathcal O(\|\dot H\|/\Delta^2)$ 之类的 adiabatic 小参数。

*证明思路.* 这是 adiabatic 定理的标准形式（Kato、Avron–Seiler–Yaffe 等）：

- 利用瞬时本征态展开，将时间演化分解为 dynamical 相位与几何相位；
- adiabatic 极限下，跃迁到激发态的幅度受 $\|\dot H\|/\Delta^2$ 控制；
- 限制到初态在基态/零模子空间的情形，激发态贡献可忽略，剩下的几何部分正是 Berry holonomy，误差给出上式中的 $\varepsilon_{\mathrm{adiab}}$。

因此，真实的绝热时间演化在编码子空间上与 Berry holonomy 等价，只差一个可控小量。

##### 3.4 引理 C：F,R‑数据与 e^{iH_P} 嵌入的一致性

**引理 C.** 在理想可积点 $H_P^{(0)}$ 上，对每一条键 $(ij)$，存在 $H_P^{(0,ij)}$ 使得
$$
R_{ij}^{(0)} := e^{iH_P^{(0,ij)}}\big|_{\mathcal H_0}
$$
在编码子空间 $\mathcal H_0$ 上实现 TQFT 给出的 R‑矩阵；同时，若 $H_P^{(ij)}=H_P^{(0,ij)}+V^{(ij)}$ 的扰动满足 $\|V^{(ij)}\|\le\epsilon$ 且沿整条路径谱隙保持，则对应的 Berry holonomy 与理想 TQFT 表示相差至多
$$
\varepsilon_{\mathrm{YBE}} \lesssim C_C\,\epsilon
$$
（可能再乘以路径总“长度”）。

*证明思路.*

- 在理想可积点，局域哈密顿量由某个量子群或模块范畴的 R‑矩阵生成，构造上保证 $e^{iH_P^{(0)}}$ 在零模子空间上等同于范畴 R‑符号的表示；
- 沿 YBE/可积流形的小扰动 $V$ 不会闭合谱隙，且其对低能有效理论只产生平滑的、局域的重整化；
- 由稳定性定理（Hastings–Wen topological order stability 等），拓扑量子数与非阿贝尔统计在小局域扰动下不变；
- 因此，$e^{iH_P}$ 在 $\mathcal H_0$ 上的表示与理想 R‑表示之间的差异由扰动范数和路径长度控制，给出 $\varepsilon_{\mathrm{YBE}}$ 项。

物理上，这一项衡量我们在 $H_P$ 空间中“对于理想 YBE/可积点的偏离程度”，也是 classical YBE / 平坦子流形偏离的定量化。



##### 3.5 引理 D：空间重构电路逼近 quasi‑adiabatic 演化

**引理 D.** 设存在一个 quasi‑adiabatic 生成元 $K(t)$，它是由局域算符指数衰减加权积分得到的（Hastings–Wen 方案），使得
$$
U_{\mathrm{QA}} = \mathcal T\exp\Bigl(-i\int_0^T K(t)\,\mathrm dt\Bigr)
$$
在编码子空间上实现沿路径 $\gamma'$ 的 adiabatic/Berry 演化。则对任意给定精度 $\delta>0$，可以构造有限轮局域 R‑门电路 $U_{\mathrm{spatial}}$，使得
$$
\bigl\|U_{\mathrm{spatial}}-U_{\mathrm{QA}}\bigr\| \le \varepsilon_{\mathrm{Trotter}} \le \delta.
$$

*证明思路.*

- quasi‑adiabatic 生成元 $K(t)$ 是局域和的积分，满足 Lieb–Robinson 速度界，因而可在空间上分块、在时间上离散；
- 用 Trotter–Suzuki 分解将时间序列上的指数拆分为各个局域块的指数乘积；
- 每一块指数 $e^{-i h_\alpha \Delta t}$ 可由有限轮局域 $R=e^{iH_P}$ 门近似实现（通过门集完备性或直接识别）；
- 误差由标准 Trotter 估计控制，随时间步长 $\Delta t$ 和分块截断尺度衰减，从而给出 $\varepsilon_{\mathrm{Trotter}}$。

这一步把“连续的 quasi‑adiabatic 演化”真正压缩为“有限深度、有限范围的空间重构电路”，具体形式可以是配对图重写、Kekulé 模式拨动、branch cut/涡旋线段的重新排布等。

##### 3.6 主不等式的拼接

最后说明如何由上述四个引理得到主不等式。分几步进行：

1. **理想拓扑门 vs 理想 Berry holonomy.**

- 在 TQFT 层，有
     $$U_{\mathrm{top}} \simeq \rho_{\mathrm{WL}}(\beta) = \rho_{\mathrm{MCG}}(\Phi(\beta)), $$
    	 其中 $\simeq$ 忽略整体相位；

- 在理想可积点，由引理 C 和 B，可在某条“世界线路径” $\gamma$ 上构造理想哈密顿量 $H^{(0)}(X)$，其 adiabatic/Berry 演化 $U[\gamma]$ 在编码子空间上等同于 $U_{\mathrm{top}}$，至多差一个 $\varepsilon_{\mathrm{YBE}}+\varepsilon_{\mathrm{adiab}}$ 的小量。

2. **世界线路径 vs 空间路径.**

- 在参数/配置空间中构造对应的空间操作路径 $\gamma'$（用 Dehn twist、配对图重写、Kekulé/Dirac 纹理变形等实现），保证 $\gamma$ 与 $\gamma'$ 在保持谱隙的区域内同伦；

- 由引理 A，二者 Berry holonomy 差异满足
    $$
    	 \bigl\|U[\gamma']-U[\gamma]\bigr\| \lesssim \kappa\,\mathcal A_{\gamma,\gamma'}.
    $$

3. **理想 Berry vs quasi‑adiabatic/电路实现.**

- 由引理 B 和 quasi‑adiabatic 构造，可用 $U_{\mathrm{QA}}$ 逼近 $U[\gamma']$，误差吸收到 $\varepsilon_{\mathrm{adiab}}$；

- 再由引理 D，用有限轮 R‑门电路 $U_{\mathrm{spatial}}$ 逼近 $U_{\mathrm{QA}}$，误差记为 $\varepsilon_{\mathrm{Trotter}}$。

4. **三角不等式合并误差.** 将上述三步的偏差用三角不等式串联：

$$
	\begin{aligned}
	 \bigl\|U_{\mathrm{spatial}}-e^{i\phi}U_{\mathrm{top}}\bigr\|
	 &\le \bigl\|U_{\mathrm{spatial}}-U_{\mathrm{QA}}\bigr\|
					 + \bigl\|U_{\mathrm{QA}}-U[\gamma']\bigr\|
					 + \bigl\|U[\gamma']-U[\gamma]\bigr\|\\
	 &\quad + \bigl\|U[\gamma]-e^{i\phi}U_{\mathrm{top}}\bigr\|.
	\end{aligned}

	下面给出每一项的可取上界（常数为模型/几何常数，已在正文各处定义）：

	- Trotter / 电路误差（引理 D）：令 quasi‑adiabatic 生成元被分解为局域项且每项范数上界为 $k_*$，总演化时间 $T$，Trotter 步数 $m$，空间截断半径 $R_k$，则存在常数 $C_D,C'_D,\mu>0$ 使得
		$$
		\bigl\|U_{\mathrm{spatial}}-U_{\mathrm{QA}}\bigr\| \le C_D\frac{T\,k_*^2}{m} + C'_D e^{-\mu R_k}.
		$$

	- 绝热 / quasi‑adiabatic 误差（引理 B）：令 $\Delta$ 为沿路径的最小谱隙，记 $L_0$ 为基态能量相关项尺度，$L_V$ 为扰动导数尺度，则存在常数 $C_B,C''_B$ 使得
		$$
		\bigl\|U_{\mathrm{QA}}-U[\gamma']\bigr\| \le C_B\frac{L_0+L_V}{\Delta^2} + C''_B\frac{\epsilon_{\mathrm{tot}}}{\Delta}.
		$$

	- 曲率 / 同伦面积项（引理 A，Stokes 精确界）：令曲面 $S$ 为 $\gamma,\gamma'$ 所围曲面、面面积为 $\mathcal A_{\gamma,\gamma'}$，并设 $\|F\|\le\kappa$，则存在几何常数 $C_A,C_A'>0$ 使得
		$$
		\bigl\|U[\gamma']-U[\gamma]\bigr\| \le C_A\,\kappa\,\mathcal A_{\gamma,\gamma'}\;\exp\Big(C_A'\,\kappa\,\mathcal A_{\gamma,\gamma'}\Big).
		$$
		在小曲率/小面积近似（$\kappa\,\mathcal A_{\gamma,\gamma'}\ll1$）下，可线性化为 $\lesssim C_A\,\kappa\,\mathcal A_{\gamma,\gamma'}$。

	- YBE / R‑嵌入偏差（引理 C）：令路径长度（或门数）为 $L$，扰动范数尺度 $\epsilon_{\mathrm{tot}}=\sum_{\alpha}\epsilon_{\alpha}$，则存在常数 $C_C,C'_C$ 使得
		$$
		\bigl\|U[\gamma]-e^{i\phi}U_{\mathrm{top}}\bigr\| \le C_C\,L\,\frac{\epsilon_{\mathrm{tot}}}{\Delta} + C'_C\frac{\epsilon_{\mathrm{tot}}^2}{\Delta^2}.
		$$

	把上述不等式代入三角不等式并合并常数，可得到一个带显式常数的主不等式：
	$$
	\begin{aligned}
	\|U_{\mathrm{spatial}}-e^{i\phi}U_{\mathrm{top}}\| &\le C_D\frac{T\,k_*^2}{m} + C'_D e^{-\mu R_k} \\
	&\quad + C_B\frac{L_0+L_V}{\Delta^2} + C''_B\frac{\epsilon_{\mathrm{tot}}}{\Delta} \\
	&\quad + C_A\,\kappa\,\mathcal A_{\gamma,\gamma'}\,e^{C_A'\,\kappa\,\mathcal A_{\gamma,\gamma'}} \\
	&\quad + C_C\,L\,\frac{\epsilon_{\mathrm{tot}}}{\Delta} + C'_C\frac{\epsilon_{\mathrm{tot}}^2}{\Delta^2}.
	\end{aligned}
	$$

	注：如果在应用情形中满足 $\kappa\,\mathcal A_{\gamma,\gamma'}\ll1$ 且 $\epsilon_{\mathrm{tot}}/\Delta\ll1$，可忽略指数项和二阶项，主不等式退化为文中先前的简洁形式（把所有模型常数合并为“≲”符号）。
$$

这就得到开头给出的主不等式形式。

把纯拓扑层的 $B_n\cong\mathrm{MCG}(D_n)$、TQFT 层的 $\rho_{\mathrm{WL}}=\rho_{\mathrm{MCG}}\circ\Phi$ 与上述 R/Berry 层的误差估计结合，即得本文物理版命题：

> 任意 braid $\beta$ 在 TQFT 中对应的拓扑门 $U_{\mathrm{top}}$ 都可以在具体晶格拓扑相中通过一族**几何/空间操作**实现：包括 Dehn twist、branch cut/Kekulé/Dirac 涡旋畸变、配对图重写以及局域 R = exp(iH\_P) 门电路。只要所选 $H_P$ 方向处在 YBE/近平坦子流形附近，且操作过程保持谱隙不闭合，这些空间操作在零模/编码子空间上与“时间中的 braid”给出同一个 SU(2) 共轭类的幺正，误差由上式中的四项可控。

从而，从严格意义上的纯数学层面（mapping class group 与 Dehn twist）到 TQFT 表示，再到具体 Hamiltonian 级别的 R = exp(iH\_P) 实现，我们完成了“braid 可以被变形成空间操作”的命题的完整链路：

- braid = mapping class（定理 1.1）；
- mapping class 由 Dehn twist 生成（定理 1.3）；
- TQFT 中世界线 braid 表示 = mapping class 表示（定理 2.1）；
- R = exp(iH\_P) + Berry 几何为这些 mapping class 提供了具体的空间操作实现，并在近平坦联络下保证 holonomy 的同伦刚性。



