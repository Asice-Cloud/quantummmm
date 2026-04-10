## 4.2 Kitaev honeycomb 在非阿贝尔相的几何/拓扑图景

这一节把前面在 1D/Kitaev 链和 2D p+ip 超导上建立的“上层空间 = 配置/模空间 + Berry 联络 + Dehn twist = $F^{-1}R^2F$”图景，翻译到更接近物理实现的 Kitaev honeycomb 模型上。在非阿贝尔相下，honeycomb 的低能有效理论正好是一个带 $\mathbb Z_2$ 规范场的 Majorana 拓扑超导；通量缺陷（vison）携带局域 Majorana 零模，其编织统计与 Ising 任何子相同。

我们关心的是：

- 如何从自旋哈密顿量出发，得到“Majorana + $\mathbb Z_2$ 背景场”的 BdG‑样描述；
- 在存在若干通量缺陷的情形下，如何定义缺陷配置空间上的 Hilbert 丛与 Berry 联络，并把涡旋/vison 的编织视作 mapping class 群元在零模子空间上的 holonomy；
- 如何利用前面 p+ip 的数值结果，把 Dehn twist = $F^{-1}R^2F$ 的解析结构当作对 honeycomb 的“低能约束”，指导对更真实格点模型的路径与参数选择。

### 4.2.1 自旋模型到 Majorana + $\mathbb Z_2$ 背景场

Kitaev honeycomb 模型的自旋哈密顿量为

$$
H_{\text{KH}} = -\sum_{\gamma\in\{x,y,z\}} J_\gamma \sum_{\langle i j\rangle_\gamma} S_i^\gamma S_j^\gamma,
$$

其中 $\langle i j\rangle_\gamma$ 表示蜂巢格子上的 $\gamma$ 型键。为进入 Majorana 表示，我们对每个自旋引入四个 Majorana 算符 $c_i, b_i^x, b_i^y, b_i^z$，满足

$$
\{c_i,c_j\}=2\delta_{ij},\qquad \{b_i^\alpha, b_j^\beta\}=2\delta_{ij}\delta^{\alpha\beta},\qquad \{c_i,b_j^\alpha\}=0,
$$

并用约束 $D_i = b_i^x b_i^y b_i^z c_i = +1$ 将物理 Hilbert 空间选出来。在这个表示下，自旋算符可以写为

$$
S_i^\alpha = \frac{i}{2} b_i^\alpha c_i.
$$

代入 $H_{\text{KH}}$，每条 $\gamma$ 键上的耦合项变为

$$
S_i^\gamma S_j^\gamma = -\frac{1}{4} (b_i^\gamma c_i)(b_j^\gamma c_j)
 = -\frac{i}{4} (i b_i^\gamma b_j^\gamma) c_i c_j.
$$

于是整个哈密顿量可写成

$$
H_{\text{KH}} = \frac{i}{4} \sum_{\langle i j\rangle_\gamma} J_\gamma u_{ij} c_i c_j,
\qquad u_{ij} = i b_i^\gamma b_j^\gamma,
$$

其中 $u_{ij} = \pm 1$ 是定义在键上的 $\mathbb Z_2$ 规范变量，满足

$$
u_{ij} = -u_{ji},\qquad [u_{ij}, H_{\text{KH}}]=0,
$$

因此在每个 $u_{ij}$ 固定的扇区内，$H_{\text{KH}}$ 退化为一个自由的 c‑Majorana 费米子在静态 $\mathbb Z_2$ 背景场中的二次哈密顿量。对每个六角 plaquette $p$，通量算符为

$$
W_p = \prod_{\langle ij\rangle\in p} u_{ij} \in \{+1,-1\},
$$

满足 $[W_p,H_{\text{KH}}]=0$。$W_p=-1$ 的 plaquette 即对应一个 $\mathbb Z_2$ 通量缺陷（vison），在非阿贝尔相中会携带 Majorana 零模。

### 4.2.2 非阿贝尔相与有效 p+ip 描述

在零磁场下，Kitaev honeycomb 的各相随 $(J_x,J_y,J_z)$ 展开，其中“B 相”是一个 gapless Dirac‑like 相；在加入弱磁场或者合适的三自旋相互作用后，Dirac 点打开能隙，系统进入一个拓扑非平凡的 gapped 相。这个 gapped 相的 c‑Majorana 模式在低能/长波极限下可以用一个带 Chern 数 $\pm1$ 的 2D p+ip 拓扑超导来描述：

- 在某个均匀 $u_{ij}$ 配置下（通常取无通量扇区 $W_p=+1$），c‑Majorana 的能带结构在布里渊区的一些点附近呈现出类似 p+ip BdG 的拓扑结构；
- 加入磁场 $h$ 或三自旋项 $K$ 后，有效哈密顿量中出现 next‑nearest‑neighbor 的 $ic_ic_j$ 项，这可以看成是在 c‑Majorana 上诱导出一种 p+ip 形式的配对；
- 将这些长程 Majorana 跳跃/配对在连续极限下重写，就得到一个形如

	$$
	H_{\text{eff}} \sim \int d^2\mathbf r\, c(\mathbf r)\, \bigl( -i v (\partial_x \gamma_y - \partial_y \gamma_x) + m \bigr)\, c(\mathbf r),
	$$

	的拓扑超导场论，其拓扑性质与 2D p+ip 模型相同。

更重要的是，在这个 gapped 非阿贝尔相内，如果在某些 plaquette 上插入 $W_p=-1$ 的通量 vison，那么在每个 vison 附近会出现局域的 c‑Majorana 零模；多个 vison 的零模空间维度与 Ising TQFT 中多个 $\sigma$ 任何子对应的 Hilbert 空间维度一致。这一点在文献中已经有详尽分析，我们在此直接采用这一标准图像：

> 非阿贝尔相下的 Kitaev honeycomb，可以在低能上等价地看作一个 2D p+ip 拓扑超导，其上的通量缺陷 vison 相当于 p+ip 中的涡旋，每个 vison 携带一支 Majorana 零模，多 vison 的零模 Hilbert 空间与 Ising 任何子的 Hilbert 空间幺正同构。

这为我们把前面在 2D p+ip 模型中数值验证过的“涡旋绕行 Berry holonomy = Ising Dehn twist = $F^{-1}R^2F$”图像，直接移植到 honeycomb 提供了理论基础。

### 4.2.3 缺陷配置空间上的 Hilbert 丛与 Berry 联络

在 honeycomb 模型的 Majorana+$\mathbb Z_2$ 描述中，通量缺陷的几何位置可以通过其所在 plaquette 的中心位置来粗略标记。假设我们在某个 genus $g$ 的曲面 $\Sigma$ 上定义了一个有限蜂巢 patch，并在若干 plaquette 上插入 $n$ 个 vison，则可以定义通量配置空间

$$
\mathcal C_{\text{vison}} = \mathrm{Conf}_n(\Sigma) = \{(p_1,\dots,p_n)\in \Sigma^n \mid p_a\neq p_b\}/S_n,
$$

其中 $p_a$ 是对应 vison 的“几何位置”（plaquette 中心），$S_n$ 把不可区分的 vison 标签模掉。对每个 $X=(p_1,\dots,p_n)\in\mathcal C_{\text{vison}}$，在给定的 $u_{ij}$ 背景场和局部扰动下，我们可以视为有一个 c‑Majorana BdG‑样哈密顿量 $H_{\text{vison}}(X)$，其最低能的若干本征态张成一个有限维零模/编码子空间 $\mathcal H_{\text{KH}}(\Sigma,X)$。

把所有这些零模子空间在配置空间上拼接，我们得到一个 Hilbert 向量丛

$$
\mathcal E_{\text{KH}} \to \mathcal C_{\text{vison}},
$$

其纤维为 $\mathcal H_{\text{KH}}(\Sigma,X)$，在非阿贝尔相内部，这个丛在拓扑上与 Ising TQFT 的丛

$$
\mathcal E_{\mathrm{Ising}} \to \mathcal C_{\text{vison}},\qquad \mathcal E_{\mathrm{Ising}}|_X \cong \mathcal H_{\mathrm{Ising}}(\Sigma,X)
$$

同构。

选定一个规范，即对每个 $X$ 选一组正交归一基 $\{\psi_i(X)\}$，我们可以像在 2D p+ip 和 4‑Majorana 模型中那样定义 Berry 联络

$$
 A_{ij}(X) = \langle\psi_i(X)|\,d_X\psi_j(X)\rangle,
$$

以及其曲率

$$
 F = dA + A\wedge A.
$$

在缺陷配置空间中选一条路径 $\gamma:[0,1]\to\mathcal C_{\text{vison}}$，可以是：

- 一条简单 braid：某两个 vison 彼此交换；
- 一条 Dehn twist 路径：某组 vison 绕着 handle 或其他 puncture 走一圈，对应 mapping class 群中的一个 $T_\gamma$。

沿着这条路径，绝热移动通量缺陷，构造出参数依赖的哈密顿量 $H_{\text{vison}}(X(t))$，对每个 $t$ 提取零模子空间基底，并用相邻基底的重叠矩阵经极分解投影到最近的幺正矩阵，如同我们在 [verify/run_pip_vortex_berry.py](verify/run_pip_vortex_berry.py) 里所做的那样。这一过程给出 Berry holonomy

$$
 U_{\text{KH}}[\gamma] = \mathcal P\exp\Bigl(-\int_\gamma A\Bigr),
$$

它在 $\mathcal H_{\text{KH}}(\Sigma,X)$ 上的作用，在非阿贝尔相内部应当与 Ising TQFT 提供的映射

$$
 \rho_{\mathrm{Ising}}: B_n,\ \mathrm{MCG}(\Sigma,X) \to U\bigl(\mathcal H_{\mathrm{Ising}}(\Sigma,X)\bigr)
$$

幺正同构（允许一个 $X$ 依赖的基变换）。

### 4.2.4 Dehn twist 与 $F^{-1}R^2F$ 在 honeycomb 中的意义

在抽象的 Ising TQFT 中，我们已经知道：

- 两个 $\sigma$ 粒子的 braid 与 half twist 由 $R^{\sigma\sigma}_c$ 的本征值给出；
- 绕包住某个内部融合通道的 Dehn twist 在适当的融合基底中表现为 $(R^{\sigma\sigma})^2$ 的对角矩阵；
- 一般曲线上的 Dehn twist 可以用有限次 F‑move 把 pants 分解重排成“围绕某条内部边”的形式，然后在该边施加 $R^2$ 或 topological spin，再 F‑move 回原基底，即

	$$
	U_{\text{Dehn}}(\gamma) \simeq F^{-1} R^2 F.
	$$

在前面的 1D/4‑Majorana toy 模型中，我们通过半解析半数值的方式，已经显示了用半扭门 $\exp(\frac{\pi}{4}\gamma_i\gamma_j)$ 和有限维基变换可以重建这一结构，并与 Berry holonomy 做 SU(2) 共轭比较；而在 2D p+ip 模型中，我们则通过 [verify/run_pip_vortex_berry.py](verify/run_pip_vortex_berry.py) 和 [verify/run_pip_vortex_scan.py](verify/run_pip_vortex_scan.py) 显示了实际涡旋绕行路径的 Berry holonomy 与 Ising 的 $(R^{\sigma\sigma})^2$ 在 SU(2) 意义下具有高 fidelity，且在一定的 $\mu$ 区间内保持近乎不变。

对 honeycomb 而言，这些结果可以理解为对下面这句话的支持：

> 在非阿贝尔相下，若我们选取合适的参数 $(J_x,J_y,J_z,h,\dots)$ 和一条在通量配置空间中对应 Dehn twist 的路径 $\gamma$，那么在低能零模子空间中，Kitaev honeycomb 的 Berry holonomy $U_{\text{KH}}[T_\gamma]$ 将与 Ising TQFT 给出的 $U_{\mathrm{Ising}}[T_\gamma]\simeq F^{-1}R^2F$ 幺正同构，只差一个局部的编码基变换和整体相位。

这为我们后续在 honeycomb 上设计具体的“通量编织/Dehn twist 实验”提供了两个重要的理论约束：

- **路径选择约束**：路径 $\gamma$ 在几何上必须对应某个 mapping class 群生成元（如 braid 或 Dehn twist），而不是任意杂乱无章的循环；在 p+ip 中我们已经看到，路径的几何形状（是否真的绕 handle 一圈）会强烈影响 $U_{\text{holo}}$ 是否接近 $(R^{\sigma\sigma})^2$；
- **参数区间约束**：即便在 nominal 的非阿贝尔相内部，$U_{\text{Dehn}}$ 的 fidelity 也会随参数（如化学势 $\mu$、涡旋间距等）变化，在某些“平坦区域”内 Berry holonomy 几乎不变，在其它区域则迅速偏离。p+ip 上的 $F_{\mathrm{Dehn}}(\mu)$ 曲线正是这种现象的一个连续模型示范，我们期望 honeycomb 在 $(J_x,J_y,J_z,h,\dots)$ 空间中也有类似的“Dehn twist 平坦子流形”。

因此，从 1D R(a,b,c,d) 模型到 2D p+ip 再到 honeycomb，我们得到了一条清晰的理论链条：

$$
	ext{自旋哈密顿量} \Rightarrow \text{Majorana}+\mathbb Z_2 \Rightarrow \text{有效 p+ip} \Rightarrow \text{Ising TQFT}(F,R) \Rightarrow \text{Berry holonomy on }\mathcal C_{\text{vison}},
$$

其中每一层都可以用 Berry 联络与曲率的几何语言来刻画，并在适当极限下导出 $U_{\text{Dehn}}\simeq F^{-1}R^2F$ 的结构性结果。后续如果在 honeycomb 上实现类似于 [verify/run_pip_vortex_berry.py](verify/run_pip_vortex_berry.py) 的数值实验，就可以把这条链条在一个更真实的格点模型里闭合起来。

