## 4. 从 1D 到 2D：2D $p+ip$ 拓扑超导与涡旋零模

这一章我们把前面在 1D/Kitaev 链和 4‑Majorana toy 模型上走通的那条链条

$$
R(a,b,c,d) \Rightarrow \text{BdG/Kitaev} \Rightarrow \text{Berry 联络}~A,~F \Rightarrow \text{YBE 子流形}~\mathcal M_R^{(\mathrm{YBE})} \Rightarrow \text{braid/Dehn twist 拓扑门} \Rightarrow \text{LQC+permutation 电路}
$$

尝试推广到一个真正的 2D 平台：2D $p+ip$ 拓扑超导，在其中插入磁通涡旋，使每个涡旋携带一个 Majorana 零模，多涡旋配置给出 Ising 任何子 Hilbert 空间。

这里的主线是：

- 选定 2D 微观模型：连续或格点形式的 2D $p+ip$ 拓扑超导；
- 明确“缺陷”是什么：超导序参量中的磁通涡旋（vortex），以及可能的边界/缺口；
- 在 BdG/Majorana 表示下，构造含有若干涡旋的 Hamiltonian，识别其零模子空间；
- 在涡旋配置空间 $\mathcal C=\mathrm{Conf}_n(\Sigma)$ 上引入 Berry 联络 $A$，计算 braid/Dehn twist 的 holonomy，并与 Ising TQFT 的 $F,R$‑数据对比；
- 探索在这类 2D 模型上是否存在“近似 YBE/平坦”的参数子流形，以及对应的 constant‑depth LQC+permutation 实现。

下面从一个最简单的 2D $p+ip$ BdG 模型与涡旋零模开始。

### 4.1 2D $p+ip$ BdG 模型与涡旋零模

在连续极限中，可以用一个 spinless 2D $p+ip$ 拓扑超导来捕捉 Ising‑like 非阿贝尔统计的本质。一个简化的 BdG 哈密顿量是

$$
H = \int d^2\mathbf r\, \Psi^\dagger(\mathbf r)
\begin{pmatrix}
\frac{\mathbf p^2}{2m} - \mu & \Delta(\mathbf r) (p_x + i p_y) \\
\Delta^*(\mathbf r) (p_x - i p_y) & -\frac{\mathbf p^2}{2m} + \mu
\end{pmatrix}
\Psi(\mathbf r),
$$

其中 $\mu>0$、$\Delta(\mathbf r)=|\Delta(\mathbf r)|e^{i\theta(\mathbf r)}$。当 $|\Delta|$ 非零且 $\mu$ 落在拓扑相区时，这个系统支持 Chern 数为 $\pm1$ 的拓扑超导态。

在这个背景上，插入一个磁通涡旋意味着让序参量的相位 $\theta(\mathbf r)$ 绕某一点 $\mathbf R_a$ 绕一圈时发生 $2\pi$ 绕转，如

$$
\oint_{|\mathbf r-\mathbf R_a|=R} \nabla \theta\cdot d\mathbf l = 2\pi.
$$

对这样的涡旋配置，BdG 方程在每个涡旋附近支持一个指数局域的零能解 $\gamma_a$，满足

$$
\{\gamma_a,\gamma_b\}=2\delta_{ab},\qquad \gamma_a^\dagger = \gamma_a,
$$

即每个涡旋携带一个 Majorana 零模。对 $2N$ 个涡旋，零模子空间张成维数为 $2^N$ 的实 Clifford 代数，在适当的总费米数约束下，给出 $2^{N-1}$ 维的编码 Hilbert 空间，这正是 Ising 任何子 $\sigma$ 粒子在 $2N$ 个 puncture 上的 Hilbert 空间维度。

从我们的视角看，这一步的作用是：给定一个固定的曲面 $\Sigma$ 和若干涡旋位置 $X=\{\mathbf R_a\}$，2D $p+ip$ 模型为每个 $(\Sigma,X)$ 提供一个零模 Hilbert 空间 $\mathcal H_{p+ip}(\Sigma,X)$，它在拓扑相区内可以与 Ising TQFT 的 $\mathcal H_{\mathrm{Ising}}(\Sigma,X)$ 幺正同构。

### 4.2 涡旋配置空间与 Berry 联络

接下来，我们把涡旋的位置当作参数，考虑配置空间

$$
\mathcal C = \mathrm{Conf}_n(\Sigma) = \{(\mathbf R_1,\dots,\mathbf R_n)\in \Sigma^n\mid \mathbf R_a\neq\mathbf R_b\}/S_n,
$$

其中 $S_n$ 把不可区分的涡旋标号模掉。对每个 $X\in\mathcal C$，我们有一个零模 Hilbert 空间 $\mathcal H_{p+ip}(\Sigma,X)$，把它们拼在一起得到一个 Hilbert 丛

$$
\mathcal E_{p+ip} \to \mathcal C,
$$

在这个丛上选定一个光滑的规范（即对每个 $X$ 选一组正交归一基 $\{\psi_i(X)\}$），定义 Berry 联络

$$
A_{ij}(X) = \langle \psi_i(X) | d_X \psi_j(X)\rangle,
$$

以及对应的曲率

$$
F = dA + A\wedge A.
$$

对任意一条涡旋配置路径 $\gamma:[0,1]\to\mathcal C$，其 Berry holonomy 为

$$
U_{p+ip}[\gamma] = \mathcal P\exp\Bigl(-\int_\gamma A\Bigr),
$$

这就是在这个 2D $p+ip$ 模型中，沿着给定涡旋编织/Dehn twist 路径所实现的拓扑门。在拓扑相区、且路径保持能隙不闭合时，$U_{p+ip}[\gamma]$ 只依赖于 $[\gamma]\in\pi_1(\mathcal C)$ 的同伦类，并与 Ising TQFT 提供的辫子群/映射类群表示

$$
\rho_{\mathrm{Ising}}: B_n,\ \mathrm{MCG}(\Sigma,X) \to U\bigl(\mathcal H_{\mathrm{Ising}}(\Sigma,X)\bigr)
$$

幺正同构。这一点在文献中已经有半解析、半数值的论证，我们的目标是：在一个有限尺寸的格点 BdG 模型上，显式地重复我们在 4‑Majorana toy 模型中做过的步骤 —— 沿具体路径 $\gamma$ 构造 $H(\lambda)$，离散化参数 $\lambda$，在零模子空间上做平行转运，得到一个具体矩阵近似 $U_{p+ip}[\gamma]$，再与 Ising TQFT 的 $F,R$‑矩阵做数值共轭比较。

### 4.3 从 2D Berry holonomy 到 Dehn twist 与 YBE

在 2D $p+ip$ 背景中，一条典型的 Dehn twist 路径可以是：

- 在 genus\,$g$ 的曲面上沿某条 handle 的非平凡闭合曲线，拖着一组涡旋绕行一圈；
- 或者让一对涡旋绕包住其他涡旋或 puncture 的圈，在配置空间中对应 mapping class 群的一个生成元 $T_\gamma$。

对这样的 $\gamma$，Berry holonomy $U_{p+ip}[T_\gamma]$ 应当与 Ising TQFT 的

$$
U_{\mathrm{Ising}}[T_\gamma] = \rho_{\mathrm{Ising}}(T_\gamma)
$$

在某个编码基底下共轭。抽象上，后者可以用 $F,R$‑矩阵写成

$$
U_{\mathrm{Ising}}[T_\gamma] \simeq F^{-1} R^2 F,
$$

而在我们的 1D 分析中，这一结构已经通过 4‑Majorana toy 模型得到了数值验证。2D 情况下，思路完全类似：

- 用 BdG 数值计算得到若干代表性路径（简单 braid、绕 handle 的 Dehn twist 等）的 $U_{p+ip}[\gamma]$；
- 在同一 Hilbert 空间中写出 Ising TQFT 的 $U_{\mathrm{Ising}}[\gamma]$（通过 $F,R$ 数据）；
- 寻找幺正矩阵 $V$，比较 $V^\dagger U_{p+ip}[\gamma] V$ 与 $U_{\mathrm{Ising}}[\gamma]$ 的差异，从而在 2D 连续模型中具体验证“Berry holonomy = F^{-1}R^2F” 的图景。

在这个基础上，我们可以进一步把 1D 中讨论的“YBE 子流形/平坦区域”和“曲率–复杂度”的思路推广到 2D：在 2D $p+ip$ 模型的参数空间（如 $\mu,|\Delta|$、外加势形状等）和几何空间（涡旋路径的形状、曲面模空间）上寻找那些使得 Berry 曲率 $F$ 在某些代表性的 2‑胞（辫子关系 2‑胞、Dehn twist 关系 2‑胞）上接近零的区域，把它们看作“2D 版的近似 YBE/平坦子流形”；然后在固定电路深度的 LQC+permutation ansatz 下，考察这些区域内外的拟合 fidelity 差异，得到一个真正意义上的“2D 上层空间几何 vs 逻辑电路复杂度”的对应关系。


