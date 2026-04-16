### R = exp(iH_P)、Kekulé 畸变与“空间诱导拓扑”的一维推导

这一份笔记的目标是：在已经建立好的 R = exp(iH_P) 框架下，

1. 系统推导“哈密顿自然分解 + Kekulé 类空间畸变”的结构；
2. 以一个简化的一维链模型为例，明确地展示：
	- 如何从均匀拓扑相出发，在 H_P 空间中选出一个“拓扑通道”（类似 Dirac 质量 / pairing）；
	- 如何对这个通道施加空间调制（类 Kekulé pattern），把原本的“时间演化路径”重写成一条“空间纹理”；
	- 怎样在这个空间纹理上得到 Berry holonomy / 拓扑门。

这里先做一维推导，二维蜂窝/Kekulé 版本可以在后续文件中平行展开。

---

#### 1. 从 R = exp(iH_P) 到自然分解

我们考虑一条由 N 个自旋‑1/2 组成的开链，对每条最近邻键 $\langle j,j+1\rangle$ 放一个两比特生成元
$$
H_P^{(j,j+1)} = \sum_{a,b\in\{0,x,y,z\}} c_{ab}^{(j)}\,\sigma_a^{(j)}\otimes\sigma_b^{(j+1)},
$$
对应的局部 R‑矩阵为
$$
R_{j,j+1} = e^{iH_P^{(j,j+1)}}.
$$
总哈密顿量取为
$$
H = \sum_{j=1}^{N-1} H_P^{(j,j+1)}.
$$

对每一项做 Jordan–Wigner / Majorana 映射，得到一组 Majorana 算符 $\gamma_{2j-1},\gamma_{2j}$。逐项检查可知：

- 形如 $\sigma^x_j\sigma^x_{j+1},\sigma^y_j\sigma^y_{j+1}$ 的双线性项在 JW 后只生成二次 Majorana $i\gamma\gamma$，对应最近邻 hopping / pairing；
- 形如 $\sigma^z_j,\sigma^z_{j+1}$ 的单点或密度项也只给出二次 Majorana（化学势）；
- 更一般的混合项（$\sigma^z_j\sigma^z_{j+1}$、含交叉 Pauli 的项等）往往生成四 Majorana 乃至更高阶的多体项。

因此总哈密顿量自然分解为
$$
H = H_{\mathrm{quad}} + H_{\mathrm{int}} + H_{\mathrm{gauge}},
$$
其中：

1. $H_{\mathrm{quad}}$：所有二次 Majorana 项之和，对应一个有效的 1D BdG/Kitaev 链；
2. $H_{\mathrm{int}}$：四 Majorana 及以上多体相互作用，描述偏离自由/可积的方向；
3. $H_{\mathrm{gauge}}$：整体能量平移、全局 U(1) 相位等，对 Berry 曲率与拓扑门只是规范自由度。

在原本的均匀拓扑 Kitaev 链中，$H_{\mathrm{int}}$ 可以忽略，$H_{\mathrm{quad}}$ 由常数参数 $(t,\Delta,\mu)$ 给出。我们现在的做法是：

> 在 R = exp(iH_P) 语言下，先在 $H_{\mathrm{quad}}$ 里选出一个“拓扑通道”（某个复序参量），再用空间依赖的模式对这个通道做 Kekulé 式调制，同时保持 $H_{\mathrm{int}}$ 尽量小/处在近平坦层，从而实现“空间纹理诱导的拓扑门”。

---

#### 2. 选出一个“拓扑通道”：复序参量形式

在线性化、小耦合极限下，$H_{\mathrm{quad}}$ 可写成
$$
H_{\mathrm{quad}} = H_0 + \sum_{j} \Bigl[ m_j\,\mathcal O_j + m_j^*\,\mathcal O_j^\dagger\Bigr],
$$
其中：

- $H_0$：作为背景的均匀拓扑 Kitaev 链（比如固定 $t,\Delta$）；
- $\mathcal O_j$：选定的某个 Pauli/Majorana 组合，对应于一个给定方向的 $H_P$ 变化；
- $m_j\in\mathbb C$：在格点 $j$ 上的复序参量（强度+相位），在 R 语言中，它是某几个 $c_{ab}^{(j)}$ 的线性组合。

对应的 R 矩阵可以近似写成
$$
R_{j,j+1} \approx e^{iH_P^{(0)}}\,\exp\bigl(i\,\delta H_P^{(j)}\bigr),
$$
其中 $H_P^{(0)}$ 给出背景拓扑链，$\delta H_P^{(j)}$ 由 $m_j$ 决定，是落在“拓扑通道”上的微扰方向。

我们要求：

1. 这一通道在 $H_P$ 空间中尽量接近 classical YBE/可积子流形（即 classical YBE 残差小，Berry 曲率小）；
2. 它只改变 $H_{\mathrm{quad}}$ 的参数（例如局域 pairing 幅度/相位），不显著增加 $H_{\mathrm{int}}$。

在这种设定下，实空间上的 $\{m_j\}$ 就像 Dirac 理论中的复质量场 $m(x)$，其相位结构 $\theta_j = \arg m_j$ 将直接决定拓扑缺陷与 Berry holonomy。

---

#### 3. 类 Kekulé 调制：从时间路径到空间纹理

传统的“时间演化路径”做法是：选一个全局参数 $m\in\mathbb C$，沿着某条闭合路径 $m(t)$ 缓慢变化，使系统始终保持拓扑谱隙不闭合，最后通过 Berry holonomy 得到一个拓扑门。

类 Kekulé 的思路是：

- 不再主要依赖 $m$ 的时间变化，而是设定一个空间依赖的模式
  $$
  m_j = |m_j|e^{i\theta_j},
  $$
  例如：
  - 简单的 bond‑density wave：$m_j = m_0\cos(Qj)$；
  - 三色 Kekulé pattern：$m_j$ 在周期 3 的子晶格上取不同相位；
  - “畸变涡旋”：在一维环上，把 $\theta_j$ 作为坐标的函数，绕一圈从 0 变到 $2\pi$。
- 随时间进行的是“relabeling / code‑deformation”：通过有限步局域 $e^{iH_P}$ 门，把编码从链的一段搬运到另一段，或绕着某个空间纹理走一圈。

在 Berry 几何的语言下：

- 原来的路径是参数空间中 $m(t)$ 的一条闭合曲线；
- 现在的路径是“序参量配置空间”中 $\{m_j(s)\}$ 的闭合曲线，其中 $s$ 是缓慢演化的参数，物理上通过修改一部分 $H_P^{(j,j+1)}$ 实现。

只要整个过程中谱隙不闭合，且路径处在 $H_P$ 的近平坦区域内，曲率 $F$ 的贡献就可控。按 kit-new2 第 6 节的定理骨架，这两类路径在同伦意义上是等价的，其 Berry holonomy 在零模/编码子空间上给出相同的拓扑门（误差受 $\varepsilon_{\mathrm{YBE}},\kappa$ 等控制）。

---

#### 4. 一维 toy 模型：两段链 + 畸变区域

为了具体化上述推导，我们考虑一个简单的一维 toy 模型：

- 链长 N，分为三段：左段 L、中间畸变区 D、右段 R；
- 左右两段 L,R 处在同一均匀拓扑 Kitaev 相（$H_{\mathrm{quad}}$ 相同，$H_{\mathrm{int}}\approx 0$）；
- 中间 D 区的 $H_P^{(j,j+1)}$ 被赋予一个类 Kekulé 的调制模式 $m_j$，例如相位在 D 区内缓慢从 0 变到 $2\pi$。

构造方式如下：

1. 选择一个背景 $H_P^{(0)}$，其 JW 后给出拓扑点 $(t,\Delta,\mu)$ 的 Kitaev 链；

2. 在 D 区，对某个拓扑通道方向的系数 $c_{ab}^{(j)}$ 加上
$$
	\delta c_{ab}^{(j)} = |\delta c|\cos(Qj+\theta_j),
$$
其中 $\theta_j$ 在 D 区从 0 到 $2\pi$ 单调变化，从而定义出一个“配对序参量涡旋”；

3. 保证调制强度 $|\delta c|$ 小到不会闭合谱隙，也不会显著引入 $H_{\mathrm{int}}$。

这时：

- 在左/右两段，零模/编码子空间与原始拓扑链相同；
- 中间畸变区 D 提供了一个“拓扑纹理”，在其两端可以出现局域化的 near‑zero modes；
- 沿着畸变区域 adiabatic 地“滑动编码”（例如通过 code‑deformation 或空间重构电路），等效于在序参量配置空间中绕着一次完整的相位 winding，给出一个非平凡的 Berry holonomy。

从 R 的角度看：

- 不同区域的键对应不同的 $R_{j,j+1} = e^{iH_P^{(j,j+1)}}$；
- 在 D 区内，$R_{j,j+1}$ 沿着某一“拓扑通道方向”在 SU(4) 中缓慢旋转一圈；
- 通过合适选择，这条旋转路径可以近似满足 classical YBE 条件，使得三体离散曲率小，从而保证 holonomy 主要由“绕圈的拓扑信息”控制。

---

#### 5. 下一步：具体数值模型与 Berry holonomy 计算

要把这个 toy 模型真正跑数值，步骤可以是：

1. 选一个具体的 $H_P^{(0)}$（例如之前 4‑Majorana half‑twist 的放大版本），写出对应的 Kitaev 链参数 $(t,\Delta,\mu)$；
2. 在中间 D 区选择一个 Pauli 通道 $\mathcal O$，构造 $m_j$ 的空间纹理（比如线性相位绕圈），并据此修改 $H_P^{(j,j+1)}$ 中相应的 $c_{ab}^{(j)}$；
3. 对这个一维模型构造出两个“编码端点”（比如两端的 Majorana 零模），
	- 数值对角化得到零模子空间；
	- 沿着一个离散化的“空间重构路径”（一串局域 R 门）计算 Berry holonomy；
4. 比较这个 holonomy 与目标拓扑门（例如 half twist 的 SU(2) 矩阵），检验它们是否在 SU(2) 中共轭，并测量随调制强度/路径形状变化的稳定性。

这样我们就用一个非常具体的一维模型证明了：

- 在 R = exp(iH_P) 框架下，哈密顿量可以自然分解出一个“拓扑通道”；
- 对这一通道施加类 Kekulé 的空间调制，本质上是在参数/序参量配置空间中画一条有 winding 的路径；
- 通过 Berry holonomy，可以把这条路径对应的“空间纹理”直接翻译成一个拓扑门，从而把“时间演化 → 空间变化”的抽象思想具体化为可运算、可数值验证的模型。

