## 4.3 从 R 到 Majorana 与 Berry 几何：一套完整流程

本节把前面在 1D R(a,b,c,d) 模型、2D p+ip 超导与 honeycomb 上零散出现的思想，整理成一套统一的、可重复使用的“从代数 R 到 Majorana 几何与参数扫描”的完整步骤。今后无论在 1D 链、2D p+ip，还是在蜂窝格 Kitaev 模型上，只要满足相同的物理假设，都可以沿用这一流程来研究 Majorana 模式及其拓扑门。

### 4.3.1 步骤 1：选定平台与参数路径

- 选定具体物理平台：
	- 1D：Kitaev 链 / 4-Majorana toy 模型（参考 kit3-5 与相应 verify 脚本）。
	- 2D：方格 p+ip 拓扑超导 + 磁通涡旋（参考 kit4-1 与 verify/run_pip_vortex_berry.py）。
	- honeycomb：Kitaev 蜂窝格自旋模型，在非阿贝尔相中用 Majorana+Z2 背景场与等效 p+ip 描述（参考 kit4-2）。
- 在所选平台上，指定一个参考参数点 $p_0$ 及其附近的一条参数路径：
	- 对 1D R(a,b,c,d) 模型：选定 $p_0=(a_0,b_0,c_0,d_0)$，例如 Ising 点附近 $(0,1,0,0)$，并选择一条简单的扫描线（如固定 $a,d$，沿 $b$ 或 $c$ 扫描）。
	- 对 2D p+ip：选定一组 $(t_0,\Delta_0,\mu_0)$ 或等价的 $(a_0,b_0,c_0,d_0)$，并指定一条或若干条参数路径（例如固定 $t_0,\Delta_0$，只扫描化学势 $\mu$，参考 kit4-1 中的 $F_{\mathrm{Dehn}}(\mu)$ 扫描）。
	- 对 honeycomb：在 $(J_x,J_y,J_z,h,\dots)$ 空间中选取处在非阿贝尔相内的参考点，并设计一条不闭合能隙的参数路径。

这一“步骤 1”明确了我们要在哪个平台、哪条参数路径上研究 Majorana 模式的稳定性与几何性质。

### 4.3.2 步骤 2：从 R 到 BdG / Majorana 哈密顿量

- 在 1D 情形，使用 [R_to_Kitaev.md](R_to_Kitaev.md) 中的映射
	$$
	t=b+c,\qquad \Delta=b-c,\qquad \mu\simeq 4d+\mu_0,
	$$
	把每个 $(a,b,c,d)$ 点转成 Kitaev 链的 $(t,\Delta,\mu)$，构造对应的 BdG 哈密顿量 $H_p$。
- 在 2D p+ip 平台，直接在格点上写 BdG 哈密顿量（参考 verify/run_pip_vortex_berry.py），必要时仍可通过 $(a,b,c,d)$→$(t,\Delta,\mu)$ 把 YBE 解家族的代数结构移植到 2D 模型的耦合常数中。
- 在 honeycomb 模型中，先用 kit4-2 所述自旋→Majorana+Z2 背景场→有效 p+ip 的低能映射，把 $(J_x,J_y,J_z,h,\dots)$ 转成等效 p+ip 的参数，再在该等效 p+ip 模型上构造 BdG 哈密顿量 $H_p(x)$，其中 $x$ 标记 vison/涡旋等缺陷的位置。

通过这一步，我们把“形式上的常数 R(a,b,c,d)”转化为可数值求解的 BdG/Majorana 哈密顿量族 $H_p(x)$，为后续的谱分析与 Berry 几何奠定基础。

### 4.3.3 步骤 3：谱与拓扑不变量（静态拓扑判据）

- 对参数路径 $p(s)$ 上的离散点 $\{p_k\}$，构造相应的 BdG 哈密顿量 $H_{p_k}$，对角化得到能谱：
	- 记录基态与第一激发态之间的最小能隙 $\Delta E(p_k)$，用于判定是否发生相变。
	- 在 1D/Kitaev 链情形，计算 Pfaffian 符号或拓扑 winding 数，判定拓扑相是否改变（参考 kit3-5-3 对 $(a,b,c,d)$ 的判据）。
- 若在一段参数区间内，能隙始终保持打开且拓扑不变量不变，则该区间内 Majorana 边界/缺陷零模的拓扑保护在静态上是稳定的；反之，若能隙闭合或不变量跳变，则该点标记为相变点或拓扑不稳定区。

这一层对应“静态拓扑相”的分析：不涉及 Berry 曲率，只关心谱间隙和拓扑不变量随参数的变化。

### 4.3.4 步骤 4：配置空间上的 Berry holonomy

- 固定一个在给定平台中具有拓扑意义的缺陷路径 $\gamma$：
	- 1D/4-Majorana toy：某一对端点的交换或 Dehn twist 路径（参考 kit3-5 及相关 verify 脚本）。
	- 2D p+ip：一个涡旋绕另一涡旋的闭合路径（矩形或圆形回路，参考 verify/run_pip_vortex_berry.py）。
	- honeycomb：某个 vison 在蜂窝格上绕另一 vison 或绕周期圈的路径（参考 kit4-2 的配置空间讨论）。
- 把这条路径离散化为 $x_0\to x_1\to\cdots\to x_N$，在每个点 $x_k$ 上：
	- 提取零模子空间 $\mathcal H_{x_k,p}$，选取一组正交归一基底 $\{\psi_a(x_k)\}$；
	- 计算相邻点间的重叠矩阵 $M_k=\langle\psi_a(x_{k+1})|\psi_b(x_k)\rangle$，通过极分解提取其酉部分 $U_k$；
	- 沿路径有序相乘得到 Berry holonomy
		$$
		U_{\text{Berry}}[\gamma;p]=\mathcal P\prod_{k} U_k.
		$$

在数值实现上，这一步与我们在 verify/run_pip_vortex_berry.py、verify/run_dehn_twist_micro_vs_berry.py 中使用的“重叠矩阵 + 极分解”方法完全一致，只是平台和几何路径不同。

### 4.3.5 步骤 5：与 TQFT F,R / Dehn twist 的对齐

- 在零模子空间的逻辑基底中，把 $U_{\text{Berry}}[\gamma;p]$ 规约到 SU(2) 或更高维的特殊酉群（去掉整体 U(1) 相位），得到“纯逻辑门”部分。
- 在 Ising TQFT 或更一般的模范畴框架中，写出相应的 $F$、$R$ 矩阵，例如 Ising 情形下的 $R^{\sigma\sigma}$ 与 $F^{\sigma\sigma\sigma}_{\sigma}$，并构造形式上的 Dehn twist
	$$
	U_{\text{Dehn}}^{\text{TQFT}}\simeq F^{-1}R^2F.
	$$
- 在同一逻辑基底下，用 SU(2)（或 SU($n$)）重合度（fidelity）比较数值 Berry holonomy 与理论 Dehn twist：
	- 若重合度 $F_{\text{Dehn}}\approx1$，则说明在该参数点及路径上，微观 Majorana 模型实现的 Berry holonomy 与理想 TQFT 层的 Dehn twist 在逻辑子空间上几乎等价；
	- 若 $F_{\text{Dehn}}$ 明显偏离 1，则该偏离量刻画了“几何/微扰对理想拓扑门的修正”。

这一层把“上层的抽象 F,R 结构”和“下层的 BdG/Majorana 数值结果”精确地扣在一起，回答“这条路径上真正实现的门在多大程度上等价于某个理想的 mapping class 群元素（如 Dehn twist）”。

### 4.3.6 步骤 6：线性稳定性与参数扫描

- 线性稳定性分析：
	- 在参考点 $p_0$ 附近取小微扰 $p=p_0+\delta p$，沿一个或多个参数方向做一阶扫描；
	- 对每个 $p$，重复步骤 3–5，记录能隙 $\Delta E(p)$、拓扑不变量以及 $F_{\text{Dehn}}(p)$ 的变化；
	- 由此得到对 $(a,b,c,d)$ 或 $(t,\Delta,\mu)$ 等参数的“线性稳定性判据”：哪些方向的微扰对 Majorana 拓扑门影响最小，哪些方向最为敏感。
- 系统参数扫描：
	- 在选定的参数子空间内（例如 1D 的 $(a,b,c,d)$ 面、2D p+ip 的 $(\mu,\Delta)$ 平面、honeycomb 的 $(J_x,J_y,J_z)$ 面）进行更大范围的扫描；
	- 对每个网格点执行步骤 3–5，绘制：
		- 能隙热图与拓扑不变量相图；
		- $F_{\text{Dehn}}(p)$ 或其他 gate fidelity 的图像，标出“Dehn twist 平坦区 vs 非拓扑区”。

在 2D p+ip 平台上，我们已经在 verify/run_pip_vortex_scan.py 中实现了沿 $\mu$ 轴的 $F_{\text{Dehn}}(\mu)$ 扫描，并在 kit4-1 中讨论了结果；同样的思想可以推广到更高维的参数空间以及 honeycomb 模型的参数族。

### 4.3.7 步骤 7：几何–代数综合与对 Majorana 的应用

- 将上述数值结果与代数先验（如 YBE 约束、Ising TQFT 的 F,R 结构）综合起来，给出关于“Majorana 模式与拓扑门稳定性”的系统结论：
	- 在 YBE 子流形 $\mathcal M_R^{(\mathrm{YBE})}$ 或其在 2D/honeycomb 平台中的对应区域内，Berry 联络可以在适当规范下视为平坦，Majorana 编织实现的门与理想 Dehn twist 高度一致；
	- 离开这些区域时，曲率增大、holonomy 对路径和参数更敏感，$F_{\text{Dehn}}$ 下降，从而需要更强的容错机制（如 LQC+permutation 电路）来补偿几何误差。
- 这为“从代数解 YBE → 微观 Majorana 模型 → 拓扑门的几何实现与容错开销评估”提供了一条原则清晰、可数值实施的路线图。

今后在具体平台（例如更真实的材料建模、含相互作用与无序的 Majorana 体系）上，只要能写出合适的 BdG/Majorana 哈密顿量并识别出缺陷配置空间，就可以按本节的 7 个步骤逐一执行，从而系统研究 Majorana 模式的存在性、稳定性以及由其实现的拓扑量子门的几何性质与资源要求。

### 4.3.8 Honeycomb vison Berry：一个最小数值实验设计

作为本节流程的具体示例，这里给出一个针对 Kitaev 蜂窝格模型的最小 vison Berry 数值实验设计方案，后续可以据此实现脚本（例如 verify/run_honeycomb_vison_berry.py）。

1. **选平台与参数点**
	- 取有限尺寸的蜂窝格（例如 $L_x\times L_y=4\times4$ 个原胞，周期或开边界均可），用 kit4-2 中的 Majorana+Z2 背景场表示。
	- 在 $(J_x,J_y,J_z,h,\dots)$ 空间中选取处在非阿贝尔相内的一组参数 $p_0$，作为参考点。

2. **构造 Majorana 哈密顿量与基态 Z2 规范场**
	- 按照 kit4-2 的做法，把自旋 Kitaev 模型映射为在给定 Z2 规范场 $\{u_{ij}\}$ 背景下的自由 Majorana 跳跃哈密顿量 $H_c[u]$。
	- 在无通量扇区（所有六边形 $W_p=+1$）内求解基态对应的 $u_{ij}^{(0)}$，作为“真空规范场”。

3. **在有限格中放入两个 vison**
	- 在有限蜂窝格中选定两个相距适中的六边形 $p_A,p_B$，希望在它们处放入通量 $W_{p_A}=W_{p_B}=-1$。
	- 从真空规范场 $u_{ij}^{(0)}$ 出发，沿一条连接 $p_A$ 与 $p_B$ 的链路（branch cut）反转其上的 $u_{ij}$ 号，从而仅在端点六边形处产生 $W_p=-1$，得到含两个 vison 的规范场 $u_{ij}^{(\text{vison})}$。
	- 在该规范场下构造 Majorana 哈密顿量 $H_c[u^{(\text{vison})}]$，对角化并确认在 $p_A,p_B$ 附近存在成对的近零能局域态（Majorana 零模）。

4. **设计 vison 绕行路径**
	- 在六边形中心的平面图上，设计一条闭合路径 $\gamma$，描述 $p_A$ 绕 $p_B$ 一圈的 Dehn twist（例如沿围绕 $p_B$ 的最小环路）。
	- 把 $\gamma$ 离散成若干步：$\gamma: p_A^{(0)}\to p_A^{(1)}\to\cdots\to p_A^{(N)}=p_A^{(0)}$，每一步相当于把一个 vison 从当前六边形移动到相邻六边形。

5. **沿路径更新规范场并计算零模子空间**
	- 对 $\gamma$ 的每一步 $k=0,\dots,N$：
	  - 通过在 vison“步进”方向上的一条短 branch cut 上反转适当的 $u_{ij}$，得到新的规范场 $u_{ij}^{(k)}$，保证在每一步只有两个 vison，且其中一个的轨迹跟随 $p_A^{(k)}$；
	  - 在 $u_{ij}^{(k)}$ 下构造 Majorana 哈密顿量 $H_c[u^{(k)}]$ 并对角化，提取最低若干个本征态，选取组成 vison 相关零模子空间的正交归一基底 $\{\psi_a^{(k)}\}$；
	  - 记录每一步的零模谱，确认在整个路径中零模与低能连续演化、能隙不闭合。

6. **计算 vison 路径上的 Berry holonomy**
	- 对相邻步 $k\to k+1$，构造零模基底之间的重叠矩阵
	  $$
	  M_{ab}^{(k)}=\langle\psi_a^{(k+1)}|\psi_b^{(k)}\rangle,
	  $$
	  通过极分解 $M^{(k)}=U^{(k)}P^{(k)}$ 提取酉部分 $U^{(k)}$，并在整个路径上做有序乘积
	  $$
	  U_{\text{Berry}}[\gamma]=\mathcal P\prod_k U^{(k)}.
	  $$
	- 将 $U_{\text{Berry}}[\gamma]$ 规范到 SU(2)（去掉整体相位），在一个固定的逻辑基底（例如以两 vison 绑定的 Majorana 零模为编码子空间）中表示。

7. **与 Ising TQFT Dehn twist 比较**
	- 在同一逻辑基底中写出 Ising TQFT 层的 $F,R$ 矩阵，构造对应的 Dehn twist
	  $$
	  U_{\text{Dehn}}^{\text{Ising}}\simeq F^{-1}R^2F.
	  $$
	- 用 SU(2) 重合度评估数值 Berry holonomy 与理论 Dehn twist 的一致性
	  $$
	  F_{\text{Dehn}}=\frac{1}{2}\big|\mathrm{tr}\big(U_{\text{Berry}}[\gamma]^{\dagger}U_{\text{Dehn}}^{\text{Ising}}\big)\big|,
	  $$
	  并考察 $F_{\text{Dehn}}$ 在参数 $p$ 附近的稳定性（例如在 $(J_x,J_y,J_z,h,\dots)$ 的微扰下）。

8. **参数扫描与相图**
	- 固定上述几何设置（格点大小、两 vison 的相对位置、路径 $\gamma$），在参数空间中选择一条或若干条扫描线（如固定 $J_x=J_y$，改变 $J_z$ 或有效场强度），重复步骤 3–7：
	  - 记录能隙与拓扑不变量的变化；
	  - 记录 vison Dehn twist 的 Berry holonomy 与 Ising Dehn twist 的重合度 $F_{\text{Dehn}}(p)$。
	- 绘制能隙、拓扑不变量与 $F_{\text{Dehn}}(p)$ 的图像，标出“honeycomb 平台上的 Dehn twist 平坦区 vs 非拓扑区”，并与 1D R(a,b,c,d) 模型和 2D p+ip 平台中的相应结果对照。

这一实验设计为在蜂窝格上实现“vison 作为 Majorana 零模载体 + Dehn twist = Berry holonomy”的数值验证提供了具体蓝图，并把 4.3.1–4.3.7 中抽象的七步流程落实到一个可实现的 honeycomb 示例之中。

9. **当前数值实现与初步结果（verify/run_honeycomb_vison_berry.py）**

在上述方案的基础上，我们已经在 [verify/run_honeycomb_vison_berry.py](verify/run_honeycomb_vison_berry.py) 中实现了一个最小的 honeycomb vison-loop 数值实验，具体设定与结果概括如下：

- **格点与谱诊断**：
	- 取 $L_x=L_y=4$ 的有限蜂窝格（brick-wall 实现，每个原胞两点，共 $N=2L_xL_y=32$ 个 Majorana 自由度），构造对应的无通量 Z2 规范场 $u_{ij}^{(0)}$ 并建立自由 Majorana 哈密顿量 $H_c[u^{(0)}]$。
	- 初始谱在 $J_x=J_y=J_z=1$ 时的最低若干本征值为
		$$
		E\approx\{-0.0144,\;0.0144,\;\pm0.0929,\;\pm0.2698,\;\pm0.5,\dots\},
		$$
		表明存在一对接近零能的态与有限能隙，为在低能子空间上定义 vison 相关零模和 Berry holonomy 提供了合理起点。

- **plaquette/vison 结构与 branch cut**：
	- 脚本自动从最近邻图中枚举 6 边形 plaquette，$L_x=L_y=4$ 时共得到 9 个六边形，并计算对应的 Z2 通量 $W_p=\prod_{(ij)\in p}u_{ij}$；在初始规范场下所有 $W_p=+1$，对应无通量扇区。
	- 通过在 plaquette 邻接图上做 BFS，寻找两给定六边形 $p_A,p_B$ 之间的最短路径，并对相邻 plaquette 对选择一条公共边 $(ij)$ 作为 branch cut 上的 bond；沿该链路反转 $u_{ij}$ 号即可在 $p_A,p_B$ 上各产生一个 vison（$W_{p_A}=W_{p_B}=-1$），其余 plaquette 维持 $W_p=+1$。

- **vison-loop 绕行路径与 Berry holonomy**：
	- 在 9 个六边形的中心构成的 $3\times3$ 网格上，选取中心 plaquette 作为固定 vison $p_A$，再在其四周 8 个 plaquette 上构造一个离散矩形回路，让第二个 vison $p_B$ 沿此回路绕 $p_A$ 一圈。
	- 在回路上的每一步，重新构造从 $p_A$ 到当前 $p_B$ 的 branch cut 并建立对应的 $u_{ij}^{(k)}$，确保在演化过程中始终只有两个 vison 存在；对每个 $u^{(k)}$，对角化 $H_c[u^{(k)}]$，抽取最低能的二维子空间作为“vison 编码子空间”，用重叠矩阵极分解法累积 Berry holonomy $U_{\text{Berry}}$。
	- 在各向同性点 $J_x=J_y=J_z=1$、$L_x=L_y=4$、子空间维数 $\dim\mathcal H=2$ 时，得到的 vison-loop Berry holonomy 数值上接近恒等：
		$$
		U_{\text{Berry}}\approx I_{2\times2},\qquad \operatorname{spec}(U_{\text{Berry}})=\{1,1\}.
		$$
		将其规范到 SU(2) 后，与 Ising TQFT 中的 $R^{\sigma\sigma}$ 和 Dehn twist $F^{-1}R^2F$ 作 SU(2) 重合度比较，得到
		$$
		F_R\approx0.707,\qquad F_{\text{Dehn}}\approx0,
		$$
		表明在当前有限尺寸与参数下，这条 vison 绕行路径尚未实现理想的 Ising Dehn twist，而是更接近于一个“R 型相位结构”或平凡门。

- **沿 $J_z$ 方向的小参数扫描**：
	- 进一步在 $J_x=J_y=1$、$L_x=L_y=4$ 下，对 $J_z\in\{0.5,0.8,1.0,1.2,1.5\}$ 重复上述 vison-loop Berry 计算，并对每个 $J_z$ 提取 SU(2) 意义下的 $F_R(J_z)$ 与 $F_{\text{Dehn}}(J_z)$：
		$$
		F_R(J_z)\approx0.707\quad(\forall J_z),\qquad F_{\text{Dehn}}(J_z)\approx0\quad(\forall J_z),
		$$
		即在这条简单的 $J_z$ 轴切片上，vison 绕行路径实现的门在 SU(2) 意义下近似保持恒定，且始终远离 Ising Dehn twist。

从几何–代数图景看，这一 honeycomb 上的初步数值扫描与我们在 2D p+ip 上的结论形成互补：

- 在 2D p+ip + 两涡旋的连续/格点 BdG 模型中，沿 $\mu$ 轴存在一段 $F_{\text{Dehn}}(\mu)\approx1$ 的“Dehn twist 平坦区”，其中 Berry holonomy 在逻辑子空间中几乎刚性地实现了 Ising Dehn twist；
- 而在当前有限尺寸的 honeycomb toy 模型中，选取的这条 $(J_x=J_y=1, J_z$ 扫描$)$ 路径尚未进入对应的 Dehn twist 区：vison-loop Berry holonomy 近似保持恒等，$F_{\text{Dehn}}$ 近零，提示我们还需要在 $(J_x,J_y,J_z,h,\dots)$ 空间中寻找其它路径或区域来逼近真正的 Ising Dehn twist 区域。

因此，verify/run_honeycomb_vison_berry.py 当前给出的是 honeycomb 平台上的一个“可运行的 vison-loop Berry 几何实验骨架 + 初步参数切片”，证明了在蜂窝格上也可以按 4.3 节的统一流程系统实现 Berry holonomy 与 Ising F,R 的对比；下一步则可以在更大的系统和更合适的参数区间内搜索可能的 honeycomb 版“Dehn twist 平坦区”，从而在 1D R(a,b,c,d)、2D p+ip 与蜂窝格三条线上同时锁定实现相同 mapping class 群元素的拓扑/几何条件。

为了便于直观比较，我们还在同一脚本中输出了两张图像：

- [honeycomb_vison_F_Dehn_vs_Jz.png](honeycomb_vison_F_Dehn_vs_Jz.png)：横轴为 $J_z$，纵轴同时给出 $F_{\mathrm{Dehn}}(J_z)$ 与 $F_R(J_z)$。数值上，$F_R(J_z)$ 在所扫的 $J_z\in\{0.5,0.8,1.0,1.2,1.5\}$ 上几乎是一条 $\approx0.707$ 的水平线，而 $F_{\mathrm{Dehn}}(J_z)$ 则基本贴在 0 轴附近，图像呈现“一条高位平坦曲线 + 一条几乎为零的曲线”的对比。这一图像形象地说明：在当前参数和系统尺寸下，vison 绕行实现的门稳定地接近某个 R 型相位结构，却在 SU(2) 意义下与 Ising Dehn twist 明显不同。
- [honeycomb_vison_gap_vs_Jz.png](honeycomb_vison_gap_vs_Jz.png)：横轴为 $J_z$，纵轴为初始 vison 配置下最低本征值的绝对值 $\min|E|(J_z)$，作为能隙的简单 proxy。图上 $\min|E|(J_z)$ 在扫描区间内保持为一条缓慢变化但始终远离 0 的曲线，表明对于这条 $J_z$ 路径 vison 配置空间的谱缝并未闭合，当前看到的“$F_{\mathrm{Dehn}}\approx0$” 并不是简单的因能隙塌陷导致拓扑量子数失效，而是 Berry holonomy 本身在这一参数切片上就没有流向 Ising Dehn twist 的迹象。

从 4.3 的统一流程视角看，[honeycomb_vison_F_Dehn_vs_Jz.png](honeycomb_vison_F_Dehn_vs_Jz.png) 与 [honeycomb_vison_gap_vs_Jz.png](honeycomb_vison_gap_vs_Jz.png) 一起提供了 honeycomb 平台上“门的拓扑形态 vs 能隙行为”的最初二维截面：一方面能隙在该截面内相对平稳，另一方面 $F_{\mathrm{Dehn}}$ 却几乎为零，这恰好说明“拓扑门是否等价于 Ising Dehn twist”这一条件是额外而精细的几何/代数约束，并不会自动由“处在非阿贝尔相内且能隙不闭合”所保证。这也为后续在更大系统与更合理的 $(J_x,J_y,J_z,h,\dots)$ 路径上寻找 honeycomb 版 Dehn twist 平坦区提供了具体的数值参照。

