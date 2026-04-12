### kit-new：基于 R = exp(i H
a) 的统一流程草案

本笔记从“指数表示”
$$
R = e^{i H_P},
\qquad
H_P = \sum_{i,j\in\{0,x,y,z\}} c_{ij}\,\sigma_i\otimes\sigma_j,\quad (\sigma_0 = I)
$$
出发，串联起整个故事：

1. **常数 YBE → 指数形式解及参数约束**：从
$$
	R_{12}R_{13}R_{23} = R_{23}R_{13}R_{12}
$$
出发，用小参数展开得到对 $H_P$ 的李代数约束（经典 YBE），再讨论在 Pauli 线性 slice 上如何还原到 $a,b,c,d$ 的多项式条件；

2. **可积子流形与“离散曲率”**：把满足（量子/经典）YBE 的 $H_P$ 看成 $c_{ij}$ 空间中的“可积子流形”，并用“三角形回路两种走法之差”来定义离散 holonomy / 曲率；

3. **1D 推广：R → Kitaev/BdG → 端点零模**：在 1D 链上用 JW/Majorana 映射把 $H_P$ 分解成自由 BdG 部分和相互作用部分，给出有效 $(t,\Delta,\mu)$ 与相互作用通道，并用这些数据连接到拓扑区间 $|\mu|<2|t|$、端点零模、以及 4‑Majorana half twist 检验；

4. **2D 推广：p+ip 与 honeycomb**：把 2D BdG/Majorana 哈密顿写成若干局域 $H_P$ 的嵌入，讨论在涡旋/vison 配置空间上的 Berry 联络与 holonomy，指出哪些参数路径对应“近平坦层”，哪些对应强曲率；

5. **配置空间、联络与 Dehn twist / half twist / holonomy**：在抽象配置空间 $\mathcal C$ 上构造 Hilbert 丛 $\mathcal E\to\mathcal C$，给出 Berry 联络 $\mathcal A$、曲率 $F$、以及 holonomy $U[\gamma]$，并说明：

- 在可积/平坦子流形上，$U[\gamma]$ 实现 braid group / mapping class group 的表示，Dehn twist 与 half twist 可以纯粹由 $H_P$ 的组合写出；

- 在偏离可积层的方向 $V$ 上，曲率与 holonomy 的偏差如何与数值上观测到的 Dehn twist fidelity、LQC 复杂度增长对应起来。

下面按这个顺序详细展开。

---

#### 1. 常数 YBE 在指数表示下的展开与参数约束

我们从一般的两比特幺正算符
$$
R = e^{iH_P},
\qquad H_P = \sum_{i,j\in\{0,x,y,z\}} c_{ij}\,\sigma_i\otimes\sigma_j
$$
出发。常数 YBE 写成
$$
R_{12}R_{13}R_{23} = R_{23}R_{13}R_{12},
$$
其中 $R_{ab}$ 是 R 作用在三比特 Hilbert 空间的 (a,b) 对上，对其他比特恒等。

为分析 YBE 对 $H_P$ 的约束，我们在 $H_P$ 前加一个形式小参数 $\lambda$：
$$
R_{ab}(\lambda) = e^{i\lambda H_{ab}},\qquad H_{ab}\equiv H_P \text{ 作用在子空间 } a,b.
$$
要求 YBE 对所有 $\lambda$ 成立：
$$
R_{12}(\lambda)R_{13}(\lambda)R_{23}(\lambda)
 = R_{23}(\lambda)R_{13}(\lambda)R_{12}(\lambda).
$$

**一阶展开（自动满足）**

对小 $\lambda$ 展开
$$
R_{ab}(\lambda) = I + i\lambda H_{ab} + O(\lambda^2),
$$
代入 YBE 后，$O(\lambda)$ 阶两边都给出 $i(H_{12}+H_{13}+H_{23})$，一阶自动满足，不产生约束。

**二阶展开：经典 YBE（classical YBE）**

继续保留到 $O(\lambda^2)$，整理得到经典结果：若写
$$
R_{ab}(\lambda) = I + \lambda r_{ab} + O(\lambda^2),
\qquad r_{ab} = iH_{ab},
$$
则 YBE 的 $\lambda^2$ 阶等价于经典 YBE：
$$
[r_{12},r_{13}] + [r_{12},r_{23}] + [r_{13},r_{23}] = 0.
$$
换回 $H_{ab}$ 记号，条件写成
$$
[H_{12},H_{13}] + [H_{12},H_{23}] + [H_{13},H_{23}] = 0.
$$

在 Pauli 线性 ansatz 下，我们选择
$$
H_P = J_x\,\sigma^x\otimes\sigma^x + J_y\,\sigma^y\otimes\sigma^y + J_z\,\sigma^z\otimes\sigma^z,
$$
对应的算符两两对易，因此上式自动成立。进一步按 Pauli 线性形式
$$
R = aI\otimes I + b\,\sigma^x\otimes\sigma^x + c\,\sigma^y\otimes\sigma^y + d\,\sigma^z\otimes\sigma^z
$$
写出 YBE 时，约束退化为一组关于 $(a,b,c,d)$ 的多项式方程，这正是 verify/ybe_eqs.txt 中的内容——其解形成了旧框架中的代数子流形 $\mathcal M_R^{(\mathrm{YBE})}\subset\mathbb C^4$。

在新的指数表示下，我们更关心的是：

- **经典 YBE** 给出 $H_P$ 空间里“最低阶曲率为零”的李代数子空间；
- 真正满足量子 YBE 的 $H_P$ 形成其中更小的一块“可积子流形”，上面 $R=e^{iH_P}$ 满足完整的 operator YBE。

在后续 1D / 2D 模型中，我们可以选取一个满足（或近似满足）经典 YBE 的 $H_P^{(0)}$，再围绕它做扰动 $H_P = H_P^{(0)} + \varepsilon V$，用 classical/quantum YBE 作为“平坦/可积”的判据。

---

#### 2. 可积子流形、三角形回路与“离散曲率”

从几何角度看，把三体系统中三对作用 $H_{12},H_{13},H_{23}$ 想象为参数空间中沿三条边的生成元。沿每条边演化一个小参数 $\lambda$，对应幺正
$$
R_{ab}(\lambda) = e^{i\lambda H_{ab}}.
$$

围绕这三个点有两条自然的走法：
$$
U_\triangle^{(L)}(\lambda) = R_{12}(\lambda)R_{13}(\lambda)R_{23}(\lambda),\\
U_\triangle^{(R)}(\lambda) = R_{23}(\lambda)R_{13}(\lambda)R_{12}(\lambda).
$$
若常数 YBE 成立，则对所有 $\lambda$ 有 $U_\triangle^{(L)}=U_\triangle^{(R)}$。在小 $\lambda$ 极限展开到二阶得到
$$
U_\triangle^{(L)}(\lambda) - U_\triangle^{(R)}(\lambda)
 = -\lambda^2\bigl([H_{12},H_{13}] + [H_{12},H_{23}] + [H_{13},H_{23}]\bigr) + O(\lambda^3).
$$

这表明：

- 经典 YBE 条件
  $$
  [H_{12},H_{13}] + [H_{12},H_{23}] + [H_{13},H_{23}] = 0
  $$
  等价于“在 $O(\lambda^2)$ 阶上，两条走法绕三角形给出的 holonomy 一致”；
- 若把这个三角形看成参数空间上的一个极小 2‑simplex，其面积 $\mathrm{area}\sim\lambda^2$，则 $U_\triangle^{(L)}-U_\triangle^{(R)}$ 正是曲率 2‑form $F$ 在该小胞上的离散化测度。

因此：

- **可积子流形**：是 $H_P$ 空间中既满足经典 YBE 又满足高阶量子 YBE 的那部分，上面所有这样的小三角形回路的离散曲率都为零；
- **近平坦层**：只满足经典 YBE（或者经典 YBE 近似成立）的区域，小回路 holonomy 的偏差很小，可以视作 Berry 联络“几乎平坦”；
- **一般扰动方向**：离开这些层的方向 V 会立刻产生非零的离散曲率，在数值上表现为 Dehn twist / braid holonomy 对路径细节和微扰高度敏感。

这一层结构，是后面讨论 1D / 2D 模型、Dehn twist 稳健性以及电路复杂度时的“几何背景板”。

---

#### 3. 1D：从 R = exp(iH_P) 到 Kitaev 链、零模与 half twist

这一节把一般的两比特生成元 $H_P$ 放在一维链上，说明：

1. 如何在 JW/Majorana 映射下把 $H_P$ 分解成“自由 BdG 部分”和“相互作用/多体部分”；
2. 在小耦合极限下，自由部分如何给出有效 $(t,\Delta,\mu)$；
3. 为什么在适当的参数点上会出现端点零模和 4‑Majorana half twist；
4. classical YBE / 可积子流形在 1D 中的意义。

**3.1 把 H\_P 放到 1D 链上：JW/Majorana 分解**

考虑一条由 N 个自旋‑1/2 组成的开链，Hilbert 空间 $\mathcal H = (\mathbb C^2)^{\otimes N}$。在每条最近邻键 $\langle j,j+1\rangle$ 上放一个局域生成元
$$
H_P^{(j,j+1)}
 = \sum_{i,k\in\{0,x,y,z\}} c_{ik}\,\sigma_i^{(j)}\otimes\sigma_k^{(j+1)}.
$$
整条链的哈密顿量取为
$$
H = \sum_{j=1}^{N-1} H_P^{(j,j+1)}.
$$

对这类自旋链做 Jordan–Wigner 映射，得到一组费米算符 $c_j,c_j^\dagger$ 或等价的 Majorana 算符 $\gamma_{2j-1},\gamma_{2j}$，满足
$$
\{\gamma_a,\gamma_b\}=2\delta_{ab}.
$$
逐项检查可知：

- 形如 $\sigma^x_j\sigma^x_{j+1},\sigma^y_j\sigma^y_{j+1}$ 的两点算符，在 JW 后给出最近邻 hopping 和 p‑wave pairing 项，即二次 Majorana $i\gamma\gamma$；
- 形如 $\sigma^z_j,\sigma^z_{j+1}$ 的单点或密度项给出化学势（对角）贡献，同样是二次 Majorana；
- 而混合型算符（如 $\sigma^z_j\sigma^z_{j+1}$、带交叉 Pauli 的 $\sigma^x_j\sigma^z_{j+1}$ 等）一般会产生四 Majorana 项 $\gamma\gamma\gamma\gamma$ 或更远程耦合，对应真正的相互作用。

因此可以把总哈密顿量自然分解为
$$
H = H_{\mathrm{quad}} + H_{\mathrm{int}},
$$
其中

- $H_{\mathrm{quad}}$ 是所有 JW 后只含 $\gamma_a\gamma_b$ 的二次项之和，对应一个 1D BdG/Kitaev 链；
- $H_{\mathrm{int}}$ 收集所有 quartic 及以上的多体项，对应相互作用/非自由部分。

在我们原来的 Pauli 线性 slice 上，取
$$
H_P^{(j,j+1)} = J_x\,\sigma^x_j\sigma^x_{j+1} + J_y\,\sigma^y_j\sigma^y_{j+1} + J_z\,\sigma^z_j\sigma^z_{j+1},
$$
并专注于 Kitaev‑样的各向异性情形（例如 $J_z$ 主要通过单点项吸收到化学势通道），可以保证 $H_{\mathrm{int}}$ 很小或为零，整条链在低能上由一个有效的 1D p‑wave BdG 模型支配。

**3.2 有效 (t,\Delta,\mu)：小耦合下的线性映射**

在二体层面，先看 2‑site 模型。把 $H_P$ 只作用在站点 1,2 上，JW 后得到一个 4×4 BdG 矩阵。由于 Pauli 算符在 2 比特空间上张成 su(4)，而 Kitaev 链的 2‑site BdG 由三个实参数 $(t,\Delta,\mu)$ 决定（忽略整体能量平移），在小耦合极限下必然存在一个线性映射
$$
(t,\Delta,\mu) = (t_0,\Delta_0,\mu_0) + L\cdot c + O(\|c\|^2),
$$
其中 $c$ 是所有非零 $c_{ij}$ 组成的向量，$L$ 是一个实矩阵。verify/mapping\_fit\_results.json 与 mapping\_from\_micro.json 数值上确定了在旧的 $(a,b,c,d)$ slice 上
$$
t = b + c,\qquad \Delta = b - c,\qquad \mu = 4d + \mu_0,
$$
对应某个特殊的 $L$。

在一般 $H_P$ 的情形下，我们可以把 $L$ 拆成两块：

- 一块作用在“自由方向” $c^{(\mathrm{free})}$ 上，只改变 $(t,\Delta,\mu)$；
- 一块作用在“相互作用方向” $c^{(\mathrm{int})}$ 上，对应 $H_{\mathrm{int}}$ 的强度，无法简单吸收到 BdG 参数中。

这给出一个自然的线性稳定性判据：

- 若扰动只沿 $c^{(\mathrm{free})}$ 方向，小幅改变 $(t,\Delta,\mu)$ 且仍然满足 $|\mu|<2|t|$，则 1D 拓扑相与端点零模保持稳定；
- 若扰动有显著 $c^{(\mathrm{int})}$ 分量，则会生成 $H_{\mathrm{int}}$，其对 half twist 与 Dehn twist 的影响需要用前面讨论的 H\_P=H\_0+V 框架量化（见 kit-exp 第 6.4 节）。

**3.3 端点零模与 4‑Majorana half twist**

对于无相互作用、齐次参数的 1D Kitaev 链
$$
H_\mathrm{K} = -\sum_j\bigl(t c_j^\dagger c_{j+1} + \Delta c_j c_{j+1} + \mathrm{h.c.}\bigr) - \mu\sum_j\Bigl(c_j^\dagger c_j - \tfrac12\Bigr),
$$
标准的动量空间求解给出能谱
$$
E(k) = \pm\sqrt{(\mu+2t\cos k)^2 + 4|\Delta|^2\sin^2 k},
$$
其能隙在 $k=0,\pi$ 处关闭，当且仅当 $|\mu|=2|t|$。因此在 $|\mu|<2|t|$ 区间内链处于拓扑相，开边界时存在指数局域在两端的 Majorana 零模，这是 [Kitaev 2001] 的标准结果，这里不再赘述证明。

在最小的 2‑site 拓扑链（或者等价的 4‑Majorana 系统）上，可以显式写出 Majorana 算符
$$
\gamma_1,\gamma_2,\gamma_3,\gamma_4
$$
并验证在特定参数点（例如 $t=\Delta=1,\mu=0$）哈密顿量简化为
$$
H_0 = \frac{i}{2}\,\gamma_2\gamma_3,
$$
对应 2 个成对耦合的 Majorana 和 2 个自由零模。对这个 H\_0 取时间 $\tau=\pi/2$，时间演化算符为
$$
U_0(\tau) = e^{-iH_0\tau} = \exp\Bigl(\frac{\pi}{4}\,\gamma_2\gamma_3\Bigr),
$$
这正是标准的 half twist 生成元。详细的矩阵表示与谱分解可以在 verify/run\_instantaneous\_braid\_crosscheck.py 中直接数值验证：在那里我们显式构造了 4×4 Majorana 表示并比较了 $U_0(\tau)$ 与 $\exp((\pi/4)\gamma_2\gamma_3)$。

若现在一般地写
$$
H_P = H_0 + V
$$
并假设 $\|V\|$ 很小，则前文在 kit-exp 6.4 节的 Dyson 展开给出
$$
\|e^{-iH_P\tau} - e^{-iH_0\tau}\| \le \tau\,\|V\| + O(\|V\|^2),
$$
在投影到零模子空间后同样成立，量化了 half twist 对一般扰动 V 的线性稳健性。这一结论在 run\_interacting\_braid\_check.py 中通过扫描相互作用强度 U,\mu 得到数值验证。

**3.4 1D 中的可积子流形**

把 1D 链看作无穷长，并在每条键上放相同的 $H_P$。若 $H_P$ 来自某个解常数 YBE 的 R‑矩阵（例如 XXZ / XYZ 等可积自旋链），则整条链对应一个可积模型，其散射矩阵满足“因子化散射”，Bethe Ansatz 可解。我们这里不重新构造完整的可积理论，只记下几点与几何 picture 直接相关的事实：

1. 在可积点上，excitation 之间的散射完全由两体 R 决定，多体散射可分解为两体散射的组合次序（这正是 operator YBE 的物理含义）；
2. 这意味着在多粒子配置空间中，“绕开哪一个粒子先后顺序”的 holonomy 只依赖于拓扑同伦类，而与具体拼接顺序无关，对应 Berry 联络的平坦性；
3. classical YBE 给出的“近平坦层”则可以理解为：在某些方向上的扰动只在更高阶（多重散射）时才积累出明显曲率，多体散射在一段能量/长度尺度上仍近似因子化。

因此，在 1D 场景中，可积子流形与 $H_P$ 空间中的“平坦/近平坦层”直接对应，前者可以通过 YBE 的代数约束刻画，后者可以通过小三角形回路的离散曲率近似刻画，两者为我们理解 half twist / Dehn twist 在 1D 中的稳健性提供了精确与近似两种层次。

---

#### 4. 2D：从 R = exp(iH_P) 到 p+ip / honeycomb 与缺陷 Berry holonomy

在 2D 中，我们不再把 $H_P$ 仅仅视为 1D 链上的邻接项，而是把它嵌入到 2D 格点上的每条键或邻居对之间。这里以两类模型为代表：

1. 2D p+ip BdG 模型（实空间上是 spinless 费米子 + 最近邻 hopping 与有相位的配对）；
2. Kitaev 蜂窝模型的 Majorana+Z\_2 描述（每条键上一个 $J_a\,\sigma^a_i\sigma^a_j$）。

**4.1 2D p+ip 模型中的 H\_P 嵌入**

在 verify/run\_pip\_vortex\_berry.py 中，我们构造了一个有限格点的 2D p+ip BdG 哈密顿量
$$
H_{\mathrm{BdG}} = \frac12\sum_{\langle ij\rangle}
\Bigl[-t\,(c_i^\dagger c_j + c_j^\dagger c_i)
  + (\Delta_{ij} c_i c_j + \mathrm{h.c.})\Bigr]
 - \mu\sum_i\Bigl(c_i^\dagger c_i - \tfrac12\Bigr),
$$
其中 $\Delta_{ij}$ 在 x,y 方向分别带有 $p_x$、$ip_y$ 的结构，并通过涡旋构造引入空间依赖的相位。若我们把每一个最近邻对 $\langle ij\rangle$ 看成一个“局域两比特 Hilbert 空间”的投影，则可以在该对上引入一个 H\_P^{(i,j)}，使得它在 JW/Majorana 映射后给出对应的 hopping + pairing 结构。

形式上可以写为
$$
H_{\mathrm{BdG}} = \sum_{\langle ij\rangle} H_P^{(i,j)} + H_{\mathrm{onsite}},
$$
其中 $H_{\mathrm{onsite}}$ 收集化学势等单点项。若我们选择每条键上的 H\_P^{(i,j)} 来自某个满足（近似）YBE 的 R‑矩阵族，那么在“拉直”到一维（例如沿某条条带）时，我们可以把这整个 2D 模型看成是一维可积链的某种二维推广，只是现在缺陷（涡旋）的位置也进入了参数空间。

**4.2 蜂窝 Majorana+Z\_2 模型中的 H\_P 嵌入**

在 Kitaev 蜂窝模型中，每条键 $\langle ij\rangle\_a$ 上的自旋耦合是
$$
H_{ij}^{(a)} = J_a\,\sigma^a_i\sigma^a_j,\quad a\in\{x,y,z\}.
$$
在 Majorana+Z\_2 表示中（见 kit2 / verify/run\_honeycomb\_vison\_berry.py），引入实 Majorana $c_j,b^a_j$ 并定义
$$
\sigma^a_j = i b^a_j c_j,\qquad u^{(a)}_{ij} = -i b^a_i b^a_j = \pm1,
$$
则
$$
H_{ij}^{(a)} = J_a\,u^{(a)}_{ij}\,(ic_i c_j).
$$
因此可以把每条键上的“生成元”看成
$$
H_P^{(i,j)} = J_a\,\sigma^a_i\sigma^a_j,
$$
其指数 $R^{(i,j)}=e^{iH_P^{(i,j)}}$ 在固定的 Z\_2 背景 $u^{(a)}_{ij}$ 下给出 Majorana 哈密顿量的局域贡献。蜂窝模型本身就是一个可积模型，其 R‑矩阵可以通过 8‑vertex / free fermion 形式导出，这里不重复完整构造，只强调：

- 在各向同性点附近（以及某些各向异性线）上，蜂窝模型的散射结构与 YBE 解直接相关；
- 在这些点附近，小的参数变化沿某些方向保持“近平坦”，对应我们在 kit4 中观察到的相图边界和平带结构。

**4.3 缺陷配置空间与 Berry holonomy：严谨定义**

不论是 2D p+ip 还是蜂窝模型，当我们在系统中插入涡旋、vison 或其他局域缺陷时，可以定义一个“缺陷配置空间”
$$
\mathcal C = \{X=(x_1,\dots,x_n)\mid x_i\in\Sigma,\ x_i\neq x_j\},
$$
其中 $\Sigma$ 是空间或流形（例如二维环面或平面上的有界区域），$n$ 是缺陷数目。对每个配置 $X\in\mathcal C$，考虑相应哈密顿量 $H(X)$，以及其零能/近零能本征空间
$$
\mathcal H_0(X) = \mathrm{span}\{\text{零模/基态}\}.
$$
这些子空间拼在一起构成一个 Hilbert 丛
$$
\pi: \mathcal E\to\mathcal C,\qquad \pi^{-1}(X)=\mathcal H_0(X).
$$

选取每个 $X$ 上的一组正交归一基 $|\psi_a(X)\rangle$（局部截面），定义投影
$$
P(X) = \sum_a |\psi_a(X)\rangle\langle\psi_a(X)|.
$$
Berry 联络 1‑form 由
$$
\bigl(A_\mu(X)\bigr)_{ab} = i\,\langle\psi_a(X)|\partial_{\mu}\psi_b(X)\rangle,
\qquad \partial_{\mu}=\frac{\partial}{\partial X^\mu},
$$
或者投影形式
$$
A_\mu(X) = i\,P(X)\,\partial_{\mu}P(X)\,P(X).
$$
曲率 2‑form 为
$$
F_{\mu\nu}(X) = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu,A_\nu]
     = P[\partial_\mu P,\partial_\nu P]P.
$$

给定一条在配置空间中的闭合路径 $\gamma:[0,1]\to\mathcal C$，Berry holonomy 定义为
$$
U[\gamma] = \mathcal P\exp\Bigl(-\int_0^1 A_\mu(\gamma(s))\,\dot X^\mu(s)\,ds\Bigr),
$$
这是在零模子空间上的一个幺正算符，在不同规范下仅相差一个共轭变换。若 $F=0$（平坦）或 $F$ 为中心元，则 $U[\gamma]$ 只依赖于 $\gamma$ 的同伦类，从而给出 braid group 或 mapping class group 的幺正表示。

在我们的具体数值脚本中（run\_pip\_vortex\_berry.py, run\_honeycomb\_vison\_berry.py），我们用离散化路径、相邻步重叠矩阵的极分解来近似计算这个 $U[\gamma]$，并在 2D p+ip 中验证了它与 Ising TQFT 的 Dehn twist 矩阵在 SU(2) 中共轭（在某个 $\mu$ 区间形成“平台”）。

---

#### 5. 配置空间、可积 H\_P 与 Dehn twist / half twist / holonomy

最后把前面的所有结构汇总到配置空间 picture 中，明确 Dehn twist、half twist 与 R = e^{iH_P} / Berry holonomy 之间的对应关系。

**5.1 braid / mapping class group 与半平坦联络**

在缺陷配置空间 $\mathcal C$ 上，纯 braid group $PB_n$ 的元素可以看作是固定终点的闭合路径类，full braid group $B_n$ 则允许整体置换。若 Berry 联络在某个参数区域内平坦（$F=0$）或仅有中心曲率，则：

1. 任意两条同伦等价的闭合路径 $\gamma_1,\gamma_2$ 给出的 holonomy 在零模子空间内相同（或仅差一个整体相位）；
2. 取一组标准生成元（如相邻粒子交换的 half twist 对应的路径），即可为 braid group 定义一个单值的幺正表示；
3. 在 genus>0 的底空间上，进一步的周期/同伦类（如环绕 handle 的 Dehn twist）给出 mapping class group 的表示。

对于 Ising 拓扑序，TQFT 中的 F,R 符号给出了一套抽象的表示：
$$
R^{\sigma\sigma} = \begin{pmatrix} e^{-i\pi/8} & 0 \\ 0 & e^{3i\pi/8} \end{pmatrix},
\qquad
U_{\mathrm{Dehn}} \sim F^{-1}(R^{\sigma\sigma})^2 F,
$$
在 2‑$\sigma$ fusion 空间中等价于一个 SU(2) 矩阵（剥去整体相位后）。数值上我们在 4‑Majorana 模型和 2D p+ip 涡旋模型中都计算了 Berry holonomy，发现其在拓扑区间内与 TQFT 的 Dehn twist 矩阵在 SU(2) 中共轭，这就是“Dehn twist plateau”的实质。

**5.2 R = e^{iH_P}、half twist 与 Berry holonomy 的统一视角**

在 1D/4‑Majorana 情形下：

- 选择特定的 $H_P^{(0)}$ 使其 JW 后等价于 2‑site Kitaev 链在拓扑点 $t=\Delta=1,\mu=0$ 的哈密顿量 $H_0=(i/2)\gamma_2\gamma_3$；
- 这时 R = e^{iH_P^{(0)}} 在适当时间步长下就是 half twist 生成元；
- 把这类 R 作为局部门放在 1D/2D 网络中，并要求它们满足（近似）YBE，就得到一个近似平坦的联络，其 Berry holonomy 与 R‑矩阵的代数组合一致。

更抽象地说：

1. **局域层面**：R = e^{iH_P} 给出“沿一条短路径段的平行移动元”；
2. **代数层面**：YBE 确保这些平行移动元在三体/多体组合时的关联方式只依赖于 braid 结构；
3. **几何层面**：在配置空间上把这些局部 R 沿路径拼在一起，就得到 Berry holonomy；
4. **拓扑层面**：在平坦/可积子流形上，该 holonomy 只依赖于路径的同伦类，从而与 TQFT 的 F,R 表示一致；
5. **数值/复杂度层面**：离开可积层的扰动 V 通过曲率 $F$ 在 holonomy 中表现出来，其对 LQC+permutation 电路可实现性的影响则可由 fidelity 函数 $F_{\max}(\varepsilon;V)$ 的“complexity 曲率” $K_{\mathrm{comp}}(V)$ 来间接测量（详见 kit-exp 第 6.5 节）。

**5.3 总结**

新的指数表示 $R=e^{iH_P}$ 把原来基于 $(a,b,c,d)$ 的 R‑矩阵故事提升到一个更统一的框架：

- 在 $H_P$ 空间中，量子 YBE / 经典 YBE 刻画出可积/近平坦子流形，它们在几何上对应 Berry 联络曲率 $F$ 的消失或弱化；
- 1D 中，这些子流形给出 Kitaev 链及其零模/half twist 的稳定区域；
- 2D 中，它们通过 p+ip 和 honeycomb 模型的嵌入，在涡旋/vison 配置空间上产生近似拓扑的 Dehn twist holonomy；
- 偏离这些子流形的扰动 V 则通过 Dyson 展开和 classical YBE 显式地体现在 holonomy 偏差和 complexity 曲率上，为数值实验和电路实现提供了一个可计算的“曲率/复杂度” 词典。

接下来要做的线性稳定性分析和数值参数扫描，可以完全在 $H_P$ 或 $c_{ij}$ 空间中表述：选一个可积点 $H_P^{(0)}$，挑选若干代表性的扰动方向 $V$，沿每个方向测量 Berry holonomy、Dehn twist fidelity 和 LQC 复杂度随 $\varepsilon$ 的变化，从而在这张统一的几何图上填充精确的数据点。

**5.4 R, half twist, Dehn twist 与链/2D 嵌入的词典**

为了便于在后续数值和讨论中快速对照，这里把上面零散出现的对应关系整理成一个简洁的“词典”：

- **局域 R 与 half twist（4‑Majorana / 2‑site）**：
  在 4‑Majorana 极限点（例如 $t=\Delta=1,\mu=0$）上，$H_0=(i/2)\gamma_2\gamma_3$，时间演化 $U_0(\tau=\pi/2)=\exp((\pi/4)\gamma_2\gamma_3)$ 就是一次 half twist。选取一个 $H_P^{(0)}$ 在 2 比特空间等价于这个 $H_0$，则
  $$
  R^{(0)} = e^{iH_P^{(0)}}
  $$
  在合适时间/标度下就是“局域 half twist 门”，作用在 2‑$\sigma$ 融合空间或 4‑Majorana 零模子空间上，与拓扑场论的 half twist 同构（差一个整体相位和基变换）。

- **R 作为 YBE 解与 braid / Dehn twist**：
  若 $R=e^{iH_P}$ 满足常数 YBE，定义多粒子空间上的
  $B_i := R_{i,i+1}$。$B_i$ 自动满足 braid 关系，从而给出 braid group 的一个表示，其中单个 $B_i$ 就是代数上的“half twist 生成元”。Dehn twist 则不是单个 R，而是若干 F‑move 与 R 组合而成的 word，例如 Ising TQFT 中的
  $U_{\mathrm{Dehn}}\sim F^{-1}(R^{\sigma\sigma})^2F$，对应“绕一个 handle 或环路做一次满绕”的同伦类。

- **1D 链中的嵌入（Hamiltonian / 电路双视角）**：
  Hamiltonian 视角下，把同一个两体生成元 $H_P$ 放到每条键上 $H=\sum_j H_P^{(j,j+1)}$，得到 Kitaev 链或自旋链；其端点零模和 4‑Majorana half twist 是整个 1D 模型的拓扑边界效应。电路视角下，把 $R=e^{iH_P}$ 当作一层局部门，构造一步演化 $U_{\text{step}}=\prod_j R_{j,j+1}$，多步叠加后，可以近似实现把端点 Majorana 沿链搬运和交换，因此“物理 half twist”是这一串 R 在零模子空间上的作用。

- **2D p+ip / honeycomb 中的嵌入与 Berry holonomy**：
  在 2D p+ip 中，每条最近邻键上的 BdG 耦合都可以看成某个 $H_P^{(i,j)}$ 的 JW 映像，指数 $e^{iH_P^{(i,j)}}$ 是该键上的局域 R。缓慢移动涡旋的过程可以离散成很多小步，每一步在零模丛上做一次平行移动；在可积/近平坦区域，这个离散平行移动和用 R 代数组合出来的 braid 表示是同一个幺正群表示。在蜂窝 Majorana+Z\_2 模型中，每条键本身就对应 $H_P^{(i,j)}=J_a\sigma_i^a\sigma_j^a$，指数给出局域 R；可积结构保证在特定参数线上缺陷绕行的 Berry holonomy 精确/近似等同于抽象的 braid/Dehn twist 表示。

- **几何与复杂度层面**：
  在 $H_P$ 空间中，满足 YBE 的点给出严格平坦或中心曲率的 Berry 联络，holonomy 只依赖于路径的同伦类；classical YBE 则描述近平坦层。离开这些区域的扰动 $V$ 在几何上通过曲率 $F$ 和小三角离散曲率显现，在数值上通过 Dehn twist fidelity 的下降和 LQC+permutation 的 complexity 曲率 $K_{\mathrm{comp}}(V)$ 显现。这样，R‑矩阵、half twist、Dehn twist、Berry holonomy 与“复杂度地形”被统一到了同一张几何图上。



