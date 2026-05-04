### kit-exp：从线性 Pauli R 到指数 R=exp(iH_P) 的 YBE / 平坦联络

本节作为独立实验笔记，专门讨论把最初的两比特 R 从线性 Pauli ansatz
$$
R = a\,I\otimes I + b\,\sigma^x\otimes\sigma^x + c\,\sigma^y\otimes\sigma^y + d\,\sigma^z\otimes\sigma^z
$$
改写为更一般的指数形式
$$
R = e^{i H_P},\qquad H_P = \sum_{i,j\in\{0,x,y,z\}} c_{ij}\,\sigma_i\otimes\sigma_j,\quad (\sigma_0=I),
$$
之后在“YBE ↔ 平坦联络”这一步会发生什么变化，尤其是小耦合极限下的李代数版 YBE（classical YBE）及其几何含义。

 

#### 1. 旧的 Pauli 线性 ansatz：代数子流形 \(\mathcal M_R^{(\mathrm{YBE})}\) 与 F=0 区域

在早期的 1D/kit3 推导中，我们选取了一个非常“刚性”的 ansatz：
$$
R(a,b,c,d) = a\,I\otimes I + b\,\sigma^x\otimes\sigma^x + c\,\sigma^y\otimes\sigma^y + d\,\sigma^z\otimes\sigma^z.
$$
这个选择有两个关键特点：

1. 只保留了四个彼此对易的算符基底：
   $$
   A=\sigma^x\otimes\sigma^x,\quad B=\sigma^y\otimes\sigma^y,\quad C=\sigma^z\otimes\sigma^z,\quad I\otimes I,
   $$
   满足 \([A,B]=[A,C]=[B,C]=0\), 且 \(A^2=B^2=C^2=I\)。
2. 因此 R 的任何多体组合（例如 YBE 中出现的 \(R_{12}R_{13}R_{23}\) 等）在计算时仍然只会落在 \(I,A,B,C\) 张成的 4 维子空间里，YBE
   $$
   R_{12}R_{13}R_{23} = R_{23}R_{13}R_{12}
   $$
   被化简成一组关于 \((a,b,c,d)\) 的代数方程（见 verify/ybe_eqs.txt）。

其结果是：

- 常数 YBE 的解集在 \(\mathbb C^4\) 中形成一个代数子流形
  $$
  \mathcal M_R^{(\mathrm{YBE})} \subset \mathbb C^4_{(a,b,c,d)},
  $$
  在 [kit3-5.md](kit3-5.md) 里我们把它解释为“R 参数空间里联络曲率 \(F=0\)（或中心元）的一块”。
- 几何上，我们可以把 \((a,b,c,d)\) 看成某个“configuration space 上的平行移动元”的坐标：
  - 对每一条局域边/路径段，赋予一个 R(a,b,c,d)；
  - YBE 确保在三点系统里，不同的局部拼接顺序给出的总平行移动是一样的；
  - 这等价于在相应的 Hilbert 丛上，Berry 联络的曲率张量 \(F=dA+A\wedge A\) 在该子流形上为零（或者只留下全局 U(1) 相因子）。

这种 picture 之所以“干净”，是因为我们把 R 限制在一个对易 Pauli 子代数里，所有组合都不会跑出这个 4 维线性空间。

 

#### 2. 新的指数 ansatz：\(R=e^{iH_P}\) 及其参数空间

现在改用更一般的两比特哈密顿量
$$
H_P = \sum_{i,j\in\{0,x,y,z\}} c_{ij}\,\sigma_i\otimes\sigma_j,\quad c_{ij}\in\mathbb R,
$$
并设
$$
R = e^{iH_P}.
$$
几点立刻发生的变化：

- Pauli 基底总体上生成 su(4)（或其实数形式），多数项彼此不对易：
  $$
  [\sigma_i\otimes\sigma_j,\,\sigma_{i'}\otimes\sigma_{j'}]\neq 0 \quad \text{(一般情形)}.
  $$
- 指数展开
  $$
  e^{iH_P} = I + iH_P - \tfrac12 H_P^2 + \cdots
  $$
  会在乘法下产生整个 16 维 Pauli 基底的混合；因此 \(R\) 不再用 4 个标量 \((a,b,c,d)\) 描述，而是“坐标”为 \(c_{ij}\) 的一块 su(4) 流形上的点。
- 原先的 Pauli 线性 ansatz 可以被视为一种特殊情形：
  - 选取只含 \(I,\sigma^x\sigma^x,\sigma^y\sigma^y,\sigma^z\sigma^z\) 的 \(H_P\)，并要求它们两两对易；
  - 此时 \(e^{iH_P}\) 的指数可以在这个 4 维子空间内被完全重新写成“线性形式” R(a,b,c,d)，从而回到旧框架。

结论：**新的指数 ansatz 把我们从一个 4 维“对易 slice” 提升到一个接近整个 su(4) 的参数空间；旧的 \(\mathcal M_R^{(\mathrm{YBE})}\) 只是这个大空间中一条非常特殊的可积子流形。**

 

#### 3. 常数 YBE 的指数化与小耦合展开：classical YBE

我们希望理解，在 \(R=e^{iH_P}\) 的情形下，常数 YBE
$$
R_{12}R_{13}R_{23} = R_{23}R_{13}R_{12}
$$
在小耦合极限（所有 \(c_{ij}\) 很小）会对 \(H_P\) 施加什么样的“李代数版”约束。

为此引入一个形式小参数 \(\lambda\)，写
$$
R_{ab}(\lambda) = e^{i\lambda H_{ab}},\quad H_{ab}\equiv H_P \text{ 作用在子空间 } a,b.
$$
要求 YBE 对所有 \(\lambda\) 成立：
$$
R_{12}(\lambda) R_{13}(\lambda) R_{23}(\lambda)
= R_{23}(\lambda) R_{13}(\lambda) R_{12}(\lambda).
$$
在 \(\lambda\to 0\) 时两边都趋近于单位算符，我们可以把两侧展开为 \(\lambda\) 的级数并逐阶比较。

**一阶 (\(\lambda\))：** 直接展开
$$
R_{ab}(\lambda) = I + i\lambda H_{ab} + O(\lambda^2)
$$
可见 YBE 左右两边在 \(O(\lambda)\) 上都给出同一个和 \(i(H_{12}+H_{13}+H_{23})\)，一阶自动满足，没有约束。

**二阶 (\(\lambda^2))：** 把三项相乘并保留 \(O(\lambda^2)\) 项，可以得到经典的结果：

- 若我们写
  $$
  R_{ab}(\lambda) = I + \lambda r_{ab} + O(\lambda^2),\quad r_{ab}\equiv iH_{ab},
  $$
  并仅保留到二阶，则 YBE 的 \(\lambda^2\) 阶等价于所谓经典 YBE（classical Yang–Baxter equation）：
  $$
  [r_{12},r_{13}] + [r_{12},r_{23}] + [r_{13},r_{23}] = 0.
  $$
- 换回 \(H_{ab}\) 的记号，则条件可以写成
  $$
  [H_{12},H_{13}] + [H_{12},H_{23}] + [H_{13},H_{23}] = 0.
  $$

这就是我们在“李代数版 YBE”中想要的约束：三体子系统中，不同成对哈密顿量 \(H_{ab}\) 之间满足一个对易关系的和为零。

再往高阶走，需要用 Baker–Campbell–Hausdorff (BCH) 展开：
$$
\log(e^A e^B) = A + B + \tfrac12[A,B] + \tfrac1{12}[A,[A,B]] - \tfrac1{12}[B,[A,B]] + \cdots,
$$
把三项指数的乘积反复合并，得到等价于某个“量子 YBE”对 \(H_{ab}\) 的无穷多李代数约束；经典 YBE 出现在最低非平凡阶 \(O(\lambda^2))\)。

**物理上**：

- 旧的 Pauli 对角 ansatz 中，由于 \(H_{ab}\) 选在彼此对易的子代数里，\([H_{12},H_{13}]\) 等天然为零，经典 YBE 自动满足，进一步的高阶约束退化为多项式代数条件（即我们在 verify/ybe_eqs.txt 里看到的那些方程）。
- 在一般的 \(H_P\) 中，经典 YBE 提供了一个“平坦极限”的一阶诊断：若二阶就已经违反，则说明我们离可积/平坦结构较远；若经典 YBE 近似成立，说明可以把对应区域视作“弱曲率”的近似平坦层。 

 

#### 4. 几何解释：李代数 YBE 与联络曲率 F

在配置空间/参数空间的 Berry 联络语言中，我们把低能 Hilbert 丛 \(\mathcal E\to\mathcal C\) 上的联络 1‑form 写成 Lie 代数值：
$$
\mathcal A = \sum_{\mu} A_\mu(\lambda)\,d\lambda^\mu,\quad A_\mu(\lambda)\in\mathfrak g,
$$
其中 \(\mathfrak g\) 在我们的 Majorana/BdG 场景下通常是 \(so(2n)\) 或其子代数。曲率为
$$
F = d\mathcal A + \mathcal A\wedge\mathcal A.
$$

下面把“联络 A、曲率 F、holonomy U[\gamma]”这几步用比较标准的 Berry‑Wilczek–Zee 形式详细推一遍，再和 \(R=e^{iH_P}\) / YBE 联系起来。

**4.1 从本征态/投影子空间构造联络 A**

设参数空间为 \(\mathcal C\)，坐标为 \(\lambda=(\lambda^1,\lambda^2,\dots)\)。在每个 \(\lambda\) 上，选定一个有限维“低能子空间” \(\mathcal H_0(\lambda)\subset\mathcal H\)，由一组正交归一的本征态
$$
\{|\psi_a(\lambda)\rangle\}_{a=1}^{N_0}
$$
张成（例如零模/基态简并空间）。定义投影算符
$$
P(\lambda) = \sum_{a=1}^{N_0} |\psi_a(\lambda)\rangle\langle\psi_a(\lambda)|.
$$

**(1) Berry 联络的矩阵元形式**

在这一基底下，非阿贝尔 Berry 联络分量定义为
$$
\bigl(A_\mu(\lambda)\bigr)_{ab}
\equiv i\,\langle\psi_a(\lambda)|\partial_{\mu}\psi_b(\lambda)\rangle,
\qquad \partial_{\mu}\equiv \frac{\partial}{\partial\lambda^\mu}.
$$
于是
$$
\mathcal A = \sum_\mu A_\mu(\lambda)\,d\lambda^\mu,\qquad
A_\mu(\lambda)\in \mathfrak u(N_0).
$$

在规范变换
$$
|\psi_a(\lambda)\rangle \mapsto \sum_b |\psi_b(\lambda)\rangle\,U_{ba}(\lambda),
\quad U(\lambda)\in U(N_0)
$$
下，联络按
$$
A_\mu \mapsto U^{-1}A_\mu U + i\,U^{-1}\partial_\mu U
$$
变换，这是标准的规范联络变换律。

**(2) 用投影算符 P 表达的形式**

Berry 联络、曲率也可以直接用投影算符来写。这对和 \(H_P\) 的算符表达打通很有用。对任意在 \(\mathcal H_0(\lambda)\) 中的向量
$$
|\phi(\lambda)\rangle \in \mathrm{Im}\,P(\lambda)
$$
施加“平行移动”条件
$$
P(\lambda)\,\frac{d}{d\lambda^\mu}|\phi(\lambda)\rangle = 0,
$$
可以推出在 \(\mathcal H_0\) 上诱导的协变导数形式
$$
\nabla_\mu = P\,\partial_\mu + i A_\mu P,
$$
从而得到一个等价的表达（在投影到子空间后）：
$$
A_\mu = i\,P\,\partial_\mu P\,P,
$$
以及曲率的投影形式
$$
F_{\mu\nu} = P\,[\partial_\mu P,\partial_\nu P]P.
$$
这两条公式只依赖于投影算符 \(P(\lambda)\)，不用显式选一组本征态基底。

在 Majorana/BdG 场景下，常见情形是：

- 在每个 \(\lambda\) 上有一组零能/近零能本征矢 \(|\psi_a(\lambda)\rangle\)；
- 用这些矢量定义 \(P(\lambda)\)，再用上式得到 \(A_\mu,F_{\mu\nu}\)。

这就是我们在数值脚本中用“重叠矩阵 → 极分解 → 步进平行移动”来近似实现的对象的连续极限版本。

**4.2 holonomy：从 A,F 到路径幺正 U[\gamma]**

给定参数空间中的一条闭合路径
$$
\gamma:[0,1]\to \mathcal C,\qquad \gamma(0)=\gamma(1)=\lambda_0,
$$
沿路径的非阿贝尔 holonomy 定义为
$$
U[\gamma]
 = \mathcal P\exp\Bigl(-\int_0^1 A_\mu(\gamma(s))\,\dot\lambda^\mu(s)\,ds\Bigr),
$$
其中 \(\mathcal P\) 为路径有序指数，\(\dot\lambda^\mu=\frac{d\lambda^\mu}{ds}\)。在选定规范下，\(U[\gamma]\) 是 \(U(N_0)\) 中的一个幺正矩阵；在不同规范之间，\(U[\gamma]\) 只差一个共轭变换。

对一个面积很小的“微元回路”（例如 \(\lambda^\mu,\lambda^\nu\) 平面上边长 \(\delta\lambda^\mu,\delta\lambda^\nu\) 的小矩形 \(\square_{\mu\nu}\)），可以把 holonomy 近似展开为
$$
U[\square_{\mu\nu}]
\approx I - F_{\mu\nu}(\lambda_0)\,\delta\lambda^\mu\,\delta\lambda^\nu + O(\delta\lambda^3).
$$
因此：

- **曲率 \(F_{\mu\nu}\) 就是“小回路 holonomy 偏离单位算符的一阶（面积阶）测度”**；
- 若 \(F_{\mu\nu}(\lambda_0)=0\)，则任意足够小、围绕 \(\lambda_0\) 的回路的 holonomy 在该点处都是 \(I+O(\text{高阶})\)，体现局部平坦；
- 积分到有限大小的回路上（例如 Dehn twist 轨道），holonomy 一般是
  $$
  U[\gamma] = \mathcal P\exp\oint_\gamma A = \exp\Bigl(-\iint_S F + \cdots\Bigr),
  $$
  其中 \(S\) 是以 \(\gamma\) 为边界的某个曲面，省略号表示非阿贝尔情形下更复杂的路径有序/面积有序项。

这套结构与我们在数值上做的事情是严格对应的：

- 在脚本中，我们离散化路径 \(\gamma\) 为若干 \(\lambda_k\)，在每一步用邻近本征态的重叠矩阵 \(W_k\) 的极分解近似 \(e^{-A_\mu\delta\lambda^\mu}\)；
- 它们的乘积 \(U_\mathrm{disc}=\prod_k U_k\) 在步长 \(\delta\lambda\to 0\) 时趋近于连续的 \(U[\gamma]\)。

**4.3 与 \(R=e^{iH_P}\) / YBE 的对应：小三角形回路与经典 YBE**

- 若 \(F=0\)，联络是平坦的，平行移动只依赖于路径的同伦类，可以很好地给出 braid group / mapping class group 的表示；
- 若 \(F\neq 0\)，不同同伦等价的路径组合之间会有差异，需要更高复杂度的电路来逼近相同 holonomy（这正是我们在 complexity_flatness 实验里看到的）。

把这一 picture 投影到三点/三体子系统上：

- 令 \(H_{12},H_{13},H_{23}\) 是对应三对边/路径段上的生成元；
- 考虑参数空间中一个极小的“二维 simplex”：三条边分别对应沿 \(H_{12},H_{13},H_{23}\) 方向演化一个很小参数 \(\lambda\)；
- 用 \(R_{ab}(\lambda)=e^{i\lambda H_{ab}}\) 表示沿对应边的“离散平行移动元”，则围绕这个三角形有两条自然的走法：
  $$
  U_\triangle^{(L)}(\lambda) = R_{12}(\lambda)R_{13}(\lambda)R_{23}(\lambda),\\
  U_\triangle^{(R)}(\lambda) = R_{23}(\lambda)R_{13}(\lambda)R_{12}(\lambda).
  $$

若常数 YBE 严格成立，则对所有 \(\lambda\) 有 \(U_\triangle^{(L)}=U_\triangle^{(R)}\)。在小 \(\lambda\) 极限，把两边展开到二阶，可得
$$
U_\triangle^{(L)}(\lambda) - U_\triangle^{(R)}(\lambda)
 = -\lambda^2\bigl([H_{12},H_{13}] + [H_{12},H_{23}] + [H_{13},H_{23}]\bigr) + O(\lambda^3).
$$
这说明：

- 经典 YBE 条件
  $$
  [H_{12},H_{13}] + [H_{12},H_{23}] + [H_{13},H_{23}] = 0
  $$
  等价于“在 \(O(\lambda^2)\) 阶上，两条走法围绕三角形给出的 holonomy 一致”；
- 把这个三角形视为一个极小曲面元 \(S\)，区域大小 \(\mathrm{area}\sim\lambda^2\)，上式正是
  $$
  U_\triangle^{(L)} - U_\triangle^{(R)} \propto F_{\mu\nu}\,\mathrm{area}_{\mu\nu} + O(\lambda^3)
  $$
  的算符版展开：**经典 YBE ↔ 嵌在三角形小胞上的曲率张量 \(F\) 在最低阶为零**。

因此：

- 旧的 \(\mathcal M_R^{(\mathrm{YBE})}\) 是“在对易 Pauli slice 上，既满足经典 YBE 又满足所有高阶量子 YBE 约束”的一块完全平坦流形；
- 新的指数 ansatz 中：
  - 经典 YBE 给出的是 **\(H_P\) 空间中“曲率最低阶近似为零”的线性/李代数子流形**；
  - 真正的量子 YBE 则进一步把这条子流形收缩为某些特别的可积模型族（如 Heisenberg/XXZ 类 R‑矩阵）。

从拓扑门/Dehn twist 的角度看：

- 在完全平坦的可积流形上，Dehn twist / braid 的 holonomy 与 YBE 解之间有严格的 SU(2) 共轭关系；
- 一旦偏离这条流形，经典 YBE 和高阶量子 YBE 逐步被破坏，几何上表现为非零曲率，直接后果就是：
  - 几何 Berry holonomy 与“某个简单 R² 组合”之间只剩下近似关系；
  - 需要更深的 LQC+permutation 电路才能逼近同一个 holonomy（与我们在 run_complexity_flatness_scan.py 和 run_interacting_braid_check.py 中观察到的现象一致）。

 

#### 5. 展望：在线性稳定性 / 微扰分析中的角色

引入 \(R=e^{iH_P}\) 之后，自然的下一步是围绕某个“理想可积点” \(H_P^{(0)}\) 做小扰动：
$$
H_P = H_P^{(0)} + \varepsilon V,
$$
其中 \(H_P^{(0)}\) 满足量子 YBE（或至少经典 YBE），\(V\) 是一般的 su(4) 扰动。

- 经典 YBE 提供了一个对 \(V\) 的线性化条件：哪些方向 \(V\) 保持 \(O(\varepsilon^2)\) 曲率为零，哪些方向立即产生非零曲率；
- 在 Majorana/BdG 语言中，这正好对应“哪些扰动只改变有效 (t,\Delta,\mu) 而不破坏拓扑/平坦结构，哪些扰动会引入相互作用/非对易结构，从而拉高实现同一个 Dehn twist 所需的电路复杂度”。

这些将是后续做“线性稳定性 / 微扰分析”和“数值参数扫描”的自然技术工具，kit-exp 就作为这部分的理论预备和记账。

 

#### 6. 沿原流程重走一遍：R=e^{iH_P} 带来的差异与细节

这里按我们在 kit3 / kit4-3 / kit_all 里已经固定下来的“主线流程”，逐步看一下把 R 从线性 Pauli 形式换成指数形式后，会出现哪些新的结构或细节。

**6.1 R 层：从 (a,b,c,d) 到 c_{ij} 的提升与分解**

- 旧流程：
  - 直接把 R 当作 4 维参数空间 \((a,b,c,d)\) 里的点；
  - 物理上，只允许 \(I,\sigma^x\sigma^x,\sigma^y\sigma^y,\sigma^z\sigma^z\) 这 4 个“互相对易的通道”。
- 新 R=e^{iH_P} 流程：
  - 参数空间升级为近似 su(4) 的 15 维 \(\{c_{ij}\}\)（除去整体相位）；
  - 物理上，可以把 \(H_P\) 按 JW/Majorana 映射拆成三块：
    1. 产生最近邻自由 BdG 的“好方向” \(H_P^{(\mathrm{free})}\)（给出有效 \(t,\Delta,\mu\)）；
    2. 产生 quartic 或远程项的“相互作用方向” \(H_P^{(\mathrm{int})}\)；
    3. 与对称性/规范选择相关的“纯规”方向 \(H_P^{(\mathrm{gauge})}\)。

在小耦合极限下，我们可以写
$$
H_P = H_P^{(\mathrm{free})} + H_P^{(\mathrm{int})} + H_P^{(\mathrm{gauge})},
$$
并把原来的 \((a,b,c,d)\) slice 看成是只打开 \(H_P^{(\mathrm{free})}\) 的一个 3~4 维截面；其他 \(c_{ij}\) 方向则对应各种“出平面”的微扰。

**额外细节**：这给 later 的“线性稳定性分析”一个非常自然的分解基底：沿着 \(H_P^{(\mathrm{free})}\) 做扫描保持自由/近自由结构，沿 \(H_P^{(\mathrm{int})}\) 做扫描测试 Dehn twist 与 braid 的稳健性，沿 \(H_P^{(\mathrm{gauge})}\) 把规范等价类整理清楚。

**6.2 R→Kitaev/BdG：有效 (t,\Delta,\mu) 与相互作用通道的区分**

- 旧流程：直接用
  $$
  t = b + c,\quad \Delta = b - c,\quad \mu = 4d + \mu_0
  $$
  把 (a,b,c,d) 映到一条 1D Kitaev 链或一个 2-site BdG 模型上；不显式考虑 quartic/多体项。
- 新流程：先通过 JW/Majorana 把 \(H_P\) 写成
  $$
  H_P = H_\mathrm{quad}(c_{ij}) + H_\mathrm{quartic}(c_{ij}),
  $$
  其中
  - \(H_\mathrm{quad}\) 是 \(\gamma_a\gamma_b\) 型的二次 Majorana；
  - \(H_\mathrm{quartic}\) 是 \(\gamma_a\gamma_b\gamma_c\gamma_d\) 等相互作用项。

在小 \(c_{ij}\) 极限下，**主导的有效 BdG 参数**来自 \(H_\mathrm{quad}\) 部分；线性响应给出
$$
(t,\Delta,\mu) = (t_0,\Delta_0,\mu_0) + L\cdot c_{ij} + O(c^2),
$$
其中 \(L\) 是一个简单的线性映射（在 verify/mapping_from_micro.json、mapping_fit_results.json 里我们已经在特定 slice 上数值确定了一部分）。

**新结论/细节**：

- 我们可以把“保持拓扑相结构不变的参数方向”定义为那些只通过 \(L\) 改变 \(t,\Delta,\mu\)，但在拓扑判据 \(|\mu|<2|t|\) 的意义下仍留在同一区间的 \(c_{ij}\) 方向；
- 真正危险的是那些主要注入 \(H_\mathrm{quartic}\) 的方向：它们在 4‑Majorana / 1D 链中已经通过 run_interacting_braid_check 显示，会逐渐破坏 braid 与 half twist 的精确一致；
- 因此，在 exp 视角下“线性稳定性”自然拆成：
  1. BdG 参数空间内部的稳定性（自由近似内）；
  2. 对相互作用通道的敏感性（超出自由近似）。

**6.3 YBE / 平坦联络：从多项式约束到 classical YBE 约束**

这一部分在上文第 3–4 节已经详细推了，这里只按“原流程节点”总结差异：

- 旧流程：
  - 直接求解 \(R(a,b,c,d)\) 的常数 YBE，得到一个代数子流形 \(\mathcal M_R^{(\mathrm{YBE})}\subset\mathbb C^4\)；
  - 解释为“在这个子流形上，联络曲率 \(F=0\)（或为中心元）”。
- 新流程：
  - 把 \(R=e^{i\lambda H_P}\) 看成小参数 \(\lambda\) 的族，YBE 在 \(O(\lambda^2)\) 阶给出经典 YBE：
    $$
    [H_{12},H_{13}] + [H_{12},H_{23}] + [H_{13},H_{23}] = 0;
    $$
  - 几何上，这是“小三角形回路上的曲率为零”的条件，即联络的最低阶平坦性；
  - 高阶 BCH 展开进一步约束 \(H_P\) 落在真正的可积子流形上。

**新细节**：我们得到一个自然的分层：

1. **严格可积层**：\(H_P\) 满足所有量子 YBE 约束（旧的 \(\mathcal M_R^{(\mathrm{YBE})}\) 就在这一层内的某个 slice）；
2. **classical YBE 层**：只在 \(O(\lambda^2)\) 意义下平坦（局部曲率最低阶为零）；
3. **一般扰动层**：上式明显非零，对应立即产生曲率的方向。

从“原流程”的角度看：这直接告诉我们 **在哪些方向上 holonomy 对路径细节不敏感（平坦/近平坦），在哪些方向上一点点扰动就会改变 braid/Dehn twist 的结果**。

**6.4 4‑Majorana half twist / Dehn twist：精确等式 vs 近似等式**

- 旧流程：
  - 在特定 (a,b,c,d) 点上，2‑site Kitaev 链的 \(H\) 简化为 \(H=(i/2)\gamma_2\gamma_3\)，从而
    $$
    e^{-iH\tau} = e^{(\pi/4)\gamma_2\gamma_3}
    $$
    在 \(\tau=\pi/2\) 时精确成立；
  - Dehn twist 则通过几何路径 H(\phi) 的 Berry holonomy，与 Ising TQFT 的 \(F^{-1}R^2F\) 在 SU(2) 中精确共轭。
- 新流程：
  - 一般的 \(H_P\) 会在 4‑Majorana 子空间里产生额外的 \(\gamma_1\gamma_2,\gamma_3\gamma_4\) 或 quartic 项；
  - 把 \(H_P\) 写成
    $$
    H_P = H_0 + V,\qquad H_0=\tfrac{i}{2}\gamma_2\gamma_3,
    $$
    其中 \(V\) 收集所有偏离“理想 half-twist 生成元”的项（包括额外的二次 Majorana 和 quartic 相互作用）；
  - 单周期时间演化算符可以写成 Dyson 级数
    $$
    U_P(\tau)=e^{-iH_P\tau}
    = e^{-iH_0\tau}\,\mathcal T\exp\Bigl(-i\int_0^\tau dt\,V_I(t)\Bigr),
    $$
    其中 \(V_I(t)=e^{iH_0 t}V e^{-iH_0 t}\) 是交互表象中的扰动；
  - 展开路径有序指数到一阶，得到
    $$
    U_P(\tau)
    = U_0(\tau) - i\int_0^\tau dt\,U_0(\tau-t)V U_0(t) + O(\|V\|^2),
    $$
    这里 \(U_0(\tau)=e^{-iH_0\tau}\) 就是理想 half twist；
  - 用算符范数估计差异：若取 \(\|\cdot\|\) 为任意一致算符范数（例如谱范数），则
    $$
    \|U_P(\tau)-U_0(\tau)\|
    \le \int_0^\tau dt\,\|V\| + O(\|V\|^2)
    = \tau\,\|V\| + O(\|V\|^2),
    $$
    因此 small-V 极限下演化偏差是线性阶的：\(\|U_P-U_0\|=O(\tau\|V\|)\)。

在 4‑Majorana 模型中，我们真正关心的是投影到零模/基态子空间后的逻辑作用：

- 记 \(P_0\) 为 \(H_0\) 在能量 0（一维或二维简并）的投影，则理想 half twist 在逻辑子空间上的作用为
  $$
  U_0^{(\mathrm{log})} = P_0 U_0(\tau) P_0,
  $$
  对于合适选择的 \(\tau\)（如 \(\pi/2\)），这是一个 2×2 的 SU(2) 矩阵；
- 在存在扰动 \(V\) 时，逻辑作用变为
  $$
  U_P^{(\mathrm{log})} = P_0 U_P(\tau) P_0,
  $$
  其与 \(U_0^{(\mathrm{log})}\) 的差异可以估计为
  $$
  \|U_P^{(\mathrm{log})}-U_0^{(\mathrm{log})}\|
  \le \|P_0\|^2\,\|U_P(\tau)-U_0(\tau)\|
  = O(\tau\|V\|).
  $$

配合我们在 verify/run_interacting_braid_check.py 中的数值实验，可以把
$$
\delta_\mathrm{braid}(V) = \max_{\tau\sim\pi/2}\|U_P^{(\mathrm{log})}-U_0^{(\mathrm{log})}\|
$$
视为“half twist 稳健性”的一个定量指标：

- 在 \(V=0\) slice 上，\(\delta_\mathrm{braid}=0\)，half twist 精确实现；
- 对小 \(V\)，有 \(\delta_\mathrm{braid}(V)\approx c\,\|V\|\)；
- 实验上，我们可以通过扫描相互作用强度（如 run_interacting_braid_check 里的 U,\mu）来拟合这一比例常数 c。

**新结论**：

- 旧框架里我们隐含地工作在 \(V=0\) slice 上；现在可以显式把 V 当成扰动，直接用 run_interacting_braid_check 之类脚本测量
  $$
  \|e^{-i(H_0+V)\tau} - e^{-iH_0\tau}\| \sim O(\|V\|),
  $$
  作为“half twist 稳健性”的定量指标；
- 对 Dehn twist 也是类似：Berry holonomy 与 \(F^{-1}R^2F\) 之间的差异可以在 \(H_P\) 空间里看成是沿“非平坦方向 V” 积累起来的曲率效应；
- 这为后面在 \(c_{ij}\) 空间里做系统的“F_\mathrm{Dehn}(c) vs 曲率/相互作用强度”的扫描提供了更精细的物理解释。

**6.5 复杂度 / LQC+permutation 层：与 H_P 空间几何的对应**

- 旧流程中，complexity_flatness 实验是通过在 4‑Majorana 模型里手工加入一个扰动 \(\varepsilon\,(i/2)\gamma_1\gamma_2\)，看 fixed-depth LQC+SWAP 电路逼近 Berry holonomy 的 fidelity \(F(\varepsilon)\) 如何下降；
- 在新视角下，这正好对应沿着某个具体的 \(V\) 方向穿出可积 slice：
  - 若 \(V\) 满足经典 YBE（或与 \(H_0\) 对易），曲率保持很小，\(F(\varepsilon)\) 在一段区间内接近 1；
  - 若 \(V\) 强烈破坏经典 YBE，则曲率迅速变大，\(F(\varepsilon)\) 很快下降，意味着需要增加 LQC 深度才能保持同样精度。

**额外细节 1：F(\varepsilon) 的连续展开与“complexity curvature”**

在 verify/run_lqc_permutation_fit.py 与 verify/run_complexity_flatness_scan.py 里，我们定义了：

- 固定一个几何路径（如 H(\phi) 环路）和一个 shallow ansatz 结构（如两层 SU(2)⊗SU(2)+SWAP）；
- 对每个扰动强度 \(\varepsilon\) 计算几何 Berry holonomy \(U_\mathrm{Berry}(\varepsilon)\)；
- 在给定浅层 ansatz 范围内，寻找最佳逻辑逼近 \(U_\mathrm{LQC}(\varepsilon)\)，并定义(fidelity function)
$$
  F_{\max}(\varepsilon;V)
  = \max_{U\in\mathcal A_\mathrm{shallow}}
    \frac{1}{2}\bigl|\mathrm{Tr}\bigl(U_\mathrm{Berry}^\dagger(\varepsilon) U\bigr)\bigr|
  \in[0,1].
$$

在 small-\(\varepsilon\) 附近，若路径和平行移动是解析的，可以把 \(F_{\max}(\varepsilon;V)\) 展开为
$$
F_{\max}(\varepsilon;V)
 = 1 - \alpha(V)\,\varepsilon^2 + O(\varepsilon^3),
$$
其中线性项通常可以通过适当 re-phase 被消去，二阶系数 \(\alpha(V)\ge0\) 则捕捉“沿 V 方向离开可积 slice 的曲率敏感度”。

这启发我们定义一个经验量：
$$
K_\mathrm{comp}(V) := 2\alpha(V)
 \approx -\left.\frac{d^2}{d\varepsilon^2}F_{\max}(\varepsilon;V)\right|_{\varepsilon=0},
$$
可以粗略理解为“固定浅层 ansatz 深度时，沿 V 方向的 complexity 曲率”：

- 若 V 沿着某个“classical YBE / 近平坦”方向，则 \(K_\mathrm{comp}(V)\) 很小，说明浅层电路在一段扰动区间内仍能高保真实现几何 holonomy；
- 若 V 显著破坏经典 YBE，则 \(K_\mathrm{comp}(V)\) 会较大，对应 F‑曲线很快从 1 下弯，需要加深 ansatz 才能恢复相同精度。

**额外细节 2：与几何曲率 F 的经验对应**

从几何端看，沿某个扰动方向 V 修改 \(H_P\) 会改变 Berry 联络 \(A_\mu\) 及其曲率 \(F_{\mu\nu}\)。在 small-\(\varepsilon\) 极限下，holonomy 的变化可以写成（示意性地）
$$
U_\mathrm{Berry}(\varepsilon)
 = U_\mathrm{Berry}(0)
   \exp\Bigl(-\varepsilon\iint_S \delta F + O(\varepsilon^2)\Bigr),
$$
其中 \(\delta F\) 是沿 V 的曲率变化，\(S\) 是 Dehn-twist 轨道围成的某个代表性曲面。

如果我们在 su(2) 逻辑子空间里把 holonomy 视为
$$
U_\mathrm{Berry}(\varepsilon)
 \sim \exp\bigl(-i\,\vec n(\varepsilon)\cdot\vec\sigma\,\theta(\varepsilon)/2\bigr),
$$
则 \(F_{\max}(\varepsilon;V)\) 的下降主要反映了 \(\vec n,\theta\) 相对 \(\varepsilon=0\) 时的偏移量。简单情形下，这种偏移的二阶系数与 \(\|\delta F\|^2\) 成正比，从而给出一个“complexity 曲率 \(K_\mathrm{comp}(V)\) 与几何曲率 \(\|\delta F\|\)” 之间的经验关系：
$$
K_\mathrm{comp}(V) \propto \|\delta F(V)\|^2 + \cdots.
$$

这并非严格定理，但提供了一张工作图景：

- 经典 YBE 保证在三点/三边单元上的 \(F\) 在最低阶为零，对应 \(K_\mathrm{comp}\) 小；
- 偏离经典 YBE 的 V 方向会产生非零 \(F\)，并在数值上表现为 \(F_\max\) 曲线的快速下弯。

因此，complexity_flatness 类型的实验可以被理解为：用固定深度的 LQC+permutation ansatz，把“几何曲率 F 随 V 的变化”间接投影到一个 operational 量 \(K_\mathrm{comp}(V)\) 上，它衡量的是“在给定资源约束下，holonomy 的可实现性有多敏感”。

**6.6 2D p+ip / honeycomb 层：R 指数视角下的统一**

- 虽然我们的 2D 数值脚本本身是直接从 BdG / Majorana 哈密顿量出发构造的，没有显式用到 R，但在新的 exp 视角下，可以把沿每条键的局域哈密顿量都看成某个 \(H_P\) 的嵌入；
- 1D 中“\(R=e^{iH_P}\) → BdG → 零模 Berry” 的链条，在 2D p+ip 和 honeycomb 中仍然存在，只是：
  - 参数空间 \(\mathcal C\) 现在不仅包含涡旋/vison 的位置，还包含一部分局域 \(H_P\)（相当于耦合强度/各向异性）的变化；
  - 我们可以把先前在 2D 中观察到的 Dehn twist 平台/非平台现象，重新解读为：
    - 平台区：等效地停留在某个“近经典 YBE / 近平坦”的参数区域；
    - 非平台区：穿出了这块区域，曲率和相互作用效应积累起来。

这样，整个“原流程”（1D R→Kitaev→Majorana→Berry / Dehn、以及 2D p+ip 与 honeycomb 的数值）在 R=e^{iH_P} 视角下被统一进一张更大的图：

- 旧的 (a,b,c,d) 只是一条对易、完全平坦的 slice；
- 新的 \(c_{ij}\) 空间给出所有可能的方向，我们可以系统地：
  1. 在其中寻找满足经典/量子 YBE 的“可积层”；
  2. 围绕这些点做小扰动，研究曲率、Dehn twist fidelity 与 LQC 复杂度的响应；
  3. 把 1D、2D（p+ip, honeycomb）的所有实验结果都嵌入到这张“\(H_P\) 空间几何”里来解释。

后续要做的“线性稳定性 / 微扰分析”和“数值参数扫描”，就可以直接在 \(H_P\) 或其投影的 \(c_{ij}\) 空间上设计：选一个理想点 \(H_P^{(0)}\)，选若干代表性的扰动方向 \(V_k\)，数值测量 F_Dehn、曲率指标和 LQC 复杂度随 \(\varepsilon_k\) 的变化，从而在这张统一图上填充更多定量细节。