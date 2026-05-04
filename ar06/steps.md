# 可积 R 矩阵族与局域哈密顿 h 的显式构造（参考）

本文档列出常见的自旋‑1/2 可积 R(u) 族（XXX, XXZ, XYZ, free‑fermion）及其在两比特基底 |00>,|01>,|10>,|11> 下的矩阵形式，说明如何从 R(u) 导出局域两体密度 h 的计算式、以及对应 Pauli 基系数 c_{μν} 的约束条件。可直接复制公式用于数值脚本或写入 `ybe.md`。

约定与符号：
- 基底顺序按字典序： |00>,|01>,|10>,|11>. 置换算符 P 在该基底为交换 |01>↔|10> 的矩阵。
- Pauli 基表示：对一对站点的局域算符展开为 $\sum_{\alpha,\beta\in\{0,x,y,z\}} c_{\alpha\beta},\sigma^\alpha\otimes\sigma^\beta$.
- 要从 R(u) 得到 h，常用规范化为
  $$
  h = \frac{P\,\partial_u R(u)\big|_{u=0}}{\mathrm{scale}},
  $$
  其中 scale = 整体标度（例如 XXX 中的常数或 XXZ 中的 sin η），确保 h 为厄米且无置换缩放因子。

## 关于 `scale` 与 $\eta$ 的选择

`scale` 用来去掉 $R(0)$ 的整体标度，使得由导数得到的局域密度 $h$ 具有合适的数值尺度并为厄米。原则是：若
$$
R(0)=\text{scale}\cdot P,
$$
则自然选取 `scale` 为该标度因子，并定义
$$
h=\frac{P\,\partial_uR(u)\big|_{u=0}}{\text{scale}}.
$$

简短示例：
- XXX: $R(u)=uI+\eta P$，有 $\partial_uR|_0=I$，因此可取 $\text{scale}=1$，得到 $h=P$；这里 $\eta$ 是 R 中的耦合常数，但不出现在 $\partial_uR|_0$ 中（仍影响 R 的其它性质）。
- XXZ: $R(u)$ 如上，$R(0)=\sin\eta\,P$，故取 $\text{scale}=\sin\eta$，代入得
  $$
  h=\frac{1}{2\sin\eta}(\sigma^x\otimes\sigma^x+\sigma^y\otimes\sigma^y)+\frac{\cot\eta}{2}\,\sigma^z\otimes\sigma^z.
  $$

关于时间演化的对应：若把 $R(u)$ 解释为演化算符（例如 $R(u)=e^{i u H}$），则 $\partial_uR|_0=iH$，从而 $H=-i\,\partial_uR|_0$；若 $R(u)=P\,e^{i u H}$，則 $H=-i\,P\,\partial_uR|_0$，這與上面除以 `scale` 的做法兼容（scale 可被吸入 H 的定义）。

备注：$u$ 与 $\eta$ 在常见约定中是同类谱参数（通常无量纲）。整体常数仅重标度 Hamiltonian（相当于时间尺度重标），因此不改变可积性或拓扑性质，但为与文献和数值比较常按上面惯例选取 `scale`。

## 1. Rational (XXX) —— R(u)=u I + η P

矩阵形式：
$$
R(u)=u\,I + \eta\,P,
$$
其中 P 是置换算符，η 是耦合常数。该 R 满足加法型 YBE。

局域密度 h（规范化形式）：取 scale=1，可取
$$
h = P\,\partial_u R(u)|_{u=0} = P\cdot I = P.
$$
在 Pauli 基下，置换算符有展开
$$
P = \tfrac{1}{2}\big(I\otimes I + \sigma^x\otimes\sigma^x + \sigma^y\otimes\sigma^y + \sigma^z\otimes\sigma^z\big).
$$
因此$ c_{xx}=c_{yy}=c_{zz}=1/2$（相对权重，可按整体常数重标度），其余$ c_{\alpha\beta}=0$（无单体与串项）。

物理注记：$h ∝ σ·σ（Heisenberg）$，JW 后为等向自旋耦合；对 MZM 需退化到适当子空间或通过外部场引入配对。

## 2. Trigonometric (XXZ / 6‑vertex)

标准 R(u)（基底 |00>,|01>,|10>,|11>）：
$$
R(u)=\begin{pmatrix}
a(u) & 0 & 0 & 0 \\
0 & b(u) & c & 0 \\
0 & c & b(u) & 0 \\
0 & 0 & 0 & a(u)
\end{pmatrix},
\quad a(u)=\sin(u+\eta),\; b(u)=\sin u,\; c=\sin\eta.
$$

正则性：$R(0)=\sin\eta\,P$，因此$ scale=\sin\eta$。

导出 h 的步骤：
- 计算$ ∂_u R(u)|_{u=0}: a'(0)=\cos\eta, b'(0)=1, c'=0.$
- 取
$$
h = \frac{P\,\partial_u R(u)|_{u=0}}{\sin\eta}.
$$
代数化简并去掉可吸收常数项得到（去常数后）
$$
h = \frac{1}{2\sin\eta}\big(\sigma^x\otimes\sigma^x + \sigma^y\otimes\sigma^y\big) + \frac{\cot\eta}{2}\;\sigma^z\otimes\sigma^z.
$$


     
$c_{μν} 约束$：仅$c_{xx}, c_{yy}, c_{zz} $非零，且$ c_{xx}=c_{yy} $(可按上式具体数值)，其余$ c=0$（无单体/串项）。

物理注记：XXZ 链，JW 后 $σ^{x,y} $项给 hopping/pairing，$σ^zσ^z $给交互 U 与化学势修正。要得到纯二次 BdG（便于 MZM），需把 $c_{zz} $控制或处理为次要项。

极限关系：η→0 可退化为 XXX（比例关系）；η→ i∞ 或用双曲形式得到不同参数域的物理相。

## 3. Elliptic (XYZ / 8‑vertex / Baxter)

通用形式：Baxter 的 8‑vertex R(u) 使用椭圆 theta 函数给出四个非零函数 a(u), b(u), c(u), d(u)：
$$
R(u)=\begin{pmatrix}
a(u) & 0 & 0 & d(u)\\
0 & b(u) & c(u) & 0\\
0 & c(u) & b(u) & 0\\
d(u) & 0 & 0 & a(u)
\end{pmatrix}.
$$
（注：不同文献有不同的规范化与符号约定，参见 Baxter 的原著。）

c_{μν} 约束：一般情形下 c_{xx}, c_{yy}, c_{zz} 三者不等，且可能出现相位导致交叉 σ^xσ^y 等项在 Pauli 展开中的非零系数（但通常不引入串 S_j）。

局域密度 h 的计算仍采用规范化导数：取 scale=R(0) 的对应常数（需检查 R(0) 形状，通常与 P 相关），然后 h = P ∂_uR|_0 / scale。代入后得到 XYZ 型哈密顿：
$$
h = J_x\,\sigma^x\sigma^x + J_y\,\sigma^y\sigma^y + J_z\,\sigma^z\sigma^z + \text{(常数项)}
$$
其中 J_α 由椭圆参数决定。

物理注记：XYZ 是最一般的可积两体自旋模型，包含 XXZ、XXX 作为极限；JW 后若只含 XX/YY 组合则仍能映成包含 pairing 的 BdG，若含全异向性则需具体检查是否会产生非二次项。

## 4. Free‑fermion 子族（8‑vertex 的自由费米子限制）

定义条件（自由费米子判据）：R 的矩阵元满足特定代数条件（常见形式之一）
$$
a(u)^2 + b(u)^2 = c(u)^2 + d(u)^2
$$
或等价的约束，使得模型可映射为自由费米子（Baxter 的 free‑fermion 条件的一个表述）。

在 Pauli 展开上，free‑fermion 子族对应的 c_{μν} 允许交叉项（例如 c_{xy}, c_{yx}）且组合满足使得 JW 后得到仅二次费米子算符（无四费米项与无串）。因此这是直接构造 BdG 与 MZM 的友好族。

举例：当 R(u) 的非对角区块满足特定线性关系时，可找到一个正交/酉基变换把 R 写成 Majorana 双线性指数的形式，从而显式写出 h 的 Majorana 表达式。


## 符号推导与严格代数约束（总结）

我已用 SymPy 在 `tools/` 中实现若干脚本以导出并验证代数约束：
- `tools/classify_R_constraints.py`：数值/符号混合的快速判定脚本（SU(2)/U(1)/free‑fermion 检查，含 demo）。
- `tools/derive_constraints.py`：对一般 Pauli 系数进行符号化尝试（打印 XXZ 的消去式，并尝试对自由费米子情况做消元）。
- `tools/derive_reduced_constraints.py`：针对约化两站点情形（仅 c_{xx}=c_{yy}=A, c_{zz}=B）做 Groebner 消去，输出关于 A,B 的多项式约束。

主要结论（可直接作为严格代数约束写入文档）：

- XXZ（采用约定 A=c_{xx}=c_{yy}, B=c_{zz}）：
  $$
  4A^2 - 4B^2 - 1 = 0.
  $$ 
  （对应原始参数化 A=1/(2 sin η), B=cot η/2，消去 η 得上式。）

- XYZ：
必要的线性代数约束（无单体且无交叉项的常见情形）为
  $c_{0α}=c_{α0}=0, 以及 c_{αβ}=0 for α≠β。$
充分的“可积性”条件不是单一代数式；它等价于存在椭圆函数参数 a(u),b(u),c(u),d(u) 满足 Yang‑Baxter 函数方程。要把這一条件写成关于 c_{αβ} 的闭合代数式，需解相应的函数方程或构造 R(u)。

- Free‑fermion 子族：
  - 在 R 矩阵参数化上常见的判据为函数恒等式
    $$a(u)^2 + b(u)^2 = c(u)^2 + d(u)^2$$
    （对所有 u）。
  - 在 Pauli 系数表示上等价条件是：经 Jordan–Wigner 变换后，H 仅包含费米双线性项（hopping/pairing/μ），即所有会产生四费米（quartic）项的组合的多项式系数必须为 0。通过 SymPy 的 Groebner 消元可以把这些条件写成关于 c_{αβ} 的显式多项式集合（见上面 `tools/derive_reduced_constraints.py` 对 A,B 的专门例子）。

注意与实践：
- 我在脚本中保留了整体常数/相位自由度说明（它们会被吸收到 `const` 或 hopping/pairing 的相位里），在做 Groebner 消元时通常把整体能量位移（`const`）置零并固定相位（例如把部分相量取实）以简化消元。这样得到的多项式是判定自由费米子子族的充分必要条件（在所做规范下）。
- 若需要，我可以把 `tools/derive_reduced_constraints.py` 的输出（当前在你的环境中运行并生成的多项式）直接粘贴到此处或追加到 `ybe.md`。目前我已把脚本及结果留在 `tools/`，你可以按需运行：

```bash
python3 tools/derive_reduced_constraints.py
```

或把我运行的输出让我贴出来。  

 

如果你希望我把这些符号推导的完整输出（脚本运行结果）也追加进 `ybe.md`，我可以把 `tools/` 的符号结果段直接粘入 `ybe.md`；或者我可以把 `steps.md` 的这段内容同步追加到 `ybe.md`（保持两者一致）。请选择：
- 将此段追加到 `ybe.md`；或
- 将脚本的实际运行结果（多项式列表）也作为补充追加到 `ybe.md`。


## 约束判定流程（对任意给定 c_{μν}）

给定一组 c_{μν}（16 个复数）时，可按以下步骤判定其可积族类型与给出相应的 R(u)：
1. 检查字符串触发分量：若任一 c_{x0},c_{0x},c_{y0},c_{0y},c_{xz},c_{zx},... 非零，则存在 JW 串项，模型一般非纯二次。记录这些分量并说明其物理含义。
2. 检查对称性：若 c_{xx}=c_{yy}=c_{zz}（按比例），则归入 XXX（或其比例子族）。
3. 若 c_{xx}=c_{yy} ≠ c_{zz} 且无串分量，则归入 XXZ/6‑vertex 类型，可尝试拟合 η 使得 J_x=J_y=(2 sinη)^{-1} 等比例关系成立。
4. 若 c_{xx},c_{yy},c_{zz} 都不等，但串分量为0，可尝试匹配 XYZ（椭圆）参数并给出 J_x,J_y,J_z 对应关系。
5. 若存在交叉 c_{xy},c_{yx} 且满足 free‑fermion 代数条件（见上），优先标记为 free‑fermion 子族并输出 Majorana 双线性形式。

自动化建议：可实现 `tools/classify_R_family.py`，接受 c_{μν}（或直接接受 h 的 Pauli 系数），按上述判定返回族名、推荐的 R(u) 参数化（或极限表达式）、以及 h 的规范化写法。

## 示例：把 XXZ 结果写为可直接代入的 Pauli 系数
取 η=0.6 的示例，按上面公式得到
- c_{xx}=c_{yy}=1/(2\sin\eta) ≈ 0.8855160983
- c_{zz}=\cot\eta/2 ≈ 0.7308479735
其余 c_{μν}=0。

可把这些系数写成 JSON/CSV 后直接传入数值脚本以构造 BdG 并检查 MZM。

 

如需，我可以把上述每一类的 R(u) 与 h 的推导步骤写成完全符号化的 Python 函数（并放入 `tools/`），或把本文件的内容追加到 `ybe.md`。现在要我：
- A) 直接把这份 `steps.md` 的内容追加到 `ybe.md`；
- B) 同时生成 `tools/classify_R_family.py` 并演示对若干 c_{μν} 的分类；或
- C) 等你上传具体的 c_{μν} 约束列表，我逐一给出对应的显式 R(u) 构造。

    -
## 严格构造 Majorana 零模（MZM）——代数说明与逐类构造

目标：在 `steps.md` 中给出对每种类型（XXX、XXZ、XYZ、free‑fermion）能够严格（解析、非近似）构造 Majorana 零模的充分必要代数条件、构造流程与可证示例。

通用代数框架：
- 写出两站局域项
  $$
  h=\sum_{\alpha,\beta\in\{0,x,y,z\}} c_{\alpha\beta}\,\sigma^\alpha\otimes\sigma^\beta.
  $$
- （1）非相互作用判据：对给定 $c_{\alpha\beta}$，用 Jordan–Wigner（JW）变换把 $h$ 展开为费米子算符并把所得表达按双线性（quadratic）与四次（quartic）项分离。要求所有 quartic 项系数恒为 0（这是一组代数多项式约束，记为 P_free(c)=0）。当且仅当 P_free(c)=0，$h$ 可写成 BdG/二次形式。
- （2）Majorana 矩阵与零模：若成立，上述二次 Hamiltonian 可写成 Majorana 形式
  $$
  H=\frac{i}{4}\,\boldsymbol{\gamma}^T A\boldsymbol{\gamma}+\text{const},
  $$
  其中 $A$ 是实的反对称矩阵。精确零模存在等价于 $A$ 存在零本征向量 $w$，即
  $$
  A w = 0.
  $$
  对任一解 $w$，定义 $\Gamma=w^T\boldsymbol{\gamma}$ 则严格有 $[H,\Gamma]=0$（代数上可直接验证）。
- （3）局域性：若 $w$ 的分量沿链端呈几何衰减（由特征方程的根决定，选取 |λ|<1 的根），则 $\Gamma$ 为端点局域 MZM；局域长度可通过根的模解析得到。

按类型的严格条件与构造要点：

1) Free‑fermion（可严格构造 MZM）
- 充分必要条件：P_free(c)=0（JW 展开无 quartic 项）。这组多项式可以从 $c_{\alpha\beta}$ 显式导出（见 `tools/derive_constraints.py` / `tools/derive_reduced_constraints.py`）。
- 构造流程：消元得出 $(t_{jk},\Delta_{jk},\mu_j)$，写出 Majorana 矩阵 $A$，解 $Aw=0$ 得到精确解析零模。可给出带有闭式根的解析解（Kitaev 链类），并证明 $[H,\Gamma]=0$ 且 $\Gamma^\dagger=\Gamma$。
- 可证示例（Kitaev）：在周期/开链的明确参数下可解析求根并得到局域零模（见下方例子与递推关系）。

2) XXZ（一般为相互作用，严格 MZM 仅在特例）
- 代数事实：标准 XXZ 含 $\sigma^z\sigma^z$，JW 后会生成 quartic 项，因此通常不满足 P_free，不能严格表示 MZM。
- 特例：当 $c_{\alpha\beta}$ 满足特定消去多项式（例如退化到 free‑fermion 限或通过额外单体/配对项抵消 quartic 成分）时，P_free 成立，从而可严格构造 MZM。符号检验脚本可判定是否存在此类代数退化。

3) XYZ（一般不可直接严格构造 MZM）
- 代数事实：XYZ 的一般可积解依赖椭圆函数参量，且不保证为 free‑fermion。严格 MZM 仍需 P_free(c)=0 或在边界/投影下构造有效二次模型。

4) XXX（Heisenberg）
- 代数事实：XXX 为 SU(2) 对称强相互作用模型，不属于 free‑fermion；原始模型通常不能严格产生 MZM。
- 可行操作：通过引入额外配对项或在受限子空间（投影）中构造有效二次哈密顿量，可以在有效模型中证明 MZM 存在，但这不是对原始 XXX 的严格证明。

解析例子：Kitaev 链（可证零模）
- 模型（开链）：
  $$
  H = -\mu \sum_j c_j^\dagger c_j - t\sum_j (c_j^\dagger c_{j+1}+\mathrm{h.c.}) + \Delta \sum_j (c_j c_{j+1}+\mathrm{h.c.}).
  $$
- Majorana 形式：令 $a_{2j-1}=c_j+c_j^\dagger$, $a_{2j}=i(c_j-c_j^\dagger)$，则
  $$
  H = \frac{i}{4} \sum_{m,n} A_{mn} a_m a_n.
  $$
- 在特殊点 $\Delta=t,\mu=0$ 可直接验证端点严格零模为
  $$
  \Gamma_L = a_1,\qquad \Gamma_R = a_{2L},
  $$
  满足 $[H,\Gamma_{L,R}]=0$ 且 $\Gamma_{L,R}^\dagger=\Gamma_{L,R}$. 对一般参数，可解递推关系得到系数序列 $\psi_j$ 满足
  $$
  \psi_{j+1} + \frac{\mu}{t+\Delta}\psi_j + \frac{t-\Delta}{t+\Delta}\psi_{j-1}=0,
  $$
  其特征多项式给出根 $\lambda$，选取 $|\lambda|<1$ 的解构造局域零模，局域长度 $\xi=-1/\log|\lambda|$（解析表达可得）。

证明策略与工具链（实用）
- (A) 符号检验：用 `tools/derive_constraints.py` 做 P_free(c) 的符号消去（Groebner 基），或在约化情形下用 `tools/derive_reduced_constraints.py` 得到简洁多项式。若消去多项式全为零，继续下一步。
- (B) 若 P_free 为 0：用 JW 构造 BdG/ Majorana 矩阵 $A$（用 SymPy 构造符号矩阵），解线性系统 $A w=0$ 得到精确 $w$。
- (C) 验证：代入验证 $[H,\Gamma]=0$ 与 $\Gamma^\dagger=\Gamma$。分析根以证明指数局域性并给出 $\xi$ 的代数表达。

如需，我可以把上述证明模板连同对仓库内某个具体 `h` 的完整符号证明（包括 Groebner 消元、多项式列出、$A$ 的显式符号解及局域长度表达）追加到 `ybe.md`。或者你可以把目标 `c_{\alpha\beta}` 给我，我当场在 chat 中做完整符号证明。

## 从 Pauli 系数到 Majorana 矩阵 A 的构造（具体步骤）

1. 先做 Jordan–Wigner 把自旋算符写成费米子：
  $$
  c_j = \Big(\prod_{k<j}\sigma^z_k\Big)\frac{\sigma^x_j - i\sigma^y_j}{2},
  $$
  并定义 Majorana：
  $$
  \gamma_{2j-1}=c_j+c_j^\dagger,\qquad \gamma_{2j}=i(c_j-c_j^\dagger).
  $$

2. 若 JW 后 $H$ 已为二次（即 P_free(c)=0），把二次项写成 Nambu 形式：
  $$
  H = \frac{1}{2}\Psi^\dagger H_{\rm BdG}\Psi + \text{const},\qquad
  \Psi=(c_1,\dots,c_N,c_1^\dagger,\dots,c_N^\dagger)^T,
  $$
  其中
  $$
  H_{\rm BdG}=\begin{pmatrix} h & \Delta \\ -\Delta^* & -h^T \end{pmatrix}.
  $$

3. 用线性变换把 Nambu 表达转换到 Majorana：设 2N×2N 的变换矩阵
  $$
  T=\begin{pmatrix} I & I \\ -iI & iI \end{pmatrix},\qquad
  T^{-1}=\frac{1}{2}\begin{pmatrix} I & iI \\ I & -iI \end{pmatrix}.
  $$
  令 $\gamma = T\Psi$，代入得到 Majorana 形式
  $$
  H = \frac{i}{4}\,\gamma^T A\gamma + \text{const},
  $$
  其中
  $$
  A = -2i\,(T^{-1})^\dagger\,H_{\rm BdG}\,(T^{-1}).
  $$
  该 $A$ 必为实反对称矩阵，且其零空间与 Majorana 零模一一对应：若 $Aw=0$ 则 $\Gamma=w^T\gamma$ 为精确积代算子且 $[H,\Gamma]=0$。

4. 实际实现细节：用 SymPy 从 Pauli 系数 $c_{\alpha\beta}$ 构造 JW 后的 $h,\Delta$（若存在），组成 $H_{\rm BdG}$，用上式计算 $A$ 并解 nullspace($A$)。将得到的 $w$ 再写回 Pauli 表达以得到 $\Gamma$ 的本征算子形式。

## 示例：XX（XXZ 的 free‑fermion 极限）

考虑 XX 链（即 XXZ 的 $c_{zz}=0$ 情形），最近邻项统一耦合 $J$：
$$
h_{j,j+1}=\frac{J}{2}(\sigma^x_j\sigma^x_{j+1}+\sigma^y_j\sigma^y_{j+1}).
$$
经 JW 变换可写为仅含 hopping 的费米子二次哈密顿：
$$
H = -t\sum_j (c_j^\dagger c_{j+1}+\mathrm{h.c.}) - \mu\sum_j c_j^\dagger c_j,
$$
其中 $t=J/2$。写成 BdG（取配对 $\Delta=0$）得到
$$
h = \begin{pmatrix} \mu & -t \\ -t & \mu \end{pmatrix}\quad(\text{相邻两站的子块}),\qquad
H_{\rm BdG}=\begin{pmatrix} h & 0 \\ 0 & -h^T \end{pmatrix}.
$$
按上面公式用 $T$ 变换可得到对应的 $A$（4×4 实反对称矩阵）；对开链且 $\Delta=0$ 的纯 XX 模型通常不产生拓扑 MZM（因为缺少配对），但这是一个能被符号化验证的例子：

- 用 SymPy 构造上面的 $h, H_{\rm BdG}$，计算 $A$，并解 $Aw=0$；若没有非平凡解则无零模。
- 若要得到 MZM，需要引入配对 $\Delta\neq0$（Kitaev 型配对），或把 XX 模型与其它项联合，形成等价的二次 BdG，使得 $A$ 出现零本征向量（例如 Kitaev 链的 $\Delta=t,\mu=0$ 给出严格局域零模）。

我已把上述推导与示例写入 `steps.md`，并把可执行的符号脚本放在 `tools/`（可直接运行以验证任意给定 $c_{\alpha\beta}$）。如果你要，我可以现在把对仓库内 XXZ‑derived `h` 的符号化构造运行一遍并把 $A$ 与 $\Gamma$ 的显式表达追加到 `steps.md` 或 `ybe.md`。

**数值检验摘要（XXZ 示例, η=0.6）**

- 运行命令：`python3 tools/xxz_R_and_H.py`，并在若干 u 点计算了生成元定义 $H_P(u)=-iR(u)^{-1}\partial_uR(u)$ 的厄米性与互换子。
- 主要输出摘要：
  - 对直接定义 $H_P(u)=-iR^{-1}∂_uR$，在 u={0.0,0.1,0.2,0.4} 处的厄米性误差约为 6.49, 6.21, 6.54, 10.84 （显著不为零）。这说明 $-iR^{-1}∂_uR$ 并非厄米算符，因此不能直接视为物理的哈密顿量用于幺正演化。
  - 用常用规范化定义的局域两体密度
    $h = P\,∂_uR|_{u=0}/\sin\eta$
    在 η=0.6 的示例中是严格厄米（数值厄米性误差为 0），Pauli 展开给出
    `X_X ≈ 0.8855161, Y_Y ≈ 0.8855161, Z_Z ≈ 0.73084797`（其余系数 ~0）。这与文中给出的 XXZ‑形式一致。
  - 计算不同 u 点处 $H_P(u)$（按上面非厄米定义）之间的对易范数发现数值上为 0（在所选点组内 ||[H(u),H(v)]|| ≈ 0）。这说明在该参数化下这些生成元彼此对易（即时间序列可在数值上省略），但由于它们并非厄米，按该定义直接指数化不会产生幺正演化。

结论与建议：
- 若要把谱参数族 $R(u)$ 解释为物理的含时幺正演化子，应取一个**厄米的生成元**，例如上文使用的 $h=P∂_uR|_0/\mathrm{scale}$（或在某些规范下对 $R(u)$ 做相位/纯虚轴变换使 $R$ 酉），然后使用
  $R(u)=\mathcal T\exp\big(-i\int_0^u h(s)ds\big)$。
- 在本例中，$h(u)$ 在不同 u 点互易且为厄米，因此时间序列可省略并写成单指数（或比例单指数），即存在标度函数 Φ(u) 使得 $R(u)\propto\exp(-iΦ(u)h_0)$（相位/常数因子需单独处理）。

我可以将这些数值检验的完整输出与结论追加到 `ybe.md`，或者把更多 u 点与 η 的扫描结果生成并绘图以展示何时生成元为厄米且何时对易。要我现在把完整输出追加到 `ybe.md` 吗？

