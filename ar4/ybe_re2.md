对角形式R

$$
R_{i,i+1} = exp(i \sum_{u=0}^3 c_{\mu \mu} \sigma^\mu_{i} \otimes \sigma^\mu_{i+1} ) \\
\sigma^0 = I, c_{\mu \mu} \in C
$$


**从指数形式直接因式分解得到的约束**

记每条键的局域生成元可写为
$$
H_{ab}=\phi_{ab}\,I\otimes I + H'_{ab},
\qquad R_{ab}=e^{iH_{ab}}=e^{i\phi_{ab}}\,\widetilde R_{ab},\qquad \widetilde R_{ab}=e^{iH'_{ab}}.
$$

下面给出两种常见 YBE 形式下相位项的严格推导与结论：

1) Braid（交换子）形式

代数等式为
$$
R_{12}R_{23}R_{12}=R_{23}R_{12}R_{23}.
$$
代入上式得
$$
e^{i(\phi_{12}+\phi_{23}+\phi_{12})}\,\widetilde R_{12}\widetilde R_{23}\widetilde R_{12}
=e^{i(\phi_{23}+\phi_{12}+\phi_{23})}\,\widetilde R_{23}\widetilde R_{12}\widetilde R_{23}.
$$
整理相位因子得
$$
e^{i(\phi_{12}-\phi_{23})}\;\widetilde R_{12}\widetilde R_{23}\widetilde R_{12}=\widetilde R_{23}\widetilde R_{12}\widetilde R_{23}.\tag{★}
$$
因此：

若所有鍵上使用相同無相位矩陣 $\widetilde R$（即 $\widetilde R_{12}=\widetilde R_{23}=\widetilde R$），(★) 要求
$$
e^{i(\phi_{12}-\phi_{23})}=1\quad\Longrightarrow\quad\phi_{12}\equiv\phi_{23}\pmod{2\pi}.
$$ 

若允許 $\widetilde R_{12}\neq\widetilde R_{23}$，則 (★) 成為一個實際的矩陣等式；要么相位因子為 1，要么矩陣產物需恰好成比例（比例常數為該相位），這是非常特殊的補償情形，需逐點檢驗。

結論（對 braid 形式）：在常見的「同一常數 $R$ 被復用在各鍵」的約定下，相鄰鍵的 $I\otimes I$ 相位必須一致；否則必須檢查是否存在特殊的 $\widetilde R_{ab}$ 補償該相位。

2) 完整三站 YBE 形式

考慮完整 YBE
$$
R_{12}R_{13}R_{23}=R_{23}R_{13}R_{12}.
$$
代入相位分解後，兩邊的相位分別為 $e^{i(\phi_{12}+\phi_{13}+\phi_{23})}$ 與 $e^{i(\phi_{23}+\phi_{13}+\phi_{12})}$，恰好相同，因而整體相位抵消，不導致對 $\phi_{ab}$ 的附加代數約束。

總結：完整 YBE 對每鍵的整體常數相位不施加限制；而在 braid 形式並且採用相同 $\widetilde R$ 的情形下，必有相鄰鍵相位相等。若需保留每鍵獨立的 $a_{00}^{(j)}$，在檢驗 YBE 時應把它們作為獨立變量並按 (★) 檢查相位條件。

我们对 $R=e^{i(c_{00}I+b_xS_x+b_yS_y+b_zS_z)}$ 直接计算 $R_{12}R_{23}R_{12}-R_{23}R_{12}R_{23}$ 并对矩阵元做因式分解，可将所有非平凡约束写为若干可识别的因子相乘为零（去掉不可约的全局相位因子）。 用记号
$$
X=e^{2i b_x},\quad Y=e^{2i b_y},\quad Z=e^{2i b_z},
$$
这些因子的主要类型为：

- $e^{2i b_x}-e^{2i b_y}=0$，等价于 $\sin(b_x-b_y)=0$（即 $b_x\equiv b_y\pmod\pi$）；
- $e^{2i(b_x+b_y)}-1=0$，等价于 $\sin(b_x+b_y)=0$（即 $b_x\equiv- b_y\pmod\\pi$）；
- $F_-(X,Y,Z)$:= $XY - XZ - YZ - 1 = 0$；
- $F_+(X,Y,Z)$:= $XY + XZ + YZ - 1 = 0$.

实际得到的 6 个非平凡约束可被因式为上述类型的乘积（并包含对 $x,y,z$ 的循环置换）。因此，指数形式的常数 YBE 的解集由这些因子的零点并集生成：对于每一个因式化出的乘积因子，至少有一个因子必须为零。

直接结论（用于文档中快速判别）：

- 离散/相等分支：任意成对满足 $b_\alpha\equiv\pm b_\beta\pmod\pi$ 的情况（来自第一类、第二类因子），在许多情形下即可使所有约束成立；其中包含非退化连续族 $b_x\equiv b_y\equiv b_z\pmod\pi$。 

- 代数簇分支：满足 $F_-(X,Y,Z)=0$ 或 $F_+(X,Y,Z)=0$（及其循环置换）的连续代数子簇，给出除简单相等关系外的其他连续解（允许复值 $b_\alpha$）。

- 全局相位 $c_{00}$ 只以不可约乘子 $e^{3 i c_{00}}$ 的形式出现在因式中，不影响“哪个因子为零”的判据（因此不限制解集结构）。

说明：上述为从矩阵代数（指数形式）直接得到的因式化摘要；如果需要，我可以把完整的符号化输出（每个约束的因子化形式）和对 $F_\pm$ 的进一步代数分析（参数化/样例）作为附录文件 `scripts/ybe_exp_factorization.txt` 并在此处插入链接。

**附录 B：原子因子与最小分支枚举（自指数形式）**

我们将 $u=e^{i b_x},\ v=e^{i b_y},\ w=e^{i b_z}$，并在代数变量 $u,v,w$ 下对约束多项式做因式分解与枚举。得到的原子因子（用下标索引）为：

- 0: $w$  （对应 $e^{i b_z}=0$，在复参数意义上表明取向使该因子为零）
- 1: $u-v$  （对应 $e^{i b_x}=e^{i b_y}$，即 $b_x\equiv b_y\pmod{2\pi}$）
- 2: $u+v$  （对应 $e^{i b_x}=-e^{i b_y}$，即 $b_x\equiv b_y+\pi\pmod{2\pi}$）
- 3: $u^2 v^2 - u^2 w^2 - v^2 w^2 - 1$ （即 $XY - XZ - YZ - 1$，与文中 $F_-$ 对应，$X=u^2$ 等）
- 4: $u^2 v^2 + u^2 w^2 + v^2 w^2 - 1$ （即 $XY + XZ + YZ - 1$，与 $F_+$ 对应）
- 5: （高次对称多项式，来自约束的高次因子，见 `scripts/exp_minimal_branches.txt`）
- 6: $u v - 1$（即 $e^{2i b_x}e^{2i b_y}=1$, 等价于 $b_x\equiv -b_y\pmod\pi$）
- 7: $u v + 1$（等价于 $b_x\equiv -b_y+\tfrac\pi2$ 等更细的相位关系）
- 8,9: 另外两类二次形式（见文件）

脚本枚举得到的“最小分支组合”列出在 `scripts/exp_minimal_branches.txt` 中；其中若干代表性组合（以索引表示）为：

- (0,3) : $w=0$ 与 $F_-=0$ 的组合；
- (0,6) : $w=0$ 与 $uv-1=0$ 的组合；
- (1,8) : $u-v=0$（即 $b_x\equiv b_y$）且第 8 类二次因子为 0；
- (3,4) : 同时满足 $F_-=0$ 与 $F_+=0$（特殊代数交点）；
- (6,7) : $uv-1=0$ 和 $uv+1=0$ 的交集（通常无实解，但允许复解）。

完整的原子因子列表与 34 个最小组合已写入：

- `scripts/exp_minimal_branches.txt`（包含所有原子因子表达式与最小组合索引）。

如何使用这些结果：

- 若要构造数值样例验证某一最小组合，把对应的因子翻回 $u=e^{ib_x}$ 形式，解出 $b_x,b_y,b_z$（取对数并注意选取对数分支），代入三体残差矩阵验证 $\|B\|\approx0$；脚本 `scripts/generate_validate_samples.py` 和 `scripts/validate_factors.py` 可用于这一步。 
- 若需我对某些组合给出显式参数化（实解或典型复解），我可以为指定组合自动求解并附上示例值与残差。

我已将完整的枚举输出存档并在仓库中可用；要我现在为若干你关心的分支生成具体的复解示例并把结果追加到本文档吗？

---

**关于因子 8 与 9（明确定义与恢复 $b_z$ 的公式）**

令 $X=e^{2i b_x},\ Y=e^{2i b_y},\ Z=e^{2i b_z}$. 因子 8 与 9 在 $u=e^{ib_x},v=e^{ib_y},w=e^{ib_z}$ 的表示为
$$
\begin{align*}
	ext{(8)}&:\;u^2v^2 - u^2 w^2 + v^2 w^2 + 1 = 0,\\
	ext{(9)}&:\;u^2v^2 + u^2 w^2 - v^2 w^2 + 1 = 0.
\end{align*}
$$
把它们转换为 $X,Y,Z$ 后分别等价于
$$
\begin{align*}
	ext{(8)}&:\;XY - XZ + YZ + 1 = 0 \quad\Longrightarrow\quad Z = \frac{XY+1}{X-Y},\\[4pt]
	ext{(9)}&:\;XY + XZ - YZ + 1 = 0 \quad\Longrightarrow\quad Z = -\frac{XY+1}{X-Y}.
\end{align*}
$$

因此，对于任意给定 $b_x,b_y$（且 $X\neq Y$），因子 8/9 给出 $Z=e^{2ib_z}$ 的代数有理表达式，进而可以写成
$$b_z = -\tfrac{i}{2}\log Z$$
（对数需固定分支；$b_z$ 在 $\pi\mathbb Z$ 上有同余性）。

要点与注意事项：
- 当 $X=Y$（即 $b_x\equiv b_y\pmod\pi$）时分母为 0，应退化到“相等分支”处理；
- 对数的取值会影响 $b_z$ 的具体实/虚部分布，数值验证时应选择合适的对数分支并检验 $\|B\|$；
- 代数上，8 与 9 对应的 $Z$ 值互为相反号，因此属于同一类代数簇的两条分支（允许复解）。

示例（可用脚本检验）：取 $b_x=0.3,\ b_y=0.6$，计算 $X,Y$ 并按上式求得 $Z_8,Z_9$，然后取 $b_{z,8},b_{z,9}=-\tfrac{i}{2}\log Z_{8,9}$，最后代入三体残差矩阵检验 $\|B\|\approx0$。

（完整原子因子与最小组合见 `scripts/exp_minimal_branches.txt`；需要我为上面示例生成并插入数值验证结果吗？）

---


**附录：映射到 Kitaev/BdG 参数（总结）**

我用脚本把二格 Pauli 系数 $c_{\mu\mu}$ 映射到 Majorana 双线性，再把 Majorana 二次型解释为常见的 Kitaev/BdG 参数。约定为：
- 动力学项记作 $H_{\rm kin}=-t(c_1^\dagger c_2+\mathrm{h.c.})$；配对项记作 $H_{\rm pair}=\Delta\,c_1 c_2 + \mathrm{h.c.}$；站内化学势用 $\mu$ 表示。

基于 `tools/derive_mapping.py` 与 `tools/map_to_kitaev.py` 的符号化输出（见 `results/kitaev_mapping.txt`），在当前对角 Pauli 子空间下得到的解析表达为：

- $t_{\rm eff} = -c_{xx} - c_{yy}$
- $\Delta_{\rm eff} = c_{xx} - c_{yy} + i\,(c_{xy}+c_{yx})$
- $\mu_{\rm eff} = c_{00}$

物理含义与拓扑判据（常见约定）：
- 若以上述约定，则周期链的简单拓扑判据为 $|\mu_{\rm eff}|<2|t_{\rm eff}|$ 且 $\Delta_{\rm eff}\neq0$，则开链应出现拓扑 MZM（端点零模，能量裂分随链长指数小）。
- 对应到 $c$ 参数：大致条件为 $|c_{00}|<2|c_{xx}+c_{yy}|$ 且 $c_{xx}-c_{yy}+i(c_{xy}+c_{yx})\neq0$。

注意事项：
- 映射的符号/因子依赖于 BdG 约定，上述不等式使用的是此处的约定（绝对值判据可移植）。
- 若加入更多非对角 Pauli 项或额外的 `I\otimes I` 常数，$\mu_{\rm eff}$ 将被修改；若有四费米交互（如显式的 c_{zz} 交互项），需要用 mean‑field 或小体系 ED 量化其对 $(t,\Delta,\mu)$ 的重整化。

脚本与输出：
- 生成脚本：`tools/derive_mapping.py`, `tools/map_to_kitaev.py`。
- 符号化输出文件：[results/kitaev_mapping.txt](results/kitaev_mapping.txt)

下一步建议：基于上述解析映射，我可以（1）符号化导出 Pfaffian/winding 的 c‑参数条件，或（2）挑若干代表性 $c$ 点做开链数值 BdG（$L$ 可变量），检验最低能随 $L$ 的缩放与端态波函数以区分 MZM 与 ABS。请选择你希望的下一步。

**Extended Mapping (c -> Kitaev/BdG)**

下面为自动生成的扩展映射摘要（详见 [results/extended_mapping.txt](results/extended_mapping.txt#L1-L100)）：

Extended mapping (conventions follow R22.md)

Symbols: c_xx, c_yy, c_xy, c_yx, c_x0, c_y0, c_z0, c_0x, c_0y, c_0z, c_00, c_zz, c_xz, c_yz, c_zx, c_zy

Kitaev/BdG parameters:
	t (real symmetric hopping) = c_xx + c_yy
	Re(Delta) = c_xx - c_yy
	Im(Delta) = c_xy + c_yx  (pairing imaginary part)
	=> Delta (complex) = (c_xx - c_yy) + i*(c_xy + c_yx)
	mu_site (linear density term per site) = 2*(c_z0 + c_0z) - 4*c_zz
	U (nearest-neighbor density-density) = 4*c_zz
	E_bond (per-bond constant) = c_00 - c_z0 - c_0z + c_zz

Chiral/antisymmetric contributions:
	Imaginary/antisymmetric hopping component ~ c_xy - c_yx (generates i*(c1d c2 - c2d c1) like terms)
	Imaginary/antisymmetric pairing component already included in Im(Delta) via c_xy+c_yx (symmetric imaginary)

String / nonlocal contributions (require Jordan-Wigner string S_j):
	These arise from single-site x/y terms multiplied by neighbor z: c_x0,c_y0,c_0x,c_0y and c_xz,c_yz,c_zx,c_zy
	They map to nonlocal operators S_j (c_j + c_j^\\dagger) (2n_{j+1}-1) etc and are NOT captured by simple local (t,Delta,mu,U) parameters.

Mapping matrix M (R22 summary) for p=[t,Delta,mu_site,U,E_bond]^T and c vector order [c_xx,c_yy,c_xy,c_yx,c_zz,c_z0,c_0z,c_00]^T:
	M = [[1,1,0,0,0,0,0,0], [1,-1,0,0,0,0,0,0], [0,0,0,0,-4,2,2,0], [0,0,0,0,4,0,0,0], [0,0,0,0,1,-1,-1,1]]
	(see R22.md for derivation and notes)

Notes:
 - Sign conventions for t depend on Hamiltonian sign: if H_kin = -t(c1^\\dagger c2 + h.c.) then the scalar t above corresponds to + (c_xx + c_yy) in the expansion; some documents absorb a minus sign into t. Be consistent.
 - Delta here is complex; both Re and Im parts come from c_xx,c_yy and c_xy,c_yx respectively.
 - mu_site excludes any global I\\otimes I constant c_00 which contributes to bond energy E_bond; include c_00 if you want a chemical potential shift.


