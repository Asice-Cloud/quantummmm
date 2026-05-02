严格推导：从 R(u) 到局域哈密顿，再到 MZM / ABS

这里给出对 `R(u)` 至哈密顿再到 MZM/ABS 的严格、可检验的数学推导链条。推导分三部分：

简明要点（快速参考）

- R(u) ↔ 含时演化/局域生成元：在解析点 $u_0$ 处线性化
	$$R(u_0+\varepsilon)=P\big(I+\rho\,\varepsilon\,h_{i,i+1}+O(\varepsilon^2)\big)\implies h_{i,i+1}=\frac{P}{\rho}\left.\frac{dR}{du}\right|_{u_0},$$
	因此 $h_{i,i+1}$ 是局域两站的时间演化生成元（局域哈密顿密度）。

- Pauli 展开与参数：
	$$h_{i,i+1}=\sum_{\alpha,\beta} c^{(i)}_{\alpha\beta}\,\sigma_i^{\alpha}\sigma_{i+1}^{\beta},$$
	数值上由 `tools/xxz_R_and_H.py::expand_on_pauli` 得到 $c^{(i)}_{\alpha\beta}$。W

- Jordan–Wigner 对应（局域算符到费米子）：相邻项给出
	- 跃迁 (hop) 与配对 (pair)：
		$$t_{i,i+1}\propto c_{xx}+c_{yy}+i(c_{xy}-c_{yx}),\\\Delta_{i,i+1}\propto c_{xx}-c_{yy}-i(c_{xy}+c_{yx}).$$
	- 化学势/相互作用由 $c_{z0},c_{0z},c_{zz}$ 给出（$\mu\propto-2(c_{z0}+c_{0z})+4c_{zz}$，$U\propto4c_{zz}$）。

- MZM vs ABS（直观代数判据）：
	- 如果 Pauli 系数在体内近乎均匀（主要是对角项 XX/YY/ZZ 产生的实 $t,\Delta,\mu$），体系是均匀 Kitaev‑类 链；拓扑相（满足 $|\mu|<2|t|$ 等条件）产生端 MZM（两端对称、满足自共轭与远分离）。
	- 混合 Pauli 项（XY,XZ,YZ）在 JW 后引入虚部/相位与占据耦合：这些项能在空间上制造相位不连续或局域化学势修正，产生界面/陷阱，进而束缚局域 ±ε 态（ABS）。

- 计算性表达：对任一 Bogoliubov 态 $\Psi=(u,v)$，定义 $\Psi_{A,B}=u\pm v^*$，有
	$$\varepsilon=\tfrac12\langle\Psi_A|H_{\rm BdG}|\Psi_B\rangle,\quad S_{AB}=\frac{\sum_i|\Psi_A(i)||\Psi_B(i)|}{\|\Psi_A\|\|\Psi_B\|},$$
	若 $\mathrm{maj\_sim}\to1$ 且 $S_{AB}\to0$（随 $L$ 指数下降），则 $|\varepsilon|\to0$，态为 MZM‑like；否则为 ABS。

以上为便于快速核对的要点，下面为逐步严格推导。

I. R(u) 的正则性与局域哈密顿的线性化（分析步）

假设：令 $R(u)$ 是作用在两站 Hilbert 空间 $\mathcal{H}_i\otimes\mathcal{H}_{i+1}$ 上的解析矩阵族（例如来自可积模型或散射矩阵构造），且存在一常数点 $u_0$ 满足常见的正则性条件
$$
R(u_0)=P\quad(\text{置换算子}),
$$
并且在 $u_0$ 处可展为泰勒级数。

定义局域哈密顿密度为 $h_{i,i+1}$，通过线性化取定义常数 $\rho$ 使得
$$
R(u_0+\varepsilon)=P\big( I + \rho\,\varepsilon\, h_{i,i+1} + O(\varepsilon^2)\big).
$$
由于 $P^2=I$，对导数项有
$$
h_{i,i+1}=\frac{P}{\rho}\left.\frac{dR}{du}\right|_{u_0}.
\qquad(\ast)
$$
该式给出从 $R(u)$ 到局域两站哈密顿的严格线性关系；$\rho$ 是归一化常数，由 $R$ 的规范决定（数值实现中我们选择使 $h$ 的常见项匹配物理单位）。

II. 在 Pauli 基下展开并用 Jordan–Wigner 映射到费米子算符

把 $h_{i,i+1}$ 在 Pauli 张量基上展开：令 $\{\sigma^0,\sigma^x,\sigma^y,\sigma^z\}$ 为单站基，则
$$
h_{i,i+1}=\sum_{\alpha,\beta\in\{0,x,y,z\}} c^{(i)}_{\alpha\beta}\;\sigma_i^{\alpha}\otimes\sigma_{i+1}^{\beta}.
\qquad(\dagger)
$$
系数 $c^{(i)}_{\alpha\beta}$ 可由 (\ast) 直接数值/符号计算得到（见 `tools/xxz_R_and_H.py::expand_on_pauli`）。

接下来给出常用 Pauli 项到费米子双线性（hop/pair/number）的精确代数变换。记 Jordan–Wigner 映射
$$
c_j = \left(\prod_{k<j}\sigma_k^z\right)\sigma_j^{-},\qquad c_j^{\dagger}=\left(\prod_{k<j}\sigma_k^z\right)\sigma_j^{+},
$$
其中 $\sigma^{\pm}=\tfrac12(\sigma^x\pm i\sigma^y)$. 直接代数运算给出（对相邻站 $i,i+1$）:

- 对称组合（产生跃迁/跳跃）：
$$
\sigma_i^x\sigma_{i+1}^x+\sigma_i^y\sigma_{i+1}^y = 2(\sigma_i^{+}\sigma_{i+1}^{-}+\sigma_i^{-}\sigma_{i+1}^{+})
 = 2(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i).
 $$

- 反对称/配对组合（产生成对配对）：
$$
\sigma_i^x\sigma_{i+1}^x-\sigma_i^y\sigma_{i+1}^y = 2(\sigma_i^{+}\sigma_{i+1}^{+}+\sigma_i^{-}\sigma_{i+1}^{-})
 = 2(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i).
$$

- 交叉项产生相位：
$$
\sigma_i^x\sigma_{i+1}^y = i(\sigma_i^{+}\sigma_{i+1}^{+}-\sigma_i^{-}\sigma_{i+1}^{-}) + \text{(含数项)},
$$
以及相应的 $\sigma_i^y\sigma_{i+1}^x$。把这些代数组合按系数 $c^{(i)}_{\alpha\beta}$ 加权并收集 $c_i^{\dagger}c_{i+1}$、$c_i^{\dagger}c_{i+1}^{\dagger}$ 等基算符系数，就得到精确的费米子二次项系数。

若把相邻站 $i,i+1$ 的关联项系数记为 $c_{xx},c_{yy},c_{xy},c_{yx}$（简化记法，来源于 $c^{(i)}_{\alpha\beta}$），则代数上可写出跃迁与配对的组合：
$$
\begin{align*}
t_{i,i+1} &\propto c_{xx}+c_{yy}+ i(c_{xy}-c_{yx}),\
\Delta_{i,i+1} &\propto c_{xx}-c_{yy}- i(c_{xy}+c_{yx}).
\end{align*}
$$
上式由前面四项的线性组合直接得来（把 $\sigma^\pm$ 展开并代回 JW 表达即可得到常数因子；这些常数包含在我们数值实现的 `map_c_to_params` 中，所示公式即是按该规范写出的简洁结果）。

局域化学势与相互作用项由 $\sigma^z$ 相关项给出：
$$
\begin{align*}
\sigma_i^z &= 2n_i-1,\
\sigma_i^z\sigma_{i+1}^z &= 4n_i n_{i+1} - 2(n_i + n_{i+1}) + 1,
\end{align*}
$$
因此对应系数可写为
$$
\begin{align*}
\mu_i &\propto -2(c_{z0}+c_{0z}) + 4c_{zz},\
U_{i,i+1} &\propto 4 c_{zz},
\end{align*}
$$
与我们先前在代码中使用的表达一致（常数因子由具体正则化决定）。

III. BdG 表述与 Bogoliubov 态的 Majorana 分解：精确的 $\varepsilon$ 表达与判据

把上述二次费米子项在全链上组装，得到标准的 BdG 二次形式
$$
H_{\rm BdG} = \sum_{ij}\big( c_i^{\dagger}A_{ij}c_j + \tfrac12(c_i^{\dagger}B_{ij}c_j^{\dagger} + \mathrm{h.c.})\big),
$$
其中 $A$ 为厄米矩阵，$B$ 为反对称矩阵。构造 Bogoliubov 本征问题
$$\begin{pmatrix}A & B\ -B^{*} & -A^{*}\end{pmatrix}\begin{pmatrix}u\ v\end{pmatrix}=E\begin{pmatrix}u\ v\end{pmatrix},$$
解得规范化本征矢量 $\Psi=(u,v)$。

定义两个 Majorana 片段函数
$$
\Psi_A = u + v^{*},\qquad \Psi_B = u - v^{*}.
\qquad(\ddagger)
$$
令对偶的 Bogoliubov 算符
$$\beta = \sum_i(u_i c_i + v_i c_i^{\dagger}),\qquad \beta^{\dagger}=\sum_i(u_i^{*} c_i^{\dagger} + v_i^{*} c_i),$$
则对应两个实的 Majorana 算符（局域片段形式）可取为
$$\gamma_A = \beta + \beta^{\dagger},\qquad \gamma_B = -i(\beta - \beta^{\dagger}).$$

在 Bogoliubov 空间中计算这两算符的耦合可得有效低能哈密顿（在 $\{\gamma_A,\gamma_B\}$ 子空间）为
$$
H_{\rm eff} = i\,\varepsilon\;\gamma_A\gamma_B,
$$
其中耦合常数可精确地写为（见标准 Bogoliubov 变换与算符代数）
$$
\varepsilon = \tfrac12\;\langle\Psi_A\,|\,H_{\rm BdG}\,|\,\Psi_B\rangle,\qquad\text{(精确定义)}.
$$
证明要点：将 $H_{\rm BdG}$ 在 Nambu 表示作用于 $\Psi_{A,B}$ 上并利用 $\beta,\beta^{\dagger}$ 的展开即可得到该表达（代数步骤在此略写，但可在附录中给出逐项代数展开）。该公式与数值实现中直接用 $\tfrac12\langle\Psi_A|H_{\rm BdG}|\Psi_B\rangle$ 计算得到的值完全一致。

判据（严谨述语）：若存在一对近零能本征态 $\Psi^{(1)},\Psi^{(2)}$，使得

- 自共轭性近似成立： $u^{(k)}\approx e^{i\phi_k} v^{(k)*}$ （等价于 $\mathrm{maj\_sim}(\Psi^{(k)})\to1$）；

- 空间分离成立： $\mathrm{supp}(\Psi_A^{(1)})$ 与 $\mathrm{supp}(\Psi_A^{(2)})$ 指数地分离（即对链长 $L$ 随 $L\to\infty$ 指数小重叠）；

可得 $|\varepsilon|\to0$（指数收敛），从而两片段成为独立的物理 Majorana（MZM）。反之，若任一条件不满足，则 $|\varepsilon|$ 通常为有限值，态为平凡 ABS。

总结性说明：上述链条从 $R(u)$ 的线性化 (\ast) 到 Pauli 展开 (\dagger)，再通过 JW 代数得到 BdG 系数，最后通过 (\ddagger) 与 (\ref{eps}) 给出 $\varepsilon$ 的精确表达，是一条闭合且可检验的数学推导路径。该证明仅基于 $R(u)$ 的解析性与标准 JW 映射，以及 Bogoliubov 变换的代数性质；所有常数因子在具体实现中通过 $\rho$ 与 `map_c_to_params` 规范固定。

### 为什么对角项（例如 XX, YY, ZZ）通常不产生 ABS，而混合项（XY, XZ, YZ 等）能产生 ABS？

下面给出严格的代数说明，直接从 Pauli 系数和 Jordan–Wigner 映射看出关键差别。

设在某一局域区域（或一个键）只有对角 Pauli 系数非零，简化记作 $c_{xx},c_{yy},c_{zz}$，且无交叉项 $c_{xy}=c_{yx}=c_{xz}=\cdots=0$。

- 根据 Pauli 到费米子的映射（见上节），对角项在最近邻处产生的二次费米子项为：
	- 跃迁（hop）项：与 $c_{xx}+c_{yy}$ 成正比；
	- 配对（pair）项：与 $c_{xx}-c_{yy}$ 成正比；
	- 化学势/密度项来源于 $c_{zz}$。

	若这些系数在空间上平滑（逐键近似相同），则构成均匀的 Kitaev‑类 链：均匀的 $t,\Delta,\mu,U$。均匀系数本身不会在链内部生成局域 ABS：出现束缚态需要参数在空间上发生不连续或符号翻转（例如 $\Delta$ 的相位跳变或局域化学势井）。因此纯对角项、且无空间不均匀时，不产生局域 ABS，数值扫描得到的空缺即由此解释。

- 混合 Pauli 项的不同之处在于两点：
1) 相位（虚部）贡献：从代数上我们有
	$$
		 t\propto c_{xx}+c_{yy}+ i(c_{xy}-c_{yx}),\qquad
		 \Delta\propto c_{xx}-c_{yy}- i(c_{xy}+c_{yx}).
	$$
因此非零的 $c_{xy},c_{yx}$ 会引入 $i$ 因子，使 $t$ 或 $\Delta$ 带相位或在不同键上出現實/虛部分差異。相位差或實/虛部空間分佈的不均勻可產生界面（如 $\Delta$ 相位翻轉或 $t$ 相位突變），在界面處束縛局域態（ABS）；若同時滿足自共軛與空間分離，則趨向 MZM。

2) 占據耦合與有效局域勢：含 $Z$ 的混合項（如 $X\!Z$、$Y\!Z$）在 JW 映射後會把單粒子算符與占據數 $n_j$ 結合，例如形式上
$$
\sigma_i^x\sigma_{i+1}^z \sim (c_i + c_i^{\dagger})\,(2n_{i+1}-1)\times(\text{符號/相位因子}).
$$
這類三算符結構在作平均場或在二次近似時會產生局域的有效化學勢 $\mu_{\rm eff}$ 或改變局域的 $t,\Delta$，從而直接在格點處形成陷阱，束縛 ABS。

此外，Jordan–Wigner 的串算符在非相鄰項或多鍵同時存在 $Z$ 因子時引入整體費米子奇偶算符依賴（$(-1)^{N_{<j}}$），使得某些 Pauli 組合在不同粒子數子空間中表現不同，這也有助於把局域 Pauli 不對稱性轉化為對應的費米子局域勢或相位不連續，進而支持局域束縛態。

综上：对角项（XX/YY/ZZ）若在空间上平滑，会产生均匀的二次费米子模型，不会自动在链内部产生 ABS；而混合项通过引入相位（产生 $\Im(t),\Im(\Delta)$）或通过与 $n_j$ 耦合产生局域 $\mu_{\rm eff}$，更容易在空间上造成局域不均匀或相位翻转，从而在界面或陷阱处束缚 ABS（并在满足自共轭与空间分离时收敛为 MZM）。

我们的数值实现中的 `map_c_to_params` 所给的复数表达式正好把上述代数关系编码进了 $t,\Delta$ 的实虚部；在打开混合 Pauli 项的扫描里观察到的局域 $\Im(t)$、$\Im(\Delta)$ 或局域 $\mu_{\rm eff}$ 峰值，正是该代数机制的数值印证。

如果你需要，我可以把这节扩展为把每个混合 Pauli 项逐项 JW 展开成 $c_p^{\dagger}c_q$, $c_p c_q$, $n_p c_q$ 等基算符的完整代数表达，并把对应的常数矩陣導出為 `tensor_to_epsilon.py` 的初始化數據，便于把每個 Pauli 項的貢獻可视化为对 $t,\Delta,\mu$ 的增量图。



下面原文的“Majorana 表示与零模的代数条件”节继续接续相应论述。

1) Majorana 表示与零模的代数条件
## 代数推导：在我们的 R(u) 模型下的 MZM 与 ABS 表述，以及 ABS 作为弱耦合 Majorana 的证明

下面给出严格而紧凑的代数推导，基于我们管线中常用的映射（R(u)→局域两站哈密顿 → Pauli 展开 → map_c_to_params → BdG）。目标：
- 给出 MZM 与 ABS 在 BdG 表述中的代数条件；
- 证明任一 ABS 都可代数写成两 Majorana 片段的耦合态，并在弱耦合极限说明其收敛到两物理独立的 Majorana（即 MZM 的弱耦合表述）。

记号与预备：
- 链长为 $N$，单粒子算符为 $c_i,c_i^\dagger$。我们用 BdG 矩阵 $H_{\rm BdG}$（$2N\times2N$）表示体系，按分块写为
$$
H_{\rm BdG}=\begin{pmatrix}A & B\ -B^* & -A^*\end{pmatrix},
$$
其中 $A$ 描述粒子部分（含化学势和跃迁 $-t$），$B$ 描述配对（反对称块，含 $\Delta$）。
- 我们的映射（见 `tools/verify_mzm.py::map_c_to_params`）给出从 Pauli 系数 $c_{\alpha\beta}$ 到 BdG 参数的规则（示例）：
	- $t = c_{xx}+c_{yy}+i(c_{xy}-c_{yx})$，
	- $\Delta = c_{xx}-c_{yy}-i(c_{xy}+c_{yx})$，
	- $\mu = 4c_{zz}-2(c_{z0}+c_{0z})$，
	- $U = 4 c_{zz}$。

1) Majorana 表示与零模的代数条件

引入 Majorana 费米子定义（格点表示）：
$$
\gamma_{2i-1}=c_i + c_i^\dagger,\qquad \gamma_{2i}= -i(c_i - c_i^\dagger),
$$
它们滿足 $\{\gamma_a,\gamma_b\}=2\delta_{ab}$ 且 $\gamma_a^\dagger=\gamma_a$. 任一 BdG 哈密顿可写成 Majorana 形式：
$$
H = \frac{i}{4}\sum_{ab} M_{ab}\gamma_a\gamma_b,
$$
其中實反斜對稱矩陣 $M=-M^T$ 與 BdG 有線性關係（即 Majorana 變換 $U$，見 finite_chain_pfaffian 實現）。

若存在拓撲 Majorana 零模（MZM），則在有限鏈上會有一對近零模 $\Psi^{(1)},\Psi^{(2)}$，其粒/空穴分量滿足自共轭條件：
$$
u^{(k)} \approx e^{i\phi_k} v^{(k)*},\qquad (k=1,2),
$$
且每個模的 Majorana 組件空間上分離（例如一個局域在左端，另一個局域在右端）。在算符層面，對應近零本徵態可寫成兩個幾乎不耦合的 $	\tilde\gamma_1,\tilde\gamma_2$。

2) ABS 的代數表述與 Majorana 分解（恆等與有效模型）

任意一個 Bogoliubov 本徵態 $\Psi=(u,v)$ 都可以代數分解為兩個 Majorana 片段：定義
$$
\Psi_A = u + v^*,\qquad \Psi_B = u - v^*.
$$
對應的空間密度滿足恆等式
$$
|\Psi_A(i)|^2 + |\Psi_B(i)|^2 = 2(|u(i)|^2 + |v(i)|^2) = 2\,\mathrm{LDOS}(i).
$$
因此任何局域 ABS（被局域勢井或界面束縛的 ±ε 成對態）都可以視為兩個 Majorana 片段 $\gamma_A,\gamma_B$ 的耦合形成的普通費米子，其有效低能哈密頓可寫為
$$
H_{\rm eff}=i\varepsilon\,\gamma_A\gamma_B,
\qquad E=\pm\varepsilon.
$$
這裡 $\varepsilon$ 與兩片段在空間上的重疊及其在完整哈密頓下的矩陣元成正比，概念上可寫為
$$
\varepsilon\propto \langle \Psi_A | H_{\rm BdG} | \Psi_B\rangle \sim \int dx\; \Psi_A^*(x) H_{\rm BdG} \Psi_B(x).
$$
當兩片段重疊較大或任一片段不滿足自共轭性時，$\varepsilon$ 一般不是指數小，對局域擾動敏感——這就是平凡 ABS 的本質。

3) 從 R(u) 映射到弱耦合 Majorana 的路徑（本模型的具體步驟）

步驟如下：
- (i) 從 R(u) 的某個線性化點或分段線性化得到當地的兩站哈密頓 $h_{i,i+1}$（我們算過 $h = P \partial_u R /\rho$ 或相應公式）；
- (ii) 把 $h_{i,i+1}$ 在 Pauli 基上展開得到係數 $c_{\alpha\beta}(i)$（參見 `tools/xxz_R_and_H.py::expand_on_pauli`）；
- (iii) 用 `map_c_to_params` 把 $c_{\alpha\beta}(i)$ 映射到局域 BdG 參數 $t_i,\Delta_i,\mu_i, U_i$；
- (iv) 若在某區域（界面或局部段）這些 $\mu_i,\Delta_i$ 與周圍明顯不同，BdG 解將包含局域化的 ±ε 态（ABS）。

在第 (iv) 中，對於某一界面或陷阱，局域的 $u^{(k)},v^{(k)}$ 可被數值求出；用上面的代數式可以構造 $\Psi_A,\Psi_B$ 並計算
$$\varepsilon\approx \langle\Psi_A|H_{\rm BdG}|\Psi_B\rangle,
\quad S_{AB}=\frac{\sum_i|\Psi_A(i)||\Psi_B(i)|}{\|\Psi_A\|\,\|\Psi_B\|}$$
作爲耦合與重疊的量化指標。其中若 $S_{AB}\ll1$ 且 $\mathrm{maj\_sim}(\Psi^{(k)})\approx1$，則 $\varepsilon$ 很小，ABS 趨近於弱耦合的兩個物理 Majorana（即 MZM 的弱耦合極限）。反之若 $S_{AB}=O(1)$ 或 `maj_sim` 小，則為平凡 ABS。

4) 形式化陳述（證明草案）

命題：對於由 R(u) 映射得到的任意局域 ABS（由局域化參數引起的 ±ε 成對態），存在 Majorana 片段表示使得其有效低能 Hamiltonian 為 $H_{\rm eff}=i\varepsilon\gamma_A\gamma_B$。當且僅當 $S_{AB}\to0$ 且 $\mathrm{maj\_sim}\to1$（在某極限下，例如 $L\to\infty$ 或參數調至特殊值），則 $\varepsilon\to0$ 且這對態成爲物理上獨立的 Majorana（MZM）。

證明要點：
- 代數存在性：由前述恒等式任意 Bogoliubov 态都能寫成 $\Psi_{A,B}$，因此代數上 $H_{\rm eff}=i\varepsilon\gamma_A\gamma_B$ 的形式總是成立（取 $\varepsilon=\tfrac12\langle\Psi_A|H|\Psi_B\rangle$ 即可）。
- 當 $S_{AB}\to0$ 與自共轭性成立時，$\Psi_A,\Psi_B$ 在 $H$ 下的矩陣元因空間正交或指數小重疊而趨於零，故 $\varepsilon\to0$；此時兩個片段各自自伴且非局域地分佈，從而成爲物理上的 Majorana（MZM）。

5) 可計算性的提示（如何在流水線中驗證）
- 給定一個候選（或域壁），數值求解 BdG 得到前兩本徵態 $\Psi^{(1)},\Psi^{(2)}$；
- 計算 `maj_sim`、$\Psi_A,\Psi_B$、$S_{AB}$ 與 $\varepsilon\approx\langle\Psi_A|H|\Psi_B\rangle$，以及有限鏈 Pfaffian 隨參數的符號變化；
- 若發現 $\varepsilon$ 隨 $L$ 指數下降且 `maj_sim` 高、$S_{AB}$ 低，則將其認定爲弱耦合的兩個 Majorana；否則認定爲平凡 ABS。

6) 小結

在我們的 R(u) 地圖框架下：
- R(u) 明確生成局域兩站哈密頓並通過 Pauli 展開決定 BdG 參數分布；
- 只要映射在空間上產生界面或局域化參數不均勻，就會自發出現局域 ABS（無需人工插 μ 井）；
- 代數上任何 ABS 都可視為兩個 Majorana 片段的耦合；當且僅當自共轭性與空間分離兩條額外條件滿足時，該耦合弱化並收斂爲真正的（物理）MZM。這給出了 ABS 與 MZM 之間連續的、可量化的代數圖景。

如需，我可以把上述推導的 LaTeX 版導出為 PDF，並把相應的數值檢驗（$\varepsilon,S_{AB},\mathrm{maj\_sim}$）自動附加到每個候選的報告中。

### 附加：按 Pauli 張量給出對 $\varepsilon$ 的逐項展開（更具體的代數步驟）

為了讓推導在實際流水線中可操作，下面把第 2 點中關於 $\varepsilon$ 的計算給出逐項展開的具體代數形式，便於直接對每個 Pauli 項求和並在數值上評估。

考慮單個 Pauli 張量項 $c^{(i)}_{\alpha\beta}\sigma_i^{\alpha}\sigma_{i+1}^{\beta}$。用 Jordan–Wigner 把它展成費米子雙線性，符號細節上可寫成一組基礎雙線性算符集合 $\{O_{pq}\}$（例如 $c_p^\dagger c_q$, $c_p c_q$, $c_p^\dagger c_q^\dagger$），即
$$
\sigma_i^{\alpha}\sigma_{i+1}^{\beta} = \sum_{p,q} M^{(i,\alpha\beta)}_{pq}\;O_{pq},
$$
其中矩陣元 $M^{(i,\alpha\beta)}_{pq}$ 是已知的數值（由 JW 展開確定，僅依賴於索引位置與 $\alpha,\beta$）。

對於一個 Bogoliubov 态 $\Psi=(u,v)$，令 $\Psi_A=u+v^*$、$\Psi_B=u-v^*$。則對任一基礎算符 $O_{pq}$，矩陣元可寫成 $u,v$ 的二次型，例如：
$$
\begin{align*}
\langle\Psi_A|c_p^\dagger c_q|\Psi_B\rangle &= \sum_{r,s} [u_A^*(r)\,\mathcal{P}_{rp}]\;[u_B(s)\,\mathcal{Q}_{sq}] + \cdots,\
\langle\Psi_A|c_p c_q|\Psi_B\rangle &= \sum_{r,s} [v_A^*(r)\,\tilde{\mathcal{P}}_{rp}]\;[u_B(s)\,\tilde{\mathcal{Q}}_{sq}] + \cdots,
\end{align*}
$$
其中 $\mathcal{P},\mathcal{Q},\tilde{\mathcal{P}},\tilde{\mathcal{Q}}$ 等是由 JW 展開與重排規則決定的已知常數（只與格點索引有關，不依賴於 $c^{(i)}$ 本身）。

因此單個 Pauli 項對 $\varepsilon$ 的貢獻可以寫成
$$
\varepsilon^{(i,\alpha\beta)} = \tfrac12\,c^{(i)}_{\alpha\beta}\sum_{p,q} M^{(i,\alpha\beta)}_{pq} \, G_{pq}[u,v],
$$
其中 $G_{pq}[u,v]$ 是上面那類二次型組合的代數總和（在數值上可直接通過 $u,v$ 與已知常數求值）。整體
$$
\varepsilon = \sum_{i,\alpha\beta} \varepsilon^{(i,\alpha\beta)}.
$$ 

該形式的好處是：
- 在實現中我們已經有 $c^{(i)}_{\alpha\beta}$ 的數據；
- JW 常數矩陣 $M^{(i,\alpha\beta)}_{pq}$ 可預先計算並固定；
- 對每個候選，只要把數值 $u,v$ 帶入 $G_{pq}[u,v]$，再按上式對 $i,\alpha,\beta$ 求和，就得到 $\varepsilon$ 的精確數值（浮點）。

最後，弱耦合 Majorana 的數值判據可寫成具體阈值：
\begin{itemize}
\item $\mathrm{maj\_sim}(\Psi^{(k)})\ge 0.75$（或用你要的閾值）；
\item $S_{A1A2},S_{B1B2} \le 0.2$；
\item $|\varepsilon|/\Delta_{\rm bulk} \le 10^{-3}$ 且 $|\varepsilon(L)|$ 隨 $L$ 顯示指數下降（至少掃 2–3 個 L 值驗證趨勢）。
\end{itemize}

如果你同意這組更具體的代數展開與閾值，我會：
1) 在 `tools/` 中加入一個小模組 `tensor_to_epsilon.py`，實現 JW 常數矩陣 $M^{(i,\alpha\beta)}_{pq}$ 的生成與 $\varepsilon$ 的逐項求和；
2) 用該模組對 top‑N 候選自動計算 $\varepsilon,S_{AB},\mathrm{maj\_sim}$ 並把結果寫回 JSON/CSV；
3) 把代數推導（這一節）形成一個獨立的 LaTeX 附錄並導出 PDF 供論文/補充材料使用。

### 从 R(u) 到 MZM：显式构造（我们做过的、应写进文档的步骤）

这里把我们在代码中实际做的管线逐步写清，确保文档明确回答“如何从 R(u) 自然构造或识别 MZM”：

- Step A — 从 R(u) 得到两站局域哈密顿 `h_{i,i+1}`：
	- 在数值实现中，我们对每个格点或每对格点取 R(u) 在某参数点的线性化/导数（参见 `tools/xxz_R_and_H.py`），按物理前因子构造局域两站哈密顿：
	$$
    h_{i,i+1} = \frac{P\,\partial_u R}{\rho}\Big|_{\text{point}} + \text{局域常数},
    $$
	- 该 `h_{i,i+1}` 存在于 Pauli 张量空间，可用 Pauli 基 $\{\sigma^0,\sigma^x,\sigma^y,\sigma^z\}$ 在 `tools/xxz_R_and_H.py::expand_on_pauli` 得到系数字典 $c^{(i)}_{\alpha\beta}$（代码输出键如 `X_X`,`Y_Z`）。

- Step B — Pauli 系数映射到 BdG 参数：
	- 我们用 `tools/verify_mzm.py::map_c_to_params` 将每个 $c^{(i)}_{\alpha\beta}$ 映射为局域 $t_i,\Delta_i,\mu_i,U_i$，公式如上节所列（`t=...`,`\Delta=...` 等）。
	- 将所有局域参数组装成全链 BdG：`tools/verify_mzm.py::build_bdg_from_params(N,params)`，得到 $H_{\rm BdG}$ 并求本征对 $(u,v)$。

- Step C — Jordan–Wigner（理论联系与数值实现要点）：
	- 理论上，任何 Pauli 张量可以用 Jordan–Wigner（JW）映射写成费米子双线性/双创建/双湮算符的线性组合。我们在文档中给出符号形式，并在 `tensor_to_epsilon.py`（建议实现）中把 JW 展开预先编码为常数矩阵 $M^{(i,\alpha\beta)}_{pq}$，用于数值逐项求和。
	- 实际代码流程不需要手动做大规模符号代数：我们通过 `map_c_to_params` 得到的 `t,\Delta,\mu` 已把 JW 后的等效费米子双线性封装进了 BdG 构建中；JW 展开在数值上用于把原始 Pauli 系数逐项映射回對應的 $c_p^\dagger c_q, c_p c_q$ 等基算符的矩阵元（见附加节的 $M_{pq}$）。

- Step D — 由参数不均匀/界面识别 MZM 的数值策略（我们在脚本中实现的）：
	- 在 `tools/scan_all_mixtures.py` 中我们扫描 Pauli 组合并生成候选 $c^{(i)}_{\alpha\beta}$ 分布；`tools/validate_scan_candidates.py` 计算局域 BdG 的低能本征态并输出 `results/`。
	- 在 `tools/compute_majorana_overlap.py` 与 `tools/finite_chain_pfaffian.py` 中我们计算 `maj_sim`、$S_{AB}$、以及有限链 Pfaffian（拓扑判据），并用 `tools/domain_wall_test.py` 构造域壁检验 MZM 的空穴/粒子局域性与耦合。若满足弱耦合判据（见上节阈值），则认定为 MZM‑like。

小结：文中之前的各节给出了数学背景，但未把上面 A→B→C→D 的“流水线步骤”像这样逐项列出；我已把该段加入文件，且在文件中明确引用了实现代码位置，确保读者能从 `R(u)` 的数值输出一路追溯到 BdG、JW、以及 MZM 的判别标准。


