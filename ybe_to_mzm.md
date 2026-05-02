(# 从 YBE R(u) 到 MZM — 椭圆型（XYZ / eight‑vertex）R 矩阵定义与 Pauli 展开)

## 引言：为什么从 R(u) 出发？YBE、含时哈密顿与 MZM/ABS 的清晰化

选择以谱参数依赖的两体算子 $R(u)$ 为出发点有三层理由：代数性（Yang–Baxter 结构）、局域生成性（谱参数导数产生局域生成元），以及把代数约束映射为物理参数（hopping/配对/化学势），从而把可积性与拓扑性质联系起来。

- 代数性：$R(u)$ 满足 Yang–Baxter 方程
$$
R_{12}(u-v)R_{13}(u)R_{23}(v)=R_{23}(v)R_{13}(u)R_{12}(u-v),
$$
 该恒等保证在多体推广时的一致性（例如 transfer matrix 的可交换性），并限制了 $R$ 在 Pauli 基下的系数，从而影响生成的物理哈密顿子类型（例如 XXZ 保持 $S^z$，不产生配对项）。

- 局域生成性（修正）：把谱参数作为时间或参数轨迹，令 $U(t)=R(u(t))$，则
$$
\frac{d}{dt}U(t)=\frac{du}{dt}\,\partial_u R\big(u(t)\big).
$$
比较薛定谔方程 $\dfrac{d}{dt}U(t)=-iH_P(t)U(t)$，得到瞬时生成元
$$
H_P(t)=i\frac{du}{dt}\,\partial_u R(u)\,R(u)^{-1}.
$$
因此定义左对数导数（左生成元）
$$
H(u):=i\partial_u R(u)\,R(u)^{-1},
$$
使得（按 $u$ 的路径序列）
$$
R(u)=\mathcal{T}_u\exp\Big(-i\int^{u} H(u')\,du'\Big),
$$
其中 $\mathcal{T}_u$ 为谱参数的路径有序算符。如果以时间 $t$ 参数化 $u(t)$，则 $H_P(t)=\dfrac{du}{dt}H(u(t))$。

注意：右对数导数 $H_{\rm right}(u)=-iR^{-1}\partial_u R$ 也可定义，且与左对数满足 $H_{\rm right}=-R^{-1}H R$；使用时应保持一致的约定。

若 $R(u)$ 带有纯标量因子（$R(u)=\rho(u)\tilde R(u)$），则
$$
\partial_u R R^{-1}=\frac{\partial_u\rho}{\rho}\,I + \partial_u\tilde R\,\tilde R^{-1},
$$
第一项仅是全局相位（可忽略），第二项提供非平凡的 Pauli 张量基生成元。

若 $\det R(u)=0$（不可逆点），左对数会发散，实际计算中需用位移/分段路径/插值等正则化策略绕过奇点。

## 正确的求解流程（严格管线：R(u) → H(u) → Pauli → BdG → MZM/ABS）

下面给出从代数到物理并可数值实现的严格流程，包含常用数值稳定化方法与验证步骤：

1) 理论关系（左对数导数、路径序列）

把 $R(u)$ 当作谱参数路径上的演化算子：
$$
R(u)=\mathcal T_u\exp\Big(-i\int^{u} H(u')\,du'\Big),
\\
H(u)=i\partial_u R(u)\,R(u)^{-1}.
$$

2) 数值稳定化与正则化

- 有限差分 + 矩阵对数：对小步长 $\delta u$ 计算
$$
H(u)\approx \frac{i}{\delta u}\log\big(R(u+\delta u)R(u)^{-1}\big)
$$
建议使用 `scipy.linalg.logm`，并对结果做 $\tfrac12(H+H^\dagger)$ 对称化以去掉数值噪声。

- 极分解取单元部分：写 $R=UP$（$U$ unitary，$P$ 正定），用 $U$ 构造厄米生成元
$$
H_U(u)=i\partial_u U\,U^{-1}.
$$

- 去除纯标量因子：若 $R=\rho(u)\tilde R$，则可以忽略 $\partial_u\rho/\rho$ 项。

- 奇点处理：若 $\det R(u)=0$，用位移/插值/分段路径避免奇点。

3) 把生成元 $H(u)$ 展开到 Pauli 张量基（得到 $h_{\alpha\beta}$）
$$
H(u)=\sum_{\alpha,\beta\in\{0,x,y,z\}} h_{\alpha\beta}(u)\,\sigma^\alpha\otimes\sigma^\beta,
\\
\qquad h_{\alpha\beta}=\tfrac14\operatorname{Tr}\big[(\sigma^\alpha\otimes\sigma^\beta)H\big].
$$

重要：使用的是 $H$ 的系数 $h_{\alpha\beta}$，而不是直接用 $R$ 的系数 $c_{\alpha\beta}$（后者仅可作为启发式筛选）。

4) 映射到 Kitaev 参数（最近邻二次项）

用仓库中 `map_c_to_params` 的约定（或以下公式）把 $h_{xx},h_{yy},h_{xy},h_{yx},h_{zz},h_{z0},h_{0z}$ 转为 $(t,\Delta,\mu)$：
$$
\boxed{\begin{aligned}
t &= h_{xx}+h_{yy}+i\big(h_{xy}-h_{yx}\big),\\[4pt]
\Delta &= h_{xx}-h_{yy}-i\big(h_{xy}+h_{yx}\big),\\[4pt]
\mu &= 4h_{zz}-2\big(h_{z0}+h_{0z}\big).
\end{aligned}}
$$
在严格管线中应使用 $h_{\alpha\beta}$，而不是 $c_{\alpha\beta}$。

5) 构造 BdG / 多体 Hamiltonian

- 若 $H$ 仅含二次项（quadratic），直接构造 BdG 并对角化；
- 若出现四阶密度–密度或带 JW 串的项，可在小 L 下做精确多体对角化，或做 mean‑field 脱耦得到有效二次 BdG（记录近似）。

6) Majorana 转换与拓扑/局域性检验

- 从 BdG 得到 Majorana 耦合矩阵 $A$（$H=(i/4)\\Gamma^T A\\Gamma$），求 $\\ker A$ 的零模；
- 计算 Pfaffian、有限链能谱、Majorana 重叠（maj_sim）、局域态密度（LDOS）；做 L‑标度与对局域扰动的稳健性测试以区分 MZM 与 ABS。

## 1. 计算基表示
（标准顺序 $|00\rangle,|01\rangle,|10\rangle,|11\rangle$）

$$
R(u)=\begin{pmatrix}
 -a(u) & 0 & 0 & d(u)\\[4pt]
 0 & b(u) & c(u) & 0\\[4pt]
 0 & c(u) & b(u) & 0\\[4pt]
 -d(u) & 0 & 0 & a(u)
\end{pmatrix}.
$$

六顶点（XXZ）为其三角极限（$d(u)\equiv0$）。

## 2. Baxter（椭圆）参数化（theta 表示）

设完整椭圆积分 $K$，共轭 $K'$，nome $q=\exp(-\pi K'/K)$。用 Jacobi theta 函数 $\theta_i(z|q)$ 写出 Baxter 权重。常用记法：

- $H(u)\propto\theta_1\big(\tfrac{\pi u}{2K}\big|q\big)$，
- $\Theta(u)\propto\theta_4\big(\tfrac{\pi u}{2K}\big|q\big)$。

一种常用权重（整体归一化为 $\rho(u)$）：
$$
\begin{aligned}
 -a(u)&=\rho(u)\,\frac{H(\eta-u)\,\Theta(0)}{H(\eta)\,\Theta(u)},\\[4pt]
 -b(u)&=\rho(u)\,\frac{\Theta(\eta-u)\,H(0)}{\Theta(\eta)\,H(u)},\\[4pt]
 -c(u)&=\rho(u)\,\frac{H(\eta)\,\Theta(u)}{H(\eta)\,\Theta(0)},\\[4pt]
 -d(u)&=\rho(u)\,\frac{\Theta(\eta)\,H(u)}{\Theta(\eta)\,H(0)}.
\end{aligned}
$$

在一般情形 $d(u)\neq0$（eight‑vertex）；三角极限 $q\to0$ 下 $d(u)\to0$ 退化为六顶点。

## 3. Pauli 张量基展开

在 Pauli 基 $\{\sigma^\alpha\otimes\sigma^\beta\}$（$\sigma^0=I$）上做线性展开
$$
R(u)=\sum_{\alpha,\beta\in\{0,x,y,z\}} c_{\alpha\beta}(u)\,\sigma^\alpha\otimes\sigma^\beta,
$$
其中（实权重情形）
$$
\begin{aligned}
 c_{00}(u)&=\tfrac{a+b}{2},\\[4pt]
 c_{zz}(u)&=\tfrac{a-b}{2},\\[4pt]
 c_{xx}(u)&=\tfrac{c+d}{2},\\[4pt]
 c_{yy}(u)&=\tfrac{c-d}{2}.
\end{aligned}
$$

## 4. 映射到 Kitaev 参数（仓库中 `map_c_to_params` 的形式）

在仓库约定下（与 `map_c_to_params` 保持一致），映射示意为
$$
\begin{aligned}
 t &= c_{xx}+c_{yy},\\[4pt]
 \Delta &= c_{xx}-c_{yy}-i(c_{xy}+c_{yx}),\\[4pt]
 \mu &= 4c_{zz}-2(c_{z0}+c_{0z}).
\end{aligned}
$$
六顶点（$d\equiv0$）时 $\Delta\equiv0$；椭圆（eight‑vertex）通常可产生 $\Delta\neq0$。

## 6. 含单体项与更一般的 JW 映射说明

若 $R$ 含有 $\sigma^{x,y}\otimes I$ 或 $I\otimes\sigma^{x,y}$ 等单体项，JW 映射会引入 parity string（非局域），导致多体项。建议在实现中：

1. 在 `expand_on_pauli` 中保留所有 $c_{\alpha\beta}$；
2. 用 `tools/jw_expand.py` 做精确 JW 展开（small‑L 精确对角化用于判别 ABS/MZM）；
3. 对大系统用 mean‑field 或收集二次项的近似生成 BdG 做快速筛选；
4. 对候选用 Pfaffian、L‑标度、局域扰动测试区分 MZM/ABS。

## 7. 增加混合项（XY, YX, XZ, ZX, YZ, ZY）的代数映射与物理后果

在更一般的 Pauli 展開中，混合项可以書寫為
$$
R=\sum_{\alpha,\beta\in\{0,x,y,z\}} c_{\alpha\beta}\,\sigma^\alpha\otimes\sigma^\beta.
$$
下面按類別討論這些混合項在 Jordan–Wigner（JW）映射下的代數影響以及物理含義。

### 7.1 XY / YX 類（雙側均為 x 或 y）——直接貢獻到二次項
對於兩側都屬於 $\{x,y\}$ 的項，可以用前面給出的恒等式把它們寫成自由費米子的 hopping 與 pairing。總結性的映射（與 `map_c_to_params` 一致）為：
$$
\boxed{\begin{aligned}
t &= c_{xx}+c_{yy}+i\big(c_{xy}-c_{yx}\big),\\[4pt]
\Delta &= c_{xx}-c_{yy}-i\big(c_{xy}+c_{yx}\big).
\end{aligned}}
$$
因此非零 $c_{xy},c_{yx}$ 會引入複數部分（相位），可破壞時間反演並產生手性/相位配對，這對於實現拓撲 p‑wave（MZM）是直接有利的。

### 7.2 XZ / ZX, YZ / ZY 類（含 z 與 x/y）——出現奇偶字符串或密度相關項（多體）
這類項在 JW 映射下通常會帶上 parity‑string（JW 串），導致非二次的 many‑body 項。記 $P_j=\prod_{l<j}\sigma^z_l=(-1)^{\sum_{l<j}n_l}$，有常用恆等式（站點 $j,j+1$）：
$$
\begin{aligned}
\sigma^z_j\sigma^x_{j+1} &= P_j\,(c_{j+1}+c_{j+1}^\dagger),\\[4pt]
\sigma^z_j\sigma^y_{j+1} &= -i\,P_j\,(c_{j+1}-c_{j+1}^\dagger),\\[4pt]
\sigma^x_j\sigma^z_{j+1} &= P_j\,(c_j+c_j^\dagger)\,(2n_{j+1}-1),\\[4pt]
\sigma^y_j\sigma^z_{j+1} &= -i\,P_j\,(c_j-c_j^\dagger)\,(2n_{j+1}-1).
\end{aligned}
$$
觀察：右端含有 $P_j$ 或 $(2n-1)$ 因子，這些因子使得對應 Hamiltonian 是交互型或包含非局域字符串；只有在特殊情況下（例如在固定全局奇偶性子空間、或對 $P_j$ 做替代取值 $\pm1$、或對 $(2n-1)$ 做 mean‑field）才能把它們近似為二次項。

物理含義：
- 這類項能產生密度依賴的 hopping/pairing（density‑assisted processes）或 parity‑string‑assisted 單體場；在某些固定子空間或平均場處理下會導致有效的 $t$ 與 $\Delta$，但一般會引入多體效應並改變零模的穩健性。

### 7.3 zz 與密度–密度項
$$\sigma^z_j\sigma^z_{j+1}=4\,(n_j-\tfrac12)(n_{j+1}-\tfrac12)$$
這是典型的交互項（四階），若要在二次自由模下處理，需要做 mean‑field 脫耦或把它重新分配為等效化學勢（見前文）。交互項能改變相圖並在某些情況下促進或抑制拓撲相。
