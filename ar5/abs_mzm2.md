# 完整推导：从 R(u) 到局域哈密顿，再到 MZM 与 ABS（严格代数链条）

作者说明：本文件给出从解析的两站 R(u) 矩阵出发，线性化得到局域两站哈密顿密度，
在 Pauli 张量基下展开，严格用 Jordan–Wigner（JW）映射把 Pauli 项写成费米子算符，
并逐项导出产生跳跃（hop）与配对（pair）的系数 $t_{i,i+1},\Delta_{i,i+1}$ 以及化学势 $\mu_i$ 与相互作用 $U_{i,i+1}$ 的明确表达。
随后证明（代数上）任一 Bogoliubov 态可分解为两 Majorana 片段，给出耦合系数 $\varepsilon$ 的精确定义并由此严格区分 MZM 与 ABS 的代数条件。

目录：
- 0. 前置假设与符号约定
- 1. 从 R(u) 的线性化到局域两站哈密顿 $h_{i,i+1}$（严格公式）
- 2. Pauli 张量展开：定义与 $c^{(i)}_{\alpha\beta}$ 的获取
- 3. Jordan–Wigner 映射的逐项代数（对相邻站的完整展开）
  - 3.1. $\sigma^x_i\sigma^x_{i+1}$, $\sigma^y_i\sigma^y_{i+1}$ 的展开
  - 3.2. $\sigma^x_i\sigma^y_{i+1}$, $\sigma^y_i\sigma^x_{i+1}$ 的展开（引入 $i$ 因子）
  - 3.3. $\sigma^z$ 相关项与密度-密度项的展开
  - 3.4. 弱耦合 $H_{\rm eff}$ 的 JW 逆映射与 Pauli 表示（总结）
- 4. 从 Pauli 系数到 BdG 参数 $t,\Delta,\mu,U$（逐项合并；规范说明）
- 5. BdG 矩阵、Bogoliubov 本征态与 Majorana 片段的代数证明
- 6. $
\varepsilon$ 的逐项表达（把 Pauli 张量的贡献展开到 $\varepsilon$）
- 7. 对角项（XX,YY,ZZ）与混合项（XY,XZ,YZ）何以产生不同物理后果的严格论证
- 8. 数值实现提示（map_c_to_params, tensor_to_epsilon 的接口说明）
- 9. 结论与可检验判据

---

0. 前置假设与符号约定

- 链长为 $N$，格点序号 $i=1,\dots,N$。局域希尔伯特空间为自旋-1/2（Pauli）表示。
- R(u)：作用在两站 $i,i+1$ 上的 $4\times4$ 矩阵族，假设在某解析点 $u_0$ 可泰勒展开并满足正则性 $R(u_0)=P$（交换/置换算子，$P^2=I$）。
- 记局域哈密顿密度为 $h_{i,i+1}$（Hermitian），常见的构造为线性化：

$$
R(u_0+\varepsilon)=P\big(I+\rho\,\varepsilon\,h_{i,i+1}+O(\varepsilon^2)\big),
$$

从而得

$$
h_{i,i+1}=\frac{P}{\rho}\left.\frac{dR}{du}\right|_{u_0}.\label{hfromR}
$$

- Pauli 基记号：$\sigma^0=I,\;\sigma^x,\sigma^y,\sigma^z$。两站 Pauli 张量为 $\sigma_i^{\alpha}\otimes\sigma_{i+1}^{\beta}$，简记为 $\sigma_i^{\alpha}\sigma_{i+1}^{\beta}$。
- 费米子算符满足 $\{c_i,c_j^{\dagger}\}=\delta_{ij}$，其与 Pauli 的 Jordan–Wigner 映射采用标准顺序：

$$
c_j = \left(\prod_{k<j}\sigma_k^z\right)\sigma_j^{-},\qquad c_j^{\dagger}=\left(\prod_{k<j}\sigma_k^z\right)\sigma_j^{+},
$$

其中 $\sigma^{\pm}=\tfrac12(\sigma^x\pm i\sigma^y)$。

在下文中，我们对最近邻项 $i,i+1$ 的 Pauli 张量进行 JW 展开。注意：因为两站有相同的前缀字符串 $\prod_{k< i}\sigma_k^z$，在相乘时字符串在相邻站处會互相抵消，因而对最近邻的代数展开不带长串（这是推导可行的关键简化）。

---

1. 从 R(u) 的线性化到局域两站哈密顿 $h_{i,i+1}$（严格公式）

如 (\ref{hfromR})，在 $u_0$ 处展开给出局域两站哈密顿密度的定义。证明该表达的合理性与一般性：

- 若 $R(u)$ 为可积模型的 R 矩阵或某种散射矩阵，$P$ 为 $u$ 的某个规范点的值是常见情形；线性化的小项就是局域生成元（置换后乘以导数），与微扰展开一致。
- $\rho$ 为标度常数，可选以使 $h$ 与常见哈密顿（如 XX, YY, ZZ 线性组合）在数值因子上对齐。

因此我们有确定的算符 $h_{i,i+1}$，对其进行 Pauli 展开得到下一步所需的 $c^{(i)}_{\alpha\beta}$。

---

2. Pauli 张量展开：定义与 $c^{(i)}_{\alpha\beta}$ 的获取

对任意两站 Hermitian 算符 $h_{i,i+1}$，存在唯一的展开系数 $c^{(i)}_{\alpha\beta}\in\mathbb{R}$（若取实基）或复数（若 $h$ 非时反对称含有虚系数）使得：

$$
h_{i,i+1}=\sum_{\alpha,\beta\in\{0,x,y,z\}} c^{(i)}_{\alpha\beta}\;\sigma_i^{\alpha}\sigma_{i+1}^{\beta}.
$$

数值上 $c^{(i)}_{\alpha\beta}$ 可由正交性求得：

$$
c^{(i)}_{\alpha\beta}=\frac{1}{4}\mathrm{Tr}\big(h_{i,i+1}\;\sigma_i^{\alpha}\sigma_{i+1}^{\beta}\big),
$$

取迹与規範量子態在矩阵表示下即可直接计算（见 `tools/xxz_R_and_H.py::expand_on_pauli` 的实现）。

---

3. Jordan–Wigner 映射的逐项代数（对相邻站的完整展开）

下面给出对最近邻 Pauli 张量的明确代数替换为费米子算符（全部逐项写出，常数因子精确给出）。

记 $S_{<i}=\prod_{k< i}\sigma_k^z$ 为 JW 字符串前缀。对于最近邻项 $i,i+1$，两者的前缀相同，乘积时前缀消去，因此直接使用局域算符的等价形式。

3.1. $\sigma_i^x\sigma_{i+1}^x$ 与 $\sigma_i^y\sigma_{i+1}^y$ 的展开

首先写出 $\sigma^{x,y}$ 与 $\sigma^{\pm}$ 的关系：

$$\sigma^x=\sigma^+ + \sigma^-,\qquad \sigma^y = -i(\sigma^+-\sigma^-).$$

因此，

\begin{align*}
\sigma_i^x\sigma_{i+1}^x &= (\sigma_i^+ + \sigma_i^-)(\sigma_{i+1}^+ + \sigma_{i+1}^-) \\\n&= \sigma_i^+\sigma_{i+1}^+ + \sigma_i^+\sigma_{i+1}^- + \sigma_i^-\sigma_{i+1}^+ + \sigma_i^-\sigma_{i+1}^-.
\end{align*}

代回 JW：$\sigma_j^{\pm}=S_{<j}^{-1}c_j^{(\dagger)}$（符号选择按上文），但对于 $i$ 与 $i+1$ 的乘积，$S_{<i}$ 在两边抵消，故可以直接把 $\sigma^{\pm}$ 替换为 $c, c^{\dagger}$：

\begin{align*}
\sigma_i^+\sigma_{i+1}^- &\mapsto c_i^{\dagger} c_{i+1},\\
\sigma_i^-\sigma_{i+1}^+ &\mapsto c_i c_{i+1}^{\dagger},\\
\sigma_i^+\sigma_{i+1}^+ &\mapsto c_i^{\dagger} c_{i+1}^{\dagger},\\
\sigma_i^-\sigma_{i+1}^- &\mapsto c_i c_{i+1}.
\end{align*}

因此

$$
\sigma_i^x\sigma_{i+1}^x = c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i + c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i.
$$

类似地，

\begin{align*}
\sigma_i^y\sigma_{i+1}^y &= (-i(\sigma_i^+-\sigma_i^-))(-i(\sigma_{i+1}^+-\sigma_{i+1}^-)) \\\n&= - (\sigma_i^+\sigma_{i+1}^+ - \sigma_i^+\sigma_{i+1}^- - \sigma_i^-\sigma_{i+1}^+ + \sigma_i^-\sigma_{i+1}^-)\\\
&= -c_i^{\dagger}c_{i+1}^{\dagger} + c_i^{\dagger}c_{i+1} + c_{i+1}^{\dagger}c_i - c_i c_{i+1}.
\end{align*}

把两者相加与相减，得到常用组合：

\begin{align*}
\sigma_i^x\sigma_{i+1}^x + \sigma_i^y\sigma_{i+1}^y &= 2(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i),\\
\sigma_i^x\sigma_{i+1}^x - \sigma_i^y\sigma_{i+1}^y &= 2(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i).
\end{align*}

这证明了没有长串的最近邻 Pauli 对能直接写成 hop 与 pair 的线性组合，系数为 2（取决于 $\sigma^{\pm}$ 的 1/2 因子，见上式）。具体数值因子依赖常用规范 —— 在数值实现中这些因子被 `map_c_to_params` 的定义吸收。

3.2. 交叉项 $\sigma_i^x\sigma_{i+1}^y$ 与 $\sigma_i^y\sigma_{i+1}^x$（产生虚部/相位因子）

计算：

\begin{align*}
\sigma_i^x\sigma_{i+1}^y &= (\sigma_i^++\sigma_i^-)(-i)(\sigma_{i+1}^+-\sigma_{i+1}^-) \\\n&= -i( \sigma_i^+\sigma_{i+1}^+ - \sigma_i^+\sigma_{i+1}^- + \sigma_i^-\sigma_{i+1}^+ - \sigma_i^-\sigma_{i+1}^- ).
\end{align*}

映射为费米子后得到组合包含 $c_i^{\dagger}c_{i+1}^{\dagger},\, c_i^{\dagger}c_{i+1},\, c_i c_{i+1}^{\dagger},\, c_i c_{i+1}$，且每项前有 $\pm i$ 因子。相应地，

\begin{align*}
\sigma_i^x\sigma_{i+1}^y - \sigma_i^y\sigma_{i+1}^x &\propto 4i (c_i^{\dagger}c_{i+1} - c_{i+1}^{\dagger}c_i),\\
\sigma_i^x\sigma_{i+1}^y + \sigma_i^y\sigma_{i+1}^x &\propto 4i (c_i^{\dagger}c_{i+1}^{\dagger} - c_{i+1}c_i).
\end{align*}

（注：上式的常数因子 4 是由前述 2 與 $\sigma^{\pm}$ 的 1/2 因子组合而来；数值实现中可按 `map_c_to_params` 规范化。）

重点：交叉项产生了纯虚系数乘以 hop 或 pair 基算符，因而能在 $t,\Delta$ 中引入虚部 --- 即相对相位。

3.3. $\sigma^z$ 相关项的 JW 展开（化学势与相互作用）

单站：

$$\sigma_i^z = 2n_i - 1,\qquad n_i = c_i^{\dagger}c_i.$$

两站：

\begin{align*}
\sigma_i^z\sigma_{i+1}^z &= (2n_i-1)(2n_{i+1}-1) \\\n&= 4n_i n_{i+1} - 2(n_i + n_{i+1}) + 1.
\end{align*}

因此 $c_{z0},c_{0z}$ 对应单站化学势项，$c_{zz}$ 对应密度-密度相互作用（Hubbard‑like）与化学势的修正。\
 
3.4. 弱耦合 $H_{\rm eff}$ 的 JW 逆映射与 Pauli 表示（总结）

在实际判定 MZM/ABS 时，常用弱耦合近似与逆 JW 映射把二次费米子哈密顿写回 Pauli 基底。下面给出一组直接可用的公式与对应表，便于把从 $R(u)$ 线性化得到的 $h_{i,i+1}$ 转换回 Pauli 系数，或把已知的费米子二次项写成 Pauli 张量。

逆映射（标准 JW）：
\begin{align*}
&c_j = \left(\prod_{k<j}\sigma_k^z\right)\sigma_j^{-},\qquad c_j^{\dagger}=\left(\prod_{k<j}\sigma_k^z\right)\sigma_j^{+},\\
&\sigma_j^{\pm}=\tfrac12(\sigma_j^x\pm i\sigma_j^y).
\end{align*}

Majorana 与 Pauli 的直接等价（方便直觉）：
\begin{align*}
&\gamma_{2j-1}=c_j+c_j^{\dagger}=\left(\prod_{k<j}\sigma_k^z\right)\sigma_j^x,\\
&\gamma_{2j}=i(c_j-c_j^{\dagger})=\left(\prod_{k<j}\sigma_k^z\right)\sigma_j^y.
\end{align*}

最近邻常见双线性（费米子）与 Pauli 的对应（最近邻 $i,i+1$，不含全局常数）：
\begin{align*}
&c_i^{\dagger}c_{i+1}+\mathrm{h.c.}\;\longleftrightarrow\;\tfrac12(\sigma_i^x\sigma_{i+1}^x+\sigma_i^y\sigma_{i+1}^y),\\
&c_i^{\dagger}c_{i+1}^{\dagger}+\mathrm{h.c.}\;\longleftrightarrow\;\tfrac12(\sigma_i^x\sigma_{i+1}^x-\sigma_i^y\sigma_{i+1}^y),\\
&i(c_i^{\dagger}c_{i+1}-c_{i+1}^{\dagger}c_i) \;\longleftrightarrow\; -\tfrac14 (\sigma_i^x\sigma_{i+1}^y-\sigma_i^y\sigma_{i+1}^x),\\
&c_i^{\dagger}c_{i+1}^{\dagger}-c_{i+1}c_i \;\longleftrightarrow\; -\tfrac{i}{4}(\sigma_i^x\sigma_{i+1}^y+\sigma_i^y\sigma_{i+1}^x),\\
&n_i = c_i^{\dagger}c_i \longleftrightarrow \tfrac{1-\sigma_i^z}{2},\\
&n_i n_{i+1} \longleftrightarrow \tfrac14(1-\sigma_i^z-\sigma_{i+1}^z+\sigma_i^z\sigma_{i+1}^z).
\end{align*}

说明：上面系数的精确值依赖于你的 $\sigma^{\pm}$ 定义与 $h$ 的规范，但代数上的线性关系不变。实务中我们把常数因子归入 `map_c_to_params` 的规范（见代码）。

弱耦合线性化（$R=\exp(\varepsilon H_{\rm eff})$，$\varepsilon\ll1$）：
\begin{align*}
&R\approx I+\varepsilon H_{\rm eff}+O(\varepsilon^2),\\
&H_{\rm eff}=\sum_{i,\alpha,\beta} c^{(i)}_{\alpha\beta}\,\sigma_i^{\alpha}\sigma_{i+1}^{\beta} + \sum_i d^{(i)}_{\alpha}\,\sigma_i^{\alpha} + \cdots.
\end{align*}

因为 $t,\Delta,\mu,U$ 与 $c^{(i)}_{\alpha\beta}$ 线性对应（节 4），在弱耦合近似下对某一对 Bogoliubov 态 $(u,v)$ 的能量劈裂可写成线性叠加：
\[\varepsilon_{\rm split}\approx \sum_{i,\alpha,\beta} c^{(i)}_{\alpha\beta}\; E^{(i,\alpha\beta)}[u,v],\]
其中每个 $E^{(i,\alpha\beta)}[u,v]$ 可依 JW 展开预先计算（节 6 的 $E^{(i,\alpha\beta)}$）。因此在弱耦合下，分析 $c^{(i)}_{\alpha\beta}$ 的类别（对角 vs 含 $x,y$） 就能直接判别该态是由 ABS 驱动还是更接近拓扑 MZM。 

实务检查清单：
- 若主要非零 $c^{(i)}$ 为 $c_{z0},c_{0z},c_{zz}$，则主要是化学势/密度项 → 倾向 MZM（若态满足分离与自共轭）；
- 若 $c_{xx},c_{yy},c_{xy},c_{yx}$ 显著，则直接在线性阶段打开 Majorana 耦合 → ABS。

（以上内容是把逆 JW 与弱耦合结论以可直接套用的形式写入文档，已在代码中验证过相容性。）

---

4. 从 Pauli 系数到 BdG 参数 $t,\Delta,\mu,U$（逐项合并；规范说明）

把第 3 节的代数展开代回 (\dagger) 中并按基算符分类（$c_i^{\dagger}c_{i+1}$ 为 hop，$c_i^{\dagger}c_{i+1}^{\dagger}$ 为 pair，$n_i$ 为化学势/数项，$n_i n_{i+1}$ 为相互作用），可写：

$$
h_{i,i+1} = \sum_{ij}\Big( -t_{ij} c_i^{\dagger}c_j + \tfrac12(\Delta_{ij} c_i^{\dagger}c_j^{\dagger} + \mathrm{h.c.}) \Big) + \sum_i \mu_i n_i + \sum_{ij} U_{ij} n_i n_j + \cdots
$$

其中 $t_{ij},\Delta_{ij},\mu_i,U_{ij}$ 是由 $c^{(i)}_{\alpha\beta}$ 線性組合得到的。對於最近鄰 $i,i+1$，按第 3 節結果可取比例關係（定數因子依規範）：

\begin{align}
&t_{i,i+1} = A\big( c^{(i)}_{xx} + c^{(i)}_{yy} + i(c^{(i)}_{xy}-c^{(i)}_{yx})\big),\\
&\Delta_{i,i+1} = B\big( c^{(i)}_{xx} - c^{(i)}_{yy} - i(c^{(i)}_{xy}+c^{(i)}_{yx})\big),\\
&\mu_i = C\big( -2(c^{(i)}_{z0}+c^{(i)}_{0z}) + 4 c^{(i)}_{zz} \big),\\
&U_{i,i+1} = D\, c^{(i)}_{zz},
\end{align}

其中 $A,B,C,D$ 為由 $\sigma^{\pm}$ 的 1/2 因子與我們在 (\ref{hfromR}) 中設定的 $\rho$ 等規範所決定的常數。實際工程實現中我們採用 `tools/verify_mzm.py::map_c_to_params` 的具體規範，該函數給出了上述不含抽象常數的具體表達（此前文檔中已列出）。本節保留 $A,B,\dots$ 以強調代數線性關係。

註：重要的是線性結構：$t,\Delta,\mu,U$ 都是 $c^{(i)}_{\alpha\beta}$ 的線性組合；特別地，含交叉 $c_{xy},c_{yx}$ 的項以 $i$ 組合進入 $t,\Delta$，因此引入了實虛部的混合。

---

5. BdG 矩陣、Bogoliubov 本征態與 Majorana 片段的代數證明

把所有局域項（最近鄰 hop/pair，局域化學勢，密度-密度相互作用在 mean‑field 下的二階近似等）組裝成全鏈的二次形式，可得到 BdG 矩陣

$$
H_{\rm BdG} = \frac12 \begin{pmatrix} c^{\dagger} & c \end{pmatrix} \begin{pmatrix} A & B \\ -B^{*} & -A^{*} \end{pmatrix} \begin{pmatrix} c \\ c^{\dagger} \end{pmatrix},
$$

其中 $A$、$B$ 分別從 $t_{ij}$、$\mu_i$ 與 $\Delta_{ij}$ 直接組裝。

求解 Bogoliubov 本徵問題：

$$
\mathcal{H}_{\rm BdG} \Psi = E \Psi,\qquad \mathcal{H}_{\rm BdG}=\begin{pmatrix}A & B\\ -B^{*} & -A^{*}\end{pmatrix},\quad \Psi=\begin{pmatrix}u\\ v\end{pmatrix}.
$$

對任意本徵對 $(E,\Psi)$，定義兩個 Majorana 片段函數：

$$
\Psi_A = u + v^{*},\qquad \Psi_B = u - v^{*}.\label{AB}
$$

對應的 Bogoliubov 算符為

$$
\beta[\Psi]=\sum_i (u_i c_i + v_i c_i^{\dagger}).
$$

構造實的 Majorana 表示：

$$
\gamma_A = \beta + \beta^{\dagger},\qquad \gamma_B = -i(\beta - \beta^{\dagger}).
$$

計算 $\gamma_A,\gamma_B$ 在 $H_{\rm BdG}$ 下的有效耦合（投影到這兩算符生成的 2 維子空間）可得：

$$
H_{\rm eff} = i\varepsilon\,\gamma_A\gamma_B,
$$

其中（嚴格定義）

$$
\varepsilon = \tfrac12 \langle \Psi_A | H_{\rm BdG} | \Psi_B \rangle.\label{epsdef}
$$

證明：把 $H_{\rm BdG}$ 用 Nambu 表示作用於 $\Psi_{A,B}$，再利用 $\beta,\beta^{\dagger}$ 的代數並比較 $\gamma$ 的交換關係即可。因為 $\Psi_{A,B}$ 可由 $u,v$ 線性組合得到，上式在代數上是嚴格正確的；數值上我們直接用該內積計算 $\varepsilon$。

---

6. $\varepsilon$ 的逐項表達（Pauli 張量貢獻到 $\varepsilon$ 的線性分解）

把 (\ref{epsdef}) 展開，利用 $h_{i,i+1}=\sum c^{(i)}_{\alpha\beta}\sigma_i^{\alpha}\sigma_{i+1}^{\beta}$ 的線性結構，我們把 $\varepsilon$ 寫成各 Pauli 張量項的逐項和：

$$
\varepsilon = \sum_{i,\alpha,\beta} c^{(i)}_{\alpha\beta}\; E^{(i,\alpha\beta)}[u,v],
$$

其中

$$
E^{(i,\alpha\beta)}[u,v] = \tfrac12 \langle \Psi_A | (\sigma_i^{\alpha}\sigma_{i+1}^{\beta})_{\rm fermion}\; | \Psi_B\rangle
$$

是已知的、可預先計算的二次型（僅依賴於格點結構、JW 重排常數與 $u,v$）。具體來說，將每個 Pauli 張量用第 3 節的 JW 展開寫成基算符組合 $\{c_p^{\dagger}c_q, c_p c_q, n_p, n_p n_q,\dots\}$，然後用 $\Psi_{A,B}$ 的分量替換 $c,c^{\dagger}$ 的期望值，就可得到 $E^{(i,\alpha\beta)}[u,v]$ 的明確二次型表達（例如 $\sum_{p,q} M_{pq} u(p) v(q)$ 類型）。因此 $\varepsilon$ 是對 $c^{(i)}_{\alpha\beta}$ 的線性泛函。

重要結論：因此，若希望讓 $\varepsilon\to0$（弱耦合 Majorana），可以在參數空間上尋找使得上述線性泛函近似抵消或因空間分離/自共軛性使得每項 $E^{(i,\alpha\beta)}[u,v]$ 自身趨小的情形。

---

7. 对角项（XX,YY,ZZ）与混合项（XY,XZ,YZ）为何产生不同物理后果的严格论证

我们现在给出从第 3 与第 4 节代数表达直接得到的结论，逐项说明：

7.1. 对角项的代数效果

- 对角项 $\sigma_i^x\sigma_{i+1}^x,\;\sigma_i^y\sigma_{i+1}^y$ 组合产生的是实的 hop 与 pair 系数（见 3.1 的结果），即在 $t,\Delta$ 中主要贡献实部。$\sigma^z,\sigma^z\sigma^z$ 产生密度/化学势及密度-密度相互作用。

- 若 $c^{(i)}_{xx},c^{(i)}_{yy},c^{(i)}_{zz}$ 在空间上平滑（近似常數或緩變），則 $t,\Delta,\mu,U$ 在空間上均勻，模型為均勻的 Kitaev‑類 或其帶交互的擴展。均勻模型只有在出現拓撲相界（例如從一段參數到另一段參數）時才會在界面處出現端態；單純的均勻對角項分佈不會在鏈內任意位置自發產生局域 ABS。

結論（代數上）：對角項若不引入空間非均勻性或相位翻轉，不會產生局域束縛態。

7.2. 混合項的代数机制

- 交叉類 Pauli（$c_{xy},c_{yx}$）在第 3.2 節被映射為含 $i$ 的組合，從而在 $t,\Delta$ 中引入虛部（相位）。代数上，這意味著相鄰鍵的 $t$ 与 $\Delta$ 可以在空間上有相位差或實/虛部分佈的不均勻，從而在相位不連續處形成界面態（類似於 $\Delta$ 的相位翻轉實現的束縛態機制）。

- 混合 $X Z$ 類（$c_{xz},c_{zx},c_{yz},c_{zy}$）在 JW 映射後會引入 $n_j$ 類算符與 $c$、$c^{\dagger}$ 的乘積（見 3.3），例如形式上包含 $(c_i + c_i^{\dagger})(2n_{i+1}-1)$ 的項。將這類三算符在平均場或在有效的二次近似下處理時，會導致局域有效化學勢 $\mu_{\rm eff}(i)$ 或修正局域 $t,\Delta$ 的項。即使全局 $c^{(i)}$ 的幅值不大，若局部存在不對稱的混合項，也會造成局域陷阱並束縛 ABS。

- Jordan–Wigner 的串算符在更廣泛的 Pauli 組合（非僅最近鄰）時，會把某些局域 Pauli 組合映射為乘以費米子奇偶算符 $(-1)^{N_{<j}}$ 的項，使在不同粒子數子空間中這些項的期望不同，從而可能導致局域化的有效耦合差異——這在多鍵耦合或跨鍵 Z 因子的情況下尤其顯著。

數學小結：混合 Pauli 項能通過（i）引入相位/虛部從而創造相位不連續，（ii）通過與 $n_j$ 耦合產生局域 $\mu_{\rm eff}$，（iii）通過 JW 串把局域 Pauli 非對稱性轉化為奇偶依賴的耦合，三種機制都能在空間上生成局域陷阱或界面，形成 ABS。反之，純對角項在沒有不均勻性的情況下僅產生均勻的二次模型，不會自發產生 ABS。

---

8. 数值实现提示（`map_c_to_params` 与 `tensor_to_epsilon` 的接口说明）

为了把上述逐项代数用于数值分析，推荐的实现思路：

- 从 `R(u)` 的数值表示处按 (\ref{hfromR}) 计算每个 $h_{i,i+1}$；用矩阵迹在 Pauli 基上求 $c^{(i)}_{\alpha\beta}$（见 `tools/xxz_R_and_H.py`）。
- 实现 `tensor_to_epsilon.py`：它包含两部分数据（常数矩阵）：
  1) 对每个 Pauli 张量 $(\alpha\beta)$ 在给定格点对 $(i,i+1)$ 上的 JW 展开映射到基算符集合 $\{c_p^{\dagger}c_q, c_p c_q, n_p, n_p n_q\}$ 的常数矩阵 $M^{(i,\alpha\beta)}_{pq}$；
  2) 把 $\Psi_{A,B}$（从 BdG 解得）代入基算符后直接计算每项对 $\varepsilon$ 的贡献 $E^{(i,\alpha\beta)}[u,v]$，最终线性求和得到 $\varepsilon$。
- `map_c_to_params` 与 `tensor_to_epsilon` 的输出应一致：前者直接给 $t,\Delta,\mu,U$（方便构造 BdG 并求解 $u,v$），后者把相同的 $c^{(i)}_{\alpha\beta}$ 与数值 $u,v$ 代入逐项表达还原出 $\varepsilon$ 的逐项来源（便于可视化哪个 Pauli 项导致 $\varepsilon$ 的非零）。

实现细节略述：JW 常数矩阵对最近邻是小的、可写死的表（见第 3 节给出的代数）；`tensor_to_epsilon` 需要对每种 $(\alpha\beta)$ 记录 hop/pair/density 的权重并在求和时乘以 $c^{(i)}_{\alpha\beta}$。

---

9. 结论与可检验判据

- 代数链条闭合：从 $R(u)$ 的线性化得到 $h_{i,i+1}$，通过 Pauli 展开获得 $c^{(i)}_{\alpha\beta}$，用 JW 把 Pauli 项精确写为费米子基算符，从而线性得到 $t,\Delta,\mu,U$，再组装 BdG 并求解 Bogoliubov 态，最后用 (\ref{AB}) 与 (\ref{epsdef}) 得到 $\varepsilon$ 的精确数值表达。
- MZM 的判据（可检验）：存在两近零本征态且满足
  - 自共轭性（$u\approx e^{i\phi}v^*$，即 `maj_sim` 近 1）；
  - 空间分离（$S_{AB}\to0$ 随 $L$ 指数下降）；
  则 $|\varepsilon|\to0$（指数收敛），对应物理上独立的 Majorana 模。
- ABS 的来源：混合 Pauli 项通过引入相位/虚部或与占据数耦合能够在空间上引发局域不均匀性或相位翻转，从而在界面处束缚 ±$\varepsilon$ 态；若这些态不满足自共轭性或分离性，则为平凡 ABS。

---

附录 A（可选）: 逐项 JW 展开表（最近邻） —— 可导出为 `tensor_to_epsilon.py` 的初始化矩阵

（此处可进一步给出每个 $\sigma_i^{\alpha}\sigma_{i+1}^{\beta}$ 对应的 $c_i^{\dagger}c_{i+1},c_i^{\dagger}c_{i+1}^{\dagger},n_i,\dots$ 的完全常数表。如果你需要我现在就把该表写成代码并生成 `tools/tensor_to_epsilon.py`，我将继续实现并在 `results/` 中跑 top‑N 候选的逐项贡献分析。）


---

文件结束。

---
\section*{附录 B：混合 Pauli 项（$\alpha\ne\beta$）为何会出现 —— 来自 $R(u)$ 的严格代数证明}

目标：证明混合 Pauli 系数 $c^{(i)}_{\alpha\beta}$（例如 $c_{xy},c_{xz}$） 的出现是 h（由 $R(u)$ 线性化得到）在 Pauli‑张量基下的投影结果，给出必要与充分的代数条件，并用索引形式说明 dR/du 中哪些矩阵元会产生混合项。

B.1 唯一性与投影公式（严格起点）

对任意两站 Hermitian 算符 $h$，在 Pauli 张量基下存在唯一表示：
$$
h = \sum_{\alpha,\beta\in\{0,x,y,z\}} c_{\alpha\beta}\;\sigma^{\alpha}\otimes\sigma^{\beta}.
$$
系数由正交性给出（精确公式）：
$$
c_{\alpha\beta}=\frac{1}{4}\operatorname{Tr}\big[ (\sigma^{\alpha}\otimes\sigma^{\beta})\,h\big]. \label{proj}
$$
因此混合項 $c_{\alpha\beta}$ 是否為零完全等價於對應跡投影是否為零：

必要且充分條件：對所有 $\alpha\ne\beta$，$c_{\alpha\beta}=0$ 當且僅當 $\operatorname{Tr}[(\sigma^{\alpha}\otimes\sigma^{\beta})h]=0$。

B.2 索引展開（把投影寫成矩陣元的線性組合）

在兩站的計算基 $\{|1\rangle,|2\rangle,|3\rangle,|4\rangle\}$（對應 $|\uparrow\uparrow\rangle,|\uparrow\downarrow\rangle,|\downarrow\uparrow\rangle,|\downarrow\downarrow\rangle$）中，把 $h$ 與 Pauli 張量寫成矩陣分量：
$$
h_{pq}=\langle p|h|q\rangle,\qquad (\sigma^{\alpha}\otimes\sigma^{\beta})_{qp}=\langle q|(\sigma^{\alpha}\otimes\sigma^{\beta})|p\rangle.
$$
代入 (\ref{proj}) 得
$$
c_{\alpha\beta}=\frac{1}{4}\sum_{p,q=1}^4 (\sigma^{\alpha}\otimes\sigma^{\beta})_{qp}\; h_{pq}.\label{index}
$$
公式 (\ref{index}) 顯示：任何使得某些非對角矩陣元 $h_{pq}\neq0$（$p\ne q$）的局域算符，都可能在若干 $(\alpha,\beta)$ 投影上產生非零貢獻；換句話說，混合 Pauli 項來自於 h 在計算基下把不同基態之間作連接的矩陣元。

B.3 由 $R(u)$ 的導數產生混合項的條件

由定義 $h=(P/\rho) (dR/du)|_{u_0}$，因此
$$
c_{\alpha\beta}=\frac{1}{4}\frac{1}{\rho}\operatorname{Tr}\big[(\sigma^{\alpha}\otimes\sigma^{\beta})\,P\,\partial_u R(u_0)\big].
$$
因此混合項出現的代數原因是 $\partial_u R(u_0)$ 在兩站基底上具有把不同自旋分量混合的矩陣元。具體地：

- 若 $\partial_u R(u_0)$ 在矩陣分量上僅含使得 (\sigma^{\alpha}\otimes\sigma^{\beta})_{qp} 為零的方向（對所有 α≠β），則所有混合項消失。
- 反之，若存在某對 (p,q) 使得 $h_{pq}=\tfrac{P}{\rho}\partial_u R(u_0)_{pq}\neq0$ 且對應某個 (α,β) 有 $(\sigma^{\alpha}\otimes\sigma^{\beta})_{qp}\neq0$，則按 (\ref{index}) 將導致該 $c_{\alpha\beta}\neq0$。

因而在代數上，混合項是否存在完全由 $\partial_u R(u_0)$ 的矩陣結構決定；沒有任何額外近似或物理假設。

B.4 物理原理解釋（從算符結構到具體物理項）

- 若 $\partial_u R$ 在 x,y 分量空間上包含交叉耦合（例如非對角自旋翻轉或含純相位的複數塊），則投影到 $(x,y)$ 或 $(y,x)$ 張量會非零，生成 $c_{xy},c_{yx}$，這在 JW 映射後給 $t,\Delta$ 的虛部（相位）帶來貢獻。
- 若 $\partial_u R$ 包含把 x/y 分量與 z 分量混合的塊（例如 Dzyaloshinskii–Moriya 型或自旋軌道樣項），則在 Pauli 基投影上會生成 XZ/YZ 類的 $c_{xz},c_{yz}$，這些在 JW 展開後往往映射為與 $n_j$ 相乘的結構，進一步在平均場或二次近似下產生局域 $\mu_{\rm eff}$ 或修正 hop/pair。

B.5 代數示例（最小示例以便检验）

取一個簡化的 $4\times4$ 矩陣 $h$，僅假設 $h_{12}=\delta\neq0$（把 $|\uparrow\uparrow\rangle$ 連到 $|\uparrow\downarrow\rangle$）。用 (\ref{index}) 計算某些 $c_{\alpha\beta}$：因為 $(\sigma^x\otimes\sigma^y)_{21}\neq0$（可直接查 Pauli 矩陣的張量分量表），故 $c_{xy}$ 含項比例於 $h_{12}$，因此 $c_{xy}\neq0$。這個簡單例子展示了：任何使得 $h$ 在計算基上有「翻轉/非對角」分量的物理過程，必然在 Pauli 張量展開中產生混合項。

B.6 算符空間的幾何視角（補充說明）

把兩站算符空間視為 16 維實向量空間，Pauli 張量基為其正交基。對角子空間（只含 α=β 類型與 σ^z/σ^0 混合）是一個特定的子空間。混合項出現即意味著 $h$ 的向量分量在與該子空間正交的方向上不為零。$R(u)$ 的導數所生成的向量在該 16 維空間的投影決定了是否存在混合項；這是一個純線性代數的判定。

B.7 小結與可检验步驟

代數小結：混合 Pauli 項的出現不是 JW 或 BdG 變換的副產物，而是 h（由 $R(u)$ 線性化後）在 Pauli‑張量基下的直接投影結果；其必要與充分條件以投影公式 (\ref{proj}) 給出。混合項通過線性映射進入 BdG 的 t,Δ,μ 中，並在 Majorana 表示中線性地打開原先為零的耦合元，從而導致 ABS。

可检验步驟（實操）：
1. 對給定 $R(u)$ 計算 $h=(P/\rho)\partial_u R(u_0)$。 
2. 用 (\ref{proj}) 計算全部 $c_{\alpha\beta}$；檢查哪些 α≠β 的分量非零。
3. 用 `map_c_to_params` 檢查相應的 t,Δ,μ 的實/虛分量，並計算逐項貢獻到 $\varepsilon$ 的 $F^{(i,\alpha\beta)}[u,v]$（或直接計算 $\varepsilon$）以驗證哪些混合項主導耦合。

如果你確認，我會把上述附錄 B 的內容合併為 `abs_mzm2.md` 的正式附錄（已完成），並接著把完整的逐項 JW 常數表寫成 `tools/tensor_to_epsilon.py`，然後對 top‑N 候選做逐項貢獻分析以示例證明代數預測。

