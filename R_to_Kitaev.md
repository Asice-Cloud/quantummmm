# 从 R 到 Kitaev 链（Jordan–Wigner 展开）
## 完整证明：从 $R_{i,i+1}$ 到 $t,\,\Delta,\,\mu$ 的逐步代数推导

目标：在严格的 Jordan–Wigner（JW）映射与链上求和约定下，证明
$$
t\propto(b+c),\qquad \Delta\propto(b-c),\qquad \mu\text{ 的线性部分由 }a,d\text{ 给出（可被重整化）}.
$$ 
并给出常数因子在常见归一化下的具体值。

一、算符
<!-- (这里把张量简写成了乘法的形式，忽略了$\otimes$符号) -->
$$
R_{i,i+1}=a\,I + b\,\sigma^x_i\sigma^x_{i+1} + c\,\sigma^y_i\sigma^y_{i+1} + d\,\sigma^z_i\sigma^z_{i+1}.
$$

### 注：Pauli 矩阵与符号说明
这里的 $\sigma^\alpha$（$\alpha=x,y,z$）是每个格点上的 Pauli 矩阵，作用在单点的二维态空间 $V$ 上。常用矩阵表示（在 $|\uparrow\rangle=(1,0)^T,\;|\downarrow\rangle=(0,1)^T$ 基下）为：
$$
\sigma^x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
\sigma^y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad
\sigma^z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$
定义升降算符
$$\sigma^\pm=\tfrac12(\sigma^x\pm i\sigma^y),\qquad
\sigma^+=\begin{pmatrix}0&1\\0&0\end{pmatrix},\;\sigma^- =\begin{pmatrix}0&0\\1&0\end{pmatrix}.
$$
局域标示 $\sigma^\alpha_i$ 表示該矩阵嵌入到全链 $V^{\otimes L}$ 的第 $i$ 位（即 $I^{\otimes(i-1)}\otimes\sigma^\alpha\otimes I^{\otimes(L-i)}$,
因此项 $\sigma_i^a \sigma_{i+1}^a$
  实际是两个格点上的张量积算符 $\sigma^a \otimes \sigma^a$
 （嵌入到全链上）。
）。它们满足代数关系
$$\{\sigma^\alpha,\sigma^\beta\}=2\delta_{\alpha\beta}I,\qquad [\sigma^\alpha,\sigma^\beta]=2i\epsilon_{\alpha\beta\gamma}\sigma^\gamma.
$$
在 Jordan–Wigner 映射中用到的关系为：
$$\sigma^+_j\mapsto c_j^{\dagger}e^{i\pi\sum_{k<j}n_k},\qquad \sigma^-_j\mapsto c_j e^{i\pi\sum_{k<j}n_k},\qquad \sigma^z_j\mapsto 2n_j-1.
$$
（对最近邻 $i,i+1$，串算符在乘积中抵消，使局域自旋兩體算符在费米子表示下仍為局域項。）

二、用升降算符展开（代数恒等）
记 $\sigma^\pm=\tfrac{1}{2}(\sigma^x\pm i\sigma^y)$，于是等价写法：
$$
\sigma^x=\sigma^++\sigma^-,\qquad \sigma^y=\frac{\sigma^+-\sigma^-}{i}.
$$ 

直接计算得（对相邻站点 $i,i+1$）：
$$
\begin{aligned}
\sigma^x_i\sigma^x_{i+1}+\sigma^y_i\sigma^y_{i+1}&=2\big(\sigma^+_i\sigma^-_{i+1}+\sigma^-_i\sigma^+_{i+1}\big),\\
\sigma^x_i\sigma^x_{i+1}-\sigma^y_i\sigma^y_{i+1}&=2\big(\sigma^+_i\sigma^+_{i+1}+\sigma^-_i\sigma^-_{i+1}\big).
\end{aligned}
$$

三、Jordan–Wigner 映射（最近邻简化）
JW 映射：
$$
\sigma^+_j=c_j^{\dagger}e^{i\pi\sum_{k<j}n_k},\qquad \sigma^-_j=c_j e^{i\pi\sum_{k<j}n_k},\qquad \sigma^z_j=2n_j-1.
$$ 
对最近邻 $i,i+1$，串算符 $e^{i\pi\sum_{k<i}n_k}$ 在乘积中相互抵消，因此有
$$
\begin{aligned}
\sigma^+_i\sigma^-_{i+1}+\sigma^-_i\sigma^+_{i+1}&\mapsto c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i,\\
\sigma^+_i\sigma^+_{i+1}+\sigma^-_i\sigma^-_{i+1}&\mapsto c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i,\\
\sigma^z_i\sigma^z_{i+1}&\mapsto (2n_i-1)(2n_{i+1}-1)=4n_in_{i+1}-2(n_i+n_{i+1})+1.
\end{aligned}
$$


四、逐项代入并识别二次项
代入上述等式到 $R_{i,i+1}$：
$$
\begin{aligned}
R_{i,i+1} &= a\,I + b\,\sigma^x_i\sigma^x_{i+1} + c\,\sigma^y_i\sigma^y_{i+1} + d\,\sigma^z_i\sigma^z_{i+1} \\
&= a\,I + (b+c)\big(c_i^{\dagger}c_{i+1}+c_{i+1}^{\dagger}c_i\big) + (b-c)\big(c_i^{\dagger}c_{i+1}^{\dagger}+c_{i+1}c_i\big)\\
&\quad + d\big(4n_in_{i+1}-2(n_i+n_{i+1})+1\big).
\end{aligned}
$$

在单个键（bond）处，二次费米子项的系数即如上所示：
- hopping (单键处)：系数 $b+c$ * $(c_i^{\dagger}c_{i+1}+h.c.)$；
- pairing (单键处)：系数 $b-c$ * $(c_i^{\dagger}c_{i+1}^{\dagger}+h.c.)$。

因此局域识别给出 $t\propto(b+c)$、$\Delta\propto(b-c)$（比例常数由你对总体哈密顿量的归一化决定）。

五、链上求和与化学势项的精确系数
取常见的做法把哈密顿量定义为所有键的和：
$$
H=\sum_{j=1}^{L-1} R_{j,j+1}\qquad(\text{开链，周期边界类似但要处理双计数}).
$$
把第四项中关于 $d$ 的线性密度部分写出：单键产生 $-2d(n_j+n_{j+1})$。把所有键求和：对于体内站点（非边界），$n_j$ 出现在键 $(j-1,j)$ 和 $(j,j+1)$ 的线性项中，各贡献 $-2d$，合计为 $-4d$。因此链上线性项为
$$
H_{\rm lin}^{(d)} = -4d\sum_{j=1}^L n_j + \text{(仅边界处的修正)}.
$$ 
同时，常数项与 $a$ 的贡献为每个键贡献 $a$；求和后产生整体能级位移 $a\times(L-1)$，该常数可以通过重定义参考能量或并入 $\mu$ 的常数偏移中吸收。

若我们采用惯例写法（常见于 Kitaev 链）
$$H= -\mu\sum_j\Big(n_j-\tfrac12\Big) - t\sum_j\big(c_j^{\dagger}c_{j+1}+h.c.\big) + \Delta\sum_j\big(c_j c_{j+1}+h.c.\big)+\cdots,$$
则将上面链上求和结果与该形式比较，可选择直接取归一化使得
$$t = b+c,\qquad \Delta = b-c,$$
并识别
$$\
mu = 4d + \mu_0,\qquad \mu_0\text{ 来自 }a\text{（常数偏移，可吸收）}.
$$ 
（若你在定义中对 $t$ 或 $\Delta$ 加了额外的负号或 $1/2$ 因子，上述等号右边应乘以相应的归一化常数。）

六、关于 $d$ 项的“可重整化”说明
四费米子项 $4d\,n_in_{i+1}$ 是相互作用项，不属于二次自由理论；但是在低能或平均场近似中，此项会带来两个重要效应：
- 平均场分解 $n_in_{i+1}\to n_i\langle n_{i+1}\rangle + \langle n_i\rangle n_{i+1}-\langle n_i\rangle\langle n_{i+1}\rangle$ 直接生成线性项，等效于修改 $\mu$（即重整化 $\mu$）；
- 在严格的重整化群分析中，四费米子在低维（1D）下可能是边缘/微扰量，会把系统流到不同的低能固定点，从而间接改变有效的二次项系数（包括对 $\mu,t,\Delta$ 的重整化）。

因此说“$\mu$ 由 $a,d$ 的线性项给出（可重整化）”是准确的：$d$ 显式产生链上线性密度项（系数 $-4d$），而 $a$ 产生可以吸收进参考能量或 $\mu$ 的常数偏移；相互作用产生的重整化效应会进一步修改 $\mu$ 的有效值。

结论（规范化约定）
- 直接把 $H=\sum_j R_{j,j+1}$ 作为 Kitaev‑型哈密顿量时，可取 $t=b+c,\,\Delta=b-c$，并且链上 $d$ 的线性贡献给出 $\mu\supset 4d$，$a$ 只产生可吸收的常数偏移。
- 若希望把 $t,\Delta$ 标准化为带负号的常见写法（例如 $-t\sum c_j^{\dagger}c_{j+1}$），只需在等号前乘上你所用的总体负号（即 $t_{paper}=- (b+c)$ 等）。

若你确认要采用哪种哈密顿量整体符号与常数因子（例如是否在定义中包括 $-1$ 前因子或 $1/2$），我可把上面的等号中的比例常数固定并将完整推导替换为与你论文完全一致的记号版本。
$$
t \propto (b+c),\qquad \Delta \propto (b-c),\qquad \mu\text{ 由 }a,d\text{ 的线性项给出（可重整化）}. 
$$

因此：
- 若要得到自由（非相互作用）Kitaev 链，需要消去四费米子项，即要求
  $d=0.$ 
  这时 $R_{i,i+1}$ 仅是二次费米子算符（加上常数项），可直接对应 $t,\Delta,\mu$。
- 若 $d\neq0$ 则得到含相互作用的扩展模型（例如 XXZ‑类的相互作用），不再是自由马约拉那链。

**5. 示例（便于对照）**

- 例 1：取 $a=0,d=0,b=1,c=1$，则
  $$t\propto 2,\quad \Delta\propto 0,$$
  即纯跳跃自由费米子链。对应 Kitaev 链的配对为 0。

- 例 2：取 $a=0,d=0,b=1,c=-1$，则
  $$t\propto 0,\quad \Delta\propto 2,$$
  即纯配对项（最大化 p‑wave 配对）。

**6. 马约拉那表示（说明）**

定义局域马约拉那算符：
$$\gamma_{2j-1}=c_j + c_j^\dagger,\qquad \gamma_{2j} = -i(c_j - c_j^\dagger).$$
二次费米子项在马约拉那表示下是相邻马约拉那的双线性 $i\gamma_a\gamma_b$，在拓扑相区会出现边界上未配对的零能马约拉那模。


**附录：YBE 代数约束与典型物理情形对照表**

下面先回顾常数 Yang–Baxter 给出的代数约束（对复参数可分实、虚部分）：

```
a d (b-c)=0,
b c (a-d)=0,
a b c - a b d - a c d + b c d = 0.
```

下表列出这些约束下若干常见物理情形及其意义：

| 物理情形 | 代数约束/参数关系 | 物理含义与备注 |
|---:|---|---|
| XXX（SU(2) 不变） | $b=c=d=t$ → 约束化为 $t^2(a-t)=0$ | 若 $t\neq0$，需 $a=t$，即 $R\propto P$（置换/交换算符）。这是 XXX 点（可 Baxter 化得到 $R(u)=uI+P$）。 |
| XXZ（平面各向同性） | $b=c\neq d$：代数约束通常要求若 $b\neq0$ 则 $a=d$，进一步退化到 $a=b=c=d$（常数情形下难以得到非平凡的常数 R） | 常见的 XXZ 可积族为带谱参数的情形（Baxter 化，双曲/三角函数），但作为常数 4 参数截面通常受限。 |
| XY / 自由 Kitaev‑样族（无相互作用） | $d=0$ 且常见解为 $a=0,d=0$（即 $R=b\,\sigma^x\sigma^x + c\,\sigma^y\sigma^y$ 为一族解） | 通过 Jordan–Wigner 映射得到二次费米子：跳跃 $t\propto(b+c)$，配对 $\Delta\propto(b-c)$。这是 free‑fermion（Kitaev/XY）类模型，适合马约拉那对角化。 |
| 纯跳跃或纯配对（单项） | 仅 b 或仅 c 非零（其余为 0） | 简单独立耦合，往往满足 YBE（方程大多退化为 0），对应纯 hopping 或纯 pairing 的局域项。 |
| 只有 $\sigma^z\sigma^z$（Ising‑型） | $b=c=0$，任意 $a,d$ 满足约束 | 对应 Ising/XXZ 中的相互作用项（四费米子项在 JW 后出现）。可包含化学势偏移与最近邻相互作用。 |
| 平凡/退化情形 | 任意两个或更多参数为 0，如 $a=0,d=0$、或 $b=0,c=0$ 等 | 多数代数方程自动满足，给出容易处理的可积或可对角化子空间。

**小结**：代数约束表明可积（满足常数 YBE）的常数四参数 $R$ 是受限的——非平凡的 SU(2)‑不变点（XXX）与某些退化/单项耦合属于解集；而一般同时包含跳跃与配对且含相互作用（任意 $b,c,d$）的自由 Kitaev 链并不总是出现在常数 R 的解集中（但可在带谱参数族或通过其他映射下出现）。

如果要，我可以把上表以更严谨的代数推导（带代数推导步骤）写入文档或把其中的若干代表性参数在 1D 链上做数值/解析对角化示例。

## YBE 约束与物理情形的对应关系（补充讨论）

目标：把上面给出的代数约束与物理上出现的相、可积性与自由/相互作用情形联系起来，提供直观判据和代数解释。

- 本质含义：常数 Yang–Baxter 方程保证两体散射矩阵的可分解性，从而构造出一族相互可交换的迁移矩阵（transfer matrices）与无穷多守恒量；满足 YBE 的 $R$ 对应可积格点模型或可积量子链的局部构件。

- 代数约束的解分支与物理限制：给定约束
  ```
  a d (b-c)=0,
  b c (a-d)=0,
  a b c - a b d - a c d + b c d = 0.
  ```
  这些方程是乘积/多项式形式——意味着若要非平凡解，至少有若干参数满足特殊关系（例如某个参数为零或几者相等）。每种代数分支对应不同的物理情形：
  - 若 $d=0$：四费米子相互作用项消失，R 映到的费米子哈密顿是二次的（自由或 mean‑field 可对角化），典型为 XY / Kitaev‑样族（自由马约拉那可对角化）。
  - 若 $b=c$ 且非零：交换對稱性增强（平面各向同性），常出现 XX/XXZ/XXX 类族；若同时 $b=c=d$ 并满足 $a=b$ 则退化为置换算符（XXX 点）。
  - 若 $b=0$ 或 $c=0$：只剩纯配对或纯跳跃分量（或纯 Ising 的 $\sigma^z_i \sigma^z_{i+1}$）；这些退化情形常使 YBE 方程大量项消失，从而容易成为常数解。

- 自由‑费米 vs 相互作用的判据：从 JW 展开看，$d\neq0$ 会产生最近邻四费米子 $4d n_i n_{i+1}$，这直接把系统从二次（free）拉入相互作用范畴。故常数 R 要同时给出可积且自由的 Kitaev 链，通常需要位于 $d=0$ 的分支（或通过特殊代数关系使四费米子项可被映射掉/退化）。

- 带谱参数的情形（Baxter 化）：常数 R 解集受限，许多重要可积模型（如 XXZ）并不是常数 R 的非平凡解，而是通过引入谱参数 $u$ 得到 $R(u)$（Baxter 化）。在这种情形下，参数间的代数关系由解析函数（双曲/三角）调节，物理上对应可调的各向异性或耦合强度。

- 物理后果与可解技术：若参数位于满足 YBE 的分支——无论是常数还是带谱参数——其格点模型/量子链可用代数 Bethe ansatz、傳輸矩陣法等可积技术求解；若不满足 YBE，则通常需依靠数值或近似手段，且守恒量减少，动力学更丰富（例如散射会产生不可分解的三体过程）。

- 实用判据（工程角度）：给定具体的 $R(a,b,c,d)$，可以按下列顺序快速判断物理情形：
  1. 检查是否 $d=0$：若是，则候选为自由‑费米/XY/Kitaev‑类；
  2. 检查是否 $b=c$：若是，考虑 XX/XXZ/XXX 对称性；
  3. 检查 $b=0$ 或 $c=0$：对应纯配对或纯 hopping（退化可积）；
  4. 若以上均不成立，查看三式能否通过特殊常数组合满足（可能仍是可积但属于更窄的常数解集），否则通常为非可积/含相互作用的模型。

结语：YBE 的多项式约束把参数空间划分为若干代数簇，每个簇对应不同的物理对称性和相互作用结构。要从代数到物理的完整映射，通常需要结合 JW 展开（把矩阵元映到费米子算符）、链上求和的计数以及在必要时对相互作用项做平均场或 RG 分析。若你想，我可以对其中 2–3 个有代表性的代数分支（例如 $d=0$ 的自由分支和 $b=c\neq d$ 的 XXZ 分支）给出更详细的代数推导和简明谱图示例。


## Majorana 耦合的预测与判据（补充猜想）

基于本文件中映射 $t=b+c,\;\Delta=b-c,\;\mu\approx4d+\mu_0$，下面给出简明且可操作的判据与耦合形式猜想：

- 拓扑相与零模存在（无相互作用）：BdG 频谱
  $$E(k)=\pm\sqrt{(-2t\cos k-\mu)^2+(2\Delta\sin k)^2}$$
  表明能隙在 $\mu=\pm2t$ 处闭合。故拓扑相判据可写为
  $$
  |\mu|<2|t|\quad\text{且}\quad\Delta\neq0.
  $$ 
  代入你的参数为
  $$|4d+\mu_0|<2|b+c|\quad\text{且}\quad b-c\neq0.$$

- Majorana 之间耦合的典型形式：
  - 有限长度导致端点 Majorana 波函数重叠，能量裂分量级为
    $$\delta\varepsilon\sim A e^{-L/\xi}\cos(k_F L),$$
    其中 $\xi\sim v/\Delta$，$v$ 为有效速度，$A$ 由归一化与边界条件决定。
  - 局域单粒子扰动（例如 onsite 化学势或 next‑nearest hopping/pairing）生成双线性项
    $$H_{\rm coup}=i\epsilon\,\gamma_a\gamma_b,$$
    其中 $\epsilon$ 与扰动强度及波函数重叠成正比；若 $\epsilon\gtrsim\delta\varepsilon$，零模将被抬升。

- 参数边界与配对/跃迁竞赛：在马约拉那表示下，邻近 Majorana 的双线性系数由 $t\pm\Delta$ 的组合给出，因此当 $t\approx\pm\Delta$ 时会改变平衡配对格局与零模的空间分布。

- 相互作用 ($d\neq0$) 的作用：$d$ 同时带来线性项（链上贡献为 $-4d\sum n_j$）和四费米子 $4d n_i n_{i+1}$，后者可通过平均场或 RG 重整化出对 $\mu,t,\Delta$ 的修正，从而间接影响零模稳定性；因此含 $d$ 情形需用平均场/DMRG 等数值方法验证。

- 实用流程（预测/检验）：
  1. 计算 $t=b+c,\;\Delta=b-c,\;\mu=4d+\mu_0$.  
  2. 若 $\Delta\approx0$ 或 $|\mu|>2|t|$：无拓扑边界零模。  
  3. 若 $|\mu|<2|t|$ 且 $\Delta\neq0$：预计存在边界 Majorana，估算 $\xi$ 并用 $\delta\varepsilon\sim e^{-L/\xi}$ 估计端点耦合。  
  4. 若 $d\neq0$：做平均场或数值 RG/DMRG，评估四费米子对零模的重整化/破坏。

我可以把上述判据做成数值脚本：给定 $(a,b,c,d)$ 自动计算 BdG 谱、边界态波函数、局域化长度 $\xi$ 与端点裂分 $\delta\varepsilon$，并绘出相图。要我现在生成该脚本示例吗？

## 把 $R$ 推广为 $R_{i,i+1}$ 的理由（补充说明）

技术上，格点体系的全 Hilbert 空间为单点态空间 $V$ 的张量积 $V^{\otimes L}$。任何作用在两点上的算符
$$R:\;V\otimes V\to V\otimes V$$
可唯一嵌入到全链为局域算符：
$$
R_{i,i+1}=I^{\otimes(i-1)}\otimes R\otimes I^{\otimes(L-i-1)}.
$$ 
这就是把 ``作用在第 i, i+1 位'' 用张量积表示的标准做法，数学上完全自然且唯一。

代数层面：若常数 $R$ 满足 Yang–Baxter 方程，则这些局域嵌入满足局域 YBE（三体交換一致性），从而给出 braid‑group 的表示和可积传输矩阵的局部生成元——这正是把代数生成元 $b_i$ 映为 $\mathrm{id}^{\otimes(i-1)}\otimes R\otimes\mathrm{id}^{\otimes(n-i-1)}$ 的代数依据。

关于 Jordan–Wigner 的技术点：自旋算符映射到费米子时通常出现串算子 $e^{i\pi\sum_{k<j}n_k}$，但对于最近邻两点 $i,i+1$，这些串因子在算符乘积中互相抵消，因此 $R_{i,i+1}$ 映射后仍是局域的二体/四体费米子算符（没有产生非局域长串）。这保证了把自旋层面的 $R_{i,i+1}$ 直接翻译为 Kitaev/XXZ‑类局域项的可行性。

物理注意事项：张量嵌入本身并不自动賦予 $R_{i,i+1}$ 為可实现的“交换”操作——若要把代数上的 braid 表示解释为物理上的拓扑交换（用于量子门），还需满足额外条件，例如 $R$ 的幺正性（或在低能简并子空间内是幺正的）、$R$ 的作用限制在基态简并子空间内，以及系统的谱允许用该局域操作在实验上实现该演化。

如果你同意，我把这段放在文档末尾，或你希望把它移到前文某节（例如在引入 $R$ 的地方）以便读者立即看到？

