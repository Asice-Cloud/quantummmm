对角形式R

$$
R_{i,i+1} = exp(i \sum_{u=0}^3 c_{\mu \mu} \sigma^\mu_{i} \otimes \sigma^\mu_{i+1} ) \\
\sigma^0 = I, c_{\mu \mu} \in C
$$

<!-- 下面给出把该对角指数族代入常数（无谱参）Yang–Baxter 方程后得到的完全代数约束与解集分类（复系数情形）。

记原始参数为
$$
c_{00}\quad(\text{整体相位}),\qquad b_x=c_{xx},\;b_y=c_{yy},\;b_z=c_{zz}\in\mathbb C.
$$
定义局域生成元
$$
H_P = c_{00}I + b_x S_x + b_y S_y + b_z S_z,
\qquad S_\alpha=\sigma^\alpha\otimes\sigma^\alpha.
$$

由于 $S_\alpha^2=I$ 且 $[S_\alpha,S_\beta]=0$（对角项两两对易），指数可写为乘积：
$$
R=e^{iH_P}=e^{i c_{00}}\prod_{\alpha\in\{x,y,z\}}(\cos b_\alpha + i\sin b_\alpha\,S_\alpha).
$$

去掉整体相位（$c_{00}$ 在常数 YBE 中可约去），把 $R$ 在基 $\{I,S_x,S_y,S_z\}$ 展开得：
$$
\widetilde R = A_0 I + A_x S_x + A_y S_y + A_z S_z,
$$
其中（记 $c_\alpha=\cos b_\alpha,\ s_\alpha=\sin b_\alpha$）
$$
\begin{align*}
A_0&=c_x c_y c_z + i\,s_x s_y s_z,\\
A_x&=i\,s_x c_y c_z + c_x s_y s_z,\\
A_y&=i\,c_x s_y c_z + s_x c_y s_z,\\
A_z&=i\,c_x c_y s_z + s_x s_y c_z.
\end{align*}
$$

常数 YBE（braid 形式）
$$
\widetilde R_{12}\widetilde R_{23}\widetilde R_{12}=\widetilde R_{23}\widetilde R_{12}\widetilde R_{23}
$$
在本 4 维基底下等价于下列三个代数多项式（记 $a=A_0,\ b=A_x,\ c=A_y,\ d=A_z$）：
$$
\begin{align}
	ext{(E1)}\;& a d (b-c)=0,\\
	ext{(E2)}\;& b c (a-d)=0,\\
	ext{(E3)}\;& a b c - a b d - a c d + b c d = 0.
\end{align}
$$

上述三式与 $a,b,c,d$ 的三角表达式（关于 $b_x,b_y,b_z$ 的 cos/sin）构成等价的方程组；下列为其完整代数分类（复数域上的全集）:

1) 非退化连续族（主要的非平凡族）
- 若 $a\neq0$ 且 $b,c,d\neq0$，由 (E1),(E2) 得 $b=c$ 且 $a=d$；代入 (E3) 再得 $a=b$，因此
	$$a=b=c=d\quad(\text{任意复数}).$$
	在指数参数上等价於
	$$b_x\equiv b_y\equiv b_z\pmod{\pi},$$
	这对应 $H_P$ 可写成常数加置换算符 $P$ 的情形，得到连续一族 $R\propto\cos(2J)I+i\sin(2J)P$。

2) 退化 / 离散族（至少有一个基系数为 0 或乘积为 0）
- 若 $b\neq c$，則從 (E1) 必有 $a d =0$；若 $a\neq d$，則從 (E2) 必有 $b c =0$。這導致以下子簇：
	- $a=0$（即 $A_0=0$）時，(E3) 要求 $b c d=0$：因此至少一個 $A_\alpha=0$，例如
		* $a=d=0$：任意 $b,c$（$\widetilde R=bS_x+cS_y$）及其循環置換；
		* $a=0,b=0$：$\widetilde R=cS_y+dS_z$ 等。
	- $d=0$ 對稱地給出 $a b c=0$ 的子簇。
	- 若某一 $A_\alpha=0$（例如 $A_x=0$），代入其顯式表達式得到關於 $b_x,b_y,b_z$ 的複三角代數方程（包括 $b_\alpha\in\pi\mathbb Z$ 或 $b_\alpha\in\tfrac\pi2+\pi\mathbb Z$ 的離散點，或滿足連續代數關係的複子集）。

3) 極退化/單項與線性組合
- 任何僅有一項非零的情況（$\widetilde R=aI$ 或 $\widetilde R=bS_x$ 等）均滿足 YBE；
- 含 $I$ 與單一 $S_\alpha$ 的任意線性組合（例如 $aI+dS_z$，及其循環置換）也在上面的退化子簇中，滿足方程。

判據（用於直接檢驗某組 $b_x,b_y,b_z$ 是否滿足 YBE）：
- 計算 $c_\alpha=\cos b_\alpha,\ s_\alpha=\sin b_\alpha$（複三角函數），構造 $A_0,A_x,A_y,A_z$；檢驗 (E1),(E2),(E3) 是否成立。
- 若希望更直接的三角條件：非退化連續族當且僅當 $\sin(b_y-b_x)=\sin(b_z-b_y)=0$（即 $b_x\equiv b_y\equiv b_z\pmod\pi$）；否則參數必須滿足上述退化簇的等式（在 $A_\alpha=0$ 的子簇中可進一步化簡為具體的 cos/sin 代数式）。

示例快速對照：
- $b_x=b_y=b_z=J$：連續非平凡族。 
- $(b_x,b_y,b_z)=(0,0,0)$：$R=I$（平凡）。
- $(\tfrac\pi2,\tfrac\pi2,0)$：$R=S_z$（單項）。

注：$c_{00}$（若視為實）僅產生整體相位 $e^{i c_{00}}$，對常數 YBE 無約束（可模 $2\pi$ 任意）。

下面可以選擇將完整的符號化代入與因式分解步驟追加為附錄——告訴我要替換原文段落還是在文末追加推導即可。 -->

**附錄：从指数形式直接因式分解得到的约束（摘要）**

我们对 $R=e^{i(c_{00}I+b_xS_x+b_yS_y+b_zS_z)}$ 直接计算 $R_{12}R_{23}R_{12}-R_{23}R_{12}R_{23}$ 并对矩阵元做因式分解，可将所有非平凡约束写为若干可识别的因子相乘为零（去掉不可约的全局相位因子）。 用记号
$$
X=e^{2i b_x},\quad Y=e^{2i b_y},\quad Z=e^{2i b_z},
$$
这些因子的主要类型为：

- $e^{2i b_x}-e^{2i b_y}=0$，等价于 $\sin(b_x-b_y)=0$（即 $b_x\equiv b_y\pmod\pi$）；
- $e^{2i(b_x+b_y)}-1=0$，等价于 $\sin(b_x+b_y)=0$（即 $b_x\equiv- b_y\pmod\\pi$）；
- F_-(X,Y,Z):= $XY - XZ - YZ - 1 = 0$；
- F_+(X,Y,Z):= $XY + XZ + YZ - 1 = 0$.

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

