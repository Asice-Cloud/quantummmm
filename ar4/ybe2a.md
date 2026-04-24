下面给出把该对角指数族代入常数（无谱参）Yang–Baxter 方程后得到的完全代数约束与解集分类（复系数情形）。

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

下面可以選擇將完整的符號化代入與因式分解步驟追加為附錄——告訴我要替換原文段落還是在文末追加推導即可。 