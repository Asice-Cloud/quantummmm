(**引理 D 详细证明：quasi‑adiabatic 生成元的局域性、范数估计与 Trotter 误差上界**)
	$$
	\Big\|\mathcal T e^{-i\int_0^T K(t)dt} - \prod_{j=1}^m e^{-iK_j\delta t}\Big\| \le C_T\frac{T}{m}\sup_{j,j'}\|[K_j,K_{j'}]\| + O\Big(\frac{T^2}{m^2}\max\|K\|^3\Big).
	$$
	交换子的范数可用局域分解与 $k_*$ 上界估计：若 $K_j=\sum_\alpha k_{\alpha}^j$，则
	$$
	\|[K_j,K_{j'}]\| \le \sum_{\alpha,\beta}\|[k_{\alpha}^j,k_{\beta}^{j'}]\| \lesssim C''\,k_*^2
	$$
	（非重叠支撑对 commutator 为零，重叠项数由局域结构给出常数因子）。因此
	$$
	常用选择使得 $\hat W(\omega)=\int e^{i\omega s}W(s)ds$ 在 $|\omega|<\Delta/2$ 附近等于 $-1/\omega$，以产生把基态子空间沿时间演化平移的生成元。

步骤 1：生成元的几何与局域性（准局域性）

1.1. 由于 $\dot P$ 为局域算符（由 $H$ 的局域性与有限导数保证），且 $e^{iH s} O e^{-iH s}$ 随 $s$ 的增长以 Lieb–Robinson 光速扩散，卷积核 $W(s)$ 的快速衰减保证了积分产生的一个准局域算符。更精确地，对任何局域算符支撑半径为 $r$ 的 $O$，存在常数 $c,\mu,v_{LR}$，使得
	$$
	\|[e^{iHs} O e^{-iHs}, X]\|\le c\,\|O\|\,\|X\|\,e^{-\mu(d-r-v_{LR}|s|))}
	$$
	（$d$ 为算符支撑间距）。把此不等式与 $W(s)$ 卷积可以导出：若截断 $K(t)$ 的支撑到半径 $R_k$（即只保留与基点距离小于 $R_k$ 的项），截断误差为
	$$
	\|K(t)-K_{R_k}(t)\| \le C' e^{-\mu' R_k},
	$$
	其中 $\mu'$ 与 $\mu$、$W$ 的衰减常数以及 $v_{LR}$ 有关。此即所谓空间截断的指数小尾项（Lieb–Robinson tail）。

步骤 2：生成元每项的范数上界 $k_*$

2.1. 从 $K(t)$ 的表达式及留数/傅里叶分析，可把 $K(t)$ 展为局域项的和：
	$$
	K(t)=\sum_{\alpha} k_{\alpha}(t),
	$$
	每个 $k_{\alpha}$ 支撑在直径 $\lesssim R_k$ 的小区域内。对单项范数上界的常规估计给出
	$$
	\|k_{\alpha}(t)\| \le C_{\alpha}\int |W(s)|\,\|e^{iHs}\dot P e^{-iHs}\|\,ds \lesssim C_{\alpha}'\frac{\|\dot H\|}{\Delta},
	$$
	其中使用了 $\dot P$ 与 $\dot H$ 的线性关系（通过投影导数的解析表示，见下），以及 $W$ 的频域特性导致一个 $1/\Delta$ 因子。于是定义
	$$
	k_*:=\sup_{t,\alpha}\|k_{\alpha}(t)\| \lesssim C_{\mathrm{qa}}\frac{\|\dot H\|}{\Delta}.
	$$
	常数 $C_{\mathrm{qa}}$ 含有 $W$ 的 L1 范数与局域性几何因子。

步骤 3：Trotter 时间离散化误差

3.1. 把时间分为 $m$ 段，每段长度 $\delta t=T/m$，令 $K_j:=K(t_j)$（取中点或左端点均可），考虑近似
	$$
	\mathcal T e^{-i\int_0^T K(t)dt} \approx \prod_{j=1}^m e^{-iK_j\delta t}.
	$$
	使用基础的时间切分误差分析或 Duhamel 展开与交换子估计，可得到二阶型（或一阶）误差上界：
	$$
	\Big\|\mathcal T e^{-i\int_0^T K(t)dt} - \prod_{j=1}^m e^{-iK_j\delta t}\Big\| \le C_T\frac{T}{m}\sup_{j,j'}\|[K_j,K_{j'}]\| + O\Big(\frac{T^2}{m^2}\max\|K\|^3\Big).
	$$
	交换子的范数可用局域分解与 $k_*$ 上界估计：若 $K_j=\sum_\alpha k_{\alpha}^j$，則
	交换子的范数可用局域分解与 $k_*$ 上界估计：若 $K_j=\sum_\alpha k_{\alpha}^j$，则
	$$
	\|[K_j,K_{j'}]\| \le \sum_{\alpha,\beta}\|[k_{\alpha}^j,k_{\beta}^{j'}]\| \lesssim C''\,k_*^2
	$$
	（非重叠支撑对 commutator 为零，重叠项数由局域结构给出常数因子）。因此
	$$
	\Big\|\mathcal T e^{-i\int_0^T K(t)dt} - \prod_{j=1}^m e^{-iK_j\delta t}\Big\| \le C_D\frac{T k_*^2}{m} + O\Big(\frac{T^2 k_*^3}{m^2}\Big).
	$$

步骤 4：把局域指数分解为基本门（R 门）与累积误差

4.1. 每个局域指数 $e^{-i k_{\alpha}\delta t}$ 由有限数目的 $R$ 类门（本研究中 $R=e^{iH_P^{(ij)}}$）近似实现，近似误差可由 Trotter‑Suzuki 或 Lie–Trotter 展开控制；若对每个局域项采用常数轮次的门分解，则总门数与局域项数成正比，累积误差为每项门分解误差乘以项数，故不改变时间离散误差的主要参数依赖（$m,k_*,T$）。

步骤 5：空间截断误差（Lieb–Robinson tail）

5.1. 如步骤 1 所述，将每个 $k_{\alpha}$ 截断到半径 $R_k$ 的局域算符会产生尾项，其范数界为
	$$
	\|k_{\alpha}-k_{\alpha}^{(R_k)}\| \le C_3 e^{-\mu R_k}.
	$$
	把所有项加总得到总截断尾项
	$$
	\|K-K^{(R_k)}\| \le C_4 e^{-\mu R_k}
	$$
	并可把该尾项对时间演化的影响用 Dyson / Duhamel 展开证明为同阶指数小量，故在最终误差中以加项的形式出现。

总结（引理 D 结论）

综上，存在常数 $C_D,C'_D,\mu>0$ 使得
	$$
	\varepsilon_{\mathrm{Trotter}}:=\Big\|\mathcal T e^{-i\int_0^T K(t)dt} - \prod_{j=1}^m \prod_{\alpha} e^{-i k_{\alpha}^j\delta t}\Big\| \le C_D\frac{T k_*^2}{m} + C'_D e^{-\mu R_k} + O\Big(\frac{T^2 k_*^3}{m^2}\Big),
	$$
	且
	$$
	k_*\lesssim C_{\mathrm{qa}}\frac{\|\dot H\|}{\Delta},
	$$
	如需更严格的常数可把上述推导中的几何因子、$W$ 的 L1 范数与 Lieb–Robinson 常数写出，并用具体模型数值计算。
