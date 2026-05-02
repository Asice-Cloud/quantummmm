**代数推导：在当前模型下的 ABS 表征**

本节给出在本仓库采用的 Pauli→BdG 映射下，Andreev bound states (ABS) 的严格代数判据，并给出可检验的算法步骤。

1) Pauli → BdG 的代数映射（引用 `tools/verify_mzm.py`）

-	$t = c_{xx} + c_{yy} + i(c_{xy}-c_{yx})$
-	$\Delta = c_{xx} - c_{yy} - i(c_{xy}+c_{yx})$
-	$U = 4 c_{zz}$
-	$\mu = 4 c_{zz} - 2(c_{z0}+c_{0z})$

2) BdG 矩阵與块结构

令 Nambu 基底为 $\Psi=(c_1,...,c_N,c_1^\dagger,...,c_N^\dagger)^T$，构造单粒子 BdG 矩阵
$$
H_{\mathrm{BdG}}=\begin{pmatrix} A & B\\ -B^* & -A^T \end{pmatrix},
$$
其中（本仓库的近邻 p‑wave 约定）
$A_{j,j}=-\mu,\;A_{j,j+1}=A_{j+1,j}=-t$，$B_{j,j+1}=\Delta,\;B_{j+1,j}=-\Delta$.

3) 均匀体色散与子隙条件

在平移不变段取 Bloch 波 $u_j=u e^{ikj},\;v_j=v e^{ikj}$，得块 2×2 方程，谱满足
$$
E(k)^2=\varepsilon(k)^2+|\Delta(k)|^2,
$$
其中 $\varepsilon(k)=-\mu-2t\cos k$，且由近邻 p‑wave 结构得 $\Delta(k)=2i\Delta\sin k$.

4) 指数衰减解与复波数代换（界面外解）

对界面束缚态考虑复波数 $k=i\kappa$（或 $\pi+i\kappa$），记 $X=\cosh\kappa,\;S=\sinh\kappa$，代入得到
$$
E^2 = (-\mu - 2 t X)^2 + 4|\Delta|^2 (X^2-1).
$$
整理为关于 $X$ 的二次方程（简洁形式）：
$$
(t^2+|\Delta|^2) X^2 + (\mu t) X + \frac{\mu^2-4|\Delta|^2-E^2}{4}=0. \quad(\star)
$$
给定 $E$（$|E|<$局部能隙）解出 $X$，并取满足 $X>1$ 的根对应指数衰减常数 $\kappa=\operatorname{arccosh} X>0$。

5) 界面匹配与 ABS 判据

设左右两段参数分别为 $(t_L,\Delta_L,\mu_L)$ 與 $(t_R,\Delta_R,\mu_R)$，在界面区域用左右衰减解构造有限自由度的线性系统（在界面相邻的若干格点寫出 BdG 差分方程，未知系数為左右衰減幅度及界面局部自由度）。存在非平凡衰減解等價於匹配矩陣 $M(E)$ 的行列式為零：
$$
\det M(E)=0. \quad(\dagger)
$$
因此 ABS 的严格代数判据為：存在 $E$（$|E|<$左右段能隙）同時滿足（$\star$）對左右兩段均有 $X>1$ 的解，且界面匹配方程（$\dagger$）成立。

6) 算法化检查（針對給定 $c_{\mu\nu}$）

- 计算 $t(c),\Delta(c),\mu(c)$。
- 计算局域能隙 $\min_k E(k)$ 并在 $|E|<$ gap 區間內搜索 $E$。
- 对每个候选 $E$ 解（$\star$）得到左右 $X_{L,R}$，選取 $X>1$ 的衰減常数並構造界面處的有限線性系統，計算 $\det M(E)$ 的值；若存在根，判定為 ABS。

7) 與 MZM 的代數區分

- MZM 對應 $E=0$ 的特例，且需滿足自共軛性與拓撲約束（例如 Kitaev 鏈的參數區使得 $E=0$ 的解在無窮遠衰減且被拓撲保護）。
- 一般 ABS 則對 $E\neq0$ 且對局部參數敏感；代數上表現為 $\det M(E_0)=0$ 但 $E_0\neq0$，微小局部擾動通常會移動 $E_0$。

8) 限制與備註

- 若原始 Pauli 哈密頓包含顯著交互項 $U=4c_{zz}\neq0$（多體四費米），則單粒子 BdG 描述不是精確的，需要先做 mean‑field 近似或更高階處理；上述嚴格判據僅在二次（free‑fermion）可表示時完全嚴格。
- 對短約瑟夫森結或短連接的特殊情形，$\det M(E)=0$ 可化簡為熟知的閉式表達（例如短結透射率公式），見 `steps.md` 中的短結討論。

---

若你確認無誤，我會接著用 SymPy 在短接口（1–2 個界面格點）上符號化地構造 $M(E)$ 並給出顯式多項式與數值根；或直接實作數值腳本做 S–N / SNS 示範。請回覆“繼續實現”以讓我開始編碼與計算。

**数值示例输出（已生成）**

我已实现并运行 `tools/compute_abs.py`，生成两个示例输出：

- **S–N 接口谱图**: [results/abs_sn_spectrum.png](results/abs_sn_spectrum.png)
- **短约瑟夫森结能量随相位的扫描**: [results/abs_josephson_phi.png](results/abs_josephson_phi.png)

脚本位置： [tools/compute_abs.py](tools/compute_abs.py)

下一步我会把选定数值图与简短诊断追加到 `abs.md` 和 `steps.md`（正在进行中）。

