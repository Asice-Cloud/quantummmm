YBE -> real physics model

$$
R(u)=\mathcal T_u\exp\Big(-i\int^{u} H(u')\,du'\Big).\\
H(u)=i\,\partial_u R(u)\,R(u)^{-1}.
$$

H 分解式：
$$
H=\sum_{\alpha,\beta\in\{0,x,y,z\}} h_{\alpha\beta}\,\sigma^\alpha\otimes\sigma^\beta.
$$
要产生 MZM 需要二次项（XX, YY, XY, YX 与 Z 类型），而不应含有四次项或带串算符的多体项。

论文中的观点：能谱 $E$ 控制动力学相位，通常可以把演化写成几何与动力相乘积
$$
U=U_{\rm geom}\,e^{-i\int E(t)dt}.
$$

物理形式上可把 $R$ 写为带连接的路径有序指数：
$$
R=\mathcal P e^{\int A(\lambda)}.
$$
因此不能把 $H(u)$ 简单等同为 Berry 连接；要区分 $H$ 本征结构中的几何相（Berry）与动力学相（能谱 $E$）。

现在模型 vs 真实 braid:
- 任意 $H$ \to Berry connection（形式上可写出）
- 任意路径 \to adiabatic manifold（是否可实现取决于谱隙与路径）
- 无能谱控制 \to $E\approx0$

能谱用来区分 MZM 与 ABS：
- MZM: $E\approx0$ 且受拓扑/局域性保护；
- ABS: $E\neq0$（或仅在有限系统中近似零能）；
- trivial: 有较大能隙，无零模。

    

下列内容把从谱参数依赖的两体算子 $R(u)$ 到物理生成元 $H(u)$、瞬时能谱 $E_n(u)$ 与 Berry 连接 $A(u)$ 的严格数学关系与数值实现要点写明，便于把代数对象对应到 BdG / Majorana 判据上。

## 1 概述与目标

目标是给出明确定义与可数值实现的管线：

$$
R(u)=\mathcal T_u\exp\Big(-i\int_{u_0}^u H(s)\,ds\Big),
$$

以及在可逆/正则化条件下如何从 $R(u)$ 恢复瞬时生成元 $H(u)$，并进一步得到瞬时本征能谱 $E_n(u)$、Berry 连接 $A(u)$，最终判定某一路径或参数组是否能生成拓扑 Majorana 零模（MZM）。

## 2 严格定义：R、H

- 路径序列定义（可逆要求）：若对区间 $[u_0,u_1]$ 上每个 $u$，$R(u)$ 可逆，则定义
$$
R(u)=\mathcal T_u\exp\Big(-i\int_{u_0}^u H(s)\,ds\Big).
$$
- 左对数导数（左生成元）定义为
$$
H(u)=i\,\partial_u R(u)\,R(u)^{-1}.
$$
- 标量因子分离：若 $R(u)=\rho(u)\tilde R(u)$，则
$$
\partial_u R\,R^{-1}=\frac{\rho'}{\rho}I+\partial_u\tilde R\,\tilde R^{-1},
$$
其中第一项为全局相位，可在物理判据中忽略。

注意：若存在 $u$ 使 $\det R(u)=0$，左对数导数会发散，需用位移/插值或重参数化路径避开奇点。

## 3 极分解与单位化处理

为保证生成元为厄米并去除非物理放大，通常做极分解
$$
R=UP\quad(U\ \text{unitary},\;P>0).
$$
取单元部分 $U$ 的生成元
$$
H_U(u)=i\,\partial_u U(u)\,U(u)^\dagger,
$$
该 $H_U$ 为厄米，适合用于将 $R$ 解释为物理演化算子的情况。

在数值上可用 `scipy.linalg.polar` 或 SVD 实现极分解。

## 4 瞬时谱、本征基与 Berry 连接

对瞬时生成元做本征分解：
$$
H(u)\,|n(u)\rangle = E_n(u)\,|n(u)\rangle,
\qquad W(u)=[\,|n(u)\rangle\,].
$$
令 $D(u)=\mathrm{diag}(E_n(u))$，则
$$
W^\dagger H W = D,\qquad A(u)=iW^\dagger \partial_u W.
$$
定义 Berry（非阿贝尔）连接 $A(u)$。

把态写作 $\psi=W\tilde\psi$，代回薛定谔方程得
$$
i\partial_u\tilde\psi=(D-A)\tilde\psi,
$$
因此演化算子为
$$
U(u)=W(u)\,\mathcal T_u\exp\Big(-i\int_{u_0}^u (D(s)-A(s))\,ds\Big)\,W(u_0)^\dagger.
$$

## 5 几何相与动力相的分离条件

严格地，若 $D(s)$ 与 $A(s')$ 不对易则不能完全分离。实用判据：

- 若 $[D(s),A(s')]=0$ 对任意 $s,s'$ 成立，则演化可分解为几何因子乘以动力因子；
- 在绝热极限下（对所有 $m\ne n$，$|A_{mn}|\ll |E_m-E_n|$），非对角耦合被抑制，可近似逐本征态分离，得到每态的 Berry 相位乘以动力学相位；
- 在退化子空间需使用 Wilczek‑Zee 非阿贝尔连接，演化仍为路径序列形式。

因此不能把 $H$ 直接等同为 Berry 连接；$A$ 是在本征基下出现的几何量，$H$ 包含动力（$D$）与基变换引入的成分。

## 6 数值稳定化与实现配方

建议的数值步骤：

1. 对小步长 $\delta u$ 计算有限差分：
$$
H(u)\approx\frac{i}{\delta u}\log\big(R(u+\delta u)R(u)^{-1}\big).
$$
- 用 `scipy.linalg.logm`，并做 $H\leftarrow\tfrac12(H+H^\dagger)$ 以消除数值引入的非厄米成分；
2. 或先做极分解 $R=UP$，对 $U$ 做差分与 `logm`，得到厄米生成元 $H_U$；
3. 若 $\det R$ 很小或出现奇异，尝试改变取点/插值，并把正则化参数记录下来以便追踪；
4. 针对一系列 $u$ 值，用 `eigh` 获得 $W(u),D(u)$，在相邻 $u$ 上进行列向量匹配和相位固定（phase continuity），再差分计算 $\partial_u W$，得到 $A=iW^\dagger\partial_u W$；
5. 计算 $\|[D,A]\|$ 与最小能隙以评估是否可以采用几何/动力分离近似。

## 7 从 H 到 BdG / Majorana 与 MZM 判据

- 把 $H$ 在 Pauli 张量基上展开得到 $h_{\alpha\beta}$：
$$
H=\sum_{\alpha,\beta\in\{0,x,y,z\}} h_{\alpha\beta}\,\sigma^\alpha\otimes\sigma^\beta,
\qquad h_{\alpha\beta}=\tfrac14\mathrm{Tr}[(\sigma^\alpha\otimes\sigma^\beta)H].
$$
- 用 Jordan–Wigner 展开并收集二次费米子项，若确认为二次项则构造 BdG，得到单粒子 Hamiltonian `h` 与配对 `\Delta`；
- 把 BdG 转为 Majorana 反对称矩阵 $A$（$H=(i/4)\\Gamma^T A\\Gamma$），计算 Pfaffian、零模与对局域扰动的敏感度。

判据（综合）：

- 零能：存在受保护的 0‑mode（$E\approx0$）且与其他能级保持间隙；
- 局域性：零模在系统两端局域化（而非局域于同一区域）；
- 尺度标度：能量分裂随系统尺寸 $L$ 指数小（MZM）或非指数（ABS）；
- 拓扑不变量：Pfaffian 符号改变或其它 Z2 指标支持拓扑相。

单凭 $H$ 写成 $i\gamma_i\gamma_j$ 的线性组合并不能区分 MZM 与 ABS；需结合以上多项检验。


