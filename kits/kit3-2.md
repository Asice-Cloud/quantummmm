### R+Z2 与 Ising TQFT / 上层空间的解析对应

本节把在讨论中反复出现的那句话用统一公式写清楚：

> 局域 R‑门和 Z2 翻转在格点上的操作，
> 在上层空间（配置空间 × spin 结构）的 TQFT 描述中，
> 等价于沿某条曲线的编织 / Dehn twist / 改变 spin 结构，
> 在纤维（零模/任意子 Hilbert 空间）上的作用就是 Ising 范畴给出的 F、R、S、T 矩阵。

这里以已识别为 Ising 相的参数区为前提，假定低能极限由 Ising MTC 唯一确定。

---

#### 1. 格点上的路径算符：由 R+Z2 构造

在 Majorana+Z2 描述中，每条边 $e=(ij)$ 有：

- Majorana 双线性
	$$K_e := \gamma_i\gamma_j,$$
- Z2 链变量
	$$u_e\in\{\pm1\}.$$

由 R_to_Kitaev 给出的局域 gate（忽略整体 U(1) 相位）取为
$$
U_e(u_e) := \exp\Big(\frac\pi4\,u_e\,K_e\Big)
				 = \exp\Big(\frac\pi4\,u_e\,\gamma_i\gamma_j\Big).
$$

对任意一条格上路径/闭曲线 $\gamma$，写成边序列 $e_1,\dots,e_m$，定义沿 $\gamma$ 的有序乘积算符
$$
U(\gamma,u) := \overrightarrow{\prod_{k=1}^m} U_{e_k}(u_{e_k})
						 = \overrightarrow{\prod_{k=1}^m}\exp\Big(\frac\pi4\,u_{e_k}\,\gamma_{i_k}\gamma_{j_k}\Big).
$$

若沿 $\gamma$ 还需要整体翻转 Z2 Wilson 线以改变 spin 结构，则再定义
$$
\Omega_\gamma(u) := \prod_{e\in\gamma} u_e \in\{\pm1\}.
$$

在整个 Hilbert 空间上的“R+Z2 操作”就是
$$
\mathcal U_\gamma(u) = U(\gamma,u)\quad\text{或}\quad\mathcal U_\gamma^{\text{(spin)}}(u) = U(\gamma,u)\,\Omega_\gamma(u),
$$
分别对应“纯编织/平移”和“编织+改变 spin 结构”的两类操作。

---

#### 2. 投影到拓扑纤维：零模/任意子子空间上的表示

固定一个曲面 $\Sigma$ 和若干涡旋/缺陷位置 $\{x_a\}$，记低能拓扑子空间（Majorana 零模/任意子空间）为
$$
\mathcal H_{\mathrm{top}}(\Sigma,\{x_a\}),
$$
其在全 Hilbert 空间中的投影算符为 $P_{\mathrm{top}}$。

定义格点构造出的路径表示：
$$
\rho_{R+Z_2}([\gamma])
 := P_{\mathrm{top}}\,U(\gamma,u)\,P_{\mathrm{top}},
$$
或在需要改变 Wilson 线/自旋结构时：
$$
\rho_{R+Z_2}^{\text{(spin)}}([\gamma])
 := P_{\mathrm{top}}\,U(\gamma,u)\,\Omega_\gamma(u)\,P_{\mathrm{top}}.
$$

这里 $[\gamma]$ 是路径/闭曲线的同伦类（在给定涡旋/端点配置和固定的 Z2 背景下）：

- 若 $\gamma$ 只在涡旋世界线之间绕行/交叉，则 $[\gamma]$ 可视为 braid 群 $B_n$ 或其曲面推广的一个元素；
- 若 $\gamma$ 绕过非平凡一维环路（如 torus 的 a‑cycle）并可能配合 Wilson 线翻转，则 $[\gamma]$ 对应 mapping class 群 $\mathrm{MCG}(\Sigma)$ 的某个元素（典型为 Dehn twist）。

到这一步为止，“所有局域 R‑门和 Z2 翻转”已经通过 $\rho_{R+Z_2}([\gamma])$ 被压缩成“一个路径类 $[\gamma]$ 在拓扑子空间上的幺正作用”。

---

#### 3. Ising TQFT/MTC 提供的上层联络与 mapping class 表示

在 Ising 相区，长波极限由 Ising 模张量范畴（MTC）唯一决定。它给出：

- 对每个曲面 $\Sigma$ 与插入点/涡旋配置 $\{x_a\}$，一个 Hilbert 空间
	$$
	\mathcal H_{\mathrm{Ising}}(\Sigma,\{x_a\}),
	$$
	等价于 CFT 中的 conformal blocks 或任意子融合空间；
- 一个 mapping class 群的单值表示
	$$
	\rho_{\mathrm{Ising}}:\;\pi_1\big(\mathrm{Conf}_n(\Sigma)\times\text{spin 结构}\big)\to U\big(\mathcal H_{\mathrm{Ising}}(\Sigma,\{x_a\})\big),
	$$
	其局部生成元由 F、R 矩阵给出（多 $\sigma$ 粒子的 braid），全局生成元由 S、T 矩阵给出（torus 上的 Dehn twist 与模变换）。

例如：

- 四个 $\sigma$ 在平面上时，$B_1,B_2,B_3\in B_4$ 的表示由
	$$
	B_1^{\text{(Ising)}}\sim R^{\sigma\sigma},\quad B_3^{\text{(Ising)}}\sim R^{\sigma\sigma},\quad B_2^{\text{(Ising)}}=F^{-1}R^{\sigma\sigma}F
	$$
	给出；
- 在 torus 上，沿 a‑cycle 的 Dehn twist $T_a$ 在粒子基 $\{1,\psi,\sigma\}$ 上的作用为
	$$
	\rho_{\mathrm{Ising}}(T_a)=T=\mathrm{diag}(\theta_1,\theta_\psi,\theta_\sigma).
	$$

这些都是由 Ising 范畴的 F、R、S、T 数据解析给出的“上层联络”的矩阵表达。

---

#### 4. 精确对应关系：存在统一的基变换 $W$

关键断言是：在 Ising 相区的每个 $(\Sigma,\{x_a\})$ 上，存在一个幺正同构
$$
W:\;\mathcal H_{\mathrm{top}}(\Sigma,\{x_a\})\;\xrightarrow{\;\simeq\;}
\mathcal H_{\mathrm{Ising}}(\Sigma,\{x_a\}),
$$
使得对所有你关心的路径/映射类 $[\gamma]$，有
$$
W\,\rho_{R+Z_2}([\gamma])\,W^{-1} = \rho_{\mathrm{Ising}}([\gamma]),
$$
或在涉及 Wilson 线/自旋结构改变时：
$$
W\,\rho_{R+Z_2}^{\text{(spin)}}([\gamma])\,W^{-1} = \rho_{\mathrm{Ising}}^{\text{(spin)}}([\gamma]).
$$

在已显式计算的例子中：

- **两个涡旋（两 $\sigma$）**：
	- $\mathcal H_{\mathrm{top}}\cong\mathrm{span}\{|0\rangle,|1\rangle\}$（偶/奇费米数），
	- $\mathcal H_{\mathrm{Ising}}\cong\mathrm{span}\{|1\rangle,|\psi\rangle\}$，
	- $\rho_{R+Z_2}(b)$ 给出 $\mathrm{diag}(1,i)$，
	- $\rho_{\mathrm{Ising}}(b)=R^{\sigma\sigma}=\mathrm{diag}(e^{-i\pi/8},e^{3i\pi/8})$，
	- 通过吸收整体相位并重定义基矢，可找到 $W$ 使两者一致。

- **四个涡旋（四 $\sigma$）**：
	- $\mathcal H_{\mathrm{top}}$ 的偶费米数子空间由 $\{|00\rangle,|11\rangle\}$ 张成；
	- $\mathcal H_{\mathrm{Ising}}$ 的四 $\sigma$ 融合空间在固定总拓扑荷的子空间内由中间通道 $\{1,\psi\}$ 张成；
	- Majorana 模型中的
		$$
		B_1,B_3\sim\begin{pmatrix}1&0\\0&i\end{pmatrix},\quad
		B_2\sim\frac{1}{\sqrt2}\begin{pmatrix}1&-i\\-i&1\end{pmatrix},
		$$
		与 Ising 范畴中的
		$$
		B_1^{\text{(Ising)}},B_3^{\text{(Ising)}}\sim R^{\sigma\sigma},\quad
		B_2^{\text{(Ising)}}=F^{-1}R^{\sigma\sigma}F\sim\frac{1}{\sqrt2}\begin{pmatrix}1&-i\\-i&1\end{pmatrix}
		$$
		完全一致（只差整体相位），即存在幺正 $W$ 把两组矩阵同时对角/对齐。

概括地说：

- R_to_Kitaev + Z2 模型在 Ising 相区给出了一个**具体的、格点级的平坦联络**，其 holonomy 由 $\rho_{R+Z_2}([\gamma])$ 编码；
- Ising TQFT/MTC 给出了同一配置空间 × spin 结构上的一个**抽象的、范畴论的平坦联络**，其 holonomy 由 F、R、S、T 矩阵构成的 $\rho_{\mathrm{Ising}}([\gamma])$ 编码；
- 两者通过一个统一的纤维基变换 $W$ 精确对应：对任意 braid、Dehn twist 或带 Wilson 线的操作 $[\gamma]$，Majorana+Z2 实现的算符与 Ising 范畴的上层空间操作是一致的，只是写在不同的基底上。

这就是“局域 R‑门和 Z2 翻转，从上层看是沿曲线的编织/Dehn twist/改变 spin 结构，在纤维上的作用就是 F、R、S、T 矩阵”的解析公式版表述。

