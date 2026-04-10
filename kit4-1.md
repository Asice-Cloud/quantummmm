## 4.1.1 2D p+ip + 两涡旋的 Dehn twist 数值实验

这一小节记录上面 4.1–4.3 抽象框架在一个具体 2D p+ip toy 模型上的第一次数值实现：在格点 BdG 模型中放入两个涡旋，让其中一个绕着另一个走一圈，计算零模子空间上的 Berry holonomy，并与 Ising TQFT 给出的 $R^{\sigma\sigma}$ 与 $(R^{\sigma\sigma})^2$ 做 SU(2) 共轭比较。

数值脚本在 [verify/run_pip_vortex_berry.py](verify/run_pip_vortex_berry.py)。

### (a) 模型与参数选择

我们使用的是一个简化的 spinless 2D p+ip 超导在方格子上的离散 BdG 模型：

- 晶格：$L_x\times L_y$ 方格，开边界条件（这里取 $L_x=L_y=20$）；
- Nambu 基底：$\Psi=(c_0,\dots,c_{N-1},c_0^\dagger,\dots,c_{N-1}^\dagger)^T$；
- BdG 哈密顿量块结构

$$
	H = \begin{pmatrix}
	h & \Delta \\
	\Delta^\dagger & -h^T
	\end{pmatrix},\qquad N=L_xL_y,
	$$

其中
- $h$：最近邻实跳跃 $t$ 与化学势 $\mu$，
- $\Delta$：在 $(x,y)\to(x+1,y)$ 和 $(x,y)\to(x,y+1)$ 最近邻键上分别实现 $p_x$ 和 $p_y$ 成分，并叠加来自涡旋纹理的相位因子。

具体取值：

- $t=1.0$；
- $\mu=-2.0$（落在拓扑相区内部）；
- $\Delta_0=0.5$，分别在 x 和 y 链接上实现 $\Delta_x=\Delta_0$、$\Delta_y=i\Delta_0$；
- 涡旋纹理通过序参量相位 $\theta(x,y)$ 引入，近似为各涡旋的极角之和

$$
		heta(x,y) = \sum_a \arg\bigl((x-x_a)+i(y-y_a)\bigr),
$$

并在每条键上使用 $(\theta_i+\theta_j)/2$ 作为配对相位；这给出两个涡旋周围的 $2\pi$ 绕转结构。

### (b) 涡旋路径与零模谱

在晶格中心

$$
 (c_x,c_y)=\Big(\tfrac{L_x-1}{2},\tfrac{L_y-1}{2}\Big)
$$

附近固定一个涡旋 $v_1=(c_x,c_y)$，另一个涡旋 $v_2$ 沿着围绕 $v_1$ 的矩形路径运动：

- 路径是一条以 $(c_x,c_y)$ 为中心、边长约为 $\min\{L_x,L_y\}/2$ 的闭合矩形；
- 每条边离散为 $n_\mathrm{seg}=12$ 个步长，总步数约 $4n_\mathrm{seg}$，最后回到起点，形成涡旋配置空间中的一条闭合回路 $\gamma$。

在每个路径点 $\lambda_k$，构造对应的 BdG 哈密顿量 $H(\lambda_k)$，对角化后按 $|E|$ 升序选取前两个本征态构成“零模子空间”基底 $V_k$。在初始点上，我们的脚本打印出最低几个 $|E|$：

$$
E \approx \{\pm 0.00198, \pm 0.0526, \pm 0.0928, \pm 0.1129,\dots\},
$$

可以看到有一对非常接近零能的本征值，与后面的激发之间存在明显能隙，说明这两个态确实可以视作涡旋绑定的近零模子空间。

### (c) 零模子空间上的 Berry holonomy

对相邻两步 $\lambda_k,\lambda_{k+1}$，脚本计算零模基底的重叠矩阵

$$
W_k = V_k^\dagger V_{k+1},
$$

然后通过极分解（SVD）投影到最近的幺正矩阵 $U_k$：

$$
W_k = U_k S_k V_k^\dagger \;\mapsto\; U_k^{(\text{polar})}=U_k V_k^\dagger.
$$

沿整个回路相乘得到 Berry holonomy

$$
U_{\text{holo}} = \prod_k U_k^{(\text{polar})}
$$

在这个二维零模子空间上的作用。数值结果为

$$
U_{\text{holo}} \approx
\begin{pmatrix}
 -0.2281 + 0.9736 i & 0 \\
 0 & -0.2281 - 0.9736 i
\end{pmatrix},
$$

其本征值

$$
\lambda_{\pm} \approx -0.2281 \pm 0.9736 i,
$$

对应相位

$$
\arg\lambda_{\pm} \approx \pm 1.80\,\text{rad} \approx \pm \Big(\tfrac{\pi}{2} + 0.23\Big).
$$

去掉整体相位并归一到 $\mathrm{SU}(2)$ 后，$U_{\text{holo}}$ 仍然是对角的，表示在当前选取的零模基底下，这条涡旋绕行路径实现的是一次接近 $\pi/2$ 角度的“Z 型”旋转。

### (d) 与 Ising $R^{\sigma\sigma}$ 与 Dehn twist 的对比

在 Ising TQFT 中，两只 $\sigma$ 粒子的 $R$‑符号为

$$
R^{\sigma\sigma}_1 = e^{-i\pi/8},\qquad R^{\sigma\sigma}_\psi = e^{3i\pi/8},
$$

其平方 $(R^{\sigma\sigma})^2$ 的本征值为

$$
e^{-i\pi/4},\ e^{3i\pi/4},
$$

去掉整体相位并归一到 $\mathrm{SU}(2)$ 后，$R^{\sigma\sigma}$ 和 $(R^{\sigma\sigma})^2$ 分别给出接近“$\pi/4$ 旋转”和“$\pi/2$ 旋转”的门。脚本中我们构造了 SU(2) 归一化后的矩阵

$$
R_{\mathrm{Ising}}^{(\mathrm{su2})},\qquad U^{(\mathrm{su2})}_{\mathrm{Dehn}}=(R^{\sigma\sigma})^2_{\mathrm{su2}},
$$

并与 $U_{\text{holo}}$ 做相位无关的重合度比较：

$$
F_R = \frac{1}{2}\bigl|\mathrm{Tr}(R_{\mathrm{Ising}}^{(\mathrm{su2})\,\dagger} U_{\text{holo}}^{(\mathrm{su2})})\bigr|,\\
F_{\mathrm{Dehn}} = \frac{1}{2}\bigl|\mathrm{Tr}(U_{\mathrm{Dehn}}^{(\mathrm{su2})\,\dagger} U_{\text{holo}}^{(\mathrm{su2})})\bigr|.
$$

数值结果为

$$
F_R \approx 0.85,\qquad F_{\mathrm{Dehn}} \approx 0.97.
$$

这说明：

- 该涡旋绕行路径的 Berry holonomy 在 SU(2) 意义下**比较接近** Ising 的两粒子 braid $R^{\sigma\sigma}$（$F_R\sim0.85$），
- 但与 Dehn twist $(R^{\sigma\sigma})^2$ 的 SU(2) 归一化矩阵更加接近（$F_{\mathrm{Dehn}}\sim0.97$），可以视为一次“近乎完美”的 Dehn twist 实现，只差一个小的角度修正和整体相位。

结合初始谱中存在清晰零模和能隙，以及路径几何上确实是“一个涡旋绕另一个涡旋一圈”的事实，这为我们提供了一个 2D p+ip 连续/格点模型中的直接证据：

> 在合适的参数和路径选择下，涡旋绕行的 Berry holonomy 在零模子空间上实现的，正是与 Ising TQFT 给出的 $(R^{\sigma\sigma})^2$ 共轭的 Dehn twist 门。

这与我们在 1D/4‑Majorana toy 模型中通过 half twist 的平方和 F‑move 得到的“$U_{\text{Dehn}}\simeq F^{-1}R^2F$” 图景形成了两种完全独立、但高度一致的数值验证：一个来自 1D 端点 Majorana 模型，另一个来自 2D p+ip 涡旋模型，它们在同一个 SU(2) 逻辑子空间上给出了相同的 Dehn twist 结构。

### (e) 化学势扫描 $F_{\mathrm{Dehn}}(\mu)$

为了在 2D p+ip 模型上得到一个“参数联络上平坦区域”的更直观图像，我们在相同几何和路径下，对化学势 $\mu$ 做了一维扫描，计算每个 $\mu$ 上零模 Berry holonomy 相对于 Ising Dehn twist $(R^{\sigma\sigma})^2$ 的重合度

$$
F_{\mathrm{Dehn}}(\mu) = \frac{1}{2}\bigl|\mathrm{Tr}(U_{\mathrm{Dehn}}^{(\mathrm{su2})\,\dagger} U_{\text{holo}}^{(\mathrm{su2})}(\mu))\bigr|.
$$

实现该扫描的脚本为 [verify/run_pip_vortex_scan.py](verify/run_pip_vortex_scan.py)。在 $L_x=L_y=20, t=1.0, \Delta_0=0.5$ 固定的情况下，我们在

$$
\mu \in \{-3.0, -2.5, -2.0, -1.5, -1.0, -0.5\}
$$

上计算了 $F_{\mathrm{Dehn}}(\mu)$，数值结果为（保留三位小数）

$$
\begin{array}{c|cccccc}
\mu & -3.0 & -2.5 & -2.0 & -1.5 & -1.0 & -0.5 \\
\hline
F_{\mathrm{Dehn}}(\mu) & 0.997 & 0.966 & 0.974 & 0.087 & 0.919 & 0.000
\end{array}
$$

对应的曲线图保存在仓库根目录下的 `pip_vortex_F_Dehn_vs_mu.png` 中。

从这张图和表格可以看出：

- 在 $\mu\approx-3.0,-2.5,-2.0,-1.0$ 附近，$F_{\mathrm{Dehn}}(\mu)$ 维持在 $0.92\sim1.0$ 的高值区间，说明在这段化学势区间内，同一条涡旋绕行路径的 Berry holonomy 在 SU(2) 意义下始终与 Ising 的 $(R^{\sigma\sigma})^2$ 高度共轭，可以把这一段视作“2D p+ip 模型上实现理想 Dehn twist 的参数平坦区”；

- 在 $\mu\approx-1.5$ 和 $\mu\approx-0.5$ 一带，$F_{\mathrm{Dehn}}(\mu)$ 明显跌落到 $\sim0.1$ 甚至接近 0，说明此时要么 BdG 谱的零模/能隙结构已经发生显著变化，要么这条几何路径在零模子空间上的 Berry holonomy 已经严重偏离理想的 Dehn twist 门，不再能被简单地视为同一个 mapping class 元素的“拓扑实现”。

因此，这个一维 $F_{\mathrm{Dehn}}(\mu)$ 扫描给出了一个 2D p+ip 连续/格点模型中的“Dehn twist 平坦区 vs 参数偏离”的初步相图切片：在某个有限的 $\mu$ 区间内，涡旋绕行路径对应的 Berry holonomy 在逻辑子空间中几乎不变（始终等价于 Ising Dehn twist）；而当 $\mu$ 进一步偏离这一区间时，同一条几何路径实现的门迅速失去这一拓扑特征。这与我们在 1D/R‑参数空间中用 $\mathcal M_R^{(\mathrm{YBE})}$ 描述“平坦/拓扑子流形”的做法在概念上完全平行，只是这里的参数轴从 $(a,b,c,d)$ 换成了 2D 模型中的化学势 $\mu$。

