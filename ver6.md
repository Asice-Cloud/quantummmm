(# 验证总结 — ver6)

**目的**
- 汇总对论文补充材料中“路径门 vs 理想 braid”以及“ABS 可视为 MZM 耦合”这两个观点的数值验证，并给出若干可视化（LDOS、端点局域性、路径 d(u) 的示例图）。

**使用的脚本（仓库路径）**
- [tools/verify_path_braid.py](tools/verify_path_braid.py) — 构造 4×4 路径门，投影到逻辑子空间，计算相位对齐残差與 braid/YBE 残差。
- [tools/scan_path_braid.py](tools/scan_path_braid.py) — 在 (steps, T) 网格搜索残差。
- [tools/test_mzm_to_abs.py](tools/test_mzm_to_abs.py) — 向 H4 中加入 Majorana‑样耦合项 ε，检查投影后 |d(u)| 随 ε 的变化。
- [tools/edge_localization.py](tools/edge_localization.py) — 把投影的有效模型映射到单粒子 BdG 链 (t,Δ,μ)，做 L 扫描并计算 E0(L) 与端点权重。
- [tools/compute_topo_invariant.py](tools/compute_topo_invariant.py) — 从 H4 得到 (t,Δ,μ)，计算体能隙并用简单判据 |μ|<2|t| 作为拓扑指示器。
- [tools/compute_ldos.py](tools/compute_ldos.py) — 计算并画出 LDOS（能量×位置热图）和零能量空间分布。

**具体模型与公式**

- 原始 4×4 哈密顿（路径生成）：对每个局部参数点我们使用四维 Pauli 展开。常用的路径映射写作
	$$
	H_4(u) = t_cigl[-g_1g_2\,Z\otimes I - g_1g_3\,Y\otimes X + g_1g_4\,Y\otimes Y - g_2g_3\,X\otimes X - g_2g_4\,X\otimes Y - g_3g_4\,I\otimes Z\bigr].
	$$
	其中每段路径的局部参数按文中约定为
	$$
	g^{(1)}(u)=(u,0,1-u,1),\qquad g^{(2)}(u)=(1-u,u,0,1),\qquad g^{(3)}(u)=(0,1-u,u,1),
	$$
	$u\in[0,1]$。

- 投影到逻辑子空间 $\mathcal S=\{\lvert01\rangle,\lvert10\rangle\}$：令
	$$
	P=\lvert01\rangle\langle01\rvert+\lvert10\rangle\langle10\rvert,
	$$
	有效两能级哈密顿
	$$
	H_{\rm eff}(u)=P H_4(u) P = d_0(u) I + d_x(u)\sigma_x + d_y(u)\sigma_y + d_z(u)\sigma_z.
	$$
	其中系数可由迹计算得
	$$
	d_i(u)=\tfrac12\mathrm{Tr}\bigl(H_{\rm eff}(u)\,\sigma_i\bigr),\quad d_0(u)=\tfrac12\mathrm{Tr}\bigl(H_{\rm eff}(u)\bigr).
	$$

- Pauli 张量到 $\vec d$ 的具体字典（在本笔记中使用的约定）:
	$$
	\begin{aligned}
	d_x(u)&=h_{xx}(u)+h_{yy}(u),\\
	d_y(u)&=-h_{xy}(u)+h_{yx}(u),\\
	d_z(u)&=h_{zI}(u)-h_{Iz}(u),
	\end{aligned}
	$$
	其中 $h_{\alpha\beta}(u)$ 是在 $\sigma^\alpha\otimes\sigma^\beta$ 基下的投影系数（$H_4=\sum h_{\alpha\beta}\,\sigma^\alpha\otimes\sigma^\beta$）。

- 投影模型的谱与几何判据：
	$$
	E_\pm(u)=\pm\lvert\vec d(u)\rvert,
	$$
	若路径 $\vec d(u)$ 围绕原点且动力学相可被化为全局相，则几何部分可产生 braid‑like 幺正变换；若 $\min_u\lvert\vec d(u)\rvert>0$ 则存在稳定能级劈裂（ABS）。

- 映射到 Kitaev 单粒子参数（用于构造 BdG 链）：对于 Pauli 系数记作
	$c_{xx}=h_{xx},\;c_{yy}=h_{yy},\;c_{xy}=h_{xy},\;c_{yx}=h_{yx},\;c_{zI}=h_{zI},\;c_{Iz}=h_{Iz},\;c_{zz}=h_{zz}$，我们使用的映射为
	$$
	\begin{aligned}
	t&=c_{xx}+c_{yy}+i(c_{xy}-c_{yx}),\\
	\Delta&=c_{xx}-c_{yy}-i(c_{xy}+c_{yx}),\\
	\mu&=4c_{zz}-2(c_{zI}+c_{Iz}).
	\end{aligned}
	$$

- 构造有限链的 BdG 单粒子哈密顿（开放边界）：对长度 $L$，定义实空间块矩阵
	$$
	H_{\rm BdG}=\begin{pmatrix}h&\Delta\\ -\Delta^*&-h^T\end{pmatrix},
	$$
	其中在格点上我们采用离散化近似
	$$
	h_{ii}=-\mu/2,\quad h_{i,i+1}=-t/2,\quad (\Delta)_{i,i+1}=\Delta/2
	$$
 （以及相应的伴随项），然后对 $H_{\rm BdG}$ 求本征值与本征向量用于 E0(L) 与端点权重分析。

- 局域态密度（LDOS）定义（用于绘图）：引入能量解析宽化 $\eta$，则第 $i$ 个格点的 LDOS
	$$
	\rho_i(E)=\sum_n\bigl(|u_n(i)|^2+|v_n(i)|^2\bigr)\frac{\eta/\pi}{(E-E_n)^2+\eta^2},
	$$
	其中 $(u_n,v_n)$ 是 BdG 本征向量，对应本征能 $E_n$。

以上公式在仓库脚本中均有对应实现（见上面脚本列表），数值结果与图像已保存于 `results/` 下相应子目录。

**关键数值结果与图**
- 路径门残差（示例）: [results/path_braid/path_braid_steps80_T1_tc1_axisx_hc1.png](results/path_braid/path_braid_steps80_T1_tc1_axisx_hc1.png) — 相位对齐残差约 1.30，braid/YBE 残差也较大，说明该路径未直接数值收敛到论文的理想 braid。
- 加耦合后的最小 |d|: [results/mzm_abs/min_d_vs_epsilon.png](results/mzm_abs/min_d_vs_epsilon.png), 示例曲线 [results/mzm_abs/d_u_examples.png](results/mzm_abs/d_u_examples.png) — 显示加入 Majorana‑样耦合 ε 会把投影后的 |d(u)| 抬高/改变，产生稳定的能级劈裂（ABS）。
- 端点局域性与尺度行为（示例）: [results/edge_localization/edge_loc_u1.571_d0.015.png](results/edge_localization/edge_loc_u1.571_d0.015.png) — 对 u≈π/2（和 3π/2）做 L={40,80,160,320,640} 扫描，最低能量 E0(L) 呈指数衰减，拟合得到衰减长度 xi≈240（表明 ABS 在这些参数下可以由两端 MZM 的重叠解释）。汇总数据保存在 [results/edge_localization/edge_loc_summary.npz](results/edge_localization/edge_loc_summary.npz)。
- 拓扑指示器: [results/topo_invariant/topo_invariant.npz](results/topo_invariant/topo_invariant.npz) — 在我们测试的 (u,δ) 网格上，映射得到的 (t,Δ,μ) 满足简单 Kitaev 判据 |μ|<2|t|（标记为 topo=True），与 E0 指标相容。
- LDOS（示例）: [results/ldos/ldos_u0_d0.015_L160.png](results/ldos/ldos_u0_d0.015_L160.png) 与 [results/ldos/ldos_u1.571_d0.015_L160.png](results/ldos/ldos_u1.571_d0.015_L160.png) — 分别给出 u=0 与 u=π/2 的能量‑位置 LDOS 热图与 E=0 空间分布，直观展示端点谱权重。

**结论（简短）**
- 我们的有效模型（4×4 → 投影到两能级 → 映射到 Kitaev 单粒子 BdG）能够同时表达“几何型的 braid‑like 演化”与“由 Majorana‑Majorana 耦合产生的 ABS”。
- 在某些路径参数（如 u≈π/2,3π/2）下，BdG 扫描显示 E0(L) 随 L 指数小（xi≈240），最低态端点局域性随 L 保持 → 支持“ABS 由两端 MZM 重叠/耦合形成”的解释。
- 但也存在路径/参数使得投影后的 d(u) 不包围原点且残差较大，此时路径门不能被直接视为论文里的理想 braid；补充材料的视角把 ABS≈耦合 MZM 作为分析框架是合理的，但要在多体层面断言存在受保护的 MZM仍需尺度与鲁棒性检验（如扰动测试、Pfaffian/winding 指标、端点态对局域扰动的敏感性）。

**建议的下一步（可选、自动化）**
- 做更密的参数网格扫描并用局部优化（T 或路径形状）最小化 braid 残差；脚本已部分就绪。
- 做鲁棒性测试：在出现“数值零模”的点加入局部扰动并观察 E0 是否立即分裂。
- 为论文级结论：把 BdG L‑scan 扩到更大 L，拟合并报告 xi 的不确定度以及端点态对扰动的响应。

—— 结束 ——

（如果你同意，我可以把上述所有图整理成一页 PDF 报告，并把 `ver6.md` 中的图像引用替换为内嵌缩略图。）

