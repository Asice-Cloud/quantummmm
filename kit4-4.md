## 4.4 Honeycomb 各向异性扫描与 Dehn twist 平坦区的搜索

本节在 4.3 节已构建的 honeycomb vison-loop Berry 实验基础上，进一步沿不同的耦合各向异性方向扫描参数，试图在蜂窝格 Kitaev 模型的参数空间中寻找与 2D p+ip 模型相似的 “Dehn twist 平坦区”。

具体而言，我们在一个 $4\times 4$ 的 brick-wall 蜂窝格上，固定 Z$_2$ 规范路径为“一个 vison 绕行另一个 vison 的最小闭合回路”（见 4.3.8），在此路径上对低能两维零模子空间计算 Berry holonomy $U_\text{holo}$，并将其归一化到 SU(2) 后，与 Ising 理论中的 $R^{\sigma\sigma}$ 与 Dehn twist $F^{-1}R^2F$ 做 SU(2) fidelity 比较：
$$
F_R = \frac{1}{2}\bigl|\operatorname{Tr}(R_\text{Ising}^\dagger U_\text{SU(2)})\bigr|,\qquad
F_\text{Dehn} = \frac{1}{2}\bigl|\operatorname{Tr}(U_\text{Dehn}^\dagger U_\text{SU(2)})\bigr|.
$$

在这一固定的 vison-loop 上，我们考察了两条最简单的一维各向异性切片：

1. $J_x=J_y=1$，沿 $J_z\in\{0.5,0.8,1.0,1.2,1.5\}$ 扫描；
2. $J_x=J_z=1$，沿 $J_y\in\{0.5,0.8,1.0,1.2,1.5\}$ 扫描。

对每个参数点，我们同时记录：

- 初始规范构型下最低若干本征值的最小模长 $\min|E|$，作为能隙尺度的 proxy；
- vison-loop Berry holonomy 的 SU(2) 形式 $U_\text{SU(2)}$；
- 与 Ising $R^{\sigma\sigma}$ 与 Dehn twist 的 fidelity $F_R,F_\text{Dehn}$。

数值结果总结如下。

**(1) $J_z$ 各向异性切片：$J_x=J_y=1$**

在这条切片上，我们得到如图 honeycomb_vison_F_Dehn_vs_Jz.png 与 honeycomb_vison_gap_vs_Jz.png 所示的行为：

- $F_R(J_z)$ 在整个扫描区间内几乎保持在 $F_R\approx 0.707$ 附近，表明 vison-loop Berry holonomy 始终与 Ising $R^{\sigma\sigma}$ 之间夹角约为 $\pi/4$；
- $F_\text{Dehn}(J_z)$ 在数值精度内始终接近于 0，没有出现类似 2D p+ip 模型中那样的 “Dehn twist 平坦区”；
- 对应的 $\min|E|(J_z)$ 始终为一有限正数，说明在这条路径上 vison 配置处于一条有能隙的相内，Berry holonomy 的变化主要是 “几何性的” 而不是由于能隙闭合。

**(2) $J_y$ 各向异性切片：$J_x=J_z=1$**

对称地，在 $J_x=J_z=1$ 固定下，对 $J_y$ 做同样的扫描并得到 honeycomb_vison_F_Dehn_vs_Jy.png 与 honeycomb_vison_gap_vs_Jy.png：

- $F_R(J_y)$ 再次几乎固定在 $0.707$ 左右；
- $F_\text{Dehn}(J_y)$ 在整个区间内依然数值上接近 0；
- $\min|E|(J_y)$ 保持有限，未见明显的能隙闭合。

从这两条最简单的各向异性切片可以看出：在当前 toy honeycomb 设置（无外场、小系统、固定 vison-loop 几何）下，“vison 绕行 vison” 所给出的 Berry holonomy 虽然在谱上是稳定且有能隙的，但在 SU(2) 意义下远离 Ising Dehn twist——对应的 $U_\text{SU(2)}$ 更接近一个“约为 $\pi/4$ 的绕轴旋转”，而非 Ising TQFT 中的 Dehn twist 门。这与 2D p+ip 模型中可以找到一整段 $\mu$ 区间使 $F_\text{Dehn}\approx 1$ 的情形形成鲜明对比。

这提示我们：在 honeycomb 上复制 p+ip 式的 “Dehn twist 平坦区”，至少还需要满足以下额外条件中的若干：

- 进入更真实的非阿贝尔区（例如加入有效磁场/三自旋项，使蜂窝格 Majorana 带结构获得非平庸 Chern 数）；
- 考察更大的系统尺寸与更远距离的 vison 对，从而弱化有限尺寸效应；
- 调整 vison-loop 的路径，使之更接近 mapping class group 中的 Dehn twist 生成元，而不仅仅是一个“局部绕行回路”。

在后续的小节中，可以结合文献中的 honeycomb 相图，在 $(J_x,J_y,J_z)$ 以及有效磁场参数空间中选取更具物理意义的路径，沿这些路径重复上述 vison-loop Berry 实验，进一步搜索可能的 honeycomb 版 Dehn twist 平坦区。

作为一个最简单的“相图导向型”尝试，我们还考察了耦合落在标准 simplex $J_x+J_y+J_z=1$ 上的一条 triangle 路径：
$$
J_z = t,\qquad J_x = J_y = \tfrac{1}{2}(1-t),\qquad t\in[0,1],
$$
这条路径从一个 $xy$ 主导的角 ($t\ll 1$) 走向 $z$ 主导的角 ($t\lesssim 1$)，在 $t\approx 0.5$ 附近满足 $J_z\approx J_x+J_y$，因此在无场 Kitaev 解中通常对应靠近 gapless 线的区域。沿这条路径，我们在每个 $t$ 点上重复同样的 vison-loop Berry 计算，并得到如 honeycomb_vison_F_Dehn_vs_triangle_t.png 与 honeycomb_vison_gap_vs_triangle_t.png 所示的结果：

- $F_R(t)$ 几乎在整个区间内锁定在 $\approx0.707$，与前两条切片的行为一致；
- $F_\text{Dehn}(t)$ 在数值精度内仍然接近 0，即使在 gap proxy $\min|E|(t)$ 最小的 $t\sim 0.1$–$0.2$ 区域附近也未出现任何向 Ising Dehn twist 靠拢的迹象；
- $\min|E|(t)$ 在小 $t$ 处非常小（提示接近 gapless 区），随后随 $t$ 单调增大，进入明显有隙的 $z$ 主导相。

因此，在当前 toy honeycomb 模型和固定 vison-loop 下，即便我们沿着更贴近标准相图的一条 simplex 路径穿过/靠近 gapless 区域，得到的 vison Berry holonomy 仍然“钉”在与 Ising Dehn twist 正交的 SU(2) 方向上，未能像 2D p+ip 那样出现一段 $F_\text{Dehn}\approx1$ 的平坦区。这一负例进一步强化了 4.3 节的观点：要在 honeycomb 上真正复现 Ising Dehn twist，需要同时满足更强的拓扑条件（例如引入非平庸 Chern 数的有效磁场项）和更严格的几何条件（更适合的映射类群路径与更大系统极限）。

### 4.4.x 与 2D p+ip 的比较与阶段性总结

把本节 honeycomb 各向异性/triangle 扫描与 2D p+ip 上的结果（见 [kit4-1.md](kit4-1.md)）并列起来，可以得到一幅较为清晰的“Dehn twist 平坦区 vs 参数空间”的阶段性图景：

1. 在 2D p+ip + 两涡旋的格点 BdG 模型中，沿化学势 $\mu$ 的一维扫描给出了典型的 Dehn twist 平坦区（图 [pip_vortex_F_Dehn_vs_mu.png](pip_vortex_F_Dehn_vs_mu.png)）：
	- 在一段有限的 $\mu$ 区间内，$F_\text{Dehn}(\mu)$ 稳定在 $\sim0.95$–$1.0$ 的高值，说明同一条涡旋绕行路径的 Berry holonomy 在逻辑子空间内几乎刚性地等价于 Ising Dehn twist；
	- 当 $\mu$ 脱离这一区间时，$F_\text{Dehn}(\mu)$ 会迅速跌落到 $\sim0.1$ 或更低，意味着同一条几何路径不再实现同一个 mapping class 群元素，对应“离开拓扑平坦子流形”的情形。

2. 在 honeycomb 的 toy vison-loop 模型中，我们在三条相对简单的切片上做了完全平行的扫描：
	- $J_x=J_y=1$，扫 $J_z$（图 [honeycomb_vison_F_Dehn_vs_Jz.png](honeycomb_vison_F_Dehn_vs_Jz.png)、[honeycomb_vison_gap_vs_Jz.png](honeycomb_vison_gap_vs_Jz.png)）；
	- $J_x=J_z=1$，扫 $J_y$（图 [honeycomb_vison_F_Dehn_vs_Jy.png](honeycomb_vison_F_Dehn_vs_Jy.png)、[honeycomb_vison_gap_vs_Jy.png](honeycomb_vison_gap_vs_Jy.png)）；
	- simplex 上的 triangle 路径 $J_x=J_y=(1-t)/2, J_z=t$（图 [honeycomb_vison_F_Dehn_vs_triangle_t.png](honeycomb_vison_F_Dehn_vs_triangle_t.png)、[honeycomb_vison_gap_vs_triangle_t.png](honeycomb_vison_gap_vs_triangle_t.png)）。
	在所有这些切片上，我们都看到：
	- gap proxy $\min|E|$ 要么保持有限、要么仅在小区间内略微减小，整体上没有系统性崩塌；
	- $F_R$ 几乎锁定在 $\approx0.707$；
	- $F_\text{Dehn}$ 在数值精度内始终接近 0，从未出现类似 2D p+ip 那样的“高平坦段”。

3. 综合 4.3 的统一七步流程，可以把上述差异理解为：
	- p+ip 平台已经满足了“合适的 BdG 拓扑结构 + 合适的缺陷路径”这两层条件，因此在某段参数区间内 Berry 联络在适当规范下近似平坦，涡旋绕行路径的 holonomy 与 Ising Dehn twist 高度吻合；
	- 当前的 honeycomb toy 模型虽然借助 Majorana+Z$_2$ 表示和 vison branch cut 复制了相似的几何流程，但由于尚未引入足够真实的非阿贝尔结构（例如有效磁场/三自旋项带来的非平庸 Chern 数）以及更大系统/更“拓扑”的路径，其 Berry holonomy 仍然停留在与 Ising Dehn twist 近乎正交的 SU(2) 方向上。

从“阶段性成果”的角度看：

- 在 1D R(a,b,c,d) 模型和 2D p+ip 平台上，我们已经通过显式的 Berry 计算和参数扫描，数值锁定了一个与 Ising Dehn twist 共轭的拓扑门，并在参数空间中找到了对应的“平坦区”；
- 在 honeycomb 平台上，我们已经建立起完整的 vison-loop Berry 数值流水线，证明了 Berry–TQFT 对比在蜂窝格上同样可行，但目前三条简单的参数切片均给出“有隙但非 Dehn twist”的负例，这反而清楚地划出了“仅有非阿贝尔相和能隙，还不足以保证实现 Ising Dehn twist”的边界。

下一步，honeycomb 方向的自然发展是：在保持当前数值框架不变的前提下，引入更贴近文献相图的非阿贝尔区（例如加入磁场项、选取已知的 B 相参数）和更大系统尺寸，沿穿越/贴近这些非阿贝尔区的路径进行新的 $F_\text{Dehn}(p)$ 扫描，继续寻找可能的 honeycomb 版 Dehn twist 平坦区，与 2D p+ip 的结果做真正“对等”的比较。