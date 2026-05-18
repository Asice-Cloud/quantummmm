# 完整结果总结

这份结果把当前仓库里的主线串成一条完整链路：

路径设计
→ eight‑vertex 的局域 `R(u)` / `H(u)`
→ 全链 BdG / Kitaev 链
→ 逻辑子空间里的 Bloch 旋转
→ 子空间与全链两套判据分别判断 MZM / ABS。

## 1. 路径设计
我们的控制路径本质上是一个随时间变化的四门控参数
$$
\mathbf g(t)=(g_1(t),g_2(t),g_3(t),g_4(t)).
$$

在仓库中的 eight‑vertex 示例里，路径常写成三段：
$$
\mathbf g^{(1)}(u)=(u,0,1-u,1),\qquad
\mathbf g^{(2)}(u)=(1-u,u,0,1),\qquad
\mathbf g^{(3)}(u)=(0,1-u,u,1).
$$

它们表示局域耦合在不同通道之间按顺序打开、关闭、再交换的过程。这个路径是后面构造 `R(u)` 和 `H(u)` 的起点。

## 2. eight‑vertex 的局域 `R(u)` 与 `H(u)`
对单个两站点局域块，我们写成 4×4 哈密顿量
$$
H_4(u,\delta)=\cos u\,X\otimes X+\tfrac{\sin u}{2}(Y\otimes X - X\otimes Y)+\tfrac{\delta}{2}(Z\otimes I - I\otimes Z).
$$

对应的局域演化门是
$$
R_4(u)=\mathcal T\exp\!\left(-i\int H_4(u,\delta)\,dt\right),
$$
如果把 `u` 当作路径参数，也可以把它理解成沿路径累积得到的局域 braid 门。

把 `H_4` 做 Pauli 展开后，非零系数是
$$
h_{xx}=\cos u,\quad
h_{xy}=-\tfrac12\sin u,\quad
h_{yx}=+\tfrac12\sin u,\quad
h_{zI}=+\tfrac12\delta,\quad
h_{Iz}=-\tfrac12\delta.
$$

## 3. 投影到逻辑子空间：Bloch 旋转
取逻辑子空间
$$
\mathcal S=\{|01\rangle,|10\rangle\}
$$
并定义投影
$$
P_L=|01\rangle\langle01|+|10\rangle\langle10|.
$$

投影后得到两能级有效哈密顿量
$$
H_{\rm eff}=P_L H_4 P_L=d_0 I+d_x\sigma_x+d_y\sigma_y+d_z\sigma_z.
$$

在我们的约定里
$$
d_x=h_{xx}+h_{yy},\qquad
d_y=-h_{xy}+h_{yx},\qquad
d_z=h_{zI}-h_{Iz},\qquad
d_0=h_{II}-h_{zz}.
$$

代入 eight‑vertex 后得到
$$
\vec d(u)=(\cos u,\sin u,\delta),\qquad d_0=0.
$$

于是逻辑子空间的瞬时本征值是
$$
E_{\rm eff}(u)=\pm|\vec d|=\pm\sqrt{1+\delta^2}.
$$

这一步描述的是**逻辑子空间里的 Bloch 旋转**：
$$
U_{\rm eff}(T)=\mathcal T\exp\!\left(-i\int_0^T H_{\rm eff}(t)\,dt\right)
=e^{-i\int d_0dt}\,
\mathcal T\exp\!\left(-i\int \vec d(t)\cdot\vec\sigma\,dt\right).
$$

- 当 `\delta=0` 时，`\vec d(u)` 在 `xy` 平面绕原点，几何相主导，呈现 MZM‑like 的 braid 轨迹。
- 当 `\delta\neq0` 时，轨迹被抬离平面，逻辑子空间会出现明显的分裂，更像 ABS‑like。

这里的判断只对**投影子空间**成立，不等同于全链的最终物理结论。

## 4. 推广到全链：两种方式
### 4.1 直接嵌入到多体空间
如果只做小链、门级别演化，可以把局域算符嵌入到全局空间：
$$
O^{(j)}_{\rm global}=I^{\otimes(j-1)}\otimes O_{\rm local}\otimes I^{\otimes(L-j-1)}.
$$

多体哈密顿量就是各局域项的和：
$$
H_{\rm total}=\sum_j O^{(j)}_{\rm global}.
$$

如果是时间序列的门操作，则全局演化算符是这些嵌入门的按序乘积。

### 4.2 我们代码里真正采用的全链路线
仓库的全链数值并不是把 4×4 局域块直接拼成一个 2^L 多体大矩阵，而是：

1. 先从局域 `H_4` 做 `pauli_expand`。
2. 再用 `map_to_kitaev_from_h` 提取 `(t,\Delta,\mu)`。
3. 然后用 `build_bdg_chain` 组装长度为 `L` 的单粒子 BdG 链（2L×2L）。

对应的映射是
$$
t=h_{xx}+h_{yy}+i(h_{xy}-h_{yx}),\qquad
\Delta=h_{xx}-h_{yy}-i(h_{xy}+h_{yx}),\qquad
\mu=4h_{zz}-2(h_{zI}+h_{Iz}).
$$

这条路是 `compute_ldos.py`、`edge_localization.py`、`compute_topo_invariant.py` 采用的主线。

## 5. 全链上的 MZM / ABS 判据
在全链 BdG 里，判断不看 `d`，而看谱、LDOS 和拓扑。

### 5.1 拓扑不变量与 bulk gap
`compute_topo_invariant.py` 先得到 `(t,\Delta,\mu)`，再检查 bulk gap，并使用简单 Kitaev 判据
$$
|\mu|<2|t|,
$$
配合 `gap>0` 作为拓扑相的筛选。

### 5.2 开链零模与链长缩放
`edge_localization.py` 对不同链长 `L` 计算
$$
E_0(L)=\min|E|,
$$
以及端点权重。若 `E_0(L)` 随 `L` 指数趋零且边权集中在端点，则支持 MZM；若有明显有限分裂或振荡，则更像 ABS。

### 5.3 LDOS 直接看局域零能峰
`compute_ldos.py` 直接对开链 BdG 对角化，画出能量—位置 LDOS 图，并提取 `E=0` 的空间剖面。端点零能峰与边缘局域是 MZM 的重要证据。

## 6. 统一结论
我们的结果可以概括成两层：

- **子空间层**：`\delta` 决定 `d_z`，因此决定逻辑子空间里是更接近理想 braid / MZM‑like，还是更像有分裂的 ABS‑like。
- **全链层**：MZM / ABS 的最终判定必须看 BdG 的 bulk gap、`E_0(L)` 和 LDOS，而不是只看 `\vec d`。

这也解释了为什么在 eight‑vertex 的均匀映射里，`\delta` 可以明显改变投影后的几何，但未必改变全链的拓扑参数：它们属于不同层次的判据。

## 7. 对应脚本
- `tools/compute_topo_invariant.py`：全链拓扑判据。
- `tools/edge_localization.py`：全链 `E_0(L)` 与边缘局域。
- `tools/compute_ldos.py`：全链 LDOS 与零能空间剖面。
- `tools/test_mzm_to_abs.py`：局域投影层的 `d` 诊断。
- `tools/generate_comparison_panels.py`：全链结果与低维代理的对照。

## 8. 最终一句话
路径设计给出局域 `R(u)` / `H(u)`；投影子空间用 `\vec d` 描述 Bloch 旋转并区分 MZM‑like / ABS‑like；推广到全链后，真正的 MZM / ABS 要由 BdG 的谱、LDOS 和拓扑不变量来判断。