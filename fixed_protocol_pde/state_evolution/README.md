# 几何状态与保真度演化图

图像由已验证的完整状态方程 \(\partial_t\mathcal B=[2A,\mathcal B]\) 直接积分得到。采用既定余弦门控、两轮六段协议，\(t_c=E_0=0.3\,\mathrm{meV}\)。积分设置沿用 DOP853、rtol \(2\times10^{-11}\)、atol \(2\times10^{-13}\)、最大局部归一化步长 0.025。图中的输出采样间隔不等于积分步长。

## 对应论文 Fig. 1(d) 的终点参数图

横轴为 \(\tau/(100/\mathrm{meV})\)，范围 0–12；纵轴为 \(\lg(|t_1|/E_1)\)，范围 −1–1。固定 \(E_1=0.01\,\mathrm{meV}\)，在 241 个时长、101 个耦合比上求终点 \(t=6\tau\) 的状态。

![终点定向与保真度参数图](endpoint_maps.png)

从左到右分别是：

1. 定向分量 \(d=\mathcal B_{12}(6\tau)\)，色标范围为 −1 到 1。
2. 我们定义的平方保真度 \(F(6\tau)\)，色标范围为 0 到 1。
3. 论文所用的重叠模长 \(f=\sqrt F\)，色标与论文的量相同。

论文的十档调色板直接从原 PDF 的矢量填色中提取，\(F\) 和 \(\sqrt F\) 使用相同的 0–1 刻度。\(\mathcal B\) 的分量具有正负号，因此色标必须明确改为 −1–1；沿用该调色板的颜色顺序，并让零值附近落入同一色带，避免极小数值残差造成伪条纹。

\(d\) 只展示完整状态的一个分量。完整的十个终点分量另见 [分量图集](B_endpoint_components.png)（[PDF](B_endpoint_components.pdf)）。定向 \(d=-1\) 对应目标态，但一般情况下不能仅由 \(d\) 确定保真度，还需其他分量：

\[
F=\frac{d^2-d}{2}+\frac14\sum_{j=1}^3(x_j^2+y_j^2).
\]

终点图的横轴扫描不同的可控设计时长，每个像素来自一次完整的六段演化；它不是单次协议中沿物理时间前进的轨迹。有限绘图网格用于显示主要条带，极细的快速结构可能需要更密的输出采样。

## 固定段时长的完整时间演化

下面取可控时长 \(\tau=100/\mathrm{meV}\)，横轴采用论文 Fig. 1(b) 的 step 形式，即 \(t/\tau\in[0,6]\)。每段之间的虚线对应 \(t=j\tau\)，第二轮继续继承第一轮终态。

![完整状态、平方保真度和重叠模长](time_evolution.png)

四列依次为理想参考 \(E_1=|t_1|=0\)，以及固定 \(E_1=0.01\) 时的 \(|t_1|/E_1=0.1,1,10\)。理想列对应论文 Fig. 1(b) 的物理设置；其余三列展示不可控耦合的影响。

- 上排展示 \(\mathcal B(t)\) 的全部十个反对称分量，顺序为 \(d,x_1,x_2,x_3,y_1,y_2,y_3,h_1,h_2,h_3\)。
- 中排为平方重叠。蓝线是目标态保真度 \(F=|\langle\psi_1^+(0)|\psi_1^-(t)\rangle|^2\)；品红线是初态存留概率 \(F_{\rm initial}=|\langle\psi_1^-(0)|\psi_1^-(t)\rangle|^2\)。
- 下排取上述两条曲线的平方根，采用论文 Fig. 1(b) 的重叠模长定义、蓝色与品红色图例。

初态存留概率也直接由几何状态得到：\(F_{\rm initial}=F+d\)。两条概率的和一般不等于 1，因为态还可以投影到其他模式。理想演化中间的模长平台约为 0.5，因此平方保真度的平台约为 0.25。

单独导出：[B 的演化图](B_time_evolution.png)（[PDF](B_time_evolution.pdf)）；[F 的演化图](fidelity_evolution.png)（[PDF](fidelity_evolution.pdf)）；[组合图 PDF](time_evolution.pdf)；[终点参数图 PDF](endpoint_maps.pdf)。

## 计算与核对

保存的数组包含所有分量及保真度：[time_history.npz](time_history.npz)、[endpoint_map.npz](endpoint_map.npz)。输出样本检查 \(\mathcal B^3=-\mathcal B\)、分量范数、保真度范围及几何恒等式，并以独立波函数积分核对选定终点。每轮第三段的 \(F\) 保持不变也单独检查。实际误差记录在 [verification.json](verification.json)。

此次输出的约束 \(\mathcal B^3+\mathcal B\) 最大分量残差为 \(5.8\times10^{-11}\)；选定样本与独立波函数的 \(F\) 差异不超过 \(4.0\times10^{-14}\)。原始数组中的极小负 \(F\) 来自浮点误差（最小约 \(-2.94\times10^{-14}\)）；图中绘制平方根前取非负值，原始 \(F\) 数组保留用于核对。

脚本：[plot_state_evolution.py](../plot_state_evolution.py)。

```sh
.venv/bin/python -B fixed_protocol_pde/plot_state_evolution.py
# 只重新绘制已保存的数组：
.venv/bin/python -B fixed_protocol_pde/plot_state_evolution.py --replot
```
