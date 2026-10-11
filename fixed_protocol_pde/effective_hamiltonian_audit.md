# 从有效哈密顿量独立核对几何状态方程

本次直接以 [model.md](../model.md) 开头的有效哈密顿量为依据，同时构造文档中的四元数旋量矩阵，以及普通费米子 Fock 空间中的算符矩阵。两条路线均严格给出当前的 \(K=2A\)。五项系数和十个状态分量方程逐项一致；当前 PDE 的系数与段间传递没有在这些检查中暴露错误。

确实发现并修正了 `model.md` 的三个旋量矩阵表头：升序指标 \(24,15,34\) 对应的系数应分别为 \(-|t_2|,+|t_1|,+|t_3|\)。原来列出的五个 \(\Sigma\) 矩阵本身正确，当前 PDE 程序也一直使用正确的系数。这次修正文档不会改变已得到的演化曲线。

## 1. 从原始算符顺序确定全部符号

记 \(E=E_1\)、\(q=|t_1|g(t)\)、\(p=|t_2(t)|\)、\(r=|t_3(t)|\)、\(e_d=E_d(t)\)。按照原始乘积，

\[
H=i e_d\gamma_a\gamma_b+iE\gamma_1\gamma_2
+ip\gamma_a\gamma_2-iq\gamma_b\gamma_1-ir\gamma_a\gamma_3.
\]

固定顺序 \((1,2,3,a,b)\)，利用不同 Majorana 之间反对易，得到

\[
\boxed{H=i\left(E\gamma_1\gamma_2+q\gamma_1\gamma_b
-p\gamma_2\gamma_a+r\gamma_3\gamma_a+e_d\gamma_a\gamma_b\right).}
\]

因此，定义 \(L_{ij}=e_i e_j^T-e_j e_i^T\) 后，

\[
A=\begin{pmatrix}
0&E&0&0&q\\
-E&0&0&-p&0\\
0&0&0&r&0\\
0&p&-r&0&e_d\\
-q&0&0&-e_d&0
\end{pmatrix}.
\tag{1}
\]

在 \(\{\gamma_i,\gamma_j\}=2\delta_{ij}\) 的约定下，\(\hbar=1\)，

\[
i[H,\gamma_i]=2\sum_j A_{ij}\gamma_j.
\tag{2}
\]

这是确定因子 2 的算符等式。它不依赖门控细节、演化时长或参数是否小。

单独检验 \(iE\gamma_1\gamma_2\) 就能看清其物理含义。令 \(c=(\gamma_1+i\gamma_2)/2\)，则

\[
iE\gamma_1\gamma_2=E(2c^\dagger c-1).
\]

两种占据态的能量为 \(-E,+E\)，激发能差为 \(2E\)，Nambu 矩阵的能级为 \(\pm2E\)。旋量哈密顿量中的 \(E\) 与单粒子矩阵中的 \(2E\) 属于同一算符的不同表示。若要让后者的正能级直接等于 \(E\)，原双线性系数必须相应写成 \(E/2\)。

## 2. 使用文档给出的四元数矩阵

选择忠实的复矩阵表示

\[
\mathbf i=-i\sigma_x,\qquad \mathbf j=-i\sigma_y,
\qquad \mathbf k=-i\sigma_z.
\]

它们满足 \(\mathbf i\mathbf j=\mathbf k\)。将其代入文档的 \(2\times2\) 四元数块，得到 \(4\times4\) 复矩阵 \(\Gamma_i\)。程序精确验证了全部 Clifford 反对易关系，以及文档列出的五个

\[
\Sigma_{ij}=\frac14[\Gamma_i,\Gamma_j]=\frac12\Gamma_i\Gamma_j.
\]

五维索引中 \(4=a,5=b\)，所以由原始哈密顿量直接得到

\[
\boxed{H_{\rm spin}=2i\left(E\Sigma_{12}+q\Sigma_{15}
-p\Sigma_{24}+r\Sigma_{34}+e_d\Sigma_{45}\right).}
\tag{3}
\]

积分 \(\dot U_4=-iH_{\rm spin}U_4\)、\(U_4(0)=I_4\)，再通过共轭作用定义

\[
R_{ij}=\frac14\operatorname{tr}(U_4^\dagger\Gamma_iU_4\Gamma_j).
\]

由式 (2) 可得 \(\dot R=2AR\)。故论文所研究的五模波函数与几何状态满足

\[
\psi(t)=R(t)\psi_0,\qquad
\mathcal B(t)=R(t)\mathcal B_0R(t)^T,
\]

\[
\boxed{\partial_t\mathcal B=[2A,\mathcal B].}
\tag{4}
\]

这里通过旋量的共轭作用取得五维旋转，再计算所选单粒子态的重叠。直接选取一个四分量旋量态并计算其重叠，会涉及另一种保真度，不能与论文的波函数定义混用。

## 3. 从普通费米子独立构造 BdG 矩阵

为使用普通费米子 Fock 空间，保留论文中与其他模式不耦合的另一端 Majorana，记为 \(\gamma_s\)，组成

\[
c_1=\frac{\gamma_1+i\gamma_2}{2},\qquad
c_2=\frac{\gamma_3+i\gamma_s}{2},\qquad
d=\frac{\gamma_a+i\gamma_b}{2}.
\]

采用 Jordan–Wigner 构造三个费米子，算符矩阵维数为 8；先验证正则反对易关系。取论文允许的相位选择 \(\phi_1=0,\phi_2=-\pi\)，则 \(t_1=q,t_2=p,t_3=-ir\)。直接代入论文式 (1)：

\[
H_E=2e_d d^\dagger d+iE\gamma_1\gamma_2
+ip(d+d^\dagger)\gamma_2
+q(d^\dagger-d)\gamma_1-ir(d+d^\dagger)\gamma_3.
\]

精确算符运算给出 \(H_E=H_{EM}+e_d I\)。最后一项只是随时间变化的整体能量平移，不影响共轭演化或波函数重叠模长。

在 Nambu 顺序

\[
\Psi=(c_1,c_2,d,c_1^\dagger,c_2^\dagger,d^\dagger)^T
\]

下，由 \([\Psi_\alpha,H_E]=\sum_\beta h_{\alpha\beta}\Psi_\beta\) 独立提取出

\[
\boxed{h=\begin{pmatrix}
2E&0&p+q&0&0&p-q\\
0&0&ir&0&0&ir\\
p+q&-ir&2e_d&-p+q&-ir&0\\
0&0&-p+q&-2E&0&-p-q\\
0&0&ir&0&0&ir\\
p-q&-ir&0&-p-q&-ir&-2e_d
\end{pmatrix}.}
\tag{5}
\]

这个矩阵来自普通费米子算符运算，没有先输入五维 \(A\)。定义酉矩阵 \(T\) 使

\[
\frac1{\sqrt2}(\gamma_1,\gamma_2,\gamma_3,\gamma_s,\gamma_a,\gamma_b)^T=T\Psi.
\]

变换后的 \(ThT^\dagger\) 中，\(\gamma_s\) 对应的行列严格为零，剩余活跃五维块精确等于 \(2iA\)。因此 Nambu 波函数变换到活跃 Majorana 坐标后满足

\[
i\dot\psi=(2iA)\psi,\qquad \dot\psi=2A\psi,
\]

再次得到式 (4)。不耦合的第六个模式可精确删去，这里没有作低能投影近似。

## 4. 当前 PDE 与接口的核对

程序按 \(E,q,p,r,e_d\) 五个独立系数，分别计算式 (4) 的十个上三角分量，并与实际使用的 `state_action` 系数逐项比较。全部差值严格为零。因为生成元对这五个系数线性，这项恒等式核对覆盖任意系数值；不只是几个数值点的偶合。

代入 \(q=\eta g(t)\) 后，\(\eta=|t_1|\) 与 \(E\) 的响应就是对式 (4) 求导：

\[
\partial_t S_E=[2A,S_E]+[2L_{12},\mathcal B],
\]

\[
\partial_t S_\eta=[2A,S_\eta]+[2gL_{1b},\mathcal B].
\]

门控在一轮三个接口的值为

| 时刻 | \(p\) | \(r\) | \(e_d\) | \(q\) |
|---|---|---|---|---|
| \(0\) | 0 | 0 | \(E_0\) | 0 |
| \(\tau\) | \(t_c\) | 0 | 0 | \(\eta\) |
| \(2\tau\) | 0 | \(t_c\) | 0 | 0 |
| \(3\tau\) | 0 | 0 | \(E_0\) | 0 |

相邻段的左右控制值一致。固定基底下，\(\mathcal B,S_E,S_\eta\) 连续传递。第二轮从第一轮的全部终态出发；控制周期性不意味着状态周期性。\(t_j=j\tau\) 不随 \(E,\eta\) 改变，因此这两个参数响应无需额外的移动接口项。

## 5. 独立演化的数值检验

三条路线是：

1. 当前十分量 PDE，连同 \(E,\eta\) 响应，按六段分别积分。
2. 式 (3) 的四维旋量传播矩阵及其响应，在整个 \(0\le t/\tau\le6\) 一次性积分；之后通过共轭作用重构 \(\mathcal B\)。
3. 式 (5) 的六维 Nambu 波函数及其响应，一次性积分后直接计算目标波函数重叠，并重构 \(\mathcal B\)。

后两条路线使用单独实现的全局时钟，没有调用 PDE 的逐段门控函数。三者保持既定余弦门控、\(t_c=E_0=0.3\)、DOP853 容差和最大归一化步长 0.025。每个样本检查全程 121 个时刻，包含初始时刻、全部接口和终点。时长仍是可控设计量；下面只是检验实例：

| \(\tau\) | \(E\) | \(\eta\) | PDE 的终点 \(F\) | Nambu 重叠的终点 \(F\) |
|---:|---:|---:|---:|---:|
| 100 | 0 | 0 | 0.999984419337031 | 0.999984419337024 |
| 100 | 0.01 | 0.001 | 0.058552332886577 | 0.058552332885562 |
| 300 | 0.01 | 0.01 | 0.102518048673699 | 0.102518048678616 |
| 1000 | 0.01 | 0.1 | 0.000921389096247 | 0.000921389096233 |

除理想极限样本外，\(E_1\) 均为已确认的 \(0.01\)。表中的 \(F\) 是平方保真度，论文色标为 \(f=\sqrt F\)。四个样本、全过程以及两种独立路线的最大差异为：

| 核对量 | 最大差异 |
|---|---:|
| \(\mathcal B\) 分量 | \(8.12\times10^{-10}\) |
| \(F\) | \(7.67\times10^{-10}\) |
| \(S_E,S_\eta\) 归一化差异 | \(3.33\times10^{-7}\) |
| 保真度参数梯度归一化差异 | \(7.19\times10^{-8}\) |

归一化差异采用 \(|a-b|/(1+|a|+|b|)\)。这些是所选样本的数值误差，不是整个参数空间的误差界。系数对应的普遍性来自前面的精确算符恒等式。

## 6. 对此前论文图像差异的含义

本次验证完成了此前缺少的一环：从原始有效哈密顿量独立推出单粒子矩阵，并检查现有 PDE 的实现。按原式系数和所声明的 Clifford 归一化，\(K=2A\) 已由两条路线确定；先前 \(K=A\) 的试验对应不同的能量归一化，不能作为同一哈密顿量的修正方案。

因此，当前没有依据修改状态 PDE 或接口条件来拟合 Fig. 1(d)。给定有效模型及既定协议的演化，可以继续用这套几何框架研究。此前与原图亮带的定量差异仍然存在；本次验证不能证明作者实际绘图使用了什么数值约定，也不能据此宣称论文图错误。

实现：[verify_effective_hamiltonian.py](verify_effective_hamiltonian.py)；精确核对项目、矩阵及全部误差：[effective_hamiltonian_verification.json](effective_hamiltonian_verification.json)。

```sh
.venv/bin/python -B fixed_protocol_pde/verify_effective_hamiltonian.py
```
