### kit3-exp：配置空间、联络与 Dehn twist / half twist 的 exp 形式

本节是对 [kit3.md](kit3.md) 的补充，专门用 exp 形式来整理：

1. 把 YBE 解族 $R=e^{iK(\lambda)}$ 看成参数空间/配置空间上的“平行移动元”，$K(\lambda)$ 看成 $so(2N)$ 值的联络 1‑form；
2. 在 Majorana 零模子空间中，这个联络如何诱导出 braid、Dehn twist、half twist 的幺正表示；
3. 改写成 exp 形式后，对比线性 R 形式，在描述 Dehn twist/half twist 时多出来的细节与几何解释。

其中 $\lambda$ 可以是：
- YBE 解族的参数（如 $(J_x,J_y,J_z)$），
- 缺陷/涡旋位置的参数（多体配置空间），
- 或者两者的拼接（既改变微观耦合，又移动缺陷）。


---

#### 1. 从 R=e^{iK(\lambda)} 到参数空间上的联络

假设有一族满足常数 YBE 的幺正算符
$$
R(\lambda)=e^{iK(\lambda)},\qquad K(\lambda)\in \mathfrak{so}(2N)
$$
作用在某个有限维 Hilbert 空间（例如一段链上的若干 Majorana 模式）。这里 $\lambda$ 是某个参数向量（可以是标量或多维）。

在 1D/2D BdG 映射后，我们通常关心某个“低能子空间” $\mathcal H_{\mathrm{low}}(\lambda)$，由一组 Majorana 零模 $\{\gamma_a(\lambda)\}$ 张成。沿参数路径
$$
\lambda:[0,1]\to \mathcal M
$$
（$\mathcal M$ 是参数/配置空间），若演化足够绝热，系统在 $\mathcal H_{\mathrm{low}}$ 中的平行移动由某个 Berry‑Wilson 线给出
$$
U_{\gamma}=\mathcal P\exp\oint_\gamma \mathcal A,
$$
其中 $\mathcal A$ 是 Berry 连接 1‑form。对于 Majorana 零模而言，可以把 $\mathcal A$ 写成 $so(2n)$ 值：
$$
\mathcal A = \frac12 \sum_{a<b}\Omega_{ab}(\lambda)\,\gamma_a\gamma_b\,d\lambda,
$$
其中 $\Omega_{ab}=-\Omega_{ba}$。

**exp 视角下的桥梁**是：

- 在微观层面，我们用一列 $R(\lambda) = e^{iK(\lambda)}$（例如分段常数的门序列）来实现参数路径上的演化；
- 在低能 Majorana 子空间中，这些门的累积作用压缩为
  $$
  U_{\gamma}^{\mathrm{(low)}} = \mathcal P\exp\oint_\gamma \frac12\Omega_{ab}(\lambda)\,\gamma_a\gamma_b\,d\lambda,
  $$
  这与“直接用 H_{\mathrm{eff}}(\lambda)$ 的 adiabatic 演化”等价；
- 而 $K(\lambda)$ 就是把“局域 R”与“低能联络”联系起来的生成元：
  - 边/局域上的 $K_{ij}(\lambda)$ 决定了 Majorana 模式的瞬时耦合矩阵；
  - 组合与投影给出 $\Omega_{ab}(\lambda)$。

YBE 在 exp 视角下可以理解为：

- 对于三粒子/三点配置，其局部平行移动的组合（$(R\otimes I)(I\otimes R)(R\otimes I)$ 与另一侧）相等；
- 这对应于“局域联络在三角形配置空间上是平坦的”，即曲率 $F=d\mathcal A+\mathcal A\wedge\mathcal A$ 在对应子空间为零（或者是中心元，对应 braid group 的表现）。


---

#### 2. Majorana 零模子空间里的联络：显式形式

设在某个拓扑相区内，我们有 $2n$ 个 Majorana 零模 $\gamma_1,\dots,\gamma_{2n}$，满足
$$
\{\gamma_a,\gamma_b\}=2\delta_{ab}.
$$
它们张成的 $2^n$ 维 Hilbert 空间可以看作一个 $\mathrm{Spin}(2n)$ 表示。任意的 $so(2n)$ 元
$$
X=\frac12\sum_{a<b}\Omega_{ab}\,\gamma_a\gamma_b
$$
指数化后得到 $\mathrm{Spin}(2n)$ 中的元 $e^{X}$。

对 exp‑YBE 方案而言，我们常见的生成元有两类：

1. **局域交换/绕行生成元（half twist 类型）**：
   $$
   X_{ab}(\theta)=\frac{\theta}{2}\,\gamma_a\gamma_b,\qquad U_{ab}(\theta)=e^{X_{ab}(\theta)}.
   $$
   - 在 Ising 任意子理论中，$\theta=\pi/2$ 对应标准的 $\sigma$‑$\sigma$ 交换（up to phase）。
   - 这类幺正正是“half twist”的 Majorana 表示：把两个端点/缺陷在配置空间中交换一次。 

2. **沿某条非平凡曲线的 Dehn twist 生成元**：
   - 把 Dehn twist $T_\gamma$ 看成沿某闭合曲线 $\gamma$ 绕行一圈的 mapping class 元；
   - 在 Majorana 表示中，$T_\gamma$ 的作用往往可以写成若干 $X_{ab}$ 的组合：
     $$
     \rho(T_\gamma)=\exp\Big(\frac12\sum_{a<b}\Theta^{(\gamma)}_{ab}\,\gamma_a\gamma_b\Big).
     $$
   - 这里 $\Theta^{(\gamma)}$ 是沿配置空间路径的联络积分：
     $$
     \Theta^{(\gamma)}_{ab}=\oint_\gamma\Omega_{ab}(\lambda)\,d\lambda.
     $$

改写成 exp 形式后，“联络”就是 $\Omega_{ab}(\lambda)$，“holonomy”就是这些 $\Theta^{(\gamma)}_{ab}$，而 half twist/Dehn twist 都是 $so(2n)$ 指数的特例。


##### 2.1 从 BdG 本征态到 Majorana 联络的具体公式

上面把联络写成了抽象的 $so(2n)$ 1‑form，下面给出从 BdG/Majorana 哈密顿量出发推导 $\Omega_{ab}(\lambda)$ 的一个比较“工程可实现”的公式。

**(a) BdG 形式与零能子空间**

对给定参数点 $\lambda$，设 Majorana 模式 $c_j$（$j=1,\dots,2N$）的哈密顿量写成标准二次型
$$
H(\lambda)=\frac i4\sum_{j,k}A_{jk}(\lambda)c_j c_k,
$$
其中 $A(\lambda)$ 是实反对称矩阵。对 $iA(\lambda)$ 做实正交对角化：
$$
O^T(\lambda)\, A(\lambda)\, O(\lambda)=\bigoplus_{m=1}^N \begin{pmatrix}0 & \epsilon_m(\lambda)\\ -\epsilon_m(\lambda) & 0\end{pmatrix}.
$$
把旋转后的 Majorana 记为
$$
	ilde c_j(\lambda)=\sum_k O_{kj}(\lambda) c_k,
$$
于是哈密顿量成为
$$
H(\lambda)=\frac i2\sum_{m=1}^N \epsilon_m(\lambda)\, \tilde c_{2m-1}(\lambda)\,\tilde c_{2m}(\lambda).
$$

若存在 $2n$ 个严格零能模（或在能隙内近似零的模），可以把对应的 $2n$ 个 Majorana 记为
$$
\gamma_a(\lambda)=\tilde c_{j_a}(\lambda),\qquad a=1,\dots,2n,
$$
它们张成零能子空间。剩余模式具有非零能隙，可在绝热条件下视为始终处于真空占据态。

**(b) Berry 联络的定义**

在零模 Hilbert 空间中选一组规一正交基 $\{|\Psi_\alpha(\lambda)\rangle\}$，Berry 联络的矩阵元为
$$
\mathcal A_{\alpha\beta}(\lambda)=i\,\langle\Psi_\alpha(\lambda)|\partial_\lambda\Psi_\beta(\lambda)\rangle.
$$
这给出一个 $u(2^n)$ 值 1‑form。对自旋/费米系统而言，$\mathcal A$ 可以始终选成 $so(2n)$ 的自旋表示，其生成元恰好是二次的 Majorana 双线性 $\gamma_a\gamma_b$。

具体地，任意的零模基变换都可以写成 $\mathrm{Spin}(2n)$ 元作用：
$$
|\Psi_\alpha'\rangle = U(\lambda)|\Psi_\alpha\rangle,\qquad U(\lambda)=\exp\Big(\frac12\sum_{a<b}\chi_{ab}(\lambda)\,\gamma_a\gamma_b\Big).
$$
在这个意义下，$\mathcal A$ 可以展开成
$$
\mathcal A(\lambda)=\frac12\sum_{a<b}\Omega_{ab}(\lambda)\,\gamma_a\gamma_b\,d\lambda
$$
的形式，其中系数
$$
\Omega_{ab}(\lambda)\in\mathbb R,\qquad \Omega_{ab}=-\Omega_{ba}
$$
就是我们在前面抽象写下的“Majorana 联络”。

**(c) 从 Majorana 基底的 \(O(\lambda)\) 变化提取 \(\Omega\)**

为了得到 $\Omega_{ab}(\lambda)$ 的具体表达，可以追踪 Majorana 基底本身随参数如何变化。设
$$
\vec c=(c_1,\dots,c_{2N})^T,\qquad \vec\gamma(\lambda)=(\gamma_1(\lambda),\dots,\gamma_{2n}(\lambda))^T.
$$
在零模子空间相关的 $2n$ 维子块中，总可以写成
$$
\vec\gamma(\lambda)=Q(\lambda)\,\vec c,
$$
其中 $Q(\lambda)$ 是一个 $2n\times 2N$ 的实矩阵，满足 $Q(\lambda)Q^T(\lambda)=I_{2n}$（在这个子空间上是正交变换）。

对 $\lambda$ 求导：
$$
\partial_\lambda\vec\gamma(\lambda)=\dot Q(\lambda)\,\vec c.
$$
另一方面，根据 Majorana 算符在 Heisenberg picture 下“随基底变化”的一般形式，可以写成某个 $so(2n)$ 元对它们的无穷小旋转：
$$
\partial_\lambda \gamma_a(\lambda)=\sum_b K_{ab}(\lambda)\,\gamma_b(\lambda),\qquad K(\lambda)\in\mathfrak{so}(2n).
$$
将两种写法比较，可得
$$
K(\lambda)=\big(\partial_\lambda Q(\lambda)\big)Q^T(\lambda).
$$
因为 $Q(\lambda)Q^T(\lambda)=I$，直接可见 $K(\lambda)$ 反对称，从而 $K(\lambda)\in\mathfrak{so}(2n)$。把指标展开就是
$$
K_{ab}(\lambda)=\sum_j \partial_\lambda Q_{aj}(\lambda)\,Q_{bj}(\lambda),\qquad K_{ab}=-K_{ba}.
$$

这时，可以把 Majorana 联络 $\Omega(\lambda)$ 与 $K(\lambda)$ 直接等同（差别只在规范选择与可能的整体 U(1) 因子上）：
$$
\Omega_{ab}(\lambda)=K_{ab}(\lambda).
$$
于是
$$
\mathcal A(\lambda)=\frac12\sum_{a<b}K_{ab}(\lambda)\,\gamma_a\gamma_b\,d\lambda.
$$
这给出了一个可操作的公式：

- 数值上，只要在每个参数点 $\lambda$ 求出零模的 Majorana 波函数（即 $Q(\lambda)$ 行向量），
- 然后沿参数方向对 $Q(\lambda)$ 做差分/导数，就能得到 $K_{ab}(\lambda)$，
- 再把它代入上式即可构造联络与路径有序指数 $\mathcal P\exp\oint\mathcal A$。

**(d) 平行移动算符的形式**

给定一条参数路径 $\lambda(t)$，$t\in[0,1]$，对零模算符 $\gamma_a(t)\equiv\gamma_a(\lambda(t))$，由上面的无穷小变换可写出平行移动的微分方程：
$$
\frac{d}{dt}\gamma_a(t)=\sum_b K_{ab}(\lambda(t))\,\gamma_b(t)\,\dot\lambda(t).
$$
其解是
$$
\vec\gamma(t)=\mathcal P\exp\Big(\int_0^t K(\lambda(t'))\,\dot\lambda(t')\,dt'\Big)\,\vec\gamma(0).
$$
对应在 Hilbert 空间上的表示就是前面写的
$$
U_{\gamma}=\mathcal P\exp\oint_\gamma \mathcal A,\qquad \mathcal A=\frac12\sum_{a<b}K_{ab}(\lambda)\,\gamma_a\gamma_b\,d\lambda,
$$
这就是 exp 形式下“联络 + 平行移动”的完整具体表达：

- $K(\lambda)=(\partial_\lambda Q)Q^T\in\mathfrak{so}(2n)$ 是在 Majorana 指标空间中的联络；
- $\mathcal A=(1/2)K_{ab}\gamma_a\gamma_b\,d\lambda$ 是其在零模 Hilbert 空间上的自旋表示；
- 对应路径有序指数给出配置空间路径的 holonomy，即 braid/Dehn/half twist 在零模 Hilbert 空间上的幺正作用。


---

#### 3. Dehn twist 与 half twist：从 exp 形式看区别

在 mapping class group 的语言中：

- **half twist** 通常指沿某条小短弧把两个端点（缺陷/任意子）交换一次；
- **Dehn twist** 则指沿某个闭合曲线做一次 $2\pi$ 扭转（在拓扑上是切开环管、扭转 360° 再粘回）。

在 Majorana/exp 形式下，它们的区别可以更清楚地写成“路径与积分区间”的区别：

1. **half twist：局域交换沿短弧**

   - 取一条在配置空间中连接两个缺陷的路径 $\gamma_{ab}$，对应把粒子 $a,b$ 交换：
     $$
     U_{\text{half}}(a,b) = \mathcal P\exp\oint_{\gamma_{ab}} \mathcal A.
     $$
   - 若 Majorana 模式主要集中在这两个缺陷附近，且沿路径变化温和，则可以近似为
     $$
     U_{\text{half}}(a,b) \approx \exp\Big(\frac{\theta_{ab}}{2}\,\gamma_a\gamma_b\Big),
     $$
     其中 $\theta_{ab}$ 是某个“交换角”，取值依联络而定（在 Ising 情形，$\theta_{ab}=\pm\pi/2$ 给出标准 braid 相位）。

2. **Dehn twist：沿闭合曲线的全局扭转**

   - 选一条拓扑上非平凡的闭合曲线 $\gamma$（如环面上的 $a$‑cycle 或 $b$‑cycle），对应 Dehn twist $T_\gamma$；
   - 其在 Majorana 子空间中的 holonomy 为
     $$
     \rho(T_\gamma)=\mathcal P\exp\oint_\gamma \mathcal A=\exp\Big(\frac12\sum_{a<b}\Theta^{(\gamma)}_{ab}\,\gamma_a\gamma_b\Big).
     $$
   - 这一般会在零模 Hilbert 空间上给出一个对角/块对角幺正，其本征值与 TQFT 中的 topological spin $\theta_x$ 相关：
     $$
     \rho(T_\gamma)\big|_{\text{sector }x}\sim e^{2\pi i (h_x - c/24)}
     $$
     （此处不展开 CFT/TQFT 细节，只强调：exp 形式允许你直接在 Majorana 语言中测出这些本征相位）。

**exp 形式带来的清晰点是**：

- 两种 twist 的区别不在于算符的“代数形状”（两者都是 $e^{X}$），而在于：
  - half twist 的路径是连接两点的“开放路径”，只积累局域的 $\Omega_{ab}$；
  - Dehn twist 的路径是“闭合的拓扑不平凡回路”，积分区域更大，可能涉及整条链/整片区域上的多个 $\gamma_a$，从而在多体空间中给出不同的全局相。
- 在 exp 视角下，这两类操作的 difference = “联络沿不同路径/不同同伦类积分的结果不同”，这与 kit3 里“上层空间结构”/mapping class 的直觉完全统一。


##### 3.1 half twist 的算符推导与几何解释

这里把 half twist 的算符形式、对 Majorana/费米算符的作用以及其在配置空间中“交换两个粒子”的几何含义写得更具体一些。

**(a) half twist 的代数定义**

选定一对 Majorana 零模 \(\gamma_a,\gamma_b\)。定义生成元
$$
X_{ab}(\theta)=\frac{\theta}{2}\,\gamma_a\gamma_b,\qquad U_{ab}(\theta)=e^{X_{ab}(\theta)}.
$$
注意 \(\gamma_a\gamma_b\) 满足
$$
(\gamma_a\gamma_b)^2=-1,\qquad (\gamma_a\gamma_b)^{2k}=(-1)^k,\ (\gamma_a\gamma_b)^{2k+1}=(-1)^k\gamma_a\gamma_b.
$$
于是指数可以显式和有限地展开：
$$
U_{ab}(\theta)=\cos\frac{\theta}{2}+\gamma_a\gamma_b\,\sin\frac{\theta}{2}.
$$

**(b) 对 Majorana 算符的作用**

用共轭变换计算 \(U_{ab}(\theta)\gamma_a U_{ab}(\theta)^{-1}\) 等。利用 Baker–Campbell–Hausdorff 或直接代数运算，可得：
$$
\begin{aligned}
U_{ab}(\theta)\,\gamma_a\,U_{ab}(\theta)^{-1} &= \gamma_a\cos\theta + \gamma_b\sin\theta,\\
U_{ab}(\theta)\,\gamma_b\,U_{ab}(\theta)^{-1} &= \gamma_b\cos\theta - \gamma_a\sin\theta,
\end{aligned}
$$
而对与 \(\gamma_a,\gamma_b\) 无关的其它零模 \(\gamma_c\)（对易/反对易但不直接耦合），有
$$
U_{ab}(\theta)\,\gamma_c\,U_{ab}(\theta)^{-1} = \gamma_c.
$$
因此，\(U_{ab}(\theta)\) 在 \(\{\gamma_a,\gamma_b\}\) 张成的二维实向量空间中实现了一个 SO(2) 旋转：
$$
\begin{pmatrix}\gamma_a'\\ \gamma_b'\end{pmatrix}=\begin{pmatrix}\cos\theta & \sin\theta\\ -\sin\theta & \cos\theta\end{pmatrix}\begin{pmatrix}\gamma_a\\ \gamma_b\end{pmatrix}.
$$
几何上这是在“零模空间”里围绕 \(\gamma_a\gamma_b\) 的一个旋转。

在 Ising 相中，取 \(\theta=\pi/2\) 得到
$$
U_{ab}\equiv U_{ab}(\pi/2)=e^{(\pi/4)\gamma_a\gamma_b},
$$
即
$$
\gamma_a' = \frac{1}{\sqrt2}(\gamma_a+\gamma_b),\qquad \gamma_b' = \frac{1}{\sqrt2}(\gamma_b-\gamma_a),
$$
与 1D 例子中计算是一致的。

**(c) 对复费米模式的作用与融合通道**

若把这两 Majorana 组合成一个复费米算符
$$
f=\frac{\gamma_a+i\gamma_b}{2},\qquad f^\dagger=\frac{\gamma_a-i\gamma_b}{2},\qquad n=f^\dagger f,
$$
则有
$$
i\gamma_a\gamma_b = 2n-1.
$$
用 \(U_{ab}(\theta)=e^{\frac{\theta}{2}\gamma_a\gamma_b}\) 作用在 \(|0\rangle,|1\rangle\) 上：
$$
\begin{aligned}
U_{ab}(\theta)|0\rangle &= e^{-i\theta/2}|0\rangle,\\
U_{ab}(\theta)|1\rangle &= e^{+i\theta/2}|1\rangle.
\end{aligned}
$$
这是一个对角幺正，其本征值只取决于占据数奇偶。若把 \(|0\rangle,|1\rangle\) 识别为两种融合通道（例如 Ising 中的 \((\sigma\sigma)\to1,\psi\)），这一幺正作用就是 braid R‑矩阵在该两维融合空间上的具体表示。

**(d) 在配置空间中的几何含义**

在配置空间图像中，选两缺陷的“标号”为 1,2，它们的位置为 \(x_1,x_2\)，并在它们附近分别有零模 \(\gamma_a,\gamma_b\)。沿某条短弧 \(\gamma_{ab}\) 把 1 绕 2 逆时针绕行一次，回到配置空间中同一无序点集的位置，对应的同伦类就是 braid 群生成元 \(\sigma_1\)。在半经典/绝热极限下，这条路径的 holonomy 即为上面的 \(U_{ab}(\theta)\)。

因此，**half twist = 沿连接两点的开放路径的 braid 操作 = 在零模空间中绕 \(\gamma_a\gamma_b\) 的 SO(2) 旋转**；是否是“拓扑的”，取决于：

- 路径是否可以连续变形（不穿过其它缺陷或闭合谱隙），
- 联络是否平坦或 projectively flat（来自 YBE/拓扑相）。


##### 3.2 Dehn twist 的算符推导与拓扑含义

Dehn twist 更“全局”，对应在含有非平凡周期的曲面上对某条闭合曲线做 2\(\pi\) 扭转。这里用环面上的简单情形说明其算符形式与拓扑含义。

**(a) Dehn twist 作为 mapping class 元素**

设基底曲面 \(\Sigma\) 为环面 \(T^2\)；其 mapping class group 由 SL(2,\(\mathbb Z\)) 生成，其中一个生成元 \(T\) 就是沿 \(a\)‑cycle 的 Dehn twist：把沿 \(a\) 向的环切开，扭转 2\(\pi\) 再粘回。对拓扑量子场论而言，这个 \(T\) 在 Hilbert 空间上的作用由所谓 T‑矩阵给出，其本征值是 topological spin 的指数，如
$$
T_{xx}=e^{2\pi i(h_x-c/24)}.
$$

**(b) 在 Majorana/BdG 模型中的实现**

对于一个在环面上实现的拓扑超导/Kitaev 蜂窝‑类模型，基态 Hilbert 空间往往是拓扑简并的（例如 Ising 相在环面上有三个基态，对应 \(1,\psi,\sigma\) 扇区）。Dehn twist 可在晶格/微观层面实现为：

- 在时空中构造一个“拉伸 + 重新贴边”的演化过程，相当于沿 \(a\)‑cycle 对哈密顿量参数进行周期性调节；
- 或者用 braid 语言，把 Dehn twist 写成一串绕整个周期的 braid word（例如在存在穿越周期的 Wilson 线/缺陷串时）。

无论哪个实现，若演化足够绝热，沿该闭合路径 \(\Gamma_T\) 的 holonomy 都可以写成
$$
U_{T} = \mathcal P\exp\oint_{\Gamma_T}\mathcal A = \exp\Big(\frac12\sum_{a<b}\Theta^{(T)}_{ab}\,\gamma_a\gamma_b\Big),
$$
其中 \(\Theta^{(T)}\) 是沿这一“全局操作”路径积分得到的角度矩阵。对于 Ising 相，这个 \(U_T\) 在适当基底下与 T‑矩阵同构，其本征值集合包括 \(\{1, e^{-i\pi/8}, e^{+3i\pi/8}\}\)（up to overall phase）。

**(c) 与 half twist 的关系：由局域 braid 组成的全局扭转**

在许多具体构造中，一个 Dehn twist 可以在配置空间中被分解为若干 half twist（braid 生成子）的组合。例如，在一条包含若干缺陷的闭合链上，沿链方向依次对相邻对施加 half twist（或其平方），可以等效于“在该闭合链上做一次扭转”。

在 Majorana 表示中，这就对应于
$$
U_T = \prod_k U_{a_k b_k}(\theta_k)= \exp\Big(\frac12\sum_{a<b}\Theta^{(T)}_{ab}\,\gamma_a\gamma_b\Big),
$$
其中每个 \(U_{a_k b_k}(\theta_k)\) 是某个 half twist（\(\theta_k\approx\pm\pi/2\)），积累的 \(\Theta^{(T)}\) 只依赖于这一串 half twist 的同伦类（即在 braid/mapping class 群中的元素），与具体路径细节无关。这种“只依赖同伦类”的性质正是拓扑性的体现。

**(d) 是否能用来表达拓扑态：路径无关性与能隙保护**

从上述推导可以看出，half twist 和 Dehn twist 能否真正描述/区分拓扑态，取决于两个关键条件：

1. **能隙保护**：在整个演化路径 \(\Gamma\) 上，体系保持 gapped，使得低能零模子空间 \(\mathcal H_{\mathrm{low}}\) 与高能励起分开，保证绝热演化和 Berry 联络的良好定义。
2. **联络的平坦/投影平坦性**：若联络 \(\mathcal A\) 的曲率 \(F=d\mathcal A+\mathcal A\wedge\mathcal A\) 在给定拓扑相中要么为零，要么仅给出中心元（相位），那么 \(U_\Gamma\) 只依赖于 \([\Gamma]\in\pi_1(\mathcal C_n)\)，而与具体代表路径无关。这时，\(\rho: \pi_1(\mathcal C_n)\to U(\mathcal H_{\mathrm{low}})\) 就给出了一个真正的“拓扑量子场论的表示”。

在 Ising 例子中，这两个条件在非阿贝尔拓扑相中得到满足：

- 能隙由拓扑超导/蜂窝模型的 bulk gap 提供；
- 联络平坦性则源自 YBE 约束和拓扑相的有效 TQFT 描述，使得 braid/Dehn 元素的表示仅由融合规则与 F、R 符号决定，与微观路径无关。

因此，half twist 与 Dehn twist 不仅能“表达”拓扑态（通过它们在基态 Hilbert 空间上的矩阵形式来区分不同拓扑扇区），而且在合适条件下，其作用就是对应该拓扑相 TQFT 的 braid/mapping class 群表示。这也说明了：

> 在 exp(iK) + Majorana 联络的框架下，只要 bulk gapped 且联络 projectively flat，沿配置空间路径的平行移动算符（特别是 half twist、Dehn twist）就可以被视为定义/操控拓扑量子态的基本“拓扑门”。


---

#### 4. 与 YBE 的平坦性条件对接

常数 YBE
$$
(R\otimes I)(I\otimes R)(R\otimes I)=(I\otimes R)(R\otimes I)(I\otimes R)
$$
在 exp 形式下可以理解为：

- 把三点配置想象成一个“小三角形”的配置空间；
- 左边和右边是两条不同的分段路径，从一个角走到另一个角；
- 若 $R=e^{iK}$，那么这两种路径对应的 holonomy 为
  $$
  U_L=\mathcal P\exp\int_{\text{path }L}\mathcal A,\quad
  U_R=\mathcal P\exp\int_{\text{path }R}\mathcal A,
  $$
  而 YBE 要求 $U_L=U_R$；
- 这在连续极限下就是“联络在该区域平坦（或中心）”：曲率 $F=d\mathcal A+\mathcal A\wedge\mathcal A$ 在该三角形内为零（或落在中心元，对应编织相位）。

对你在 kit3 中的计划而言，这意味着：

- 从代数可积的 YBE 族出发，本身就自带一个“平坦联络”的结构；
- 把 R 写成 exp(iK) 后，这个结构可以直接解释为“在缺陷/参数配置空间上的平坦（或 projectively flat）联络”；
- 进一步把这联络下的 holonomy 投影到 Majorana 零模子空间，就得到 braid、Dehn twist、half twist 的具体幺正表示。


---

#### 5. 小结：kit3 视角下 exp 形式的作用

相对于线性 R 形式，在配置空间/holonomy 这一层，exp 形式提供了：

- 一个显式的联络描述：$K(\lambda)$ 和投影后的 $\Omega_{ab}(\lambda)$ 是 $so(2N)$ 值的 1‑form，holonomy 即为其路径有序指数；
- Dehn twist / half twist 都统一为“在某个路径/闭合回路上对联络积分再指数化”的特例，只是路径拓扑类型不同；
- YBE 的内容可重述为“局域联络在某些基本三角/胞腔上的平坦性条件”，与 mapping class group/TQFT 中的 F,R 符号平坦性条件呼应。

这使得你在 kit3 中构想的“从 YBE → R → Majorana → mapping class / TQFT 上层结构”的路线，在数学结构上更加统一：从一开始就工作在“exp(iK) + 平坦联络 + holonomy”的语言，而不是先用线性参数再事后取对数。这样，在后续数值实验里，一旦你从 BdG/Majorana 哈密顿量中抽取出低能耦合矩阵 $A_{ab}(\lambda)$，就可以直接构造对应的 $\mathcal A,\Theta^{(\gamma)}$，并与 TQFT 里的 Dehn/half twist 数据对比。


---

#### 6. Ising 例子：从 exp(iK) 到标准 braid 矩阵

下面给出一个尽量贴近前文的简单 Ising 例子，说明：

1. 如何从具体的 $K$（由 YBE 解/自旋哈密顿量得到）构造出一个交换两 Majorana 零模的算符 $U_{ab}(\theta)=\exp(\tfrac{\theta}{2}\gamma_a\gamma_b)$；
2. 这个 $U_{ab}(\theta)$ 在适当选择的 $\theta$ 下，其本征相位与 Ising 任意子理论中 $\sigma$‑$\sigma$ 交换的 R 矩阵本征值一致（up to overall phase）。


**6.1 选取一个具体的 exp(iK) 局域生成元**

在 [R_to_Kitaev.md](R_to_Kitaev.md) 和 [R_to_kitaev2.md](R_to_kitaev2.md) 中，我们已经展示过这样一个简单的局域选择：

- 取单键上的参数 $a=0,d=0,b=1,c=0$，则
  $$
  K_{i,i+1}=J_x\sigma^x_i\sigma^x_{i+1},\qquad J_x=1,J_y=J_z=0,
  $$
  对应的 1D Kitaev 链参数满足 $t=\Delta=1$、$\mu=0$；
- 在 Jordan–Wigner + Majorana 表示中，这一键上的主耦合为
  $$
  K_{i,i+1}\;\leadsto\; i\gamma_2\gamma_3,
  $$
  其中 $\gamma_2,\gamma_3$ 是该键中间的两 Majorana（具体标号见原文的 Majorana 排列约定）。

把这个 $K$ 当作局域哈密顿密度，短时间演化算符为
$$
U(\tau)=e^{-iK\tau}=\exp\big(-i\tau\,i\gamma_2\gamma_3\big)=\exp\big(\tau\,\gamma_2\gamma_3\big).
$$
选取 $\tau=\tfrac{\pi}{4}$，得到
$$
U\equiv U(\tau=\tfrac{\pi}{4})=\exp\Big(\frac{\pi}{4}\gamma_2\gamma_3\Big),
$$
这正是前文所写的标准局域 half twist 生成元 $U_{23}(\theta)$，其中 $\theta=\pi/2$（因为 $X_{23}(\theta)=\tfrac{\theta}{2}\gamma_2\gamma_3$）。

因此，这个具体的 exp(iK) 选择在低能 Majorana 子空间中诱导出
$$
U_{23}=\exp\Big(\frac{\pi}{4}\gamma_2\gamma_3\Big),
$$
我们现在只需把它的本征值与 Ising 任意子的 R 矩阵比较。


**6.2 在二 Majorana 子空间中对角化 U=exp((\pi/4)\gamma_2\gamma_3)**

定义复费米模式
$$
f=\frac{\gamma_2+i\gamma_3}{2},\qquad f^\dagger=\frac{\gamma_2-i\gamma_3}{2},\qquad n=f^\dagger f.
$$
可以检验
$$
i\gamma_2\gamma_3 = 2n-1,
$$
因此 $i\gamma_2\gamma_3$ 的本征值在基 $\{|0\rangle,|1\rangle\}$（$n|0\rangle=0,\;n|1\rangle=|1\rangle$）中为
$$
i\gamma_2\gamma_3|0\rangle=+|0\rangle,\qquad i\gamma_2\gamma_3|1\rangle=-|1\rangle.
$$
而
$$
\gamma_2\gamma_3 = -i(2n-1),
$$
其本征值为 $\lambda_0=-i$（在 $|0\rangle$ 上）与 $\lambda_1=+i$（在 $|1\rangle$ 上）。

于是
$$
U=\exp\Big(\frac{\pi}{4}\gamma_2\gamma_3\Big)
$$
在这两个本征态上的本征值为
$$
\begin{aligned}
U|0\rangle &= e^{(\pi/4)\lambda_0}|0\rangle = e^{(\pi/4)(-i)}|0\rangle = e^{-i\pi/4}|0\rangle,\\
U|1\rangle &= e^{(\pi/4)\lambda_1}|1\rangle = e^{(\pi/4)(+i)}|1\rangle = e^{+i\pi/4}|1\rangle.
\end{aligned}
$$
因此 $U$ 的本征相位集合为
$$
\{e^{-i\pi/4},\;e^{+i\pi/4}\}.
$$


**6.3 与 Ising 任意子的 R^{\sigma\sigma} 矩阵比较**

在 Ising 任意子理论中，两条 $\sigma$ 线的 braid R 矩阵在融合基 $\{|(\sigma\sigma)\to 1\rangle, |(\sigma\sigma)\to \psi\rangle\}$ 上为对角形式
$$
R^{\sigma\sigma}=\begin{pmatrix}R^{\sigma\sigma}_1 & 0\\ 0 & R^{\sigma\sigma}_\psi\end{pmatrix},
$$
其中（取一组常见约定）
$$
R^{\sigma\sigma}_1=e^{-i\pi/8},\qquad R^{\sigma\sigma}_\psi=e^{3i\pi/8}.
$$
把我们得到的本征值
$$
e^{-i\pi/4},\quad e^{+i\pi/4}
$$
同时乘以一个整体相位 $e^{+i\pi/8}$（物理上不可观测），得到
$$
e^{+i\pi/8}\cdot e^{-i\pi/4}=e^{-i\pi/8},\qquad e^{+i\pi/8}\cdot e^{+i\pi/4}=e^{3i\pi/8},
$$
这正好与 $\{R^{\sigma\sigma}_1, R^{\sigma\sigma}_\psi\}$ 一致。

也就是说：

- 把 $|0\rangle$（偶费米数）识别为 $(\sigma\sigma)\to 1$ 融合通道，$|1\rangle$（奇费米数）识别为 $(\sigma\sigma)\to \psi$；
- 则由局域 exp(iK) 生成的 Majorana 交换算符
  $$
  U_{23}=\exp\Big(\frac{\pi}{4}\gamma_2\gamma_3\Big)
  $$
  在这两个通道上的本征相位，up to 整体相位，正好等于 Ising TQFT 中的 $R^{\sigma\sigma}_1$ 和 $R^{\sigma\sigma}_\psi$。

从联络/holonomy 的角度看，这意味着：

- $\Omega_{23}(\lambda)$ 沿一条“half twist 路径” $\gamma_{23}$ 的积分
  $$
  	heta_{23}=\oint_{\gamma_{23}}\Omega_{23}(\lambda)\,d\lambda
  $$
  在这个例子里等于 $\pi/2$；
- 相应 holonomy $U_{23}=\exp(\tfrac12\theta_{23}\gamma_2\gamma_3)$ 的本征值与 Ising 任意子的 braid R 矩阵完全匹配；
- 若再考虑 Dehn twist（例如把包含两 $\sigma$ 的环整圈扭转），则可以由适当的 half twist 组合构造，对应于某个 $\Theta^{(\gamma)}_{ab}$，其本征相位给出 topological spin $\theta_\sigma=e^{i\pi/8}$，与上面整体相位的选择一致。

这给出了一个具体的、完全在 exp(iK) 与 Majorana 语言内部的 Ising 例子，验证了：

> 通过把 YBE 解写成 $R=e^{iK}$ 并投影到零能 Majorana 子空间，可以直接得到 Ising 任意子理论中的 half twist（braid）相位；Dehn twist 则是这些 half twist 沿闭合曲线的组合，其本征相位与 $\theta_\sigma=e^{i\pi/8}$ 等 topological spin 一致。


---

#### 7. 配置空间视角：缺陷位置空间、braid 群与 Dehn/half twist

为把 kit3 系列里“配置空间 + mapping class”的直观与上面的 exp(iK)–Majorana 形式严格对上，这里整理一下：

1. 多体配置空间 $\mathcal C_n$；
2. 其基本群与 braid 群/映射类群的关系；
3. exp(iK) 在 Majorana 子空间中给出 $\pi_1(\mathcal C_n)$ 的一个幺正表示；
4. 其中哪些生成元对应 half twist，哪些对应 Dehn twist。


**7.1 配置空间 $\mathcal C_n$ 与 braid 群**

考虑一个 2D 基底流形 $\Sigma$（例如平面、圆盘、环面等），在其上放 $n$ 个不可区分的拓扑缺陷（如 vison、涡旋或端点）。配置空间为
$$
\mathcal C_n(\Sigma)=\frac{\{(x_1,\dots,x_n)\in \Sigma^n\mid x_i\ne x_j\}}{S_n},
$$
即“$n$ 个无序点在 $\Sigma$ 上的位置集合”（排除重合）。其基本群
$$
\pi_1\big(\mathcal C_n(\Sigma)\big)
$$
就是拓扑意义上的 braid 群/广义 braid 群（具体取决于 $\Sigma$ 与边界条件）。

在物理上，一条闭合路径 $\Gamma$ 在配置空间中是“把 $n$ 个缺陷在 $\Sigma$ 上移动一圈后回到原始位置（作为集合）的一个过程”，而这正是拓扑量子计算里所谓的“编织”（braiding）。


**7.2 exp(iK) 给出的 $\pi_1(\mathcal C_n)$ 的幺正表示**

给定一套微观 BdG/Majorana 模型（如 kit2‑exp 的二维 exp(iK)+Z2 模型），在每个缺陷配置下，有一个相应的 BdG 哈密顿量 $H(\lambda)$，其中 $\lambda\in\mathcal C_n$ 表示缺陷位置。对每个 $\lambda$，低能零模子空间
$$
\mathcal H_{\mathrm{low}}(\lambda)\cong\mathbb C^{2^m}
$$
由一组 Majorana 零模 $\{\gamma_a(\lambda)\}$ 张成。

沿着配置空间中的一条闭合路径 $\Gamma: [0,1]\to\mathcal C_n$，若演化足够绝热且能隙不闭合，系统在低能子空间内的演化算符为
$$
U_{\Gamma}=\mathcal P\exp\oint_{\Gamma}\mathcal A(\lambda),
$$
其中 $\mathcal A(\lambda)$ 是前文的 $so(2m)$ 值联络 1‑form。于是我们得到一个映射
$$
\rho: \pi_1(\mathcal C_n)\to U(\mathcal H_{\mathrm{low}}),\qquad [\Gamma]\mapsto U_{\Gamma},
$$
这就是“任意子编织在零模 Hilbert 空间上的幺正表示”。

在 exp(iK) 的具体实现下：

- 路径 $\Gamma$ 被离散为若干局域时空步骤，每一步由某些 $R_{ij}(\lambda)=e^{iK_{ij}(\lambda)}$ 的序列近似；
- 把这些序列在 BdG Hilbert 空间上相乘，再投影到零模子空间，即得到 $U_{\Gamma}$；
- 由于 $K(\lambda)$ 直接是 Majorana 双线性，其投影自然落在 $so(2m)$ 中，$U_{\Gamma}$ 落在 $\mathrm{Spin}(2m)$ 表示中。


**7.3 half twist 与 braid 群生成元**

对于平面/圆盘上的 $n$ 个缺陷，$\pi_1(\mathcal C_n)$ 就是标准的 Artin braid 群 $B_n$，其生成元通常记为 $\sigma_1,\dots,\sigma_{n-1}$，几何操作是：

- $\sigma_i$：把第 $i$ 个与第 $i+1$ 个缺陷绕行交换一次（通常约定为逆时针）。

在 Majorana 描述中：

- 若第 $i,i+1$ 个缺陷分别携带零模 $\gamma_{2i-1},\gamma_{2i}$（具体编号可按链顺序或几何顺序），
- 则与 $\sigma_i$ 对应的幺正算符可写为
  $$
  U(\sigma_i)=\exp\Big(\frac{\theta}{2}\,\gamma_{a(i)}\gamma_{b(i)}\Big),
  $$
  其中 $\theta$ 由联络的积分决定，而在 Ising 相中 $\theta=\pm\pi/2$ 给出标准 braid 相；
- 这些 $U(\sigma_i)$ 满足 braid 关系
  $$
  U(\sigma_i)U(\sigma_{i+1})U(\sigma_i)=U(\sigma_{i+1})U(\sigma_i)U(\sigma_{i+1}),\quad U(\sigma_i)U(\sigma_j)=U(\sigma_j)U(\sigma_i)~(|i-j|>1),
  $$
  这在 exp 形式下可以追溯到 YBE 的平坦性条件与 $so(2m)$ 代数的 Baker–Campbell–Hausdorff 结构。

因此：**half twist = braid 群生成元 $\sigma_i$ 的 Majorana 实现 = 某个 $so(2m)$ 元指数的幺正**。


**7.4 Dehn twist 与 mapping class group 元素**

当基底流形 $\Sigma$ 本身有非平凡拓扑（如环面、带柄的曲面）时，配置空间的基本群不仅包含 braid 群部分，还包含与基底映射类群（mapping class group）相关的额外生成元。典型的例子是在环面上有一条简并轨道/边界上的 Majorana 零模时，沿一个周期做 Dehn twist 会在零模 Hilbert 空间里引入一个非平凡相位。

在 exp(iK)+Majorana 框架下：

- 选一条闭合曲线 $\gamma\subset\Sigma$（如环面上的 $a$‑cycle），考虑把所有缺陷/边界条件一并沿着 $\gamma$ 做一个 Dehn twist；
- 对应的配置空间路径 $\Gamma_\gamma$ 是“把系统切开、绕 $\gamma$ 扭转 $2\pi$ 再粘回”的同伦类；
- 其 holonomy
  $$
  U_{T_\gamma}=U_{\Gamma_\gamma}=\exp\Big(\frac12\sum_{a<b}\Theta^{(\gamma)}_{ab}\,\gamma_a\gamma_b\Big)
  $$
  在每个拓扑电荷 sector 上给出一个 topological spin：例如 Ising 中的 $\theta_\sigma=e^{i\pi/8}$；
- 具体来说，在 Ising 例子中，我们已经看到 half twist 的 eigenphase 组合出 $R^{\sigma\sigma}_1,R^{\sigma\sigma}_\psi$，而适当组合这些 half twist（例如 $\sigma_i^2$ 或更长的 braid word）就可以实现沿包含线段的闭合曲线的 Dehn twist，其本征值正是 $e^{2\pi i h_\sigma}=e^{i\pi/8}$ 等。

在配置空间语言中，这意味着：

- braid 群部分的生成元（half twist）控制“局部交换两点”的 holonomy；
- mapping class group 部分的生成元（Dehn twist 等）控制“沿基底非平凡回路的整体扭转”的 holonomy；
- exp(iK) + Majorana 联络 $\mathcal A$ 提供了一个统一的方式来计算这两类 holonomy：它们都只是 $\pi_1(\mathcal C_n(\Sigma))$ 中不同同伦类的代表路径，区别只在于路径类型和所包围的拓扑特征。

这样，kit3 系列提出的“用配置空间和 holonomy 来理解 Dehn twist 和 half twist”的想法，就在 exp(iK) 语言下具体化为：

> 选好 BdG/Majorana 模型与 YBE‑型 R 矩阵族，视 $K(\lambda)$ 为 $so(2N)$ 联络，沿配置空间中的不同路径积分得到 $U_{\Gamma}$，从而得到 braid 群与 mapping class 群在零模 Hilbert 空间中的具体幺正表示；Ising 例子则说明了这些表示如何与已知的任意子 R/T 数据精确对上。

