# 模型综述：从 Yang–Baxter R(u) 到两能级有效 Hamiltonian（MZM 与 ABS）

作者：整理自 `ybe2.md`, `ybe22.md`, `ybe223.md`, `ybe224.md`。

## 摘要

本文把四份笔记中的核心内容整理为一条清晰的管线：从参数化的 Yang–Baxter 算子 $R(u)$ 数值导出瞬时生成元 $H(u)$；将 $H(u)$ 投影为低能两能级有效 Hamiltonian $H_{\rm eff}=d_0 I+\vec d\cdot\vec\sigma$；用轨迹几何（Bloch 球上的曲线）把 MZM 与 ABS 的物理判据统一化；并给出把 Pauli 张量基嵌入到费米子 BdG/Kitaev 链的 Jordan–Wigner 映射要点与数值实现建议。

## 1. 基本对象与定义

- 路径算子（参数 $u$）：
  $$
  R(u)=\mathcal T_u\exp\Big(-i\int_{u_0}^u H(s)\,ds\Big).
  $$
- 左生成元（左对数导数）：
  $$
  H(u)=i\,\partial_u R(u)\,R(u)^{-1}.
  $$

说明：若 $R(u)$ 含整体标量因子 $\rho(u)$，该因子只产生全局相位，可在物理判据中忽略。若 $R(u)$ 在某点不可逆或含奇异值，应当对路径重参数化或用正则化（伪逆 / 极分解）处理。

## 2. 从 $R(u)$ 到数值生成元的实现要点

- 差分与求逆：在离散网格上用中心差分（端点用单边差分）估计 $\partial_u R$，再计算 $H=i\,dR\,R^{-1}$。若 $R$ 近奇异，降级为伪逆。
- 极分解法：代替直接求逆，可先做极分解 $R=U P$（$U$ unitary），用 $U$ 计算生成元 $H_U=i\partial_u U U^\dagger$，保证厄米性。
- 数值厄米化：对数值结果做 $H\leftarrow\tfrac12(H+H^\dagger)$ 以消除残留非厄米部分。
- 本征向量相位连续性：做本征分解时按相邻重叠最大化排列本征向量并固定相位，方便累加几何相位。

## 3. 投影到两能级、Pauli 映射字典

- 选定低能子空间（常取 $\mathcal S=\{|01\rangle,|10\rangle\}$），投影得到
  $$
  H_{\rm eff}(u)=P H(u) P,
  $$
  任意 $2\times2$ 厄米矩阵可表示为
  $$
  H_{\rm eff}=d_0 I + d_x\sigma_x + d_y\sigma_y + d_z\sigma_z.
  $$
- 完整投影字典（基于标准基序与相位约定）：
  $$
  \begin{aligned}
  X\otimes X &\mapsto +\sigma_x,\\
  Y\otimes Y &\mapsto +\sigma_x,\\
  X\otimes Y &\mapsto -\sigma_y,\\
  Y\otimes X &\mapsto +\sigma_y,\\
  Z\otimes I &\mapsto +\sigma_z,\\
  I\otimes Z &\mapsto -\sigma_z,\\
  Z\otimes Z &\mapsto -I.
  \end{aligned}
  $$

- 因此对任意 Pauli 系数 $h_{\alpha\beta}(u)$，有效分量为（常用封装）：
  $$
  \boxed{\begin{aligned}
  d_x(u) &= h_{xx}(u) + h_{yy}(u) \\
  d_y(u) &= -h_{xy}(u) + h_{yx}(u) \\
  d_z(u) &= h_{zI}(u) - h_{Iz}(u) \\
  d_0(u) &= h_{II}(u) - h_{zz}(u)
  \end{aligned}}
  $$

## 4. Eight‑vertex 示例（直观案例）

- 设
  $$
  R(u)=\begin{pmatrix} \cos u & 0 & 0 & i\delta\sin u\\
  0 & \sin u & \cos u & 0\\
  0 & \cos u & \sin u & 0\\
  i\delta\sin u & 0 & 0 & \cos u\end{pmatrix}.
  $$
- 投影后得到
  $$
  H_{\rm eff}(u)=\begin{pmatrix} \delta & e^{-iu}\\ e^{iu} & -\delta\end{pmatrix}
  $$
  对应
  $$
  \vec d(u)=(\cos u,\ \sin u,\ \delta).
  $$

物理观测：当 $\delta=0$ 时 $\vec d(u)$ 为穿过原点的单位圆 —— 对应理想 MZM（几何主导）；当 $\delta\neq0$ 时轨迹抬高、绕原点失败 —— 对应 ABS（动力学主导）。

## 5. MZM 与 ABS 的统一判据（轨迹几何）

- 能谱：瞬时本征值为 $E=\pm|\vec d(u)|$。最小能隙决定绝热性与几何/动力学分离可行性。
- 路径是否包围原点为判据之一：
  $$
  \text{MZM} \Longleftrightarrow \vec d(u)\ \text{在 Bloch 球原点处有包围（或穿过）}.
  $$
- 其它判据（补充）：Pfaffian / Z2 指标、零能态的局域化与随体系尺度 $L$ 的退化行为（MZM 的分裂随 $L$ 指数小）。

## 6. 几何相 / 动力学相 的分离与 α 缩放

- 在本征基下演化方程为 $i\partial_u\tilde\psi=(D-A)\tilde\psi$，其中 $D=\mathrm{diag}(E_n)$，$A=iV^\dagger\partial_u V$ 是 Berry 连接。
- 数值实现：对每步求本征值/本征向量并按相邻重叠修正排序与相位，可累加动力学相 $\theta_n=\sum E_n\,dt$ 与几何相由重叠积给出 $\gamma_n=-\arg\prod\langle n(u_i)|n(u_{i+1})\rangle$。
- 缩放思想：令 $H\to\alpha H$，定义
  $$
  I_0=\int |\vec d(u)|\,du\approx\sum |\vec d|\,\Delta u,
  \qquad \alpha_{\rm target}=\frac{\pi}{I_0}.
  $$
  在理想情形下取 $\alpha=\alpha_{\rm target}$ 可使得动力学因子在单圈后退化为全局相位（$U_{\rm dyn}\approx\pm I$），从而 $U\approx$ 几何门 $U_{\rm geom}$。

## 7. Jordan–Wigner 映射与 BdG 嵌入

- 常用替换关系（两站点示例）：
  - $X_1X_2+Y_1Y_2\mapsto 2(c_1^\dagger c_2 + c_2^\dagger c_1)$（hopping），
  - $X_1X_2-Y_1Y_2\mapsto 2(c_1^\dagger c_2^\dagger + c_2 c_1)$（pairing），
  - $Z_j\mapsto 2n_j-1$（化学势）。
- 对应 BdG 矩阵（Nambu 基 $\Psi=(c_1,c_2,c_1^\dagger,c_2^\dagger)^T$）：
  $$
  H_{\rm BdG}=\frac12\Psi^\dagger\begin{pmatrix}
   -\tilde\mu_1 & -t & 0 & \Delta \\
   -t & -\tilde\mu_2 & -\Delta & 0 \\
   0 & -\Delta^* & \tilde\mu_1 & t \\
   \Delta^* & 0 & t & \tilde\mu_2
  \end{pmatrix}\Psi.
  $$

## 非阿贝尔性、Berry 连接、Bloch 旋转与 ABS 的详细说明

下面把关键概念与判据补充到本文中，以便把代数对象、数值实现与直观几何统一起来。

### 非阿贝尔性的来源与判据

- 定义（本征基形式）：对瞬时本征基矩阵 $V(u)=[|1(u)\rangle,|2(u)\rangle,\dots]$，定义 Berry 连接
$$
  A(u)=i\,V(u)^\dagger\partial_u V(u),
$$

这是一个矩阵值的一维连接。当作用在 $m$ 维简并子空间时，$A(u)$ 为 $m\times m$ 的非阿贝尔连接（Wilczek–Zee）。

- 非阿贝尔性的判据：若对路径上的不同 $u_1,u_2$，矩阵 $A(u_1)$ 与 $A(u_2)$ 不对易，或等价地子空间内的路径序乘积
$$
  U_{\rm geom}=\mathcal P\exp\Big(\int A(u)\,du\Big)
$$
依赖于顺序，则产生非阿贝尔几何门；必要条件通常包括
  1) 存在简并或近简并子空间（使得 $A$ 为矩阵而非标量）；
 2) 子空间内的混合矩阵元不为零（$\langle m|\partial_u H|n\rangle\neq0$）；
 3) 路径导致的 $H(u)$ 在不同 $u$ 下不两两对易（$[H(u_1),H(u_2)]\neq0$）。

- 近似表达（非对角元，微扰）：对 $m\ne n$，有
  $$
  A_{mn}(u) \approx i\frac{\langle m(u)|\partial_u H(u)|n(u)\rangle}{E_n(u)-E_m(u)},
  $$
  因此当能隙小或矩阵元大时，$A_{mn}$ 将变大，产生显著的态间耦合與非阿贝尔效应。

### Berry 相 / Berry 连接的引入与动力学‑几何分离

- 从参数化的薛定谔方程出发：$i\partial_u|\psi(u)\rangle = H(u)|\psi(u)\rangle$。写作 $|\psi\rangle=V(u)\tilde\psi$ 并左乘 $V^\dagger$ 得
  $$
  i\partial_u\tilde\psi=(D(u)-A(u))\tilde\psi,
  $$
  其中 $D(u)=\mathrm{diag}(E_n(u))$ 是动力学能谱，$A(u)$ 是几何连接。

- 因此总演化為
  $$
  U(u)=V(u)\,\mathcal P\exp\Big(-i\int_{u_0}^u (D(s)-A(s))\,ds\Big)\,V(u_0)^\dagger.
  $$

- 分离条件：只有当 $[D(s),A(s')]=0$（例如在本征基下 $A$ 为对角）或在严格绝热极限（跨本征的 $A_{mn}\ll |E_m-E_n|$）下，动力学因子 $\exp(-i\int D)$ 与几何因子 $\mathcal P\exp(\int A)$ 可近似分离为产物形式 $U\approx U_{\rm geom}U_{\rm dyn}$。

### 等价为 Bloch 旋转（两能级案例）

- 对于有效两能级模型
  $$
  H_{\rm eff}(t)=\vec d(t)\cdot\vec\sigma = |\vec d|\,\hat n\cdot\vec\sigma,
  $$
  其演化算子
  $$
  U=\mathcal T\exp\Big(-i\int H_{\rm eff}(t)\,dt\Big)=\exp\Big(-i\int |\vec d|(t)\,dt\;\hat n\cdot\vec\sigma\Big).
  $$

- 记 $\Phi=2\int |\vec d|(t)\,dt$（注意 SU(2) ↔ SO(3) 的因子 2），则
  $$
  U=\exp\Big(-i\frac{\Phi}{2}\,\hat n\cdot\vec\sigma\Big)
  $$
  对应 Bloch 球上绕轴 $\hat n$ 的旋转角为 $\Phi$。因此时间依赖的 $\vec d(t)$ 轨迹直接决定旋转轴与角度累积。

### ABS（Andreev bound state）对几何相与拓扑判据的影响

- 本质：ABS 对应于 $\vec d$ 的位移项（例如 $d_z\ne0$），使得轨迹抬高到不包围原点，从而破坏围绕原点所需的拓扑缠绕。结果是能级有显著劈裂，动力学相占主导。

- 直接后果：
  - 能谱 $E=\pm|\vec d|$ 不会在循环中接近零，动力学相 $\int E\,dt$ 大幅累积；
  - 若两个逻辑态的动力学相不同，几何门将被动力学因子混淆，导致门保真度下降；
  - ABS 通常不是由体系全局拓扑保护，其能量和局域化性质对局域扰动敏感（与 MZM 的鲁棒性不同）。

- 诊断：检查 $\vec d(u)$ 是否包围原点、最小能隙随系统尺度 $L$ 的缩放（MZM 指数小，ABS 则通常非指数）、以及局域态密度（LDOS）在边界的空间分布。

### 恢复或近似纯几何门的条件

- 缩放策略：令 $H\to\alpha H$，若能够调节 $\alpha$ 或演化时长使得动力学相满足
  $$
  \int \alpha |\vec d|\,dt = k\pi\qquad(k\in\mathbb Z),
  $$
  则动力学因子对逻辑子空间可退化为全局相（$\pm1$），从而恢复近似的纯几何门；实现上见 `tools/verify_from_R.py` 中的 `alpha_scan` 与 `alpha_target` 计算。

- 限制：真正在实验或数值上做到精确的相位抵消需要高精度控制，且在能隙极小或存在非绝热跃迁通道时方法失效。

### 在代码中的对应实现

- 参考代码：`tools/verify_from_R.py`：
  - `compute_Hs_from_R`：从 $R(u)$ 用差分计算 $H(u)$ 并做厄米化；
  - `build_H2_list_from_H4s`：投影到两能级并构造 $\vec d(t)$；
  - `adiabatic_decomposition`：按邻近重叠修正相位并累加几何/动力学相；
  - `alpha_scan`：对缩放因子做扫描以寻找使动力学因子为全局相的近似点。

  ---

  ## 附录（完整代数推导与可检验约束）

  下面把从 Pauli 系数 $h_{\alpha\beta}(u)$ 到两能级 Bloch 向量 $\mathbf d(u)$ 的闭式映射、关于 $\mathfrak{su}(2)$ 生成的代数判据、动力学相退化条件，以及几何/动力学可分离的精确代数约束写清楚，便于在数值上严格验证。

  定义三个便捷组合：
  $$
  f_x(u)=h_{XX}(u)+h_{YY}(u),\qquad f_y(u)=-h_{XY}(u)+h_{YX}(u),\qquad f_z(u)=h_{ZI}(u)-h_{IZ}(u).
  $$

  投影到子空间 $\{|01\rangle,|10\rangle\}$ 后的有效 Bloch 向量为
  $$
  \mathbf d(u)=(d_x,d_y,d_z)=(f_x(u),\;f_y(u),\;f_z(u)).
  $$

  一、生成 $\mathfrak{su}(2)$ 的代数充要条件

  - 等价表述：集合 $\mathcal D=\{\mathbf d(u)\}$ 的线性包为 $\mathbb R^3$ 当且仅当存在两点 $u_1,u_2$ 使得
    $$\mathbf d(u_1)\times\mathbf d(u_2)\neq0.$$
  - 可检验的三点判据：存在 $u_1,u_2,u_3$ 使得
    $$\det[\mathbf d(u_1),\mathbf d(u_2),\mathbf d(u_3)]\neq0.$$

  证明要点：对任意两向量 $\mathbf d_1,\mathbf d_2$，有
  $$[i\mathbf d_1\cdot\boldsymbol\sigma,\;i\mathbf d_2\cdot\boldsymbol\sigma] = -2i(\mathbf d_1\times\mathbf d_2)\cdot\boldsymbol\sigma,$$
  因此两不共线方向通过对易生成第三个方向，李代数闭包为 $\mathfrak{su}(2)$。

  二、动力学相退化（$U_{\rm dyn}$ 为全局相）的精确条件

  - 两能级动能为 $E_\pm(t)=\pm\|\mathbf d(t)\|$。定义
    $$\Theta=\int_0^T\|\mathbf d(t)\|\,dt.$$
    则 $U_{\rm dyn}\propto I$ 当且仅当 $\Theta\in\pi\mathbb Z$。若允许缩放系数 $\alpha$，则存在形式解 $\alpha=\pi k/\Theta$，但需检验该 $\alpha$ 是否在可行控制域内。

  三、几何/动力学可分离的代数条件（$[D,A]=0$）

  - 在本征基下，$D=\mathrm{diag}(E_+,E_-)$，严格可分离要求 $A_{+-}=0$（在能量基上 $A$ 为对角）。微扰近似给出
    $$A_{mn}(u)\approx i\frac{\langle m|\partial_u H|n\rangle}{E_n-E_m}\quad(m\ne n).$$
  - 等价方向条件：要求 $\mathbf d(u)$ 的方向不随 $u$ 变化（仅改变幅度），即
    $$\mathbf d(u)\times\mathbf d'(u)=0\quad\forall u.$$
    展开为 $h$ 的导数形式可写成三条恒等式（见下），在数值上按阈值检验：
$$
    \begin{align*}
    (h'_{XX}+h'_{YY})(-h_{XY}+h_{YX})-( -h'_{XY}+h'_{YX})(h_{XX}+h_{YY})=0,\\
    ( -h'_{XY}+h'_{YX})(h_{ZI}-h_{IZ})-(h'_{ZI}-h'_{IZ})(-h_{XY}+h_{YX})=0,\\
    (h'_{ZI}-h'_{IZ})(h_{XX}+h_{YY})-(h'_{XX}+h'_{YY})(h_{ZI}-h_{IZ})=0.
    \end{align*}
$$
  四、显式多项式约束（用于直接代数检查）

  - 两点叉积的三个分量（以 $h$ 写出）示例：
$$
    \begin{align*}
    c_x(u_1,u_2)&=f_y(u_1)f_z(u_2)-f_z(u_1)f_y(u_2),\\
    c_y(u_1,u_2)&=f_z(u_1)f_x(u_2)-f_x(u_1)f_z(u_2),\\
    c_z(u_1,u_2)&=f_x(u_1)f_y(u_2)-f_y(u_1)f_x(u_2).
  \end{align*}
$$
  若任一 $c_{x,y,z}\neq0$ 則兩點不共線、可生成 $\mathfrak{su}(2)$。

  五、数值检验步骤（可复制）

  1. 用 `tools/compute_Hs_from_R`（或 `tools/verify_from_R.py` 的 `compute_Hs_from_R`）生成 $H(u_i)$ 的离散样本；用投影或矩阵元计算得到 $H_{\rm eff}^{(2)}(u_i)$. 
  2. 通过 $d_k(u_i)=\tfrac12\operatorname{Tr}[H_{\rm eff}^{(2)}(u_i)\sigma_k]$ 提取 $\mathbf d(u_i)$。
  3. 计算 $\max_{i<j}\|\mathbf d(u_i)\times\mathbf d(u_j)\|$ 或 $\operatorname{rank}(D)$ 判定 su(2) 生成性。
  4. 计算離散积分 $\Theta\approx\sum_i\|\mathbf d(u_i)\|\Delta u$，检查是否存在可行 $\alpha$ 使 $\alpha\Theta\in\pi\mathbb Z$。
  5. 计算本征基 $V(u_i)$、估算 $A(u_i)=iV^\dagger\partial_u V$（中心差分），计算 $C(u_i)=\|[E,A]\|/(\|E\|\|A\|)$ 作分离性指标。

  以上公式与数值步骤可以直接映射到仓库脚本（`tools/verify_from_R.py`、`tools/diagnose_decomposition.py`、`tools/majorana_lift.py`），并可将对应的数值检验结果归入本摘要作为“数值示例”。


