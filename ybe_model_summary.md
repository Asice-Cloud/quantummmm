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

物理观测：在此模型下，当 $\delta=0$ 时 $\vec d(u)$ 为穿过原点的单位圆 —— 对应理想 MZM（几何主导）；当 $\delta\neq0$ 时轨迹抬高、绕原点失败 —— 对应 ABS（动力学主导）。“x-y平面 + z偏移”

## 5. MZM 与 ABS 的统一判据（轨迹几何）

- 能谱：瞬时本征值为 $E=\pm|\vec d(u)|$。最小能隙决定绝热性与几何/动力学分离可行性。$\sigma_z$会改变对角元素，因此决定能谱打开
- 路径是否包围原点为判据之一：MZM 对应 $\vec d(u)$ 在 Bloch 球原点附近形成包围；理想极限下可以近似看成穿过原点，但严格 braid 计算应理解为“绕开闭合点而保持能隙打开”。
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

**证明与示例**

- 证明（等价性）: 令
$$
  f_x=h_{XX}+h_{YY},\quad f_y=-h_{XY}+h_{YX},\quad f_z=h_{ZI}-h_{IZ}.
$$
  直接计算向量叉积的导数分量得
$$
f\times f'=(f_y f'_z-f_z f'_y,\;f_z f'_x-f_x f'_z,\;f_x f'_y-f_y f'_x).
$$
  将 $f$ 的定义代回，上述三标量分量恰好与上面的三条标量恒等式等价（至一侧乘 $-1$ 问题不影响零等式）。因此原方程组等价于
$$
  f(u)\times f'(u)=0,
$$
  即 $f$ 与 $f'$ 在每点共线。

-  通解与构造：在任一连通区间若 $f\not\equiv0$，存在标量函数 $\alpha(u)$ 使
  $$
  f'(u)=\alpha(u)f(u).
  $$
  积分得
  $$
  f(u)=s(u)\,v,\qquad s(u)=\exp\Big(\int^u\alpha(s)\,ds\Big),
  $$
  其中 $v\in\mathbb R^3$ 为常向量。对于 $h$ 的分量我们可引入任意 "gauge" 函数 $g_1(u),g_2(u),g_3(u)$，得到一般表示

$$
  \begin{align*}
  h_{XX}&=\tfrac12 s(u)v_x+g_1(u),\\
  h_{YY}&=\tfrac12 s(u)v_x-g_1(u),\\
  h_{XY}&=-\tfrac12 s(u)v_y+g_2(u),\\
  h_{YX}&=\tfrac12 s(u)v_y+g_2(u),\\
  h_{ZI}&=\tfrac12 s(u)v_z+g_3(u),\\
  h_{IZ}&=-\tfrac12 s(u)v_z+g_3(u).
  \end{align*}
$$

 这正是此前在正文中以 $f_x,f_y,f_z$ 给出的分量结构的显式化。

-  特殊与奇异情形：若 $f\equiv0$（即 $s\equiv0$），方程恒成立，对应约束为 $h_{YY}=-h_{XX}$、$h_{YX}=h_{XY}$、$h_{IZ}=h_{ZI}$，其余自由函数任取；若 $f$ 在孤立点为零，则可在不含零点的分段区间上按上述形式构造并在零点处作相容性处理。

-  具体示例（可直接验证）: 取 $v=(1,2,3)$, $s(u)=e^{2u}$, $g_1(u)=\sin u$, $g_2(u)=0$, $g_3(u)=\cos u$，则
$$
  \begin{align*}
  h_{XX}(u)&=\tfrac12 e^{2u}\cdot1 +\sin u,\\
  h_{YY}(u)&=\tfrac12 e^{2u}\cdot1 -\sin u,\\
  h_{XY}(u)&=-\tfrac12 e^{2u}\cdot2,\\
  h_{YX}(u)&=\tfrac12 e^{2u}\cdot2,\\
  h_{ZI}(u)&=\tfrac12 e^{2u}\cdot3 +\cos u,\\
  h_{IZ}(u)&=-\tfrac12 e^{2u}\cdot3 +\cos u.
  \end{align*}
$$

  简单的符号验证脚本已加入仓库：[tools/solve_h_constraints.py](tools/solve_h_constraints.py)，可以用来做符号化简与数值检查（脚本示例中已运行并验证上面示例令三式为零）。

## 五、Majorana 二次项到 Pauli 张量积的严格对应

下面把“杂化操作对应哪个 Pauli 张量积、理想 braid 对应什么、XX/YY/XY/YX 是否就是理想操作”这几个容易混淆的问题，严格写成统一的代数对应。

### 1. 统一的 JW 约定

对链上第 $j$ 个费米模，取标准 Jordan--Wigner 约定
$$
\gamma_{2j-1}=\Big(\prod_{m<j}\sigma_m^z\Big)\sigma_j^x,\qquad
\gamma_{2j}=\Big(\prod_{m<j}\sigma_m^z\Big)\sigma_j^y.
$$

这意味着 Majorana 二次项首先对应为“端点 Pauli + 中间 $Z$ 串”的 Pauli 串，而不是一开始就等于单个两体张量积。

### 2. 同一费米模内的杂化项

同一模上的局域杂化是
$$
H_{\rm hyb}=i\,\varepsilon_j\gamma_{2j-1}\gamma_{2j}.
$$
由 JW 约定直接得到
$$
i\gamma_{2j-1}\gamma_{2j}=-\sigma_j^z.
$$
因此杂化在逻辑上对应的是 $Z$ 轴能级分裂，而不是 XX、YY 或 XY。

### 3. 不同站点的四类 Majorana 二次项

若 $j<k$，则有
$$
\gamma_{2j-1}\gamma_{2k-1}
=
\sigma_j^x\Big(\prod_{m=j}^{k-1}\sigma_m^z\Big)\sigma_k^x,
$$
$$
\gamma_{2j}\gamma_{2k}
=
\sigma_j^y\Big(\prod_{m=j}^{k-1}\sigma_m^z\Big)\sigma_k^y,
$$
$$
\gamma_{2j-1}\gamma_{2k}
=
\sigma_j^x\Big(\prod_{m=j}^{k-1}\sigma_m^z\Big)\sigma_k^y,
$$
$$
\gamma_{2j}\gamma_{2k-1}
=
\sigma_j^y\Big(\prod_{m=j}^{k-1}\sigma_m^z\Big)\sigma_k^x.
$$

所以严格地说，Majorana 二次项对应的是一条 Pauli 串；只有在最近邻或在特定投影后，才会退化成纯两体的 $\sigma^\alpha\otimes\sigma^\beta$。

### 4. 最近邻时的纯两体形式

在最近邻情况下，串会简化成纯两体 Pauli。典型例子是
$$
i\gamma_{2j-1}\gamma_{2j+1}=\sigma_j^y\sigma_{j+1}^x,
$$
$$
i\gamma_{2j-1}\gamma_{2j+2}=\sigma_j^y\sigma_{j+1}^y,
$$
$$
i\gamma_{2j}\gamma_{2j+1}=-\sigma_j^x\sigma_{j+1}^x,
$$
$$
i\gamma_{2j}\gamma_{2j+2}=-\sigma_j^x\sigma_{j+1}^y.
$$

因此 XX、YY、XY、YX 这些形式只是 Majorana 二次项在 JW 展开后的具体表示，不是“理想”的唯一类别。

### 5. tetron 逻辑子空间中的 Pauli 对照

在 tetron 的偶宇称逻辑子空间中，常用生成元取为
$$
\Sigma_x=i\gamma_2\gamma_3,\qquad
\Sigma_y=i\gamma_1\gamma_3,\qquad
\Sigma_z=i\gamma_1\gamma_2.
$$
它们满足
$$
\Sigma_a\Sigma_b=\delta_{ab}I+i\epsilon_{abc}\Sigma_c,
$$
因此存在固定逻辑基变换 $W$ 使得
$$
W\Sigma_aW^\dagger=\sigma_a\qquad(a=x,y,z).
$$

这一步说明：在有效单比特上，Majorana 双线性本身就是 Pauli 轴，但具体是哪个轴取决于你选取的逻辑编码。

### 6. 理想 braid 的严格对应

Majorana 交换算符定义为
$$
B(\gamma_i,\gamma_j)=\exp\!\left(\frac{\pi}{4}\gamma_i\gamma_j\right).
$$
在逻辑比特上，它等价于一个固定角度的 Pauli 旋转
$$
B\sim \exp\!\left(-i\frac{\pi}{4}\sigma_x\right),\qquad
\exp\!\left(-i\frac{\pi}{4}\sigma_y\right),\qquad
\exp\!\left(-i\frac{\pi}{4}\sigma_z\right),
$$
具体轴由 braid 的 Majorana 对和所选编码共同决定。

因此，理想 braid 不是“XX/YY/XY 的某一个”，而是由某个 Majorana 二次项生成的固定角度幺正门；XX/YY/XY/YX 只是其在 JW/Pauli 基下的具体展开。

### 7. 一句话总结

- 杂化对应：
  $$
  i\gamma_{2j-1}\gamma_{2j}=-\sigma_j^z.
  $$
- braid 对应：
  $$
  \exp\!\left(\frac{\pi}{4}\gamma_a\gamma_b\right)
  \leftrightarrow
  \exp\!\left(-i\frac{\pi}{4}\sigma_{x,y,z}\right).
  $$
- XX/YY/XY/YX：只是某些 Majorana 二次项的 Pauli 张量积展开，不是唯一的“理想结构”；Z 型杂化同样是理想生成元，只是它对应局域能级分裂。

### 8. 完整对照表：Majorana 双线性、Pauli 张量积、编码空间 sigma

下面给出最直接的严格对应。这里分成三层：

- 第 1 层：微观 Majorana 双线性 $i\gamma_a\gamma_b$；
- 第 2 层：JW 展开后的 Pauli 串或两体张量积；
- 第 3 层：投影到 tetron 逻辑子空间后的单比特 $\sigma_x,\sigma_y,\sigma_z$。

#### 8.1 杂化项

同一费米模的杂化：
$$
i\gamma_{2j-1}\gamma_{2j}=-\sigma_j^z.
$$

因此在编码空间中，杂化对应的就是 $\sigma_z$ 轴分裂。

#### 8.2 braid 生成元

在 tetron 逻辑子空间里，常用三组生成元取为
$$
\Sigma_x=i\gamma_2\gamma_3,\qquad
\Sigma_y=i\gamma_1\gamma_3,\qquad
\Sigma_z=i\gamma_1\gamma_2,
$$
并选取固定逻辑基使
$$
W\Sigma_aW^\dagger=\sigma_a\qquad(a=x,y,z).
$$

因此 braid 算符
$$
B(\gamma_i,\gamma_j)=\exp\!\left(\frac{\pi}{4}\gamma_i\gamma_j\right)
$$
在编码空间中对应为
$$
\exp\!\left(-i\frac{\pi}{4}\sigma_x\right),\qquad
\exp\!\left(-i\frac{\pi}{4}\sigma_y\right),\qquad
\exp\!\left(-i\frac{\pi}{4}\sigma_z\right),
$$
具体是哪一个轴由你 braid 的那一对 Majorana 和所选编码共同决定。

#### 8.3 最近邻两体 Pauli 的完整对应

在最近邻情况下，Majorana 二次项退化为纯两体 Pauli：
$$
i\gamma_{2j-1}\gamma_{2j+1}=\sigma_j^y\sigma_{j+1}^x,
$$
$$
i\gamma_{2j-1}\gamma_{2j+2}=\sigma_j^y\sigma_{j+1}^y,
$$
$$
i\gamma_{2j}\gamma_{2j+1}=-\sigma_j^x\sigma_{j+1}^x,
$$
$$
i\gamma_{2j}\gamma_{2j+2}=-\sigma_j^x\sigma_{j+1}^y.
$$

因此下面这几类两体项只是同一 Majorana 生成元在 JW/Pauli 基下的具体展开：

| Majorana 二次项 | JW / Pauli 张量积 | 编码空间里的有效轴 |
| --- | --- | --- |
| $i\gamma_{2j-1}\gamma_{2j}$ | $-\sigma_j^z$ | $\sigma_z$ |
| $i\gamma_{2j-1}\gamma_{2j+1}$ | $\sigma_j^y\sigma_{j+1}^x$ | 投影到逻辑子空间后对应单比特生成元 |
| $i\gamma_{2j-1}\gamma_{2j+2}$ | $\sigma_j^y\sigma_{j+1}^y$ | 投影到逻辑子空间后对应单比特生成元 |
| $i\gamma_{2j}\gamma_{2j+1}$ | $-\sigma_j^x\sigma_{j+1}^x$ | 投影到逻辑子空间后对应单比特生成元 |
| $i\gamma_{2j}\gamma_{2j+2}$ | $-\sigma_j^x\sigma_{j+1}^y$ | 投影到逻辑子空间后对应单比特生成元 |

#### 8.4 结论

所以，严格对应关系不是“XX/YY/XY/YX 直接等于 braid”，而是：

1. Majorana 双线性先经 JW 变成 Pauli 串或两体张量积；
2. 再投影到逻辑子空间，变成单比特的 $\sigma_x,\sigma_y,\sigma_z$；
3. braid 则是这些生成元的固定角度指数门。

换句话说：

杂化对应 $\sigma_z$；braid 对应 $\exp\!\left(-i\frac{\pi}{4}\sigma_{x,y,z}\right)$；XX/YY/XY/YX 对应微观展开形式。

### 9. 证明：为什么 $\sigma_x,\sigma_y,\sigma_z$ 都可以对应理想 braid

这里要证明的不是“同一个 braid 同时等于三个轴”，而是：在四 Majorana 的固定逻辑子空间里，存在三类独立的 braid 生成元，它们分别在逻辑比特上实现 $\sigma_x,\sigma_y,\sigma_z$ 方向的 $\pi/2$ 旋转。

#### 定理

设四个 Majorana 满足 Clifford 代数
$$
\{\gamma_a,\gamma_b\}=2\delta_{ab},
$$
并把偶宇称子空间作为逻辑子空间。定义三组双线性算符
$$
\Sigma_x=i\gamma_2\gamma_3,\qquad
\Sigma_y=i\gamma_1\gamma_3,\qquad
\Sigma_z=i\gamma_1\gamma_2.
$$
则：

1. $\Sigma_x,\Sigma_y,\Sigma_z$ 都是厄米且满足 Pauli 代数；
2. 存在固定酉变换 $W$，使得 $W\Sigma_aW^\dagger=\sigma_a$（$a=x,y,z$）；
3. 对应的 braid 算符
$$
B_{23}=\exp\!\left(\frac{\pi}{4}\gamma_2\gamma_3\right),\qquad
B_{13}=\exp\!\left(\frac{\pi}{4}\gamma_1\gamma_3\right),\qquad
B_{12}=\exp\!\left(\frac{\pi}{4}\gamma_1\gamma_2\right)
$$
分别等价于逻辑比特上的
$$
\exp\!\left(-i\frac{\pi}{4}\sigma_x\right),\qquad
\exp\!\left(-i\frac{\pi}{4}\sigma_y\right),\qquad
\exp\!\left(-i\frac{\pi}{4}\sigma_z\right),
$$
即三种理想 braid 轴。

#### 证明

第一步，验证 Pauli 代数。由 Majorana 反对易关系可得
$$
\Sigma_a^\dagger=\Sigma_a,\qquad \Sigma_a^2=I,
$$
并且
$$
\Sigma_a\Sigma_b=\delta_{ab}I+i\epsilon_{abc}\Sigma_c.
$$
这正是三维 Pauli 代数，因此三者生成一个 $\mathfrak{su}(2)$ 表示。

第二步，证明与标准 Pauli 矩阵等价。任意一组满足上式的三矩阵都与 $\sigma_x,\sigma_y,\sigma_z$ 酉等价，因此存在固定逻辑基变换 $W$ 使
$$
W\Sigma_aW^\dagger=\sigma_a.
$$
这说明在逻辑子空间里，$\Sigma_x,\Sigma_y,\Sigma_z$ 就是三条正交旋转轴。

第三步，证明 braid 生成元的固定角度性质。对任意 $a\neq b$，对应 Majorana 交换算符定义为
$$
B_{ab}=\exp\!\left(\frac{\pi}{4}\gamma_a\gamma_b\right)
=\exp\!\left(-i\frac{\pi}{4}\,i\gamma_a\gamma_b\right).
$$
若取 $i\gamma_a\gamma_b=\Sigma_c$，则
$$
B_{ab}=\exp\!\left(-i\frac{\pi}{4}\Sigma_c\right).
$$
在逻辑基下再做一次酉变换，就得到
$$
W B_{ab} W^\dagger = \exp\!\left(-i\frac{\pi}{4}\sigma_c\right).
$$
这正是 Bloch 球上绕 $c$ 轴旋转 $\pi/2$ 的理想 braid 门。

第四步，给出三轴对应关系。取上面这组标准编码时，有
$$
i\gamma_2\gamma_3=\Sigma_x\mapsto\sigma_x,
$$
$$
i\gamma_1\gamma_3=\Sigma_y\mapsto\sigma_y,
$$
$$
i\gamma_1\gamma_2=\Sigma_z\mapsto\sigma_z.
$$
因此三条独立 braid 轴分别实现 $\sigma_x,\sigma_y,\sigma_z$ 方向的理想交换门。

#### 结论

所以，“$\sigma_x,\sigma_y,\sigma_z$ 都能对应理想 braid”是对的，但含义必须理解为：

- 不是同一个 braid 同时等于三者；
- 而是存在三种独立的 Majorana 交换生成元，它们在固定逻辑子空间里分别对应 $\sigma_x,\sigma_y,\sigma_z$；
- 每一种都是理想 braid，因为它们都是固定角度 $\pi/2$ 的纯几何旋转。

### 10. 从路径 $u(t)$ 到哈密顿量，再到杂化 / braid 的判别

这里把“选取路径 $u(t)$”这件事说成一个严格的三步链：

1. 先给出参数路径 $R(u)$ 或器件路径 $p(u)$；
2. 再由左对数导数得到瞬时生成元；
3. 最后看这条路径在逻辑子空间里对应的是局域杂化还是 braid 交换。

#### 10.1 由路径得到哈密顿量

若给定的是谱参数路径 $R(u)$，则定义
$$
H(u)=i\,\partial_u R(u)\,R(u)^{-1}.
$$
如果再把参数 $u$ 用物理时间参数化为 $u=u(t)$，那么总演化算符可以写成
$$
U(t)=R(u(t)),
$$
于是物理时间上的有效哈密顿量为
$$
H(t)=i\,\dot U(t)U(t)^{-1}=\dot u(t)\,H(u(t)).
$$
这说明：路径决定哈密顿量，时间重参数化只会给出一个速度因子 $\dot u(t)$。

如果给定的是论文中的器件路径 $p(t)$，则先得到微观 Hamiltonian $H_M(t)=H_M[p(t)]$，再投影到低能子空间
$$
H_{\rm eff}(t)=P\,H_M(t)\,P=d_0(t)I+\vec d(t)\cdot\vec\sigma.
$$

#### 10.2 怎么判断是杂化

杂化的特征是：哈密顿量对应某一对 Majorana 的局域双线性，或者等价地，在编码空间里是一个固定 Pauli 轴上的能级劈裂。

最标准的情况是
$$
H_{\rm hyb}=i\varepsilon\gamma_{2j-1}\gamma_{2j}\quad\Longrightarrow\quad i\gamma_{2j-1}\gamma_{2j}=-\sigma_j^z.
$$
因此在我们的 tetron 编码里，杂化最直接地表现为 $\sigma_z$ 方向的局域分裂。

更一般地说，若在一段演化中
$$
[H(t_1),H(t_2)]=0\quad\forall\,t_1,t_2,
$$
而且 $H(t)$ 始终只沿一个固定逻辑轴变化，那么这段演化就是一个单轴杂化段；它不依赖 Majorana 交换顺序，只是同一能级差的积累。

#### 10.3 怎么判断是理想 braid

理想 braid 不是“某一个瞬时 Hamiltonian 长什么样”，而是整条路径生成的总演化是否等于 braid 算符（差一个全局相位）。

在 Majorana 语言里，理想 braid 的定义是
$$
U_{\rm braid}=\exp\!\left(\frac{\pi}{4}\gamma_a\gamma_b\right),
$$
等价地，在逻辑子空间里是
$$
U_{\rm braid}\sim \exp\!\left(-i\frac{\pi}{4}\sigma_c\right).
$$

因此判别理想 braid 的最直接标准是：

1. 演化过程中能隙始终打开，路径不穿过闭合点；
2. 末态幺正门与 braid 门一致（允许一个全局相位）；
3. 若把动力学相调到全局相，则剩下的几何部分正好是 $\pi/2$ 的交换旋转。

对一个两能级有效模型来说，若
$$
H_{\rm eff}(t)=d_0(t)I+\vec d(t)\cdot\vec\sigma,
$$
并且 $\vec d(t)$ 在参数空间中围绕原点作闭合绕行，那么该路径的几何部分就会产生 braid 型的 $SU(2)$ 旋转；若同时满足动力学相可被调成全局相，则就是“理想 braid”。

#### 10.4 进入我们模型后的对应

在 eight-vertex / 两能级模型里，我们已经有
$$
H_{\rm eff}(u)= -\delta I+\cos u\,\sigma_x+\sin u\,\sigma_y.
$$
因此

- $\delta=0$ 时，$\vec d(u)$ 落在赤道平面并围绕原点，属于 MZM 主导的几何轨迹；
- $\delta\neq0$ 时，$\vec d(u)$ 被抬离原点，动力学相增强，表现为 ABS 偏移。

若把 $u=u(t)$ 看成实际控制路径，那么
$$
H_{\rm eff}(t)=\dot u(t)\big[-\delta I+\cos u(t)\,\sigma_x+\sin u(t)\,\sigma_y\big].
$$
这时：

- 当 $u(t)$ 只在固定区间停留时，得到单轴的杂化式演化；
- 当 $u(t)$ 绕一圈并且动力学相可被压成全局相时，得到 braid 型几何旋转。

#### 10.5 结论

所以，完整链条是
$$
u(t)\;\Longrightarrow\;R(u(t))\;\Longrightarrow\;H(t)=\dot u(t)H(u(t))\;\Longrightarrow\;H_{\rm eff}(t)=d_0I+\vec d\cdot\vec\sigma\;\Longrightarrow\;U.
$$

其中：

- **杂化** = 固定 Majorana 双线性 / 固定逻辑轴上的能级分裂；
- **理想 braid** = 路径生成的总演化等于 braid 算符，且动力学相可化为全局相；
- **Bloch 轨迹** = $\vec d(t)$ 在编码空间中的几何投影。
