# Pauli 张量基、Hamiltonian 与非阿贝尔结构的关系（系统总结）

 

## 一、核心结论

> **Pauli 张量基是“坐标系”，不是“非阿贝尔性的来源”。**
> 非阿贝尔性来源于**本征态随参数变化的几何结构（Berry connection）**，
> 而 Pauli 基提供了一个**可计算、可控制的表达框架**。

 

## 二、基本结构链条

从你的模型出发：

$$
R(u)=\mathcal{T}\exp\Big(-i\int_0^u H(s)ds\Big)
$$

定义：

$$
H(u)=\sum_{\alpha,\beta} h_{\alpha\beta}(u)\sigma^\alpha\otimes\sigma^\beta
$$

 

### 层级结构：

1. **Pauli 基（参数空间）**

随着 u 变化
$$
h_{\alpha\beta}(u)
$$
定义了一条路径
这条路径就是：在算符空间里的“演化轨道”

2. **Hamiltonian（生成元）**
$$
H(u)
$$

3. **本征问题**
$$
   H(u)|n(u)\rangle = E_n(u)|n(u)\rangle
$$

4. **Berry connection（几何结构）**

非阿贝尔性不在 $h_αβ$本身,而在本征态随 u 的变化结构

$$
A_{mn}(u)= i\langle m(u)|\partial_u n(u)\rangle
$$

这意味着：

同一组 $h_αβ(u)$
- 如果本征态不混合 → 阿贝尔
- 如果子空间内发生旋转 → 非阿贝尔


5. **非阿贝尔演化**
$$
U_{\text{geom}}=\mathcal{P}\exp\left(\int A(u)du\right)
$$



## 三、Hamiltonian 与 Berry connection 的严格关系

定义本征态矩阵：

$$
V(u) = (|1(u)\rangle, |2(u)\rangle, \dots)\\
V(u)V\dagger = I
$$

瞬时本征问题:
$$
H(u)|n(u)\rangle = E_n(u)n(u)\rangle
$$

通过:
$$
H(u)=i\partial_uR(u)R^{-1}(u)
$$

则有瞬时本征分解：

$$
H(u)= V(u)E(u)V^\dagger(u),\qquad E(u)=\mathrm{diag}(E_1(u),E_2(u),\dots).
$$

把任意态写成本征基的展开：
$$
|\psi(u)\rangle = V(u)\,\tilde\psi(u),
$$
并代入参量化的薛定谔方程
$$
i\partial_u |\psi(u)\rangle = H(u)\,|\psi(u)\rangle.
$$
左乘 $V(u)^\dagger$ 并利用 $V^\dagger V = I$，得
$$
i\partial_u \tilde\psi + iV^\dagger(\partial_u V)\,\tilde\psi = E\,\tilde\psi.
$$
移项可写为本征基下的演化方程：
$$
i\partial_u\tilde\psi = \bigl(E - A\bigr)\tilde\psi,\qquad A(u)=iV^\dagger(u)\partial_u V(u),
$$
其中 $A(u)$ 为（非阿贝尔）Berry 连接。由此可见，在本征基下演化生成元为 $E-A$，而不能简单将 $H$ 写成 “$V E V^\dagger - i(\partial V)V^\dagger$”。



### 在简并子空间 $\mathcal{S}$ 上投影：

$$
P(u)=\sum_{a \in S}|a(u)\rangle \langle a(u)| \\
H_{\mathcal{S}}(u)=P(u)H(u)P(u)
$$

得到：

$$
\boxed{
H_{\mathcal{S}}(u)=E_{\mathcal{S}}(u)-A(u)
}
$$

其中：

* $(E_{\mathcal{S}})$：能谱（对角）
* $(A(u))$：Berry connection

 

## 四、Berry connection 的来源（关键公式）
pauli基如何参与非阿贝尔性：Pauli 基的作用是决定本征态如何变化
非阿贝尔Berry Connection:
$$
A(u)=iP(u)\partial_uP(u)\\
A_{ab}(u)=i\langle a(u)|\partial_u |b(u)\rangle
$$

利用微扰论：

$$
A_{mn}(u)=i\frac{\langle m|\partial_u H|n\rangle}{E_n-E_m}
\quad (m\neq n)
$$

 

代入 Pauli 展开：

$$
\partial_u H(u)=\sum_{\alpha\beta} \dot h_{\alpha\beta}(u)\sigma^\alpha\otimes\sigma^\beta
$$

得到：

$$
\boxed{
A_{mn}(u)=i\sum_{\alpha\beta}
\dot h_{\alpha\beta}(u)
\frac{\langle m|\sigma^\alpha\otimes\sigma^\beta|n\rangle}{E_n-E_m}
}
$$

 

### 重要结论：

> **Berry connection 完全由 Pauli 系数的变化决定**

 

## 五、非阿贝尔性的产生条件

要出现非阿贝尔结构，必须满足：

 

### 条件 1：存在（准）简并

$$
E_a(u)\approx E_b(u)
$$
否则
$$
A_{ab} \approx \frac{1}{E_a-E_b} \rightarrow 0
$$

 

### 条件 2：子空间内有 mixing

$$
\langle a|\sigma^\alpha\otimes\sigma^\beta|b\rangle \neq 0
$$
如果所有项都对角：

- 不发生 mixing
- 只有 Abelian phase


### 条件 3：路径非对易

$$
[H(u_1), H(u_2)] \neq 0
$$
否则：
- 演化可交换
- 不可能非阿贝尔


## 六、直接从 Pauli 基判断非阿贝尔性

 

### 判据 1：生成元不对易

$$
[H(u_1), H(u_2)] \neq 0
$$
因此必然产生 path-ordering
⇒ 非阿贝尔结构可能存在




### 判据 2：子空间内非对角

$$
H_{\mathcal{S}}(u) \neq \text{对角}
$$

如果：H_S(u)不是对角的
⇒ 非阿贝尔 mixing

 

### 判据 3：YBE 偏离

$$
\Delta = R_{12}R_{23}R_{12} - R_{23}R_{12}R_{23}
$$

Pauli 系数的变化 → 决定 Δ

 

## 七、Pauli 基的真正作用

$$
H(u) \rightarrow {|n(u)\rangle} \rightarrow A(u)
$$

### ✔ 作用 1：完备参数化

任意两比特算符都可展开

 

### ✔ 作用 2：物理可读性

| 项      | 物理意义              |
|    |      -- |
| XX, YY | hopping / pairing |
| Z      | 化学势               |
| ZZ     | 相互作用              |

 

### ✔ 作用 3：控制路径

$$
h_{\alpha\beta}(u)
$$

定义演化轨迹

 

## 八、局限性

当出现：

* 强相互作用
* 非局域项（JW string）
* 多体效应

 

则：

> Pauli 展开仍成立，但不再等价于 BdG / Majorana 结构

 

## 九、统一理解

 

| 层级 | 对象               | 作用    |
| -- |      - |  -- |
| 1  | Pauli 系数         | 坐标    |
| 2  | Hamiltonian      | 生成元   |
| 3  | 本征态              | 动力学结构 |
| 4  | Berry connection | 几何结构  |
| 5  | braid            | 非阿贝尔性 |

 

## 十、最终总结

> **Pauli 张量基不是用来定义非阿贝尔性，而是提供一个可计算的坐标框架，使得非阿贝尔结构可以通过 Hamiltonian → 本征态 → Berry connection 的路径被精确提取。**

 

## 十一、建议的研究方向

构造最小模型：

$$
H(u)= a(u),XX + b(u),YY + c(u),XY
$$

步骤：

1. 保持简并
2. 设计路径 (a(u), b(u), c(u))
3. 计算 Berry connection
4. 得到 braid

 

> 该模型可以直接展示：
> **Pauli 系数如何驱动非阿贝尔结构的产生**

## 附录：严格推导与代数判据

下面给出从 $H_{\rm spin}$ 的 Pauli 展开到两能级有效 Hamiltonian $H_{\rm eff}=\mathbf d\cdot\boldsymbol\sigma$ 的完整代数推导，并列出可用于数值检验的充要代数条件。

假设总体自旋 Hamiltonian 在 Pauli 张量基下展开为
$$
H_{\rm spin}=\sum_{\alpha,\beta\in\{I,X,Y,Z\}} h_{\alpha\beta}(u)\;\sigma^{\alpha}\otimes\sigma^{\beta},
$$
系数有唯一的正交投影表示
$$
h_{\alpha\beta}(u)=\frac{1}{4}\operatorname{Tr}\big[(\sigma^{\alpha}\otimes\sigma^{\beta})\,H_{\rm spin}(u)\big],
$$
利用 $\operatorname{Tr}(\sigma^a\sigma^b)=2\delta_{ab}$。

取逻辑子空间为 $\{|01\rangle,|10\rangle\}$（按基序 $|00\rangle,|01\rangle,|10\rangle,|11\rangle$）。投影后得到 2×2 有效矩阵
$$
H_{\rm eff}^{(2)}(u)_{ij}=\langle\psi_i|H_{\rm spin}(u)|\psi_j\rangle,\qquad i,j\in\{1,2\},
$$
其中 $|\psi_1\rangle=|01\rangle,\;|\psi_2\rangle=|10\rangle$。由于张量算子的矩阵元可分解为单次的乘积（例如 $\langle0|A|0\rangle\langle1|B|1\rangle$），直接计算得：

令
$$
H_{\rm eff}^{(2)}(u)=\begin{pmatrix}a(u) & b(u)\\[4pt] \overline{b(u)} & d(u)\end{pmatrix},
$$
则
$$
\begin{align*}
a(u)&=h_{II}-h_{IZ}+h_{ZI}-h_{ZZ},\\
d(u)&=h_{II}+h_{IZ}-h_{ZI}-h_{ZZ},\\
b(u)&=h_{XX}+h_{YY}+i\bigl(h_{XY}-h_{YX}\bigr),
\end{align*}
$$
其中省略了显式依赖 $u$ 以简化书写。

任取 2×2 厄米矩阵都可用 Pauli 展开：
$$
H_{\rm eff}^{(2)}=d_0I + d_x\sigma_x + d_y\sigma_y + d_z\sigma_z.
$$
用线性代数直接解出 $d$ 的分量：
$$
\begin{align*}
d_0(u) &= \tfrac{1}{2}\bigl(a(u)+d(u)\bigr) = h_{II}(u)-h_{ZZ}(u),\\
d_x(u) &= \operatorname{Re} b(u) = h_{XX}(u)+h_{YY}(u),\\
d_y(u) &= -\operatorname{Im} b(u) = -h_{XY}(u)+h_{YX}(u),\\
d_z(u) &= \tfrac{1}{2}\bigl(a(u)-d(u)\bigr) = h_{ZI}(u)-h_{IZ}(u).
\end{align*}
$$

因此关于 Pauli 系数的闭式映射为
$$
\boxed{\;\mathbf d(u)=\bigl\langle\,h_{XX}+h_{YY},\; -h_{XY}+h_{YX},\; h_{ZI}-h_{IZ}\,\bigr\rangle\; }.
$$

---

### su(2) 生成的代数判据（必要且充分）

把两能级 Hamiltonian 写成 $H(u)=\mathbf d(u)\cdot\boldsymbol\sigma$。令集合 $\mathcal D=\{\mathbf d(u):u\in\Gamma\}$ 为参数轨迹的像。

命题（等价陈述）：李代数 $\operatorname{Lie}\{i\mathbf d(u)\cdot\boldsymbol\sigma\}_{u\in\Gamma}=\mathfrak{su}(2)$ 当且仅当 $\operatorname{span}_{\mathbb R}\mathcal D=\mathbb R^3$。

简明可检验条件：存在 $u_1,u_2\in\Gamma$ 使得
$$
\mathbf d(u_1)\times\mathbf d(u_2)\neq 0.
$$
更强的三点判据（等价）：存在 $u_1,u_2,u_3$ 使得
$$
\det\bigl[\mathbf d(u_1),\mathbf d(u_2),\mathbf d(u_3)\bigr]\;=\;\mathbf d(u_1)\cdot\bigl(\mathbf d(u_2)\times\mathbf d(u_3)\bigr)\neq 0.
$$

证明要点：对任意两个方向 $\mathbf d_1,\mathbf d_2$，有 Lie 对易子
$$
[i\mathbf d_1\cdot\boldsymbol\sigma,\; i\mathbf d_2\cdot\boldsymbol\sigma]=-2i(\mathbf d_1\times\mathbf d_2)\cdot\boldsymbol\sigma.
$$
因此只要存在两不共线向量，便可通过一次对易得到第三个方向；三线性独立则生成全 $\mathfrak{su}(2)$。

实用数值检验（离散采样）:
\begin{itemize}
\item 计算采样矩阵 $D=[\mathbf d(u_1),\dots,\mathbf d(u_N)]^T$，检查 $\operatorname{rank}(D)$；若 $\operatorname{rank}(D)=3$ 则生成 $\mathfrak{su}(2)$。
\item 或直接计算最大叉积模长
$$\max_{i<j}\|\mathbf d(u_i)\times\mathbf d(u_j)\|>\varepsilon$$
为判据（数值阈值 $\varepsilon$ 由机器精度与物理尺度决定）。
\end{itemize}

---

### 动力学相退化（$U_{\rm dyn}$ 为全局相）的精确条件

令瞬时本征能量为 $E_{\pm}(u)=\pm\|\mathbf d(u)\|$（两能级情形）。动力学相为
$$
	heta_{\pm}=\int_0^T E_{\pm}(t)\,dt=\pm\int_0^T\|\mathbf d(t)\|\,dt.
$$
则 $U_{\rm dyn}\propto I$（模整体相）当且仅当
$$
	heta_+ - \theta_- = 2\int_0^T\|\mathbf d(t)\|\,dt \in 2\pi\mathbb Z
$$
即等价于
$$
\boxed{\;\int_0^T\|\mathbf d(t)\|\,dt \in \pi\mathbb Z\; }.
$$
数值上常用縮放參數 $\alpha$（即把 $H\mapsto\alpha H$）來嘗試達成上述等式：存在 \(\alpha\) 使得 \(\alpha\int\|\mathbf d\|dt=\pi k\)。但物理可行性取決於 \(\alpha\) 的大小與實驗可控域。

---

### 幾何/動力學可分離的對易條件

在本征基下，有有效演化生成元 $E(u)-A(u)$，若欲嚴格分解為幾何與動力學乘積需要
$$
[E(u),A(s)]=0\quad\text{對所有 }u,s
$$
尤其的弱化指標（數值可檢驗）為歸一化對易度量
$$
C(u)=\frac{\|[E(u),A(u)]\|_F}{\|E(u)\|_F\,\|A(u)\|_F}.
$$
若 $C(u)\ll 1$ 且能隨網格收斂到 0，則可近似把 $U\approx U_{\rm geom}U_{\rm dyn}$；反之不可。

實用近似表達（非對角元微擾）：對 $m\neq n$，
$$
A_{mn}(u)\approx i\frac{\langle m(u)|\partial_u H(u)|n(u)\rangle}{E_n(u)-E_m(u)}
$$
並代入 $\partial_u H = \sum_{\alpha\beta}\dot h_{\alpha\beta}\,\sigma^\alpha\otimes\sigma^\beta$ 可得到在 Pauli 系數上的明確表達。

---

### 数值核对与算法摘要（便于复现）
\begin{enumerate}
\item 从 $H_{\rm spin}(u)$ 通过投影或直接矩阵元计算得到 $H_{\rm eff}^{(2)}(u)$，用上面的闭式关系提取 $\mathbf d(u)$。\
\item 计算采样矩阵 $D$ 的秩，或计算最大叉积模长验证是否能生成 $\mathfrak{su}(2)$。\
\item 计算 $\theta_{\pm}$ 并判定是否满足 $\int\|\mathbf d\|\,dt\in\pi\mathbb Z$。\
\item 计算本征基 $V(u)$，用中心差分估计 $\partial_u V$ 并求 $A(u)=iV^\dagger\partial_u V$，随后计算 $C(u)=\|[E,A]\|/(\|E\|\|A\|)$ 作分离性判据。\
\end{enumerate}

以上公式与步骤可直接对应你仓库中的实现（`tools/verify_from_R.py`、`tools/diagnose_decomposition.py`、`tools/majorana_lift.py`），并可作为文档中数值检验的“规范性附录”。

 
