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
H(u)=\sum_{\alpha,\beta} h_{\alpha\beta}(u),\sigma^\alpha\otimes\sigma^\beta
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
H(u)|n(u)\rangle = E_n(u)_n(u)\rangle
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
\partial_u H(u)=\sum_{\alpha\beta} \dot h_{\alpha\beta}(u),\sigma^\alpha\otimes\sigma^\beta
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

 
