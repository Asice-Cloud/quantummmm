# 描述非阿贝尔 SO(5) 演化的数学工具

> 针对 PRB111 5-Majorana 体系（SO(5) / Sp(2) 四元数 Riccati 框架）的数学描述工具
>
> 核心问题：满 10 维紧半单李代数 so(5) 没有非平凡理想，不能被分解为更小子代数的直和，
> 导致非线性耦合不可消除。

---

## 一、SO(5) 的 Cartan 分解 —— 区分可分解与不可分解部分

### KAK 分解

紧李群 G 有 Cartan 分解：**g = k ⊕ p**，k 为极大紧子代数，p 为 k 在 g 中的正交补。

```
so(5) = so(4) ⊕ p
       = (su(2) ⊕ su(2)) ⊕ p    (p 含 4 个"真正非阿贝尔"方向)
```

### 对本系统的应用

| 编织段 | 活跃生成元 | 李代数闭包 | Cartan 分量 |
|--------|-----------|-----------|-------------|
| Step 1 | E₁, t₂, t₁ | so(4) (6D) | **k** — 可分解 |
| Step 2 | 全 5 项 | so(5) (10D) | **k ⊕ p** — 满，不可分解 |
| Step 3 | E₁, t₃, E_d | u(1)⊕su(2) (4D) | k 的子代数 |

### Riccati 对应

```
A = D  →  F_{μν} 在 k = so(4) 内 → 非阿贝尔曲率仅来自 [A_μ, A_ν]
A ≠ D  →  p 分量被激活 → 联络中包含威尔切克-泽的"带电标量场" ψ 的贡献
```

**参考**：Helgason, *Differential Geometry, Lie Groups, and Symmetric Spaces* (1978), Ch. V-VII

---

## 二、SO(5) 的 B₂ 根系统 —— 解释"轴锁定"

### so(5) ≅ sp(2) 的根系统

so(5) 对应 B₂ = C₂ 型（秩 2），共 8 个根：

```
         α₁+2α₂  (短根，最高根，对应 σ_z)
              ↑
 α₁ ←──●──→ α₁+α₂ ←──●──→ 2α₁+α₂
              ↑
             α₂  (长根，对应 σ_x, σ_y 平面内的旋转)
```

- **2 个短根** (±α₁, ±(α₁+2α₂))：对应 iγ₁γ₂（σ_z）方向
- **6 个长根** (±α₂, ±α₁+α₂, ±2α₁+α₂)：对应 iγ₂γ₃, iγ₃γ₁ 等 xy 平面方向

### Pauli 矩阵在根空间中的位置

| Pauli | Majorana 生成元 | 对应根 |
|-------|----------------|--------|
| σ_x | iγ₂γ₃ | 长根 α₂ |
| σ_y | iγ₃γ₁ | 长根 α₂ 的 Weyl 群反射 |
| σ_z | iγ₁γ₂ | 短根 α₁ |

### 解释 E₁=0 的轴锁定

```
E₁=0 → 短根 α₁ 方向的生成元无法被激活
A=D → Weyl 群反射仅在长根平面内
→ 旋转轴锁在 xy 平面（σ_x, σ_y 方向）
```

> **E₁ 打破 A=D 对称性 → 短根被激活 → 满 B₂ 根系统 → 三 Pauli 全部释放**

### Cartan 矩阵和 Dynkin 图

```
B₂ Dynkin:  ○=<=○
α₂(长)    α₁(短)

Cartan:  [ 2 -2 ]
         [-1  2 ]
```

**参考**：Humphreys, *Introduction to Lie Algebras and Representation Theory* (1972), Ch. 12; Fulton & Harris, *Representation Theory* (1991), Ch. 18

---

## 三、余伴随轨道方法（Kirillov-Kostant-Souriau）—— 从李代数到辛几何

### 核心思想

李代数 **g** 的对偶空间 **g*** 被自然划分为**余伴随轨道**（coadjoint orbits）：

```
O_x = {Ad*_g(x) : g ∈ SO(5)},   x ∈ so(5)*
```

每个轨道是**自然的辛流形**（Kirillov-Kostant-Souriau 辛形式）。动力学在轨道上是哈密顿的。

### SO(5) 的余伴随轨道分层

SO(5) 的轨道按**秩**（半单元素的独立特征值数）分类：

| 轨道类型 | 维数 | 辛结构 | 物理对应 |
|---------|------|--------|---------|
| 零轨道 (秩 0) | 0 | 平凡 | 平凡态 |
| 最小轨道 (秩 1, 短根) | 4 | ℂℙ³ | **我们的 q 变量空间** |
| 最小轨道 (秩 1, 长根) | 6 | 旗流形 SO(5)/U(2) | Step 1 so(4) |
| 一般轨道 (秩 2, 正则) | 8 | SO(5)/T² | **Step 2 满 so(5)** |

### 与本系统的精确对应

**Riccati 变量 q ∈ ℍ 的 4 个实分量 → 最小余伴随轨道 (4 维 ℂℙ³)**

推导链：
```
q = ZX⁻¹
|X|² + |Z|² = 1 (S⁷)
U(1) 商 → ℂℙ³ (4 维)
ℂℙ³ ≅ SO(5)/U(2) — 一个余伴随轨道
```

**工具价值**：
1. 轨道上的 Kirillov-Kostant 辛形式 ω → q(t) 的泊松括号结构
2. Casimir 函数在轨道上为常数 → 演化约束
3. 轨道间的跃迁对应对称破缺（E₁=0 → E₁≠0 的轨道分岔）

### 辛约化

```
so(5)* (10 维)
    ↓ 辛约化 (除以未扰动对称性)
ℂℙ³ (4 维)  →  q(t) 的相空间
```

这解释了为什么从 10 个生成元到 4 分量 Riccati ODE 的约化是辛的。

**参考**：Kirillov, *Lectures on the Orbit Method* (2004); Kostant, *Quantization and Unitary Representations* (1970); Arnold, *Mathematical Methods of Classical Mechanics* (1989), Appendix 2

---

## 四、参数空间上的规范理论 —— Wilczek-Zee 联络

### 在 (E₁, t₁, t₂, t₃, E_d, τ) 参数空间上

每步编织 = 参数空间中的一条路径。三步编织 = 闭合回路。

在 MZM 子空间（简并基态扇区）上，参数空间承载一个**非阿贝尔 U(2) 主丛**：

```
Wilczek-Zee 联络： A_{ab}(λ) = ⟨ψ_a(λ)|∂_λ|ψ_b(λ)⟩
```

### 曲率的非阿贝尔分量

```
F_{μν} = ∂_μ A_ν - ∂_ν A_μ - [A_μ, A_ν]
         \_________ ________/   \_______ ________/
            阿贝尔部分           非阿贝尔部分
            (Dirac 磁单极)      (杨-米尔斯型)
```

**非阿贝尔部分 [A_μ, A_ν] ≠ 0 当且仅当不同参数方向的生成元不对易。**

### 对应本系统

| 参数 | 在联络中的角色 | 物理起源 |
|------|--------------|---------|
| E₁ | 对角分量 A_z（短根方向） | γ₁-γ₂ 杂化 |
| t₁ | 横截分量 A⊥_y | γ₁-ancilla 动态耦合 |
| t₂, t₃ | 横截分量 A⊥_x | 编织耦合 |
| E_d | 对角分量 | ancilla 自身动力学 |

### 与 Dai et al. (2512.24798) 的精确定量映射

| Dai 论文 | 本系统 |
|---------|--------|
| 形状球面 S²_K (θ, φ) | (E₁, t₁) 参数平面 |
| 对角 U(1) 联络 A·σ_z | E₁ 主导的 Abelian Berry 相 |
| 带电标量 ψ = ρ + iσ | A≠D 的横截非阿贝尔混合 |
| ∮A = ½ Ω_Γ | ∮ 编织 = π/2 (几何角度) |
| Tr W_Γ = 2 cos(qΩ/2) | Fidelity = |⟨ψ₁⁻|U|ψ₁⁺⟩|² |

**参考**：Wilczek & Zee, *PRL* 52, 2111 (1984); Dai et al., arXiv:2512.24798 (2025); Nakahara, *Geometry, Topology and Physics* (2003), Ch. 10

---

## 五、量子度规与测地线 —— 最优绝热演化

### 定义

在简并子空间上，量子度规：

```
g_{μν} = Tr(∂_μ P ∂_ν P)
```

其中 P(λ) = Σ_i |D_i(λ)⟩⟨D_i(λ)| 是基态投影算符。

### 测地线方程

```
δ ∫ √(g_{μν} ∂_z λ^μ ∂_z λ^ν) dz = 0
```

→ 最小非绝热误差的演化路径。

### 对本系统的应用

在 (E₁, t₁, τ) 参数空间中：
- **小 τ 区**：演化偏离测地线 → 大非绝热误差 → 快速振荡
- **大 τ 区**：演化趋近测地线 → 绝热极限 → 保真度收敛
- **保真度峰** = 参数空间中的测地线交点（相长干涉）

### 实验优化

```
目标：在 (E₁, t₁, τ) 空间中沿量子度规的测地线驱动参数
实现：最小非绝热误差，保真度最优
```

这与 Kremer et al. (arXiv:1902.02559) 的"基于量子度规的最优非阿贝尔几何相设计"完全对应。

**参考**：Kremer et al., arXiv:1902.02559 (2019, 2024更新); Provost & Vallee, *Commun. Math. Phys.* 76, 289 (1980)

---

## 六、Casimir 不变量 —— 演化约束

### so(5) 的两个独立 Casimir

索引约定：X_{ij} = iγ_i γ_j (i < j)，生成元 adj 表示。

```
C₂ = ½ Σ_{i<j} X_{ij}²                  （二次 Casimir, 10 项求和）
C₄ = Σ ε^{ijklmn...} X_{ij}X_{kl}...     （四次 Casimir, Pfaffian 型）
```

### 轨道上的定值

每个余伴随轨道上 Casimir 取定值 → **运动常数**：

| 情形 | C₂ 值 | 含义 |
|------|-------|------|
| 纯几何 (E₁=0, t₁=0) | 特定值 | 约束 q 在赤道面 |q|=const |
| E₁≠0 | 不同值 | q 被释放到 ℂℙ³ 的一般位置 |

### 数值验证用途

- Casimir 沿着数值解是否保持恒定 → 代码正确性的强判决
- 不同 C₂ 值对应不同轨道 → 不同物理相的区分

**参考**：Racah, *Group Theory and Spectroscopy* (1949); Wybourne, *Classical Groups for Physicists* (1974), Ch. 17

---

## 七、几何控制理论 —— 在群流形上的可达性

### Agrachev 的"好李括号"

群流形 G 上的控制问题：给定漂移项和可控生成元，能到达群流形上的哪些点？

```
U̇(t) = (K_drift + Σ u_i(t) K_i) U(t),   K_i ∈ g
```

小时间可控性条件（Agrachev et al., 2025）：
- 李括号生成整个李代数（Lie algebra rank condition）
- 存在"好李括号"——可以在任意小时间内被近似

### 对本系统的对应

我们的 Riccati 对此问题给出了一个具体的解法：

```
问题                          对应
─────────────────────────────────────────
K_drift = A+D 的部分（漂移）   → A 和 D 的四元数表示
K_i：时变生成元               → t₂, t₃, E_d 的时变系数
q(t) 的可达集                 → so(5) 中某个余伴随轨道
"好李括号"条件                → A≠D 释放的非线性耦合
```

**关键问题**：在你的参数约束下（固定门控包络 f_±），可达旋转的集合是什么？

**参考**：Agrachev, Kazandjian, Pozzoli, *Systems & Control Letters* 205, 106233 (2025); Jurdjevic, *Geometric Control Theory* (1997)

---

## 八、旗流形 —— 演化态的几何空间

### 定义

李群 G 模去抛物子群 P 的商：

```
G/P = { 0 ⊂ V₁ ⊂ V₂ ⊂ ... ⊂ V_k = V }
```

对于 so(5)，满旗流形是 SO(5)/T²（T² 为极大环面，12-2=10 维，减去满旗的多余维度得 8 维正则轨道）。

### 对应本系统

你的旋量表示中：
- 总希尔伯特空间 ℍ² ≅ ℂ⁴ (4 维)
- MZM 子空间 ≈ 2 维子空间
- ancilla 子空间 ≈ 正交补

→ **演化在 Gr(2,4) ≅ SO(4)/S(O(2)×O(2))** 上，嵌入在 SO(5) 的旗流形中。

### 为什么是 4 分量

Gr(2,4) 的实维数 = 2 × (4-2) × 2 = 8。加上 U(1) 规范自由度和全局约束 → **恰好 4 个独立实分量**。

**参考**：Griffiths & Harris, *Principles of Algebraic Geometry* (1978); Borel, *Linear Algebraic Groups* (1991)

---

## 九、总结：工具之间的逻辑关系

```
                    SO(5) Cartan 分解
                 so(5) = so(4) ⊕ p
                 /                  \
        根系统 B₂                  Casimir 不变量
       (轴锁定机制)              (轨道分类约束)
       /                              \
  余伴随轨道方法                      几何控制理论
  (辛结构, 轨道分层)                 (可达集, 好李括号)
       \                              /
        参数空间上的规范理论 (WWZ)
         (联络, 曲率, Holonomy)
                |
          量子度规与测地线
         (最优绝热路径设计)
                |
          旗流形 G/P
       (整体几何图景)
```

### 从最具体到最抽象的层次

| 层次 | 工具 | 回答的问题 |
|------|------|-----------|
| **代数** | B₂ 根系统 | 为什么 E₁=0 轴锁 xy？ |
| **分解** | Cartan KAK | 为什么 Step 2 不能分解？ |
| **轨道** | 余伴随轨道 | q 在什么空间上运动？ |
| **规范** | WWZ 联络 | 参数如何产生非阿贝尔混合？ |
| **度量** | 量子度规 | 最优演化路径？ |
| **控制** | 好李括号 | 给定协议能到达哪些态？ |
| **整体** | 旗流形 | 演化在几何上是什么？ |

---

## 十、推荐文献

1. **Wilczek & Zee**, "Appearance of Gauge Structure in Simple Dynamical Systems", *Phys. Rev. Lett.* **52**, 2111 (1984) — 非阿贝尔几何相的起源

2. **Dai, Molochkov, Niemi, Westerholm**, "Non-Abelian Geometric Phases in Triangular Structures and Universal SU(2) Control in Shape Space", arXiv:2512.24798 (2025) — Wilczek-Zee 联络在形状空间上的完整构造

3. **Kremer, Teuber, Szameit, Scheel**, "Optimal Design Strategy for Non-Abelian Geometric Phases Using Abelian Gauge Fields Based on Quantum Metric", arXiv:1902.02559 (2019, 2024更新) — 量子度规定义最优绝热路径

4. **Agrachev, Kazandjian, Pozzoli**, "Good Lie Brackets for Classical and Quantum Harmonic Oscillators", *Systems & Control Letters* **205**, 106233 (2025) — 李群上的几何控制理论

5. **Kirillov**, *Lectures on the Orbit Method*, AMS (2004) — 余伴随轨道方法的标准参考

6. **Humphreys**, *Introduction to Lie Algebras and Representation Theory*, Springer (1972) — 根系统、Cartan 矩阵、Dynkin 图

7. **Fulton & Harris**, *Representation Theory: A First Course*, Springer (1991) — SO(5) 的旋量表示和 B₂ 根系统

8. **Jurdjevic**, *Geometric Control Theory*, Cambridge (1997) — 李群上的可达性和可控性

9. **Nakahara**, *Geometry, Topology and Physics*, IoP (2003), Ch. 10 — 规范场、联络、曲率、Wilson 环

10. **Bärligea, Sims-Goh, Kottmann**, "Enabling Lie-Algebraic Classical Simulation beyond Free Fermions", arXiv:2604.16701 (2026) — 多项式维数 DLA 的轨道结构
