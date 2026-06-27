# PRB111 动力学抵消条件：修正后的严格证明

---

## §0. 问题定义

PRB111 系统 $\dot U = KU$，$U(0)=I$，$U \in Sp(2)$。

演化矩阵 $U(T)$ 决定 SO(5) 旋转 $R_{ij} = \frac12\operatorname{Tr}(\Gamma_i U \Gamma_j U^\dagger)$。
MZM 编织由 $3\times3$ 子块 $R_{123} := R|_{i,j\in\{1,2,3\}} \in SO(3)$ 描述。

**抵消目标**：对 $E_1 \neq 0$，选择控制波形使 $R_{123}(T) = R_{123}^{\text{ideal}}$，
其中 $R_{123}^{\text{ideal}}$ 是 $E_1=t_1=0$ 的编织输出。

---

## §1. 定理一：抵消的充要条件（三方程）

**定理 1** $R_{123}(T) = R_{123}^{\text{ideal}}$ 等价于以下三个独立条件同时成立：

$$\boxed{X_0(T) + W_0(T) = S_1^{\text{ideal}}} \tag{C1}$$

$$\boxed{(X^2)_0(T) + (W^2)_0(T) - (\bar qWqX)_0(T) - (qX\bar qW)_0(T) = S_2^{\text{ideal}}} \tag{C2}$$

$$\boxed{\frac12\operatorname{Tr}_{\mathbb{C}}\!\big[\mu(\Gamma_1)\,\hat U(T)\,\mu(\Gamma_1)\,\hat U(T)^\dagger\big] = R_{11}^{\text{ideal}} = 1} \tag{C3}$$

其中 $S_1^{\text{ideal}} = \cos\theta_1^{\text{ideal}}+\cos\theta_2^{\text{ideal}}$，
$S_2^{\text{ideal}} = \cos 2\theta_1^{\text{ideal}}+\cos 2\theta_2^{\text{ideal}}$，
$\theta_{1,2}^{\text{ideal}}$ 是 $U_{\text{ideal}}$ 的本征角。

**证明分两步：**

### 1.1 (C1)+(C2) ⇔ 本征值匹配

引理 1.1–1.3 和引理 2.1–2.2 同前文，给出

$$\operatorname{Tr}[\hat U] = 2(X_0+W_0),\qquad \operatorname{Tr}[\hat U^2] = 2\big[(X^2)_0+(W^2)_0 - (\bar qWqX)_0 - (qX\bar qW)_0\big]$$

和

$$\operatorname{Tr}[\hat U] = 2(\cos\theta_1+\cos\theta_2),\qquad \operatorname{Tr}[\hat U^2] = 2(\cos 2\theta_1+\cos 2\theta_2)$$

**推论**：(C1)+(C2) ⇔ $\theta_1 = \theta_1^{\text{ideal}}$ 且 $\theta_2 = \theta_2^{\text{ideal}}$ ⇔
$U(T) \sim U_{\text{ideal}}$（共轭等价）。

但共轭等价 $\neq$ $R_{123}$ 相等：$U = S U_{\text{ideal}} S^{-1}$ 时，
$R[U] = \operatorname{Ad}_{\pi(S)}(R[U_{\text{ideal}}])$。$\pi(S)$ 是 $S \in Sp(2)$ 到 $SO(5)$ 的投影。
$\pi(S)$ 可非平凡地旋转 MZM 子空间。

### 1.2 第三条条件固定 MZM 子空间方向

$SO(3)$ 群有 3 个自由度。(C1)+(C2) 固定了 $R_{123}$ 的 2 个本征值 → 还剩 1 个自由度
（等价于旋转轴在固定本征值约束下的一个角度）。需一条额外条件来固定它。

自然选择：$\gamma_1$ 在编织中应映射到自己。
$$R_{11} = \frac12\operatorname{Tr}_{\mathbb{C}}[\mu(\Gamma_1)\hat U\mu(\Gamma_1)\hat U^\dagger] = 1$$

$R_{11}^{\text{ideal}} = 1$（$\gamma_1 \to \gamma_1$ 在理想编织中严格成立）。

因此 (C3) 等价于固定 $R_{123}$ 的最后一个自由度。

(C1)+(C2)+(C3) 共同等价于 $R_{123} = R_{123}^{\text{ideal}}$。$\square$

**注**：(C3) 在 Riccati 变量下的显式展开较冗长，但原则上可由 $q, X, W$ 表达。
在 Newton 迭代中，直接通过数值传播计算 $\hat U(T)$ 然后计算 $R_{11}$ 即可——
无需解析展开。

---

## §2. 定理二：解的存在性

**定理 2**（Chow）若 $\operatorname{Lie}\{K_{\text{drift}}, K_1,\dots,K_m\} = \mathfrak{g}$，
则右不变系统 $\dot U = (K_{\text{drift}}+\sum u_i K_i)U$ 的可达集是整个连通李群 $G$。

对本系统：
$$K_{\text{drift}} = \frac{E_1}{2}\Sigma_{12},\quad \{K_i\} = \{\Sigma_{24}, \Sigma_{34}, \Sigma_{45}, \Sigma_{15}\}$$

**引理**：$\operatorname{Lie}\{K_{\text{drift}}, K_1,\dots,K_4\} = \mathfrak{sp}(2)$（10 维满秩，report §2.4 已证）。

**推论**：$Sp(2) \cong Spin(5)$ 紧连通，可达集 $= Sp(2)$。$\exists$ 控制波形使 $U(T) = U_{\text{ideal}}$ 精确成立。

$U(T) = U_{\text{ideal}} \Rightarrow R_{123} = R_{123}^{\text{ideal}}$（精确，非近似）。$\square$

---

## §3. 有限参数化下的求解

将控制参数化为 $\mathbf{p} \in \mathbb{R}^n$（$n \ge 6$：三步独立 $\tau, t_c, E_0, t_1$ 等）。

定义残差 $\mathbf{R}: \mathbb{R}^n \to \mathbb{R}^3$：
$$R_1(\mathbf{p}) = X_0+W_0 - S_1^{\text{ideal}}$$
$$R_2(\mathbf{p}) = (X^2)_0+(W^2)_0 - (\bar qWqX)_0 - (qX\bar qW)_0 - S_2^{\text{ideal}}$$
$$R_3(\mathbf{p}) = \frac12\operatorname{Tr}[\mu(\Gamma_1)\hat U\mu(\Gamma_1)\hat U^\dagger] - 1$$

求解 $\mathbf{R}(\mathbf{p}) = \mathbf{0}$。$n \ge 6 > 3$，解空间 $\ge n-3$ 维。

**Jacobi 条件**：在解 $\mathbf{p}_0$（由定理 2 保证存在）处，需 $\operatorname{rank}(d\mathbf{R}|_{\mathbf{p}_0}) = 3$。
若成立，则由隐函数定理，解集是 $n-3$ 维光滑流形。
若秩不足，则需增加控制参数或放宽参数化约束。

---

## §4. 汇总

| 条件 | 物理含义 | 方程数 |
|------|---------|--------|
| (C1) $\operatorname{Tr}U$ 匹配 | 本征值和 | 1 |
| (C2) $\operatorname{Tr}U^2$ 匹配 | 本征值平方和 | 1 |
| (C3) $R_{11}=1$ | $\gamma_1\to\gamma_1$，固定 MZM 子空间方向 | 1 |
| **合计** | | **3 方程，$\ge 6$ 参数** |

PRB113 的 $\varepsilon+\sigma\lambda=0$ 是瞬时代数方程（1 方程，1 未知数）。
PRB111 的 $\mathbf{R}(\mathbf{p})=\mathbf{0}$ 是终态泛函方程（3 方程，$\ge 6$ 未知数）。
差异根源于 $\mathfrak{so}(5)$ 为单李代数，无不变子空间可供逐点对角化。
