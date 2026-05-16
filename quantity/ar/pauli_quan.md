# Pauli基下的算符展开与非阿贝尔性量化（完整模型整理）

## 1. 局域Hilbert空间与Pauli基

考虑一维格点系统，每个格点 \(i\) 的局域Hilbert空间为
\[
\mathcal{H}_i \simeq \mathbb{C}^2
\]

全局Hilbert空间为：
\[
\mathcal{H} = \bigotimes_{i=1}^L \mathcal{H}_i
\]

每个格点引入 Pauli 基：
\[
\sigma_i^\mu,\quad \mu \in \{0,x,y,z\}, \quad \sigma_i^0 = \mathbb{I}
\]

满足代数关系：
\[
[\sigma_i^a, \sigma_j^b] = 2i\,\delta_{ij}\,\epsilon^{abc}\sigma_i^c
\]

---

## 2. Pauli string 基展开

任意算符 \(O\) 可唯一展开为 Pauli string 线性组合：

\[
O = \sum_{\vec{\mu}} c_{\vec{\mu}} \bigotimes_{i=1}^L \sigma_i^{\mu_i}
\]

定义：
\[
P_{\vec{\mu}} := \bigotimes_{i=1}^L \sigma_i^{\mu_i}
\]

则：
\[
O = \sum_{\vec{\mu}} c_{\vec{\mu}} P_{\vec{\mu}}
\]

其中：
- \(\vec{\mu} = (\mu_1,\dots,\mu_L)\)
- \(c_{\vec{\mu}} \in \mathbb{C}\)

---

## 3. Pauli代数与乘法结构

单点Pauli乘法：
\[
\sigma^a \sigma^b = \delta^{ab}\mathbb{I} + i\epsilon^{abc}\sigma^c
\]

因此 Pauli string 满足封闭代数：
\[
P_{\alpha} P_{\beta} = e^{i\phi(\alpha,\beta)} P_{\alpha \circ \beta}
\]

其中：
- 相位来自局域反交换结构
- \(\alpha \circ \beta\) 表示逐点组合规则

---

## 4. 非阿贝尔性来源

考虑两个一般算符：
\[
A = \sum_{\alpha} a_\alpha P_\alpha,\quad
B = \sum_{\beta} b_\beta P_\beta
\]

对易子为：
\[
[A,B] = \sum_{\alpha,\beta} a_\alpha b_\beta [P_\alpha, P_\beta]
\]

关键事实：

- 若Pauli string 在同一格点包含不同方向（\(x,y,z\)混合）
- 则局域产生非零对易子
- 从而全局产生非阿贝尔结构

---

## 5. 非阿贝尔性量化定义

### 5.1 对易子范数定义

\[
\mathcal{N}(A,B) := \|[A,B]\|
\]

（可取 Hilbert-Schmidt 范数）

---

### 5.2 局域非阿贝尔密度

设局域算符：
\[
A_i = \sum_{a=x,y,z} a_i^a \sigma_i^a
\]

定义：
\[
\mathcal{N}_{\text{loc}} =
\frac{1}{L}\sum_i \sum_{a,b}
\|[A_i^a, A_i^b]\|^2
\]

---

### 5.3 显式Pauli形式

局域对易子：
\[
[A_i, A_i]
= 2i \sum_{a,b,c}
a_i^a a_i^b \epsilon^{abc}\sigma_i^c
\]

因此非阿贝尔强度：
\[
\mathcal{N}_i \sim \sum_{a,b,c} |a_i^a a_i^b|^2
\]

---

## 6. Abelian极限

当满足：
\[
[A_i, A_j] = 0,\quad \forall i,j
\]

或局域取定方向：
\[
A_i \propto \sigma_i^z
\]

则：
\[
\epsilon^{abc} \to 0 \Rightarrow \mathcal{N} = 0
\]

系统退化为 Abelian 可对角结构。

---

## 7. 非阿贝尔结构生成机制

一般局域算符写为：
\[
A_i = \vec{n}_i \cdot \vec{\sigma}_i
\]

其中：
\[
\vec{n}_i \in \mathbb{R}^3
\]

对易结构由：
\[
[A_i, A_j] = 2i\delta_{ij}(\vec{n}_i \times \vec{n}_i)\cdot \vec{\sigma}_i
\]

结论：
- \(\vec{n}_i \parallel \vec{n}_j\)：Abelian
- 方向变化：Non-Abelian

---

## 8. 完整模型哈密顿量

### 8.1 一般Pauli展开形式

\[
H = \sum_{\vec{\mu}} J_{\vec{\mu}} P_{\vec{\mu}}
\]

---

### 8.2 局域场项

\[
H_{\text{loc}} = \sum_i \vec{h}_i \cdot \vec{\sigma}_i
\]

---

### 8.3 交换相互作用

\[
H_{\text{int}} =
\sum_{\langle i,j \rangle}
\sum_{a,b}
J_{ij}^{ab}\sigma_i^a \sigma_j^b
\]

---

## 9. 非阿贝尔代数结构总结

定义算符代数：
\[
\mathfrak{A} = \mathrm{span}\{P_{\vec{\mu}}\}
\]

满足：
\[
[\mathfrak{A}, \mathfrak{A}] \subseteq \mathfrak{A}
\]

非阿贝尔性测度：
\[
\mathcal{N} \sim \|[\mathfrak{A},\mathfrak{A}]\|
\]

---

## 10. 核心结论

在 Pauli 张量基中：

> 非阿贝尔性本质上来自 SU(2) 局域生成元在张量积空间中的非交换结构累积。

核心来源：
\[
\sigma^x \sigma^y \neq \sigma^y \sigma^x
\]

并推广为：
\[
[P_\alpha, P_\beta] \neq 0
\Rightarrow \text{全局非阿贝尔代数}
\]