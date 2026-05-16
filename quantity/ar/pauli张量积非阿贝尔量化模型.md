# Pauli张量积非阿贝尔量化模型（完整整理版）

## 1. 局域Hilbert空间与Pauli基

考虑一维格点系统，每个格点 $i$ 的局域Hilbert空间为：

\[
\mathcal{H}_i \simeq \mathbb{C}^2
\]

全局Hilbert空间：

\[
\mathcal{H} = \bigotimes_{i=1}^L \mathcal{H}_i
\]

每个格点引入Pauli基：

\[
\sigma_i^\mu,\quad \mu \in \{0,x,y,z\},\quad \sigma_i^0 = \mathbb{I}
\]

满足代数：

\[
[\sigma_i^a, \sigma_j^b] = 2i\,\delta_{ij}\,\epsilon^{abc}\sigma_i^c
\]

---

## 2. Pauli张量积展开（算符基）

任意全局算符 $O$ 可唯一展开为：

\[
O = \sum_{\vec{\mu}} c_{\vec{\mu}} \bigotimes_{i=1}^L \sigma_i^{\mu_i}
\]

定义Pauli string：

\[
P_{\vec{\mu}} := \bigotimes_i \sigma_i^{\mu_i}
\]

因此：

\[
O = \sum_{\vec{\mu}} c_{\vec{\mu}} P_{\vec{\mu}}
\]

---

## 3. Pauli代数与非交换结构

单点乘法：

\[
\sigma^a \sigma^b = \delta^{ab}\mathbb{I} + i\epsilon^{abc}\sigma^c
\]

从而：

\[
P_\alpha P_\beta = e^{i\phi(\alpha,\beta)} P_{\alpha\circ\beta}
\]

结论：Pauli string代数封闭但非交换。

---

## 4. 非阿贝尔性来源

考虑：

\[
A = \sum_\alpha a_\alpha P_\alpha,\quad B = \sum_\beta b_\beta P_\beta
\]

对易子：

\[
[A,B] = \sum_{\alpha,\beta} a_\alpha b_\beta [P_\alpha,P_\beta]
\]

若同一格点存在 $x,y,z$ 混合，则：

\[
[P_\alpha,P_\beta] \neq 0
\]

因此系统天然非阿贝尔。

---

## 5. 非阿贝尔性量化定义

### 5.1 对易子范数

\[
\mathcal{N}(A,B) := \|[A,B]\|
\]

---

### 5.2 局域非阿贝尔密度

\[
A_i = \sum_{a=x,y,z} a_i^a \sigma_i^a
\]

\[
\mathcal{N}_i = \sum_{a,b,c} |a_i^a a_i^b|^2
\]

---

### 5.3 Abelian极限

\[
[A_i,A_j]=0 \Rightarrow \mathcal{N}=0
\]

或：

\[
A_i \propto \sigma_i^z
\]

---

## 6. 完整哈密顿量结构

### 6.1 Pauli展开形式

\[
H = \sum_{\vec{\mu}} J_{\vec{\mu}} P_{\vec{\mu}}
\]

---

### 6.2 局域场项

\[
H_{loc} = \sum_i \vec{h}_i \cdot \vec{\sigma}_i
\]

---

### 6.3 相互作用项

\[
H_{int} = \sum_{\langle i,j\rangle}\sum_{a,b} J_{ij}^{ab}\sigma_i^a \sigma_j^b
\]

---

## 7. 非阿贝尔代数结构

定义代数：

\[
\mathfrak{A} = \mathrm{span}\{P_{\vec{\mu}}\}
\]

非交换结构：

\[
[\mathfrak{A},\mathfrak{A}] \subseteq \mathfrak{A}
\]

非阿贝尔强度：

\[
\mathcal{N} \sim \|[\mathfrak{A},\mathfrak{A}]\|
\]

---

## 8. 非阿贝尔量化（泛函积分形式）

### 8.1 非阿贝尔作用量

\[
\mathcal{N}[c] = \sum_{\alpha,\beta} |c_\alpha c_\beta|^2 \|[P_\alpha,P_\beta]\|^2
\]

---

### 8.2 泛函积分定义

\[
Z_{NA} = \int Dc\; e^{-\lambda \mathcal{N}[c] - \beta \mathcal{E}[c]}
\]

---

## 9. 谱流定义（Spectral Flow）

### 9.1 流算符

\[
O(s)=e^{-s\mathcal{D}} O e^{s\mathcal{D}}
\]

---

### 9.2 非阿贝尔强度流

\[
\mathcal{N}(s)=\|[O(s),O(s)]\|
\]

---

### 9.3 固定点

- Abelian：$\mathcal{N}=0$
- 非阿贝尔：$\mathcal{N}>0$
- 临界：$d\mathcal{N}/ds=0$

---

## 10. PRP闭合结构

### 10.1 投影定义

\[
\Pi_{PRP}:\mathfrak{A}\to\mathcal{S}
\]

---

### 10.2 闭合条件

\[
[\mathcal{S},\mathcal{S}]\subseteq \mathcal{S}
\]

---

### 10.3 泄漏算子

\[
\mathcal{L} = \|(1-\Pi_{PRP})[A,B]\|^2
\]

PRP闭合 ⇔ $\mathcal{L}=0$

---

## 11. 零模（Zero Mode / MZM）条件

### 11.1 基本定义

\[
[H,\gamma]=0
\]

\[
\gamma^2=1,\quad \gamma^\dagger=\gamma
\]

---

### 11.2 PRP投影形式

\[
\Pi_{PRP}([H,\gamma])=0
\]

\[
(1-\Pi_{PRP})([H,\gamma])=0
\]

---

### 11.3 统一判据

\[
\mathcal{F}[\gamma]
= \|(1-\Pi_{PRP})[H,\gamma]\|^2 + \|\Pi_{PRP}[H,\gamma]\|^2
\]

零模存在 ⇔ $\mathcal{F}[\gamma]=0$

---

## 12. 总结

该模型统一描述：

- Pauli张量积代数结构
- 非阿贝尔性来源（SU(2)结构）
- 非阿贝尔量化（泛函积分）
- 谱流动力学
- PRP闭合结构
- 零模/MZM存在条件

核心统一图像：

> 非阿贝尔性 = Pauli代数非交换结构
>
> PRP闭合 = 非阿贝尔泄漏为零
>
> 零模 = PRP闭合子空间中的中心算符

