# 离域费米子模的推导与物理含义

## 目标

推导由两端 Majorana 组合得到的“离域费米子模”算符的代数性质，证明其占据数为 0 或 1、在理想极限下为何与哈密顿量对易，并解释这些结果的物理意义（基态简并、拓扑保护与有限体系的能量裂分）。

## 一、定义与代数推导

设系统在一维链的两端存在局域 Majorana 算符 \(\gamma_L\) 与 \(\gamma_R\)，它们满足
\[
\{\gamma_i,\gamma_j\}=2\delta_{ij},\qquad \gamma_i^\dagger=\gamma_i.
\]

定义复费米子（离域费米子）
\[
d=\frac{\gamma_L + i\gamma_R}{2},\qquad d^\dagger=\frac{\gamma_L - i\gamma_R}{2}.
\]

验证反对易关系：
\[
\{d,d^\dagger\}=1,\qquad d^2=(d^\dagger)^2=0.
\]

占据数算符定义为
\[
n_d = d^\dagger d.
\]

计算其平方以证明投影性质：
\[
n_d^2=(d^\dagger d)(d^\dagger d)=d^\dagger(1-d d^\dagger)d=d^\dagger d=n_d.
\]

因此 \(n_d\) 是投影算符，其本征值只能为 0 或 1，说明该离域费米子模式要么为空，要么被占据一个费米子。

等价地，用 Majorana 表达 \(n_d\)：
\[
n_d=\frac{1+i\gamma_L\gamma_R}{2}.
\]
由于算符 \(i\gamma_L\gamma_R\) 是厄米且满足 \((i\gamma_L\gamma_R)^2=1\)，其本征值仅为 \(\pm1\)，代入后得到 \(n_d=(1\pm1)/2\in\{0,1\}\)。

## 二、为何 \([d,H]=0\)（或 \([n_d,H]=0\)）在理想极限成立

在 Kitaev 链或文中讨论的极限（例如 \(t=|\Delta_p|\) 且 \(\mu=0\)）下，哈密顿可重写成由链内部成对费米子占据项构成的形式，例如
\[
H = 2t\sum_{j=1}^{L-1}\left(\tilde f_j^\dagger\tilde f_j - \tfrac12\right),
\]
此处所有出现的算符都与链两端的 Majorana （\(\gamma_L,\gamma_R\)）不相耦合：两端的 Majorana 被哈密顿“解耦”出来，不出现在 H 的项中。故对构成离域费米子的算符 \(d\) 有
\[
[d,H]=0,\qquad [d^\dagger d,H]=0.
\]

物理含义：哈密顿不包含能改变 \(d\) 占据数的项，因此 \(n_d\) 在该模型极限下是守恒量。于是基态分为两类——\(n_d=0\) 與 \(n_d=1\)——它们在能量上简并（若没有其它约束），从而产生拓扑基态的二重简并。

注意：在更一般的参数处（非理想极限或链為有限长度）可能存在微小耦合项将两端 Majorana 间接耦合，导致 \([d,H]\neq0\) 并造成能量裂分（参見下节）。

## 三、有限长度与能量裂分

在有限长度链中，左右 Majorana 的波函数有指数尾巴並會重叠，重叠導致離域模不再严格為零能，而是產生能量劈裂 \(\Delta E\)。常見估計為
\[
\Delta E \propto e^{-L/\xi},
\]
其中 \(L\) 为链长，\(\xi\) 为相干长度（与谱隙成反比）。因此要獲得近似零能且拓扑保护良好的離域模，需要保证 \(L\gg\xi\)。

## 四、物理意義總結

- 離域費米子模式将量子信息编码在空间上分離的两个 Majorana 之间，信息对局域干扰有天然鲁棒性（局域扰动通常不能翻转 \(n_d\)）。
- 在理想極限下（兩端 Majorana 完全從哈密頓解耦）占據數 \(n_d\) 是守恒量，導致基態在費米子占據上成對簡并，這是拓扑相的一大特徵。
- 在有限體系或非理想參數下，左右 Majorana 波函數重疊會引入小的能量裂分，限制拓扑保护的時长与质量；裂分隨系统尺寸指數衰减。
- 讀出与操控需要將兩端模式耦合或進行非局域測量，因此存在讀寫與保護之間的折衷。

## 五、参考與延伸

- 可參考 Kitaev 鏈原始推導與綜述（例如 A. Yu. Kitaev, Physics-Uspekhi 2001；以及近期綜述與實驗報告）以獲得更多背景和數值示例。

---
（已在項目中保存為 `delocalized_fermion.md`）

## 附：Kitaev 链的显式推导（端点 Majorana 解耦）

考虑标准的 Kitaev 链哈密顿（自旋极化费米子）：
\[
H=-\mu\sum_{j=1}^L c_j^\dagger c_j - t\sum_{j=1}^{L-1}(c_j^\dagger c_{j+1}+\mathrm{h.c.}) + \Delta\sum_{j=1}^{L-1}(c_j c_{j+1}+\mathrm{h.c.}).
\]
引入 Majorana 算符（每个格点两个）：
\[
\gamma_{2j-1}=c_j + c_j^\dagger,\qquad \gamma_{2j}=\frac{c_j-c_j^\dagger}{i},
\]
它们满足 \(\{\gamma_a,\gamma_b\}=2\delta_{ab}\)。把 H 写成 Majorana 形式（省略常数项），可得（经代数化简）
\[
H=\frac{i}{2}\sum_{j=1}^L\left[-\mu\,\gamma_{2j-1}\gamma_{2j} + (t+\Delta)\,\gamma_{2j}\gamma_{2j+1} + ( -t+\Delta)\,\gamma_{2j-1}\gamma_{2j+2}\right].
\]

在常用的理想化参数点取 \(\mu=0\) 且 \(t=\Delta\)，上式簡化為
\[
H = i t\sum_{j=1}^{L-1}\gamma_{2j}\gamma_{2j+1}.
\]
注意这里的和从 \(j=1\) 到 \(L-1\)，因此算符 \(\gamma_1\)（即第一个 Majorana）和 \(\gamma_{2L}\)（最后一个 Majorana）並未出现在 H 的任何项中——它们与哈密顿解耦。

于是定义离域费米子
\[
d=\frac{\gamma_1 + i\gamma_{2L}}{2},\quad d^\dagger=\frac{\gamma_1 - i\gamma_{2L}}{2},
\]
有
\[
[d,H]=0,
\]
因为 \(H\) 不含 \(\gamma_1,\gamma_{2L}\)。同理 \([n_d,H]=0\)，其中 \(n_d=d^\dagger d=(1+i\gamma_1\gamma_{2L})/2\)。

这就从具体哈密顿层面证明了在该理想点端点 Majorana 是严格零模并与 H 对易，从而产生守恒的离域占据数与基态二重简并。

在偏离该理想点（例如 \(t\neq\Delta\) 或 \(\mu\neq0\)）或有限链长时，上述解耦被打破，端点 Majorana 会以指数方式耦合进而导致能级劈裂（见正文第 III 节）。
