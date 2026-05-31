# PRB111 推导续篇：Step 3 的 `Sp(2)` 主线与单位四元数化

这份短文只承接前面的 full-space 推导，不再重复 `so(5)` 的整体搭建。这里的目标很明确：**先保留 `Sp(2)` 的精确时间有序指数，再只在真正闭合到 `su(2)` 的部分写单位四元数闭式**。

## 1. 处理原则

对三段门控协议，主线始终是

$$
U(t)=\mathcal T\exp\left(\int_0^t K(t')\,dt'\right),
\qquad K(t)\in\mathfrak{sp}(2).
$$

这里不预设全局 `SU(2)`，也不把整段协议先验压缩成单位四元数。只有在某一段真的闭合到 `su(2)` 时，才把那一段改写成 `Sp(1)\cong SU(2)` 的四元数形式。

## 2. Step 3 的严格分解

从前文的 `5\times 5` 反对称矩阵可读出，Step 3 的非平凡生成元可以写成

$$
K^{(3)}(t)=\kappa_0\,X_{12}+\alpha(t)\,X_{3a}+\beta(t)\,X_{ab}.
$$

在当前的归一化下，`\kappa_0` 对应 `E_1` 那条不随时间变化的通道，而 `\alpha(t),\beta(t)` 分别对应 `t_3^{(3)}(t)` 和 `E_d^{(3)}(t)` 的系数；若沿用 `5\times 5` 矩阵的约定，它们只差一个固定的数值归一化因子，不影响李代数闭包结论。

这里的关键是：

$$
[X_{12},X_{3a}]=[X_{12},X_{ab}]=0.
$$

所以 `X_{12}` 是一个纯粹对易的 `u(1)` 因子；真正决定可解析化结构的是剩下的三维子代数。

## 3. 真正闭合的 `su(2)` 子块

取

$$
J_1:=X_{3a},\qquad J_2:=X_{3b},\qquad J_3:=X_{ab}.
$$

则由 `so(5)` 的标准 commutator 公式可直接得到

$$
[J_1,J_2]=2iJ_3,\qquad [J_2,J_3]=2iJ_1,\qquad [J_3,J_1]=2iJ_2.
$$

因此

$$
\mathfrak s:=\operatorname{span}\{J_1,J_2,J_3\}\cong su(2).
$$

于是 Step 3 的李代数结构是

$$
K^{(3)}(t)=\kappa_0\,X_{12}+K^{(3)}_{su(2)}(t),
\qquad K^{(3)}_{su(2)}(t)\in\mathfrak s,
$$

也就是

$$
\mathfrak{u}(1)\oplus su(2).
$$

## 4. 精确因子分解

因为 `X_{12}` 与 `\mathfrak s` 中的所有生成元都对易，所以时间演化算符可以严格分解为

$$
U^{(3)}(t)=e^{\Phi(t)X_{12}}\,V(t),
$$

其中

$$
\Phi(t)=\int_0^t \kappa_0\,dt' = \kappa_0 t
$$

在 `\kappa_0` 常数时成立，而 `V(t)` 满足

$$
\dot V(t)=K^{(3)}_{su(2)}(t)\,V(t).
$$

这一步是严格等式，不是近似。

## 5. 单位四元数表示

现在把 `\mathfrak s\cong su(2)` 识别成纯虚四元数空间 `\operatorname{Im}\mathbb H`。选定四元数单位 `\mathbf i,\mathbf j,\mathbf k`，并把 `J_1,J_2,J_3` 视作一组标准的 `su(2)` 生成元。于是 `V(t)` 可写成一个单位四元数

$$
q(t)=q_0(t)+q_1(t)\mathbf i+q_2(t)\mathbf j+q_3(t)\mathbf k,
\qquad q_0^2+q_1^2+q_2^2+q_3^2=1.
$$

在左乘约定下，可令

$$
\Omega(t)=\omega_1(t)\mathbf i+\omega_3(t)\mathbf k,
$$

其中 `\omega_1,\omega_3` 是 `\alpha,\beta` 在所选归一化下的对应系数。则

$$
\dot q(t)=\frac12\,\Omega(t)\,q(t).
$$

展开成分量就是

$$
\dot q_0=-\frac12\big(\omega_1 q_1+\omega_3 q_3\big),
$$

$$
\dot q_1=\frac12\big(\omega_1 q_0-\omega_3 q_2\big),
$$

$$
\dot q_2=\frac12\big(-\omega_1 q_3+\omega_3 q_1\big),
$$

$$
\dot q_3=\frac12\big(\omega_1 q_2+\omega_3 q_0\big).
$$

因此，Step 3 的 reduced 演化可以写成精确的单位四元数时间有序指数

$$
q(t)=\mathcal T\exp\left(\frac12\int_0^t \Omega(t')\,dt'\right).
$$

## 6. 这一步到底得到什么

这说明 Step 3 的结构并不是“整段都要硬压成单位四元数”，而是准确地分成两层：

- 一个与 `su(2)` 对易的 `u(1)` 相位因子 `e^{\Phi(t)X_{12}}`；
- 一个真正落在 `Sp(1)\cong SU(2)` 里的单位四元数 `q(t)`。

所以 Step 3 是严格可四元数化的，但这个结论只对它自己的 `su(2)` 子块成立，不应外推到 Step 1 或 Step 2。

## 7. 后续如果要继续

如果还想把 `q(t)` 进一步化成更显式的解析表达，那么下一步不是再强行做代数压缩，而是对这组 `su(2)` 方程再用 `Wei–Norman` 分解或 `Magnus` 展开。

对 Step 1 和 Step 2，则继续保留完整的 `Sp(2)` 时间有序指数，不做单位四元数假设。

## 8. 再往前一步：把 Step 3 写成可解的 Euler / Wei-Norman 系统

上面的单位四元数方程已经把 Step 3 降到了 `su(2)`。如果还想继续得到更具体的解析表达，就要区分“固定方向”与“真正时变方向”两种情况。

### 8.1 固定方向时的显式四元数闭式

设在某一时间区间内，`su(2)` 部分可以写成

$$
K_{su(2)}^{(3)}(t)=\lambda(t)\,\hat n,
\qquad \hat n^2=-1,
\qquad \hat n\ \text{为固定的纯虚单位四元数}。
$$

则有严格闭式

$$
q(t)=\exp\left(\frac12\Theta(t)\hat n\right),
\qquad
\Theta(t)=\int_0^t \lambda(t')\,dt'.
$$

展开以后就是

$$
q(t)=\cos\frac{\Theta(t)}{2}+\hat n\,\sin\frac{\Theta(t)}{2}.
$$

这就是单位四元数的标准旋转公式。它是精确等式，不是近似。

### 8.2 当前 Step 3 的真正情况：Euler / Wei-Norman 参数化

对当前这个协议，更自然的做法是把 `SU(2)` 变量写成 Euler 形式

$$
q(t)=e^{a(t)\mathbf i/2}\,e^{b(t)\mathbf j/2}\,e^{c(t)\mathbf k/2}.
$$

若

$$
\dot q(t)=\frac12\big(\omega_1(t)\mathbf i+\omega_2(t)\mathbf j+\omega_3(t)\mathbf k\big)\,q(t),
$$

则三条标量方程等价于

$$
\omega_1=\dot a+\dot c\sin b,
$$

$$
\omega_2=\dot b\cos a-\dot c\sin a\cos b,
$$

$$
\omega_3=\dot b\sin a+\dot c\cos a\cos b.
$$

反过来，解出参数导数就是

$$
\dot b=\omega_2\cos a+\omega_3\sin a,
$$

$$
\dot c=\frac{-\omega_2\sin a+\omega_3\cos a}{\cos b},
$$

$$
\dot a=\omega_1-\dot c\sin b.
$$

这就是 Step 3 的精确 `Wei-Norman` / Euler-angle 版本。它把一个 `SU(2)` 的时间有序指数，降成了三个耦合的实标量方程。

### 8.3 对当前协议的简化

在前文采用的 Step 3 识别里，`su(2)` 部分只有两条独立通道，因此通常可取 `\omega_2(t)=0`。此时上面的系统进一步化成

$$
\dot b=\omega_3\sin a,
$$

$$
\dot c=\omega_3\frac{\cos a}{\cos b},
$$

$$
\dot a=\omega_1-\omega_3\cos a\tan b.
$$

这三个方程已经是当前 Step 3 能得到的最直接解析化结果：它比单纯写 `\mathcal T\exp` 更具体，但又没有额外假设方向恒定。

### 8.4 这一层的意义

所以，后续如果要继续推进，路线就很清楚了：

- 若某一段方向固定，就直接用单位四元数闭式；
- 若方向仍随时间变化，就保留这组 `Wei-Norman` 方程，或再转入 `Magnus` 展开；
- 对 Step 1 和 Step 2，则仍然维持完整的 `Sp(2)` 时间有序指数。

## 9. 论文情形下的具体演化形式

如果你要的是“论文这套门控协议到底长什么样”，那么目前能写到的最具体层次就是下面这个分段演化。

### 9.1 三段总演化

设三个切换时刻为 `T_1<T_2<T_3`，并记初值为 `\Psi(0)`。则全协议的演化算符可写成

$$
U(t)=
\begin{cases}
U^{(1)}(t,0), & 0\le t\le T_1,\\[4pt]
U^{(2)}(t,T_1)\,U^{(1)}(T_1,0), & T_1\le t\le T_2,\\[4pt]
U^{(3)}(t,T_2)\,U^{(2)}(T_2,T_1)\,U^{(1)}(T_1,0), & T_2\le t\le T_3.
\end{cases}
$$

其中每一段都是对应的时间有序指数

$$
U^{(s)}(t_b,t_a)=\mathcal T\exp\left(\int_{t_a}^{t_b}K^{(s)}(t')\,dt'\right),
\qquad s=1,2,3.
$$

因此状态矢量的具体形式就是

$$
\Psi(t)=U(t)\Psi(0).
$$

这就是 full-space 下最直接、也最严格的“具体演化形式”。

### 9.2 Step 3 的显式因子分解

对 Step 3，前文已经得到

$$
K^{(3)}(t)=\kappa_0X_{12}+K^{(3)}_{su(2)}(t),
\qquad [X_{12},K^{(3)}_{su(2)}(t)]=0.
$$

所以 Step 3 的传播子可以严格写成

$$
U^{(3)}(t,T_2)=e^{\Phi(t,T_2)X_{12}}\,\widetilde U^{(3)}(t,T_2),
$$

其中

$$
\Phi(t,T_2)=\int_{T_2}^{t}\kappa_0\,dt'.
$$

若 `\kappa_0` 是常数，则直接化成 `\Phi(t,T_2)=\kappa_0(t-T_2)`。

而 `\widetilde U^{(3)}` 就是那个真正的 `SU(2)` / 单位四元数部分。

### 9.3 Step 3 的四元数演化

把 `\widetilde U^{(3)}` 识别成单位四元数 `q(t)`，则有

$$
q(t)=\mathcal T\exp\left(\frac12\int_{T_2}^{t}\Omega(t')\,dt'\right),
\qquad
\Omega(t)=\omega_1(t)\mathbf i+\omega_3(t)\mathbf k.
$$

对应的分量方程是

$$
\dot q_0=-\frac12\big(\omega_1 q_1+\omega_3 q_3\big),
$$

$$
\dot q_1=\frac12\big(\omega_1 q_0-\omega_3 q_2\big),
$$

$$
\dot q_2=\frac12\big(-\omega_1 q_3+\omega_3 q_1\big),
$$

$$
\dot q_3=\frac12\big(\omega_1 q_2+\omega_3 q_0\big).
$$

如果在某个子区间里方向固定，那么还可以进一步直接写成

$$
q(t)=\cos\frac{\Theta(t)}{2}+\hat n\sin\frac{\Theta(t)}{2},
\qquad
\Theta(t)=\int_{T_2}^{t}\lambda(t')\,dt'.
$$

### 9.4 现在已经具体到什么程度

所以，对论文的情形来说，“具体的演化形式”目前已经明确到了下面这个层次：

- Step 1 和 Step 2：保留为分段 `Sp(2)` 时间有序指数；
- Step 3：严格分解为一个对易的 `u(1)` 相位因子和一个 `SU(2)` 单位四元数演化；
- 若 Step 3 的方向进一步固定，就能直接写成余弦-正弦的闭式。

如果你下一步要我继续，我建议直接把 Step 3 的 `\omega_1(t),\omega_3(t)` 依据论文的具体门控函数代入，然后把这组方程再往可积分形式整理一次。