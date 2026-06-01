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

## 9. 把论文的分段门控代入（符号化的最终 ODE）

根据论文中给出的门控 schedule，令

$$
x:=\frac{\pi t}{\tau},\qquad f_+(x):=\frac{1+\cos x}{2},\qquad f_-(x):=\frac{1-\cos x}{2}.
$$

Step 3 的原始物理量为

$$
|t_3^{(3)}(t)|=t_c\,f_+(x),\qquad E_d^{(3)}(t)=E_0\,f_-(x).
$$

在之前的约定里我们用

$$
K^{(3)}(t)=\kappa_0 X_{12}+\alpha(t)X_{3a}+\beta(t)X_{ab}.
$$

因此可取符号化替换（归一化常数以所选基底为准）：

$$
\alpha(t)=C_{t3}\,t_c\,f_+(x),\qquad
\beta(t)=C_{Ed}\,E_0\,f_-(x),\qquad
\kappa_0=C_{E1}\,E_1,
$$

其中 $C_{t3},C_{Ed},C_{E1}$ 是基底归一化常数（例如在 4.2 节的矩阵记号下常为 2 或 −2），可以根据你具体采用的 `so(5)\leftrightarrow sp(2)` 同构固定。

把 $\alpha,\beta$ 映到四元数角速度分量，我们写

$$
\Omega(t)=\omega_1(t)\mathbf i+\omega_3(t)\mathbf k,\qquad \omega_1(t)=\gamma_1\,\alpha(t),\quad \omega_3(t)=\gamma_3\,\beta(t),
$$

其中 $\gamma_1,\gamma_3$ 也是由基底選擇導出的常数（可取 1 如果已在 $C_{t3},C_{Ed}$ 中包含比例因子）。

代入第 8 节的 Euler / Wei–Norman 公式（取 $\omega_2\equiv0$），得到 Step 3 的最终三标量 ODE：

$$
\dot b(t)=\omega_3(t)\sin a(t),
$$

$$
\dot c(t)=\frac{\omega_3(t)\cos a(t)}{\cos b(t)},
$$

$$
\dot a(t)=\omega_1(t)-\omega_3(t)\cos a(t)\tan b(t).
$$

把 $\omega_{1,3}(t)$ 展开回论文变量就是：

$$
\omega_1(t)=\gamma_1 C_{t3} t_c f_+(\tfrac{\pi t}{\tau}),\qquad
\omega_3(t)=\gamma_3 C_{Ed} E_0 f_-(\tfrac{\pi t}{\tau}).
$$

初始条件一般取为 $a(0)=b(0)=c(0)=0$（即 $q(0)=1$），而 $\Phi(t)=C_{E1}E_1 t$ 給出可分离的 $u(1)$ 相位因子。

备注：若你愿意我可以把上面常数取为具体值（例如 $C_{t3}=2,C_{Ed}=2,C_{E1}=2,\gamma_1=\gamma_3=1$），并把完整的 ODE 以符号或数值形式写入 `derivation_1.md`，或直接生成用于数值积分的 Python 脚本。