# 论文流程的准确复现：从 modified tetron 到 Fig. 1(c) 的五段序列

这份文档只做一件事：按照论文的真实步骤，把“器件路径设计 → Majorana 杂化 → braid → Bloch 球旋转”严格串起来，并把它完整复现到我们的有效模型里。

核心结论是：论文中的流程不是先写 Pauli 门再倒推器件，而是先通过门电压、相位差和耦合强度构造一条时间路径，再把这条路径投影到低能 Majorana 子空间，最后才得到单比特上的 $SU(2)$ 演化。

## 1. 论文要复现的物理对象

论文考虑的是一个 modified tetron 结构。总 Hamiltonian 写成

$$
H_T(t)=\sum_{i=1}^4 H_i(t)+H_S+H_{Tc}(t),
$$

其中：

- $H_i(t)$ 是第 $i$ 条纳米线的 BdG Hamiltonian；
- $H_S$ 是中间 trivial superconductor；
- $H_{Tc}(t)$ 是纳米线与中心区的时变耦合；
- $t_{ic}(t)=g_i(t)t_0$ 由门 $G_i$ 控制。

论文的关键点在于：

1. 通过门电压把某些 Majorana 对先杂化；
2. 通过依次开关 $G_1,G_2,G_3,G_4$ 实现空间交换；
3. 保持有限相位差 $\phi_0\neq 0$，让拓扑能隙在交换期间不闭合；
4. 在 ABS 存在时，把低能态看成一对弱耦合 Majorana，再把杂化和 braid 组合成一个可控的单比特旋转。

这就是论文想证明的完整流程。

## 2. modified tetron 的交换步骤

论文在 braiding 前取如下初态：四条臂 $W_1\sim W_4$ 都处于拓扑非平庸相；$G_1,G_2$ 打开，$G_3,G_4$ 关闭，因此有三对 Majorana

$$
(\gamma_{2j-1},\gamma_{2j}),\qquad j=1,2,3,
$$

分别局域在三个分段的端点上。

随后 braid $\gamma_2$ 和 $\gamma_3$ 的三步过程是：

1. 关闭 $G_1$，打开 $G_3$，把 $\gamma_2$ 传送到 $W_3$；
2. 关闭 $G_2$，打开 $G_1$，把 $\gamma_3$ 传送到 $\gamma_2$ 的原位置；
3. 关闭 $G_3$，打开 $G_2$，完成 $\gamma_2$ 与 $\gamma_3$ 的空间交换。

这一过程的时间演化由

$$
U(t)=\hat T\exp\!\left(i\int_0^t H(\tau)\,d\tau\right)
$$

给出。为了让交换结果稳定有效，论文要求 $\phi_0\neq 0$，这样在交换过程中能隙始终打开，不会因为 domain wall 再生成一对零模而破坏 braid。

## 3. ABS 在论文里的精确角色

论文把 ABS 看成一对弱耦合 Majorana。也就是说，ABS 并不是一个抽象的“杂散低能态”，而是可以在低能有效理论中写成两个 Majorana 生成元的杂化。

把这对 Majorana 记作 $\gamma_1,\gamma_2$，那么 ABS 的局域杂化就是

$$
H_{\rm ABS}=i\,\varepsilon_{12}(t)\,\gamma_1\gamma_2.
$$

而另外两条仍然是真正用于 braid 的拓扑 Majorana，记为 $\gamma_3,\gamma_4$。

论文图 1(c) 的意思就是：先让 $\gamma_1,\gamma_2$ 杂化，再 braid $\gamma_2,\gamma_3$，再让 $\gamma_1,\gamma_3$ 杂化，再做一次反向 braid，最后再回到 $\gamma_1,\gamma_2$ 的杂化。这个顺序不是随意挑的，而是论文用来展示“hybridization + braiding 可以合成任意 Bloch 旋转”的具体协议。

## 4. Fig. 1(c) 的五段操作，按论文原顺序写出

图中的五段操作可以直接写成五个时间演化算符的乘积。按时间顺序从 section 1 到 section 5：

$$
U_1=\exp\!\left(\frac{\theta_1}{2}\,\gamma_1\gamma_2\right),
$$

$$
U_2=\exp\!\left(-\frac{\pi}{4}\,\gamma_2\gamma_3\right),
$$

$$
U_3=\exp\!\left(\frac{\theta_2}{2}\,\gamma_1\gamma_3\right),
$$

$$
U_4=\exp\!\left(-\frac{\pi}{4}\,\gamma_3\gamma_2\right),
$$

$$
U_5=\exp\!\left(-\frac{\theta_3}{2}\,\gamma_1\gamma_2\right).
$$

因此总演化是

$$
U_{\rm tot}=U_5U_4U_3U_2U_1.
$$

这就是论文图 1(c) 的严格数学表达。

需要强调的是：

- $U_1,U_3,U_5$ 是 hybridization 段，对应可调角度 $\theta_1,\theta_2,\theta_3$；
- $U_2,U_4$ 是 braid 段，对应固定的 $\pi/4$ 交换角；
- $U_2$ 和 $U_4$ 的方向相反，因此它们不是同一个门的重复，而是两个相反向的 braid pulse。

这五段拼起来，才是论文要复现的“准确流程”。

## 5. 低能逻辑比特与 Majorana 生成元

为了把上面的五段序列变成单比特门，需要把四个 Majorana 的偶宇称子空间投影成一个逻辑 qubit。

取逻辑子空间后，定义三组双线性生成元：

$$
\Sigma_x=i\gamma_2\gamma_3,\qquad
\Sigma_y=i\gamma_1\gamma_3,\qquad
\Sigma_z=i\gamma_1\gamma_2.
$$

它们满足 Pauli 代数

$$
\Sigma_a\Sigma_b=\delta_{ab}I+i\epsilon_{abc}\Sigma_c.
$$

因此存在一个固定的逻辑基变换 $W$，使得

$$
W\Sigma_aW^\dagger=\sigma_a,\qquad a=x,y,z.
$$

在这个逻辑基下，每个 Majorana 杂化或 braid 都变成一个 Pauli 旋转。

具体地：

$$
i\gamma_1\gamma_2 \mapsto \sigma_z,
$$

$$
i\gamma_2\gamma_3 \mapsto \sigma_x,
$$

$$
i\gamma_1\gamma_3 \mapsto \sigma_y.
$$

所以图 1(c) 的五段流程在逻辑比特上就是一串绕不同轴的旋转。它的本质不是某个单独 braid，而是“hybridization + braid + hybridization + braid + hybridization”的组合门。

## 6. 把 Fig. 1(c) 变成有效单比特 Hamiltonian

现在把论文的器件控制写成时间路径

$$
p(t)=\bigl(g_1(t),g_2(t),g_3(t),g_4(t),V_x(t),\phi_0(t),V_D(t)\bigr).
$$

这些参数在低能 Majorana 描述下对应耦合系数

$$
\varepsilon_{12}(t),\qquad \varepsilon_{23}(t),\qquad \varepsilon_{13}(t),
$$

以及标量偏置 $\mu(t)$。于是 Majorana Hamiltonian 可以写成

$$
H_M(t)=i\varepsilon_{12}(t)\gamma_1\gamma_2+i\varepsilon_{23}(t)\gamma_2\gamma_3+i\varepsilon_{13}(t)\gamma_1\gamma_3+\mu(t)I.
$$

投影到逻辑 qubit 后得到

$$
H_{\rm eff}(t)=d_0(t)I+d_x(t)\sigma_x+d_y(t)\sigma_y+d_z(t)\sigma_z,
$$

其中

$$
d_x(t)=-\varepsilon_{23}(t),\qquad d_y(t)=-\varepsilon_{13}(t),\qquad d_z(t)=-\varepsilon_{12}(t),\qquad d_0(t)=\mu(t).
$$

这一步是论文流程中最关键的“桥”：器件参数不是直接给出门，而是先给出 Majorana 双线性，再投影成 Pauli 生成元。

## 7. 按论文步骤实现五段控制

为了把图 1(c) 的五段操作写成显式的时间路径，可以把时间分成五个窗口：

$$
f_j(t)=\chi_{[T_{j-1},T_j]}(t),\qquad j=1,\dots,5,
$$

其中 $T_0=0$，$T_5=T_{\rm tot}$。

然后令

$$
\varepsilon_{12}(t)=\varepsilon_{12}^{(1)}f_1(t)+\varepsilon_{12}^{(5)}f_5(t),
$$

$$
\varepsilon_{23}(t)=\varepsilon_{23}^{(2)}f_2(t),
$$

$$
\varepsilon_{13}(t)=\varepsilon_{13}^{(3)}f_3(t),
$$

并在第 4 段取与第 2 段相反方向的 braid：

$$
\varepsilon_{23}^{(4)}=-\varepsilon_{23}^{(2)}.
$$

这样就得到与图 1(c) 完全对应的时变 Hamiltonian：

$$
\begin{aligned}
H_{\rm eff}^{(1)} &= d_0 I + d_z\sigma_z,\\
H_{\rm eff}^{(2)} &= d_0 I + d_x\sigma_x,\\
H_{\rm eff}^{(3)} &= d_0 I + d_y\sigma_y,\\
H_{\rm eff}^{(4)} &= d_0 I - d_x\sigma_x,\\
H_{\rm eff}^{(5)} &= d_0 I - d_z\sigma_z.
\end{aligned}
$$

于是总演化变成

$$
U_{\rm tot}=e^{-i\phi_5}e^{-i\frac{\theta_3}{2}\sigma_z}
e^{-i\frac{\pi}{4}\sigma_x}
e^{-i\phi_3}e^{-i\frac{\theta_2}{2}\sigma_y}
e^{+i\frac{\pi}{4}\sigma_x}
e^{-i\phi_1}e^{-i\frac{\theta_1}{2}\sigma_z},
$$

其中 $\phi_j$ 是各段的全局相位，不影响 Bloch 球上的旋转。

这就是论文图 1(c) 在我们有效模型中的精确复现。

## 8. 为什么这一定是一个 Bloch 旋转

因为 $\sigma_x,\sigma_y,\sigma_z$ 生成 $\mathfrak{su}(2)$。只要路径中允许在至少两个不共线方向之间切换，李代数闭包就已经足够生成整个单比特代数。

最直接的证明是：

$$
[\sigma_x,\sigma_y]=2i\sigma_z,
$$

所以 $i\sigma_x,i\sigma_y,i\sigma_z$ 张成三维空间，正好是 $\mathfrak{su}(2)$。

因此，图 1(c) 里的五段序列虽然写成 Majorana 双线性的乘积，但在逻辑 qubit 上就是一个标准的 $SU(2)$ Euler 型控制协议。

换句话说：

$$
U_{\rm tot}\in SU(2),
$$

而其 Bloch 球意义就是“绕某个轴旋转某个角度”。

## 9. 论文里“动态演化可被消去”的严格含义

论文强调：如果把各段参数分别调好，就可以让动力学相退化成全局相，只留下 braid 的几何部分。

在我们的写法里，这意味着可以通过控制 $\theta_1,\theta_2,\theta_3$ 和每段的驻留时间，使得

$$
\int d_0(t)\,dt \in \pi\mathbb Z,
$$

或者更一般地，让动力学因子在逻辑子空间上化成 $\pm I$。

一旦做到这一点，总演化就只剩下由 $\sigma_x,\sigma_y,\sigma_z$ 组成的几何门，而这正是论文所说的“non-Abelian braiding statistics, independent of the braiding time, retrieves”。

## 10. 论文图 1(c) 与我们模型的逐段对应表

可以把图 1(c) 直接翻译成下面的对照：

| 论文段落 | 物理操作 | Majorana 形式 | 逻辑 qubit 形式 |
| --- | --- | --- | --- |
| section 1 | hybridization | $\exp(\frac{\theta_1}{2}\gamma_1\gamma_2)$ | $z$ 轴旋转 |
| section 2 | braid | $\exp(-\frac{\pi}{4}\gamma_2\gamma_3)$ | $x$ 轴旋转 |
| section 3 | hybridization | $\exp(\frac{\theta_2}{2}\gamma_1\gamma_3)$ | $y$ 轴旋转 |
| section 4 | braid | $\exp(-\frac{\pi}{4}\gamma_3\gamma_2)$ | 反向 $x$ 轴旋转 |
| section 5 | hybridization | $\exp(-\frac{\theta_3}{2}\gamma_1\gamma_2)$ | 反向 $z$ 轴旋转 |

这个表就是论文步骤的最简洁版本。

## 11. 结论：我们到底复现了什么

这份复现不是只复现了“某个结论”，而是把论文的流程顺序完整保留下来：

1. 先由 modified tetron 和门控 $G_1\sim G_4$ 定义可控的 Majorana 网络；
2. 再把 ABS 识别为一对弱耦合 Majorana；
3. 然后按 Fig. 1(c) 的五段顺序进行 hybridization / braiding / hybridization / braiding / hybridization；
4. 再把这些 Majorana 双线性投影到逻辑 qubit；
5. 最终得到单比特上的 $SU(2)$ Bloch 旋转。

这就是论文流程下的准确复现。

## 12. 和我们的有效模型的最终闭环

如果只保留我们模型里最抽象但也最关键的一层，那么整条链可以压成：

$$
p(t)
\longrightarrow
\{\varepsilon_{12}(t),\varepsilon_{23}(t),\varepsilon_{13}(t),\mu(t)\}
\longrightarrow
H_M(t)
\longrightarrow
H_{\rm eff}(t)
\longrightarrow
U_{\rm tot}.
$$

其中：

- $H_M(t)$ 负责把器件路径翻译成 Majorana 双线性；
- $H_{\rm eff}(t)$ 负责把 Majorana 双线性翻译成 Pauli 生成元；
- $U_{\rm tot}$ 负责把整条路径翻译成单比特门。

所以，论文的“杂化 + braid”不是两套独立的现象，而是同一条控制链在不同层级上的表现。

这也是我们模型能够准确复现论文的原因。