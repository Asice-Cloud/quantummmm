解析推导：从二次 Majorana 到缠绕算符

1) Majorana 双线性与幺正变换

若哈密顿仅含二次 Majorana 项，可写为
$$
H_{\rm quad}=\frac{i}{2}\sum_{a<b}A_{ab}\,\gamma_a\gamma_b,
$$
其中 $A$ 为实反对称矩阵。对应的幺正算符
$$
U=e^{iH_{\rm quad}}
$$
在 Majorana 算符上的共轭作用是线性的：
$$
U^{\dagger}\,\gamma\,U = O\,\gamma,\qquad O=\exp(A)\in SO(2N).
$$

2) 理想两 Majorana 的缠绕算符

交换一对 Majorana 的理想算符为
$$
B(\gamma_p,\gamma_q)=\exp\Big(\frac{\pi}{4}\,\gamma_p\gamma_q\Big).
$$
它在 (p,q) 平面上产生角度 $\pi/2$ 的正交旋转：
$$
B^{\dagger}\gamma_p B=\gamma_q,\qquad B^{\dagger}\gamma_q B=-\gamma_p.
$$ 
因此在二次情形中，充要条件是 $A$ 在 Majorana 基下只在 p–q 平面有非零塊，且该塊指数化后给出角 $\pi/2$。

3) 从自旋/Pauli 系数到 A 的线性关系（示例）

文档中常用映射（以相邻两站点为例）给出：
$$
\sigma^x\otimes\sigma^x\mapsto \tfrac12 i(\gamma_{2i-1}\gamma_{2(i+1)}+\gamma_{2i}\gamma_{2(i+1)-1}),
$$
$$
\sigma^y\otimes\sigma^y\mapsto \tfrac12 i(-\gamma_{2i-1}\gamma_{2(i+1)}+\gamma_{2i}\gamma_{2(i+1)-1}).
$$

因此任意二次部分的系数 $c_{\mu\nu}$ 是 A 的线性组合，具体对偶索引的系数可直接读出并填入 $A_{ab}$。

γ1…γ6 的具体例子（3 个站点）

约定：站点 1→(γ1,γ2)，站点 2→(γ3,γ4)，站点 3→(γ5,γ6)。

(1) 对于仅含 XX/YY 最近邻耦合的格点模型

（键 1–2 的系数记为 $J_{x1},J_{y1}$；键 2–3 为 $J_{x2},J_{y2}$），最近邻贡献在 Majorana 表示为：
$$
  \begin{aligned}
  H_{1-2} &= \tfrac{i}{2}\big[ (J_{x1}-J_{y1})\,\gamma_1\gamma_4 + (J_{x1}+J_{y1})\,\gamma_2\gamma_3\big],\\
  H_{2-3} &= \tfrac{i}{2}\big[ (J_{x2}-J_{y2})\,\gamma_3\gamma_6 + (J_{x2}+J_{y2})\,\gamma_4\gamma_5\big].
  \end{aligned}
$$

(2) 对应 6×6 的 A 矩阵（只列非零上三角分量）：
  
  - $A_{1,4}=J_{x1}-J_{y1}$
  - $A_{2,3}=J_{x1}+J_{y1}$
  - $A_{3,6}=J_{x2}-J_{y2}$
  - $A_{4,5}=J_{x2}+J_{y2}$

其余 $A_{ab}=0$，并且 $A_{ba}=-A_{ab}$。

(3) 实现理想缠绕的操作条件（以实现 $B(\gamma_1,\gamma_4)$ 为例）：

要使 $U=e^{iH_{\rm quad}}$ 等价于 $B(\gamma_1,\gamma_4)$，需要在 A 中只保留 1–4 平面的块并把该块的指数化角调整为 $\pi/2$。具体操作为：

1. 令 $A_{2,3}=J_{x1}+J_{y1}=0$，即 $J_{x1}=-J_{y1}$，以关闭 γ2–γ3 平面旋转；

2. 令键 2–3 的耦合关闭或为 0（$J_{x2}=J_{y2}=0$），以清除其它块；

3. 调整 $A_{1,4}=J_{x1}-J_{y1}$ 与操作时间 $t$ 使得该块的有效角满足（按 $H=(i/2)A\gamma\gamma$ 的约定）指数化结果等于旋转角 $\pi/2$，例如满足 $(J_{x1}-J_{y1})\,t=\pi/2$（具体因子依你的时间尺度与归一化可微调）。

(4)
若 A 有多块（比如同时存在 $A_{1,4}$ 与 $A_{2,3}$），则 $O=\exp(A)$ 为这些平面旋转的复合，通常不是单一对的理想交换，但可以通过时间分段或把不想要的块置零来近似分解为理想的单次交换。

备注（交互与串项）
- 本推导假定无串算符与无四费米交互（即 JW 后只保留二次项）。若模型含 $\sigma^z\sigma^z$ 等产生四费米或含单侧 x/y 导致串项，则需把 $R$ 投影到低能子空间或做数值验证：计算 $U_{\rm eff}=P e^{iH_P} P$（P 为低能子空间投影），再与理想 $B$ 做 operator‑fidelity / state‑overlap 比较