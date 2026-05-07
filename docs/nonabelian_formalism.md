# 严格的非阿贝尔性与几何/动力学分离：形式化量化方法

整理自项目笔记（ybe2, ybe22, ybe223, ybe224）——给出可直接数值实现的数学公式与处方，用以量化非阿贝尔 Berry 连接、曲率、Wilson holonomy 以及区分 MZM 与 ABS 的可计算指标。

---

## 1. 对象与基本定义

- 给定参数化的 Yang–Baxter 算子 $R(u)$（可推广为 $R(u,\lambda)$），定义生成元：

$$
H(u)=i\,\partial_u R(u)\,R(u)^{-1}.
$$

- 在 $N$ 维希尔伯特空间中选取目标子空间的正交列基矩阵 $W(u)\in\mathbb C^{N\times m}$，满足 $W^\dagger W = I_m$，其中 $m$ 为子空间维数（通常 $m=2$）。

## 2. Berry 连接与曲率（非阿贝尔）

- 定义参数簇 $\{\lambda^\mu\}$（例如 $\lambda^1=u,\ \lambda^2=\delta$），矩阵值 Berry 连接为：

$$
A_\mu(u,\lambda)=i\,W(u,\lambda)^\dagger\partial_\mu W(u,\lambda),
\qquad A_\mu\in\mathfrak{u}(m).
$$

- 曲率（场强）：

$$
F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu + [A_\mu,A_\nu].
$$

- Gauge 变换：若 $W\mapsto W G$（$G\in U(m)$），则

$$
A_\mu\mapsto G^\dagger A_\mu G + i G^\dagger \partial_\mu G,
\qquad F_{\mu\nu}\mapsto G^\dagger F_{\mu\nu} G.
$$

因此 $\mathrm{Tr}(F^\dagger F)$ 与 Wilson 特征值为 gauge‑不变诊断量。

## 3. 可计算的标量指标

- 局域非阿贝尔强度（Frobenius 范数）：

$$
\mathcal F(u,\lambda)=\|F_{\mu\nu}(u,\lambda)\|_F=\sqrt{\mathrm{Tr}[F^\dagger F]}.
$$

- 参数面 $S$ 上的总体量（归一化）：

$$
\mathrm{NA}=\frac{1}{\mathrm{Area}(S)}\int_S \mathrm{Tr}[F^\dagger F] \, dS.
$$

- 离散 Wilson 小胞近似：对网格小胞面积 $\Delta\lambda^1\Delta\lambda^2$，用相邻点重叠矩阵 $M_{ij}=W(i)^\dagger W(j)$ 计算小胞 holonomy

$$
U_{\square}=M_{12}M_{23}M_{34}M_{41}\approx\exp(F_{u\lambda}\,\Delta u\Delta\lambda).
$$

诊断量：

$$
D_{\square}=\frac{\|U_{\square}-I_m\|_F}{\Delta u\Delta\lambda}.
$$

以及 Wilson 的特征值集合 $\{e^{i\varphi_j}\}$（gauge‑invariant）。

## 4. 两能级情形：Bloch 旋转與几何相

- 若投影后 $m=2$ 且有效 Hamiltonian 写为

$$
H_{\rm eff}(t)=\vec d(t)\cdot\vec\sigma,
$$

令 $\hat d=\vec d/|\vec d|$，则总演化可近似视为 SU(2) 旋转：

$$
U=\mathcal T\exp\Big(-i\int H_{\rm eff}(t)\,dt\Big)=\exp\Big(-i\frac{\Phi}{2}\,\hat n\cdot\vec\sigma\Big),
$$

其中累积动力学角 $\Phi = 2\int |\vec d(t)|\,dt$（SU(2)↔SO(3) 的因子）。

- 单圈几何相（上能带）：若 $\hat d=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$，则

$$
\gamma_+ = -\tfrac{1}{2}\Omega,\qquad \Omega=\oint(1-\cos\theta)\,d\phi.
$$

- 绕数（（$d_x,d_y$ 平面））：

$$
W=\frac{1}{2\pi}\oint\frac{d_x\partial_u d_y - d_y\partial_u d_x}{d_x^2+d_y^2}\,du.
$$

## 5. ABS 与 MZM 的可计算判据

- 最小能隙：

$$
\Delta_{\min}=\min_{u\in S} |\vec d(u)|.
$$

若 $\Delta_{\min}\approx 0$ 且轨迹包围原点（$W\ne0$ 或相关 Chern 数非平凡），倾向 MZM；若 $\Delta_{\min}\gg0$ 则为 ABS‑like（动力学主导）。

- 体系尺度检验（链长 $L$）：测定零模分裂 $E_{\rm split}(L)$ 并拟合

$$
E_{\rm split}(L)\sim A e^{-L/\xi} \quad(\text{MZM})
$$

若衰减非指数或收敛慢，倾向 ABS。另可计算 BdG 的 Pfaffian 指标（若有动量表示）：

$$
\nu=\mathrm{sgn}\big(\mathrm{Pf}[B(0)]\,\mathrm{Pf}[B(\pi)]\big).
$$

## 6. 绝热性与几何/动力学分离

- 态间绝热参量：

$$
\epsilon_{mn}(u)=\frac{|A_{mn}(u)|}{|E_n(u)-E_m(u)|}.
$$

若 $\epsilon_{\max}=\max_{u,m\ne n}\epsilon_{mn}(u)\ll1$，则可把演化近似分解为几何因子乘以动力学因子。

- α‑缩放策略（使动力学相退相）：定义

$$
I_0=\int |\vec d(u)|\,du,\qquad \alpha_{\rm target}=\frac{\pi}{I_0}.
$$
若可控制缩放 $\alpha$ 使得满足上式，动力学相在单圈后退化为全局相（近似），从而 $U\approx U_{\rm geom}$。

## 7. 数值处方（伪代码）

```
# Inputs: discrete grid in parameters (u, lambda), target subspace dimension m
for each grid point (u,lambda):
  compute H(u,lambda) (or via finite-difference from R(u,lambda))
  compute eigenpairs (E_k, v_k)
  choose m eigenvectors (e.g. smallest |E| or chosen band) -> form W(u,lambda)
  align gauge by parallel-transport (SVD-based) along grid lines
compute finite-difference derivatives dW/dlambda
compute A_mu = i W^dagger dW
compute F_{mu nu} = d_mu A_nu - d_nu A_mu + [A_mu, A_nu]
compute scalar maps: ||F||_F, D_square, trace(F^dagger F)
integrate over parameter region to get NA
compute Wilson loops over loops of interest and their eigenvalues
compute epsilon_mn = |A_mn|/|E_n - E_m|
end
```

## 8. 实现参考

- 示例脚本： [tools/quantify_nonabelian.py](tools/quantify_nonabelian.py) 与 [tools/verify_from_R.py](tools/verify_from_R.py)。
- 输出建议：`results/nonabelian_grid.npz`（包含 `Fnorm`, `Udev` 等），以及 `results/nonabelian_heatmaps.png`、布洛赫轨迹图和 `results/deriveR_report_delta*.pkl` 等。

## 9. 结论

给出了一套数学上自洽且易于数值实现的指标（$F_{\mu\nu}$、Wilson eigenvalues、winding、$\epsilon_{mn}$、$\Delta_{\min}$、$E_{\rm split}(L)$ 拟合），可用于在你的模型管线中严格量化非阿贝尔性并区分 MZM/ABS。

---

如果你需要，我可以把该 Markdown 转为 HTML 报告并嵌入热图与关键数值表；或直接把数值扩展到包含 $\epsilon_{mn}$ 分布与局部 Wilson 特征值表。 
