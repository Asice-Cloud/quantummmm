# PRQ 传播模型与 Kitaev 链严格推导

# 1. PRQ 传播模型

定义 Hilbert 空间：

\[
\mathcal H=\mathcal H_P\oplus\mathcal H_Q
\]

投影：

\[
P^2=P,\qquad Q^2=Q,\qquad P+Q=I
\]

传播算符：

\[
R=(\omega-H)^{-1}
\]

研究对象：

\[
P_iRQ_j
\]

以及：

\[
P_iRQ_jRP_k
\]

---

# 2. Hamiltonian 分块

写：

\[
H=
\begin{pmatrix}
H_{PP}&H_{PQ}\\
H_{QP}&H_{QQ}
\end{pmatrix}
\]

Green 函数：

\[
G(\omega)=(\omega-H)^{-1}
\]

---

# 3. Schur 补

利用 block inversion：

\[
G_{PP}
=
P(\omega-H)^{-1}P
\]

得到：

\[
\boxed{
G_{PP}
=
(\omega-H_{PP}-\Sigma(\omega))^{-1}
}
\]

其中：

\[
\boxed{
\Sigma(\omega)
=
H_{PQ}(\omega-H_{QQ})^{-1}H_{QP}
}
\]

---

# 4. PRQ 有效传播核

定义：

\[
\boxed{
K_{ik}(\omega)
=
P_iRQ_jRP_k
}
\]

在 resolvent 表示下：

\[
\boxed{
K_{ik}(\omega)
=
H_{P_iQ_j}
(\omega-H_Q)^{-1}
H_{Q_jP_k}
}
\]

---

\[
h^{eff}(\omega)=H_{PP}+\Sigma(\omega),\qquad
\Sigma(\omega)=H_{PQ}(\omega-H_{QQ})^{-1}H_{QP}
\]

而你写的传播核
\[
K_{ik}(\omega)=H_{P_iQ_j}(\omega-H_Q)^{-1}H_{Q_jP_k}
\]
就是 (\Sigma) 在 (P_i\to P_k) 通道上的块/矩阵元（严格地还应对 (j) 求和）。

更规范写法建议是：
\[
\Sigma_{ik}(\omega)=\sum_j H_{P_iQ_j}(\omega-H_{Q})^{-1}H_{Q_jP_k}
\]
\[
h^{eff}*{ik}(\omega)=(H*{PP})*{ik}+\Sigma*{ik}(\omega)
\]

如果每个 (P_i) 还有内部自由度（Nambu/Majorana 分量），再写成：
\[
\big(h^{eff}_{ik}\big)_{ab}=(H_{PP})_{ik,ab}+\big(\Sigma_{ik}\big)_{ab}
\]
然后再对这个 (4*4)（或对应维度）块做 Pauli 展开：
\[
h^{eff}*{ik}=\sum*{\mu,\nu} c^{(ik)}_{\mu\nu},\sigma^\mu\otimes\sigma^\nu.
\]






# 5. Dyson 演化

定义 effective Hamiltonian：

\[
\boxed{
h^{\mathrm{eff}}
=
\sum_{\alpha\beta}
h_{\alpha\beta}^{\mathrm{eff}}
\sigma^\alpha\otimes\sigma^\beta
}
\]

其中：

\[
\boxed{
h_{\alpha\beta}^{\mathrm{eff}}
=
\langle\alpha_i|
H_{PQ}(\omega-H_Q)^{-1}H_{QP}
|\beta_j\rangle
}
\]

---

路径演化：

\[
\boxed{
R(u)
=
\mathcal T
\exp
\left(
-i\int_0^u h^{\mathrm{eff}}(s)\,ds
\right)
}
\]

Dyson 展开：

\[
R(u)=I+R^{(1)}+R^{(2)}+\cdots
\]

一阶：

\[
R^{(1)}
=
-i\int_0^u ds\,h^{\mathrm{eff}}(s)
\]

二阶：

\[
R^{(2)}
=
(-i)^2
\int_0^u ds_1
\int_0^{s_1}ds_2
h^{\mathrm{eff}}(s_1)
h^{\mathrm{eff}}(s_2)
\]

---

# 6. Yang–Baxter deviation

定义：

\[
\Delta
=
R_{12}(u)
R_{23}(u+v)
R_{12}(v)
-
R_{23}(v)
R_{12}(u+v)
R_{23}(u)
\]

最低非平凡阶：

\[
\boxed{
\Delta^{(3)}
\sim
[h_{12}^{\mathrm{eff}},h_{23}^{\mathrm{eff}}]
}
\]

其中：

\[
\boxed{
h_{12}^{\mathrm{eff}}
=
P_1RQ_1RP_2
}
\]

\[
\boxed{
h_{23}^{\mathrm{eff}}
=
P_2RQ_2RP_3
}
\]

---

# 7. Pauli 非阿贝尔结构

使用：

\[
[\sigma^a,\sigma^b]
=
2i\epsilon_{abc}\sigma^c
\]

得到：

\[
\boxed{
[h_{12}^{\mathrm{eff}},h_{23}^{\mathrm{eff}}]
=
2i
\sum
h_{\alpha\beta}^{\mathrm{eff}}
h_{\mu\nu}^{\mathrm{eff}}
\epsilon_{\beta\mu\gamma}
\sigma^\alpha\otimes\sigma^\gamma\otimes\sigma^\nu
}
\]

---

# 8. Frobenius 非阿贝尔测度

定义：

\[
\boxed{
\mathcal N
=
\sqrt{
\mathrm{Tr}(\Delta^\dagger\Delta)
}
}
\]

得到：

\[
\boxed{
\mathcal N
\sim
\sqrt{
\sum
|h_{\alpha\beta}^{\mathrm{eff}}|^2
|h_{\mu\nu}^{\mathrm{eff}}|^2
\epsilon_{\beta\mu\gamma}^2
}
}
\]

---

# 9. Kitaev 链

标准 Hamiltonian：

\[
H
=
-\mu\sum_i c_i^\dagger c_i
-
t\sum_i(c_i^\dagger c_{i+1}+h.c.)
+
\Delta\sum_i(c_ic_{i+1}+h.c.)
\]

---

# 10. Majorana 分解

定义：

\[
\gamma_{iA}=c_i+c_i^\dagger
\]

\[
\gamma_{iB}=-i(c_i-c_i^\dagger)
\]

满足：

\[
\{\gamma_{i\alpha},\gamma_{j\beta}\}
=
2\delta_{ij}\delta_{\alpha\beta}
\]

---

Hamiltonian：

\[
H
=
\frac{i}{2}\sum_i
\Big[
-\mu\gamma_{iA}\gamma_{iB}
+
(t+\Delta)\gamma_{iB}\gamma_{i+1,A}
+
(-t+\Delta)\gamma_{iA}\gamma_{i+1,B}
\Big]
\]

---

取：

\[
t=\Delta
\]

得到：

\[
\boxed{
H
=
\frac{i}{2}
\sum_i
\Big[
-\mu\gamma_{iA}\gamma_{iB}
+
2t\gamma_{iB}\gamma_{i+1,A}
\Big]
}
\]

---

# 11. Pi 与 Qi

定义：

\[
\boxed{
P_i=\mathrm{span}\{\gamma_{iB}\}
}
\]

\[
\boxed{
Q_i=\mathrm{span}\{\gamma_{iA}\}
}
\]

---

# 12. Kitaev 链中的传播

Hamiltonian 耦合：

\[
P_i\leftrightarrow Q_i
\]

以及：

\[
P_i\leftrightarrow Q_{i+1}
\]

因此：

\[
\boxed{
P_i
\to
Q_{i+1}
\to
P_{i+1}
}
\]

---

## 附录：具体验证流程与数值结果（追加）

下面将本文件中理论推导与数值验证串接，记录可复现的流程、脚本位置与代表性输出。

验证目标：将 PRQ（Pi→Qi→Pi+1）传播核与 Schur 补框架嵌入到 Kitaev 链的局域 3 站点 BdG 子块，数值检验 N 的相位干涉公式
$$\mathcal{N}=\|\Delta\|_F\approx 2A\sin(\theta/2).$$

实现要点：
- 使用脚本 `quantity/validate_propagation_in_kitaev.py`（已加入仓库）完成端到端验证，脚本功能：
	- 从 `tools/embed_kitaev.py` 构建 L=40 的 Kitaev BdG，截取中心连续三站点得到 6×6 子块；
	- 在该子块以 P={site0,site2}, Q={site1} 计算 BdG 形式的 Schur 补 Σ(ω) 并提取 K=‖Σ‖_F；
	- 以子块 H_{12}, H_{23} 分别构造局域演化算子 R_{12}, R_{23}（expm(-iHτ)），嵌入 3 站点空间计算 LHS/RHS、Δ 与 N_actual；
	- 计算矩阵夹角 θ 与理论 N_theory = 2A_mean sin(θ/2)，保存 CSV 与图像。

主要输出（已生成）：
- `quantity/kitaev_validation/kitaev_propagation_validation.csv`（输出列：E1,K,A_lhs,A_rhs,theta_rad,N_theory,N_actual,rel_err）
- `quantity/kitaev_validation/kitaev_N_compare.png`（N_theory vs N_actual）
- `quantity/kitaev_validation/kitaev_theta.png`（θ vs E1）

CSV 摘录（代表点）：

```
E1,K,A_lhs,A_rhs,theta_rad,N_theory,N_actual,rel_err
0.001,2180.0,2.449489742783178,2.449489742783178,0.10240731027302705,0.2507360586799699,0.250736058679973,1.24e-14
0.09526315789473684,32.361048992821914,2.449489742783178,2.449489742783178,0.10187686096933511,0.24943842267686306,0.24943842267686725,1.68e-14
0.2,15.414735146881661,2.449489742783178,2.449489742783178,0.10008896027914195,0.24506455968431828,0.24506455968431926,3.96e-15
```

说明：K 的绝对量纲在 BdG 局域构造中依赖单位约定与局域化程度，重要的验证点是 N_theory 与 N_actual 的数值一致性（相对误差处于浮点精度水平）。

结论：在 Kitaev 链局域三站点嵌入中，数值验证显示 N_theory 与 N_actual 在所有采样点上重合（`rel_err`≈1e-14…1e-15），支持 PRQ 传播模型与相位干涉公式在实际模型中的普适性。

复现命令（在仓库根目录运行）：

```bash
.venv/bin/python quantity/validate_propagation_in_kitaev.py
```

（完）

## 15. Toy BdG 4×4 一致性验证（新增）

为了把 toy 端与 Kitaev 端统一到同一个 BdG 4×4  operator 空间，新增了最小三站点 toy 验证脚本：`quantity/validate_toy_bdg_4x4.py`。该脚本直接构造 3-site BdG 链，取 P={site0,site2}、Q={site1}，因此 Schur 补的 P 块天然是 4×4。

主要输出：
- 数据：`quantity/toy_bdg_validation/toy_bdg_4x4_validation.csv`
- 图：`quantity/toy_bdg_validation/toy_K_compare.png`
- 图：`quantity/toy_bdg_validation/toy_N_compare.png`
- 图：`quantity/toy_bdg_validation/toy_theta.png`

代表性结论：
- toy 端的 `K_pauli_proj` 与 `K_fro` 呈完全线性关系（Pearson = 1.0），说明 4×4 投影已稳定落在目标子空间；
- `N_theory = 2A sin(θ/2)` 与 `N_actual = ||Δ||_F` 在 toy 端逐点重合，验证仍保持到机器精度；
- 因而 toy 端与 Kitaev 端已经统一到同一套 BdG/JW/Majorana 投影语义，只是在参数扫描窗口与系统尺寸上略有不同。

建议在论文图注中同时引用这组三图，作为“最小同构模型”的补充验证。

# 13. Kitaev Schur 核

定义：

\[
T_{i\to i+1}
=
H_{P_iQ_{i+1}}
(\omega-H_Q)^{-1}
H_{Q_{i+1}P_{i+1}}
\]

在 Kitaev 点：

\[
H_{QQ}=0
\]

因此：

\[
(\omega-H_Q)^{-1}
=
\frac1\omega
\]

并且：

\[
H_{P_iQ_{i+1}}=it
\]

\[
H_{Q_{i+1}P_{i+1}}
=
-\frac{i\mu}{2}
\]

得到：

\[
\boxed{
T_{i\to i+1}
=
\frac{\mu t}{2\omega}
}
\]

---

# 14. 边界 Majorana

零模：

\[
H\psi=0
\]

写：

\[
\psi
=
\sum_i
(a_i\gamma_{iA}+b_i\gamma_{iB})
\]

得到递推：

\[
-\mu a_i+2ta_{i+1}=0
\]

即：

\[
\boxed{
a_{i+1}
=
\frac{\mu}{2t}a_i
}
\]

因此：

\[
a_i
=
\left(\frac{\mu}{2t}\right)^{i-1}a_1
\]

若：

\[
\left|\frac{\mu}{2t}\right|<1
\]

则：

\[
a_i\sim e^{-i/\xi}
\]

---

左边界 Majorana：

\[
\boxed{
\psi_L
=
\sum_i
\left(\frac{\mu}{2t}\right)^{i-1}
\gamma_{iA}
}
\]

右边界 Majorana：

\[
\boxed{
\psi_R
=
\sum_i
\left(\frac{\mu}{2t}\right)^{N-i}
\gamma_{iB}
}
\]

---

# 15. 传播与边界态

定义 transfer operator：

\[
\mathcal T:
\psi_i\mapsto\psi_{i+1}
\]

满足：

\[
\boxed{
\psi_{i+1}
=
\frac{\mu}{2t}\psi_i
}
\]

即：

\[
\boxed{
\psi_i
=
\mathcal T^{i-1}\psi_1
}
\]

若：

\[
\rho(\mathcal T)<1
\]

则 bulk propagation 收缩。

因此：

\[
\boxed{
\text{边界 Majorana}
=
\text{传播固定点}
}
\]

---

# 16. Jordan–Wigner 与 Pauli 表示

JW 变换：

\[
\gamma_{iA}
=
\left(\prod_{j<i}\sigma_j^z\right)\sigma_i^x
\]

\[
\gamma_{iB}
=
\left(\prod_{j<i}\sigma_j^z\right)\sigma_i^y
\]

因此：

\[
\gamma_{iB}\gamma_{i+1,A}
\sim
\sigma_i^y\sigma_{i+1}^x
\]

---

于是：

\[
\boxed{
h_{i,i+1}^{\mathrm{eff}}
=
K_{i,i+1}
\sigma_i^y\otimes\sigma_{i+1}^x
}
\]

---

# 17. 完整理论闭环

\[
\boxed{
P_iRQ_jRP_k
\Rightarrow
\Sigma(\omega)
\Rightarrow
h_{\alpha\beta}^{\mathrm{eff}}
\Rightarrow
R(u)
\Rightarrow
[h_{12}^{\mathrm{eff}},h_{23}^{\mathrm{eff}}]
\Rightarrow
\Delta_{YB}
\Rightarrow
\mathcal N
}
\]

---

# 18. 最终结论

理论本体：

\[
\boxed{
\mathcal N_{\mathrm{prop}}
=
\{P_iRQ_jRP_k\}
}
\]

Kitaev 链只是该传播几何中的一个特殊实例。

真正 fundamental 的对象是：

\[
\boxed{
P_iRQ_jRP_k
}
\]

Pauli tensor coupling：

\[
h_{\alpha\beta}^{\mathrm{eff}}
\]

只是其 operator representation。
