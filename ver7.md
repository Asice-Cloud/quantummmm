# ver7 — Eight‑vertex: 从 $d$ 到全局 BdG 能谱的完整推导

目的：把仓库中 eight‑vertex 示例的代数推导写清楚，解释为什么文档中说 $\delta=0$ 看起来像 MZM，$\delta\neq0$ 看起来像 ABS，并给出用你们映射得到“全局能量谱”的明确公式与数值检验步骤。

## 1. 模型与 Pauli 展开
我们采用仓库中定义的 eight‑vertex 4×4 哈密顿（参见实现 `tools/*`）：
$$
H_4(u,\delta)=\cos u\,X\otimes X+\tfrac{\sin u}{2}(Y\otimes X - X\otimes Y)+\tfrac{\delta}{2}(Z\otimes I - I\otimes Z).
$$

按 $\tfrac14\operatorname{tr}[P H_4]$ 展开可读出非零 Pauli 系数：
$$
\begin{aligned}
h_{xx}&=\cos u,\\
h_{xy}&=-\tfrac12\sin u,\\
h_{yx}&=+\tfrac12\sin u,\\
h_{zI}&=+\tfrac12\delta,\\
h_{Iz}&=-\tfrac12\delta.
\end{aligned}
$$

其余 $h_{\alpha\beta}=0$。

## 2. 投影到逻辑子空间（$H_{\rm eff}$）
定义子空间 $\{|01\rangle,|10\rangle\}$ 的投影后得到
$$
H_{\rm eff}=d_0 I + d_x\sigma_x + d_y\sigma_y + d_z\sigma_z
$$
其中（见 `ybe223.md` 的映射）
\begin{aligned}
d_x&=h_{xx}+h_{yy},\\
d_y&=-h_{xy}+h_{yx},\\
d_z&=h_{zI}-h_{Iz},\\
d_0&=h_{II}-h_{zz}.
\end{aligned}

代入 eight‑vertex 的 $h$ 得到
$$
\vec d(u)=(d_x,d_y,d_z)=(\cos u,\;\sin u,\;\delta),\qquad d_0=0.
$$
因此投影子空间的瞬时本征值（2×2）為
$$
E_{\rm eff}(u)=d_0\pm|\vec d|=\pm\sqrt{1+\delta^2}.
$$
注意：这是**逻辑子空间的即时本征能级**，并不是全局 BdG 单粒子谱的 $E=0$ 判据。

## 3. Pauli → Kitaev 的映射（仓库约定）
仓库中用于构造 BdG 的线性组合为：
$$
\begin{aligned}
t &= h_{xx}+h_{yy}+i(h_{xy}-h_{yx}),\\
\Delta &= h_{xx}-h_{yy}-i(h_{xy}+h_{yx}),\\
\mu &= 4h_{zz}-2(h_{zI}+h_{Iz}).
\end{aligned}
$$

把 eight‑vertex 的 $h$ 代入得：
$$
\begin{aligned}
t &= \cos u + i\big(-\tfrac12\sin u-\tfrac12\sin u\big)=\cos u - i\sin u = e^{-iu},\\
\Delta &= \cos u - i\big(-\tfrac12\sin u+\tfrac12\sin u\big)=\cos u,\\
\mu &= 4\cdot0 - 2\big(\tfrac12\delta + (-\tfrac12\delta)\big)=0.
\end{aligned}
$$

**重要观察**：在此约定下 $\delta$（即 $h_{zI}=-h_{Iz}$ 的反对称组合）恰好在构造 $\mu$ 的线性组合中相互抵消，因此映射出的 $(t,\Delta,\mu)$ **不依赖于** $\delta$。

## 4. 构造全局 BdG 与能谱
采用仓库對周期链的約定，2×2 Bloch‑BdG 為
$$
\begin{aligned}
h_k &= -\frac{\mu}{2} - \frac12\big(t e^{ik} + t^* e^{-ik}\big),\\
\Delta_k &= \tfrac12(\Delta e^{ik} - \Delta e^{-ik}) = i\Delta\sin k.
\end{aligned}
$$
代入上面得到的 $(t,\Delta,\mu)$：
$$
\begin{aligned}
h_k &= -\frac{0}{2} - \frac12\big(e^{-iu} e^{ik} + e^{iu} e^{-ik}\big) = -\cos(k-u),\\
\Delta_k &= i\cos u\;\sin k.
\end{aligned}
$$
因此全局單粒子能谱為
$$
E(k)=\pm\sqrt{h_k^2+|\Delta_k|^2}=\pm\sqrt{\cos^2(k-u)+\cos^2u\;\sin^2k }.
$$
bulk gap 即 $\min_k |E(k)|$；對本模型常見的拓扑判據 $|\mu|<2|t|$ 在這裡退化為 $0<2|e^{-iu}|$，即 $\mu=0$ 滿足拓撲不變量條件（只要 $\Delta\neq0$ 且 gap>0，則開鏈會有零模）。

## 5. 解析矛盾：為什麼文檔說 $\delta=0$ 是 MZM，$\delta\neq0$ 是 ABS？
- 在 **邏輯子空間（門/演化層面）** 看：$\delta$ 是 $d_z$，改變瞬時本征值；$\delta\neq0$ 會在投影後的 2×2 能級上產生劈裂，使局域測量/門操作上表現出明顯的動力學相與能級分裂（看起來像 ABS）。當 $\delta=0$ 時軌跡在 $(d_x,d_y)$ 平面繞原點，幾何相主導，門更接近理想 braid（文中稱的 MZM‑like 行為）。
- 在 **全局 BdG（拓撲/開鏈譜）** 看：對於 eight‑vertex 的這一嵌入映射，$\delta$ 不進入 $(t,\Delta,\mu)$，因此 bulk 譜與拓撲不變量不受 $\delta$ 影響；若 $\Delta\neq0$ 且 gap>0，開鏈應仍有拓撲零模（MZM）。

因此兩者並不矛盾：文中所說的 ABS/MZM 區分更多是從“投影後的態與門操作的可觀測行為”來描述的，而判定是否存在真正的拓撲零模必須參照**全局 BdG 的體譜與開鏈邊態**。

## 6. 建議的數值檢驗步驟（可複製命令）
1. 計算映射與拓撲指示（批量）
```bash
python tools/compute_topo_invariant.py --u-list '0,1.5708,3.1416,4.7124' --delta-list '0,0.015,0.1'
```
2. 開鏈低能掃描（檢查 $E_0(L)$ 與端點權）
```bash
python tools/edge_localization.py --u 0.0 --delta 0.0 --Ls 40,80,160,320
```
3. 單點 LDOS 可視化（查看 $E\approx0$ 的空間分佈）
```bash
python tools/compute_ldos.py --u 0.0 --delta 0.0 --L 160
```

## 7. 結論（一句話）
在你們的 eight‑vertex 嵌入約定下，$\delta$ 改變的是投影後 $H_{\rm eff}$ 的 $d_z$（影響門操作和瞬時能級），但對映射到的 Kitaev 參數 $(t,\Delta,\mu)$ 無影響；因此判定“是否為拓撲 MZM”必須基於全局 BdG（bulk gap + 開鏈邊態），而不能僅憑 $\vec d$ 的形式。

---
（文檔初稿，如需我可把若干典型 $(u,\delta)$ 數值點跑一遍並把結果附在此文件下方。）

## 8. 脚本层级：哪些在看全链，哪些只是在看局域 $d$

为了避免把“局域投影诊断”和“全链谱判断”混在一起，这里把仓库中的相关脚本分成三层。

### 8.1 全链 BdG 层
这类脚本先走 `H4 → pauli_expand → (t,Δ,μ) → build_bdg_chain`，再看谱、LDOS 或 `E0(L)`。

- `tools/compute_ldos.py`：输出完整 LDOS 和零能空间剖面，判断依据是全链 BdG 本征谱。
- `tools/edge_localization.py`：扫描 `E0(L)` 和边权重，判断零模是否随链长趋零并端点局域。
- `tools/compute_topo_invariant.py`：先映射到 `(t,Δ,μ)`，再算 bulk gap 和简单拓扑判据。

### 8.2 局域投影层
这类脚本只是在逻辑子空间里看 `d`，适合诊断路径几何或局域劈裂，但不能单独当全链结论。

- `tools/test_mzm_to_abs.py`：直接投影到 `|01>, |10>`，看 `|d(u)|` 是否被耦合抬离原点。

### 8.3 过渡层
这类脚本会同时比较全链结果和低维代理量。

- `tools/generate_comparison_panels.py`：一边算全链 `min|E|` / LDOS，一边算低维拟合能量 `E_pred`，明确区分“全链”与“代理模型”。

### 8.4 结论
- 核心脚本没有把 `d` 误当成全链谱。
- 真正判断 MZM / ABS，优先看 `compute_ldos.py`、`edge_localization.py`、`compute_topo_invariant.py`。
- `test_mzm_to_abs.py` 这类脚本只能说明局域有效模型的几何变化，不能替代全链 BdG 结论。


