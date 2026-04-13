holonomy（回转/回环单连）是指在带联络的纤维丛或流形上把纤维沿某条闭合路径做平行移动后所得到的纤维自同构（变换）。所有闭路的这些变换构成的群称为 holonomy 群

##### 补充（平直联络与拓扑表示）

当联络的曲率 $F$ 在考虑的区域上严格为零（平直），则沿闭合曲线的平行搬运（holonomy）只依赖于曲线的同伦类；记联络为 $A$，闭合曲线 $\gamma$ 的 holonomy 写为
$$
U[\gamma]=\mathcal P\exp\Big(\oint_\gamma A\Big).
$$
利用非阿贝尔 Stokes（面次序指数），可以把它写成所围曲面 $S$ 上曲率的面次序指数：
$$
U[\gamma]=\mathcal P_S\exp\Big(\iint_S F\Big).
$$
若 $F\equiv0$，右侧对任意两個以同一同伦类为边界的填充面均相同，故 $U[\gamma]$ 仅取决于 $[\gamma]$（同伦类）。在参数空间或带标记点的配置空间情形，基礎群 $\pi_1$ 即为编织群（braid group）；在带洞的曲面情形，相关的是映射类群（mapping class group）。因此平直联络直接给出这些群在零模子空间上的表示：每个同伦类对应一个幺正算符，并且该算符由拓扑资料确定，成为拓扑不变量。

若曲率不是严格为零但在该子空间上为标量倍的单位算子（即 $F(x)=\lambda(x)I$，称为中心化曲率），则面次序指数退化为一个标量因子乘以单位算符，holonomy 仍由同伦类决定，但只给出群的投影表示（projective representation），这在带整体相位不可观测时很常见。

若 $F$ 既不为零也不在中心，則 holonomy 会依赖所选填充面的具体几何，从而不再是仅由同伦类决定的拓扑量。

直观举例：在 Ising‑型或任意子系统中，若在有效基态子空间上联接近似平直，则绕粒子交换的闭环给出的 $U[\gamma]$ 正是 braid group 的矩阵表示；在曲面上对基圈做 Dehn twist，得到的矩阵即映射类群的表示（half‑twist / Dehn twist 矩阵），这些矩阵在联接平直或曲率中心化时只由拓扑数据决定。

简明示例（Ising 任何子）:

- 交换两个 $\sigma$ 任何子（按融合基底 $\\{1,\psi\\}$）的 Braid 矩阵为对角形式：
	$$
	B_{\sigma\sigma}=\mathrm{diag}(R^{\sigma\sigma}_1,\;R^{\sigma\sigma}_\psi)=\mathrm{diag}(e^{-i\pi/8},\;e^{3i\pi/8}).
	$$
	这里 $R^{\sigma\sigma}_a$ 是交换两個 $\sigma$ 后若它们融合到通道 $a\in\{1,\psi\}$ 得到的相因子。

- Dehn twist / half‑twist（绕单个任何子做 $2\pi$ 自旋）由拓扑自旋 $\theta_a$ 给出：例如
	$$
	\theta_\sigma=e^{i\pi/8},\qquad \theta_\psi=-1,\qquad \theta_1=1.
	$$ 
	因此在对应的零模子空间上，绕包含某个任何子的圈做 Dehn twist 会乘以相应的 $\theta$（若作用在局域标记的子空间上则为标量因子），这就是映射类群元素在该子空间上的具体实现。

该示例展示了平直或中心化曲率下 holonomy 如何退化为由拓扑自旋和 $R$‑相位组成的群表示（或投影表示）。

##### 关于 $H_P$ 的精确定义与分解必要性（补充说明）

为便于后续的引理和误差估计，并与最初的三项分解保持一致，我们把每个局域生成元写成三部分的严格分解：
$$
H_P^{(ij)} = H_{\mathrm{quad}}^{(ij)} + H_{\mathrm{int}}^{(ij)} + H_{\mathrm{string}}^{(ij)},
$$
其中：
- $H_{\mathrm{quad}}$ 是二次（BdG）主项，易于 JW 映射与对角化，决定基态子空间与基准谱隙 $\Delta_0$；
- $H_{\mathrm{int}}$ 是局域的多体相互作用项（例如四费米），通常以算符范数 $\epsilon_{\mathrm{int}}:=\|\sum_{\langle ij\rangle}H_{\mathrm{int}}^{(ij)}\|$ 作为量化；
- $H_{\mathrm{string}}$（含此前称为 $H_{\mathrm{gauge}}$ 的串/规范项）包含 JW‑string 或非局域串项，其整體范數以 $\epsilon_{\mathrm{str}}:=\|\sum_{\langle ij\rangle}H_{\mathrm{string}}^{(ij)}\|$ 量化。

定義總微擾範數
$$
\epsilon_{\mathrm{tot}} := \epsilon_{\mathrm{int}} + \epsilon_{\mathrm{str}}.
$$

下面按原來四個引理說明各項的作用與如何把之前的 $\epsilon,\Delta_0$ 替換為三項範數：

- 引理 B（绝热性）和绝热误差估计：
	- 需要基态间隙下界 $\Delta\gtrsim\Delta_0-\epsilon_{\mathrm{tot}}>0$；
	- 絕熱誤差中 $\varepsilon_{\mathrm{adiab}}$ 的主要項仍由 $H_{\mathrm{quad}}$ 的時間導數控制（$\|\dot H_{\mathrm{quad}}\|/\Delta^2$），但 $H_{\mathrm{int}}$ 和 $H_{\mathrm{string}}$ 的存在會引入 $O(\epsilon_{\mathrm{tot}}/\Delta)$ 的修正。

- 引理 C（$e^{iH_P}$ 在 $\mathcal H_0$ 上逼近 R‑符号）：
	- 要求 $e^{iH_{\mathrm{quad}}^{(ij)}}$ 在編碼子空間上提供理想的 R‑表示；
	- 小範數的 $H_{\mathrm{int}},H_{\mathrm{string}}$ 會把這一表示擾動，量級為 $O(\epsilon_{\mathrm{tot}}/\Delta_0)$，因此引理 C 的技術假設可寫為 $\epsilon_{\mathrm{tot}}/\Delta_0\ll1$。

- 引理 D（quasi‑adiabatic 生成元與電路逼近）：
	- quasi‑adiabatic 生成元的局域性、衰減速率與截斷半徑主要由 $H_{\mathrm{quad}}$ 決定；
	- 在把連續生成元用有限輪 $R=e^{iH_P}$ 門逼近時，$H_{\mathrm{int}}$ 和 $H_{\mathrm{string}}$ 對 Trotter 誤差和局域截斷的影響可通過 $\epsilon_{\mathrm{tot}}$ 估算並作為次要修正。

- 引理 A（Berry 曲率与平坦性估计）：
	- 曲率 $F$ 的主贡献來自 $H_{\mathrm{quad}}$ 對投影的導數；若在 YBE/可积子流形附近 $H_{\mathrm{quad}}$ 產生小 $\partial P$，那麼
		$$\|F\|\lesssim O\big(\|\partial H_{\mathrm{quad}}\|/\Delta_0^2\big) + O\big(\epsilon_{\mathrm{tot}}/\Delta_0^2\big),$$
		因而 $\kappa=\kappa_0+O(\epsilon_{\mathrm{tot}}/\Delta_0^2)$。

此外，对各项的物理影响可以总结为：

- $H_{\mathrm{quad}}$：二次（BdG）主项；决定低能 A 矩阵、零模存在与体间隙 $\Delta_0$；在投影后直接產生 $so(2N)$ 聯絡的主贡献，易于解析与数值对角化；
- $H_{\mathrm{int}}$：纯四费米或更高次相互作用；会产生能谱重整化、零模能级劈裂或有效交換（非平凡的多体修正），并在微扰理论中以 $O(\epsilon)$ 改变 Berry 联络与 holonomy；
- $H_{\mathrm{gauge}}$：常数项、化学势修正和含 JW‑string 的非局域項（在费米子基底上表现为串项）；其主要影响是改变局域占据数与边界条件，并可能引入非局域耦合，需在设计门序列时尽量抑制或補償。

数值实现方面，你已有的脚本正是按此思想把 $H_0$ 作为二次可解部分、把 $V$ 当作扰动计算 $\Delta$、$\|V\|$ 与 $\eta=\|V\|/\Delta$ 并据此评估 BdG 近似的有效性——这在证明中作为量化假设出现。

因此在证明语境下，建议在每个引用上述引理的陈述处明确添加一条技术假设：存在三项分解
$$
H=H_{\mathrm{quad}} + H_{\mathrm{int}} + H_{\mathrm{string}}
$$
且 $\epsilon_{\mathrm{tot}}=\|H_{\mathrm{int}}\|+\|H_{\mathrm{string}}\| \ll \Delta_0$。有了该假设，所有引理的常数与误差项可被写成 $O(\epsilon_{\mathrm{tot}})$ 或 $O(\epsilon_{\mathrm{tot}}/\Delta_0)$ 的形式，从而完成定量的估计链。




*證明（較完整的技術推導草稿）.*

設 $H_0:=\sum_{\langle ij\rangle}H_P^{(0,ij)}$ 為理想可積點的累積局域生成元，實際哈密頓為 $H=H_0+V$，其中 $V=\sum_{\langle ij\rangle}V^{(ij)}$。令 $\Delta_0$ 為 $H_0$ 的基態-激發譜隙，假定 $\|V\|=\epsilon<\Delta_0/2$，則下述推導成立。

步驟 1（Duhamel 展開）—— 對單個局域門的指數差採用 Duhamel 展開：
$$
e^{i(H_0+V)}-e^{iH_0}=i\int_0^1 e^{i(1-s)H_0} V e^{is(H_0+V)}\,ds.
$$
將兩邊在基態投影 $P_0$ 上夾住，並用算符範數不等式得到粗略界 $\|P_0(e^{i(H_0+V)}-e^{iH_0})P_0\|\le\|V\|$。為獲得 $\Delta_0$ 的縮放，我們要利用 $V$ 與基態-激發間的弱耦合性。

步驟 2（Schrieffer–Wolff / 解析重整化）—— 構造準局域反對角生成元 $S$ 使得
$$
e^{S}(H_0+V)e^{-S}=H_0^{\mathrm{eff}}+R,
$$
其中 $H_0^{\mathrm{eff}}$ 作用在基態子空間上的有效哈密頓與 $H_0$ 等價（高階修正被吸收），殘餘項 $R=O(\epsilon^2/\Delta_0)$。典型地有 $\|S\|\lesssim O(\epsilon/\Delta_0)$ 且 $S$ 準局域。由此可得
$$
e^{i(H_0+V)}\big|_{\mathcal H_0}=U_{\mathrm{SW}}\,e^{iH_0^{\mathrm{eff}}}\,U_{\mathrm{SW}}^{-1} + O\Big(\frac{\epsilon^2}{\Delta_0^2}\Big),
$$
其中 $U_{\mathrm{SW}}=e^{S}$ 為準局域酉，且 $\|U_{\mathrm{SW}}-I\|\lesssim O(\epsilon/\Delta_0)$。

步驟 3（投影穩定性與拓撲序穩定性定理）—— Hastings–Wen 及 Bravyi–Hastings–Michalakis 的穩定性結果給出更嚴格的算符界：存在一個準局域酉 $U_{\mathrm{qa}}$ 使得
$$
\|U_{\mathrm{qa}}P_0U_{\mathrm{qa}}^{-1}-P\| \lesssim C\frac{\epsilon}{\Delta_0},
$$
因此基於此變換，$e^{i(H_0+V)}$ 在實際基態子空間與 $e^{iH_0}$ 在理想基態子空間上的表示相差 $O(\epsilon/\Delta_0)$。

步驟 4（路徑累積）—— 若一段空間重構路徑需要施加 $L$ 個局域門，保守地把每個局域門的基態子空間誤差線性累加，得到
$$
\varepsilon_{\mathrm{YBE}} \lesssim C_Y\,L\,\frac{\epsilon}{\Delta_0} + O\Big(\frac{\epsilon^2}{\Delta_0^2}\Big).
$$

結論：在 $\epsilon/\Delta_0\ll1$ 的物理範圍內，上述步驟給出了從局域微擾到 R‑表示偏差的可量化界（常數 $C_Y$ 與局域性尺度、格點幾何與路徑長度有關）。這就導出命題所述的比例關係，並可作為數值判據與穩定性檢驗的依據。





---

##### 3.7 量化推导：从三项分解到各误差项的显式界

下面给出完整且一致的推导链，展示如何把 $H=H_{\mathrm{quad}}+H_{\mathrm{int}}+H_{\mathrm{string}}$ 的分解转化为引理 A–D 中的量化界，并据此判断哪一项主导、哪一项可视为微扰。

1) 投影与导数的谱表示与界

令 $P$ 为 $H$ 的基态（或基态简并子空间）投影。用围绕基态谱的复曲线 $\Gamma$ 有
$$
P=-\frac{1}{2\pi i}\oint_{\Gamma}(z-H)^{-1}\,dz.
$$
对参数求导得到
$$
\partial_\mu P = -\frac{1}{2\pi i}\oint_{\Gamma}(z-H)^{-1}(\partial_\mu H)(z-H)^{-1}\,dz.
$$
若曲线与谱的距离下界为 $\Delta/2$，则 $\|(z-H)^{-1}\|\le 2/\Delta$，从而（忽略常数因子）可得常用界
$$
\|\partial_\mu P\| \lesssim \frac{\|\partial_\mu H\|}{\Delta}.
$$

2) 曲率的尺度

由于 $F_{\mu\nu}=P[\partial_\mu P,\partial_\nu P]P$，有
$$
\|F_{\mu\nu}\| \le 2\,\|\partial_\mu P\|\,\|\partial_\nu P\| \lesssim C\frac{\|\partial H\|^2}{\Delta^2}.
$$
因此曲率上界的典型量级为 $\kappa\sim O(\|\partial H\|^2/\Delta^2)$，这给出引理 A 中 $\kappa\,\mathcal A$ 项的自然尺度。

3) 投影与曲率的微扰修正（Davis–Kahan 类界）

把 $H=H_{\mathrm{quad}}+V$，$V=H_{\mathrm{int}}+H_{\mathrm{string}}$，记 $\Delta_0$ 为 $H_{\mathrm{quad}}$ 的基准谱隙，令 $\epsilon_{\mathrm{tot}}=\|V\|$，则 Davis–Kahan 型估计给出
$$
\|P-P_0\| \lesssim \frac{\epsilon_{\mathrm{tot}}}{\Delta_0}.
$$
进而投影导数與曲率的修正为低阶项，例如
$$
\|\partial P-\partial P_0\| \lesssim O\Big(\frac{\|\partial V\|}{\Delta_0}\Big) + O\Big(\frac{\epsilon_{\mathrm{tot}}\,\|\partial H\|}{\Delta_0^2}\Big).
$$
因此可以把曲率分解为 $\kappa=\kappa_0+O(\epsilon_{\mathrm{tot}}/\Delta_0^2)$。

4) 绝热误差的明确界（引理 B 的量化）

标准绝热估计表明，沿总时长 $T$ 的演化，跃迁幅度由 $\sup_t\|\dot H(t)\|/\Delta^2$ 控制：
$$
\varepsilon_{\mathrm{adiab}} \lesssim C_{\mathrm ad}\frac{\sup_t\|\dot H(t)\|}{\Delta^2}.
$$
将 $H=H_{\mathrm{quad}}+V$ 代入并分离主项與微扰，得到
$$
\varepsilon_{\mathrm{adiab}} \lesssim C_{\mathrm ad}\frac{\sup_t\|\dot H_{\mathrm{quad}}\|}{\Delta^2} + C'_{\mathrm ad}\frac{\sup_t\|\dot V\|}{\Delta^2}.
$$
若 $\sup_t\|\dot V\|\ll\sup_t\|\dot H_{\mathrm{quad}}\|$ 且 $\epsilon_{\mathrm{tot}}\ll\Delta_0$，則第二項為次要修正。

5) R‑表示偏差（引理 C）

對局域門 $R_{ij}=e^{iH_P^{(ij)}}$，用譜投影穩定性與 BCH/Dyson 展開可得在基態子空間上的表象差異受 $\epsilon_{\mathrm{tot}}/\Delta_0$ 控制，且沿路徑門數 $L$ 會乘上該誤差因子：
$$
\varepsilon_{\mathrm{YBE}} \lesssim C_Y\,L\,\frac{\epsilon_{\mathrm{tot}}}{\Delta_0} + O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^2}\Big).
$$

6) quasi‑adiabatic 生成元與 Trotter 誤差（引理 D）

quasi‑adiabatic 生成元通常尺度為 $k_*\sim O(\|\dot H\|/\Delta)$，對其採用 $m$ 步 Trotter 分解可得常見估計
$$
\varepsilon_{\mathrm{Trotter}} \lesssim C_T\frac{T\,k_*^2}{m} + C'_T e^{-\mu R_k} \,\sim\, C_T\frac{T}{m}\Big(\frac{\|\dot H\|}{\Delta}\Big)^2 + \text{(LR 截斷)}.
$$
因此增大 $m$ 或選擇更平滑的生成元可減小此項；微擾 $V$ 的影響主要通過改變 $\|\dot H\|$ 與 $k_*$，若 $\|\dot V\|\ll\|\dot H_{\mathrm{quad}}\|$，則為次要項。

7) 实用判据與数值流程

综上，判据与实践步骤为：

- 计算 $H_{\mathrm{quad}}$ 的谱并取基准间隙 $\Delta_0$；
- 计算 $\epsilon_{\mathrm{int}}=\|\sum H_{\mathrm{int}}^{(ij)}\|$ 与 $\epsilon_{\mathrm{str}}=\|\sum H_{\mathrm{string}}^{(ij)}\|$，得 $\epsilon_{\mathrm{tot}}$；
- 若 $\eta:=\epsilon_{\mathrm{tot}}/\Delta_0 \ll 1$，则 $H_{\mathrm{quad}}$ 为主导项，其他项可视为微扰；
- 用参数微小改变量估计 $\|\partial H\|,\|\dot H\|$，把它们代入上式估计 $\kappa,\varepsilon_{\mathrm{adiab}},\varepsilon_{\mathrm{Trotter}},\varepsilon_{\mathrm{YBE}}$，以确定 $T,m,L$ 等实际参数；
- 在数值上，`kits/compute_jw_mapping.py` 已包含计算 $\Delta,\|V\|,\eta$ 的模块，可直接把这些输出对照上式以判断是否满足 $\eta\ll1$ 及其他绝热/离散化约束。

最终结论：若满足
$$
\eta=\frac{\epsilon_{\mathrm{tot}}}{\Delta_0}\ll1,\qquad \frac{\sup_t\|\dot V\|}{\Delta_0^2}\ll1,
$$
则 $H_{\mathrm{quad}}$ 决定低能拓扑结构，$H_{\mathrm{int}}$ 与 $H_{\mathrm{string}}$ 是可控的微扰；否则必须把这些项纳入主项，重新评估谱结构与联络曲率。

（以上推导把常数依赖显式化为可数值估算的形式；若需，我可以把常数 $C_\ast$ 用具体曲线/格点参数给出并用你现有的矩阵规模做示例估算。）



#### 附录：对引理 A–D 的量化命題汇总

為方便讀者和數值評估，這裡把引理 A–D 的量化版本列成命題，明確寫出它們對 $H_0,V,\Delta$ 等參數的依賴關係。

- 命題 A.1 (量化非阿貝爾 Stokes): 若 $\sup_{x\in S}\|F(x)\|\le\kappa$，則存在常數 $C_A,C'_A$ 使得
$$
	\bigl\|U[\gamma']-U[\gamma]\bigr\| \le C_A\,\kappa\,\mathcal A + C'_A\,\kappa^2\,\mathcal A^2.
$$
在 $\kappa\mathcal A\ll1$ 的實用情形下，線性項主導。

- 命題 B.1 (量化絕熱誤差): 若 $H(t)=H_0(t)+V(t)$，整程譜隙下界為 $\Delta>0$，定義 $L_0=\sup_t\|\dot H_0\|,L_V=\sup_t\|\dot V\|,\epsilon=\sup_t\|V\|$，則
$$
	\varepsilon_{\mathrm{adiab}} \lesssim C_B\frac{L_0+L_V}{\Delta^2} + C''_B\frac{\epsilon}{\Delta}.
$$
第一項為傳統絕熱項，第二項為微擾對投影的直接修正。

- 命題 C.1 (R‑embedding 穩定性): 若 $H_P^{(ij)}=H_P^{(0,ij)}+V^{(ij)}$ 且整體微擾範數 $\epsilon_{\mathrm{tot}}:=\sup_\lambda\|\sum_{\langle ij\rangle}V^{(ij)}(\lambda)\|$，在理想點譜隙 $\Delta_0>0$ 下有
$$
	\varepsilon_{\mathrm{YBE}} \lesssim C_C\frac{\epsilon_{\mathrm{tot}}}{\Delta_0} + O\Big(\frac{\epsilon_{\mathrm{tot}}^2}{\Delta_0^2}\Big).
$$

- 命題 D.1 (Trotter/電路誤差): 若 quasi‑adiabatic 生成元分解為局域項且每項范數上界為 $k_*$，採用 $m$ 步 Trotter 分解與局域截斷半徑 $R_k$，則
$$
	\varepsilon_{\mathrm{Trotter}} \le C_D\frac{T\,k_*^2}{m} + C'_D e^{-\mu R_k}.
$$

這些命題把主不等式中的四個誤差項具體化為可數值估算的形式。建議在文本中引用這些命題位置同時標註用 `kits/compute_jw_mapping.py` 計算得到的 $\Delta,\|V\|,\eta$ 數值，以便把解析界和數值檢驗直接關聯起來。

---

#### 附錄 B：引理 C 与引理 D 的完整技術證明（詳盡版）

以下內容把上文中對引理 C、D 的要點推導展開為更詳細的技術證明草稿，便於檢查常數依賴、局域性假設與數值對接。

B.1 引理 C 的詳細證明（回顧與步驟）

- 假設與符號：令 $H_0$ 為理想可積/可解點的哈密頓，$V$ 為小局域擾動，$H=H_0+V$。令 $\Delta_0=\mathrm{gap}(H_0)$，$\epsilon=\|V\|$ 且假設 $\epsilon<\Delta_0/2$。

- 目標：在基態子空間上比較 $e^{iH}$ 與 $e^{iH_0}$ 的表示，並在 $\epsilon/\Delta_0\ll1$ 下給出以 $\epsilon/\Delta_0$ 為主的誤差界。

- Step 1（解析工具）: 用 Duhamel 展開
$$
	e^{iH}-e^{iH_0}=i\int_0^1 e^{i(1-s)H_0} V e^{isH}\,ds.
$$
把該等式夾在 $P_0$（$H_0$ 的基態投影）之間得到初步界，但要得到 $\Delta_0$ 的縮放需進一步分析 $V$ 在基態-激發子空間間的耦合。

- Step 2（Schrieffer–Wolff / 近似對角化）: 構造反對角生成元 $S$ 解決 $[H_0,S]+V_{\mathrm{off}}=0$（$V_{\mathrm{off}}$ 為 $V$ 的反對角部分），常規擴展得到
$$
\|S\|\lesssim \frac{\|V_{\mathrm{off}}\|}{\Delta_0}=O\Big(\frac{\epsilon}{\Delta_0}\Big).
$$
由此可把 $H$ 近似對角化為塊對角形式，並得到
$$
e^{iH}\big|_{\mathcal H_0}=U_{\mathrm{SW}}\,e^{iH_0^{\mathrm{eff}}}\,U_{\mathrm{SW}}^{-1}+O\Big(\frac{\epsilon^2}{\Delta_0^2}\Big),
$$
其中 $\|U_{\mathrm{SW}}-I\|=O(\epsilon/\Delta_0)$，$H_0^{\mathrm{eff}}$ 在編碼子空間上的有效作用接近 $H_0$。

- Step 3（穩定性定理補強）: 引用 Hastings–Wen 與 Bravyi–Hastings–Michalakis，可保證存在準局域酉 $U_{\mathrm{qa}}$（由 quasi‑adiabatic continuation 構造）滿足
	$$\|U_{\mathrm{qa}}P_0U_{\mathrm{qa}}^{-1}-P\|\le C\frac{\epsilon}{\Delta_0}.$$
	因此基態子空間與其投影的改變受 $\epsilon/\Delta_0$ 控制，並把此界傳遞到 $e^{iH}$ 的基態表示上。

- Step 4（路徑累積與最終界）: 把單個局域門在基態子空間上造成的誤差 $O(\epsilon/\Delta_0)$ 按門數 $L$ 保守累加，得到
	$$\varepsilon_{\mathrm{YBE}} \le C_Y\,L\,\frac{\epsilon}{\Delta_0} + O\Big(\frac{\epsilon^2}{\Delta_0^2}\Big).$$

B.2 引理 D 的詳細證明（回顧與步驟）

- 假設與符號：對投影族 $P(t)$（沿路徑參數化），採用平滑濾波函數 $W$ 構造 quasi‑adiabatic 生成元 $K(t)$，其分解為局域項 $k_\alpha(t)$，令 $k_*:=\sup_{t,\alpha}\|k_\alpha(t)\|$ 且每項支撐直徑不超過 $R_k$（截斷半徑）。

- Step 1（生成元構造與范數估計）: 取
$$
K(t)=i\int_{-\infty}^\infty W(\tau) e^{iH t}(\partial_t P) e^{-iH t}\,d\tau,
$$
對於合適的 $W$ 可證明 $K(t)$ 為準局域算符，且
$$
	k_*\lesssim C_{\mathrm{qa}}\frac{\|\partial_t HW\|}{\Delta}.
$$

- Step 2（時間離散化）: 將 $[0,T]$ 劃分為 $m$ 段，令 $K_j=K(t_j)$，採用二階 Trotter 展開，有
$$\Big\|\mathcal T e^{-i\int_0^T K(t)dt} - \prod_{j=1}^m e^{-iK_j\Delta t}\Big\| \lesssim C_T\frac{T}{m} \max_{j,j'}\|[K_j,K_{j'}]\|.
$$
	
用局域和結構和 $k_*$ 可得時間離散誤差上界約為 $\sim C_T T k_*^2/m$。

- Step 3（空間截斷）: 截斷每個 $k_\alpha$ 到半徑 $R_k$ 的局域算符會帶來 $\sim e^{-\mu R_k}$ 的尾項（由 Lieb–Robinson），因此得到
$$
	\varepsilon_{\mathrm{Trotter}} \lesssim C_T\frac{T k_*^2}{m} + C'_T e^{-\mu R_k}.
$$ 

- Step 4（門分解）: 每個局域指數 $e^{-i k_\alpha\Delta t}$ 可由有限輪 $R=e^{iH_P^{(ij)}}$ 門近似實現；實際門數與局域項的支撐、欲達到的精度共同決定常數因子，但不改變上述誤差的參數依賴（$m,R_k,k_*,T$）。

註：以上均為在局域性與譜隙假設下的標準技術推導；若需要，我可以把每一步中的常數用格點大小、耦合常數與哈密頓量矩陣範數做具體數值示例，並把建議的 $m,R_k,T,L$ 自動化為 `kits/compute_jw_mapping.py` 的輸出。

