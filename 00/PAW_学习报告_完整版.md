# 常压冷等离子体活化水（PAW）技术：从等离子体物理到液体化学的跨尺度学习报告

> **摘要**：等离子体活化水（Plasma Activated Water, PAW）是常压冷等离子体（Cold Atmospheric Plasma, CAP）领域中最具工程转化潜力的技术分支之一。本报告从常压冷等离子体的非平衡物理本质出发，系统阐述等离子体-液体界面化学的核心机理与反应动力学方程，在此基础上综述PAW在农业、医疗、环境及绿色化工四大领域的研究现状，剖析标准化缺失、活性机理不明、能效瓶颈和规模化困难等关键挑战，最后展望精准调控、微气泡强化、连续流工业化及AI辅助优化等未来方向。全文以定义明确、方程完整、逻辑连贯的方式构建从基础到应用的知识框架。

---

## 一、研究背景

### 1.1 常压冷等离子体的物理本质

等离子体被称为物质的"第四态"，是由电子、离子、中性粒子和激发态粒子组成的部分或完全电离气体。**常压冷等离子体（Cold Atmospheric Plasma, CAP）** 是一类特殊的非平衡等离子体体系，其核心特征在于热力学非平衡性：电子温度 $T_e$ 远高于离子温度 $T_i$ 和气体温度 $T_g$，即满足

$$
T_e \gg T_i \approx T_g \approx 300\text{--}400\,\text{K}
$$

这一非平衡特征使得CAP能够在接近室温的宏观条件下，通过高能电子（$T_e$ 通常为 $1\text{--}10\,\text{eV}$，对应约 $10^4\text{--}10^5\,\text{K}$）的碰撞过程产生大量化学活性物种，而整体气体温度却保持在可安全接触的水平。这是CAP区别于热等离子体（如电弧等离子体，$T_g\sim 10^4\,\text{K}$）的根本优势。

从放电物理角度看，CAP的产生依赖于外部电场对自由电子的加速。电子在电场中获得能量，通过以下碰撞过程将能量传递给重粒子：

$$
e^- + \text{M} \xrightarrow{\text{弹性碰撞}} e^- + \text{M} \quad (\text{能量转移效率} \propto 2m_e/m_M \ll 1)
$$

$$
e^- + \text{M} \xrightarrow{\text{非弹性碰撞}} e^- + \text{M}^* \quad (\text{激发})
$$

$$
e^- + \text{M} \xrightarrow{\text{非弹性碰撞}} e^- + \text{M}^+ + e^- \quad (\text{电离})
$$

由于电子质量 $m_e$ 远小于分子质量 $m_M$，弹性碰撞中每次碰撞仅转移 $\sim 2m_e/m_M$ 的能量，因此电子与重粒子之间的能量耦合效率极低。这种"选择性加热"机制正是CAP能够在室温下维持高活性粒子密度的物理根源。

### 1.2 活性氧氮物种（RONS）的产生

在空气放电条件下，CAP中高能电子与 $\text{N}_2$、$\text{O}_2$、$\text{H}_2\text{O}$ 分子碰撞，引发生成一系列**活性氧氮物种（Reactive Oxygen and Nitrogen Species, RONS）**。主要的初级反应包括：

**电子碰撞解离与激发：**

$$
\begin{aligned}
\text{O}_2 + e^- &\rightarrow \text{O}(^3\text{P}) + \text{O}(^3\text{P}) + e^- \quad &(\text{解离能 } 5.12\,\text{eV}) \\[4pt]
\text{N}_2 + e^- &\rightarrow \text{N}(^4\text{S}) + \text{N}(^4\text{S}) + e^- \quad &(\text{解离能 } 9.79\,\text{eV}) \\[4pt]
\text{H}_2\text{O} + e^- &\rightarrow \text{OH} + \text{H} + e^- \quad &(\text{解离能 } 5.1\,\text{eV}) \\[4pt]
\text{O}_2 + e^- &\rightarrow \text{O}_2(a^1\Delta_g) + e^- \quad &(\text{激发态，}\sim 0.98\,\text{eV})
\end{aligned}
$$

**气相次级反应（生成NO_x）：**

$$
\begin{aligned}
\text{N} + \text{O}_2 &\rightarrow \text{NO} + \text{O} \\[4pt]
\text{N} + \text{O} &\rightarrow \text{NO} \\[4pt]
\text{NO} + \text{O} &\rightarrow \text{NO}_2 \\[4pt]
\text{NO} + \text{O}_3 &\rightarrow \text{NO}_2 + \text{O}_2
\end{aligned}
$$

反应速率常数受电子能量分布函数（EEDF）的显著影响。电子能量分布函数 $f(\varepsilon)$ 满足玻尔兹曼方程：

$$
\frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_{\mathbf{r}} f - \frac{e\mathbf{E}}{m_e} \cdot \nabla_{\mathbf{v}} f = \left(\frac{\partial f}{\partial t}\right)_{\text{coll}}
$$

其中右侧碰撞项包括弹性、激发、电离和附着等所有碰撞过程。通过求解玻尔兹曼方程或使用BOLSIG+等工具，可以获得不同放电条件下的电子能量分布，进而确定各反应通道的速率。

### 1.3 CAP直接应用的局限性

尽管CAP在医学灭菌、创面修复、表面改性、污染物降解等领域展现出显著效果，但其直接应用面临三个根本性限制：

1. **活性物种寿命极短**：许多关键自由基的寿命仅在微秒到毫秒量级。例如，气相OH自由基的半衰期约 $10^{-5}\,\text{s}$，$\text{O}(^3\text{P})$ 的反应时间尺度约 $10^{-6}\,\text{s}$。这意味着只能在等离子体产生点附近发挥作用。

2. **作用距离有限**：等离子体中的活性粒子受限于扩散长度和气体流速，其有效作用距离通常不超过数厘米，难以深入多孔介质、复杂几何结构或深层组织内部。

3. **装置现场依赖**：等离子体激励电源、气体供应管路、放电模块等装置需要现场部署，极大限制了其在农业大面积喷施、远程医疗等场景中的应用。

### 1.4 PAW的技术思想

在上述限制下，**等离子体活化水（PAW）** 的核心思想应运而生：利用水作为活性物种的**储存介质**和**传递载体**。当CAP与水（或水溶液）接触时，气相中产生的RONS通过气-液界面传质进入液相，并在液相中进一步发生级联化学反应，最终形成含有 $\text{H}_2\text{O}_2$、$\text{NO}_2^-$、$\text{NO}_3^-$、$\text{O}_3$、过氧亚硝酸盐（ONOO$^-$）等稳定组分的水溶液。这一过程实现了"等离子体效应液体化"，使等离子体的化学活性突破了时间和空间的限制——类似于将闪电的化学能量"封存"在一杯水中。

---

## 二、核心机理：等离子体-液体界面化学

### 2.1 气-液传质过程

RONS从气相进入液相的过程受**亨利定律（Henry's Law）** 支配：

$$
c_i = k_{H,i} \cdot p_i
$$

其中 $c_i$ 为物种 $i$ 在液相中的平衡浓度（mol/L），$p_i$ 为气相分压（atm），$k_{H,i}$ 为亨利常数（mol·L$^{-1}$·atm$^{-1}$）。对于不同物种，亨利常数差异很大。例如，$\text{H}_2\text{O}_2$ 的 $k_H \approx 10^5\,\text{mol·L}^{-1}\text{·atm}^{-1}$（极易溶于水），而 $\text{NO}$ 的 $k_H \approx 1.9\times 10^{-3}\,\text{mol·L}^{-1}\text{·atm}^{-1}$（难溶），这决定了不同物种在液相中的富集程度。

传质通量可由双膜理论描述：

$$
J_i = k_{L,i}(c_i^* - c_i) = K_{G,i}(p_i - p_i^*)
$$

其中 $k_{L,i}$ 为液膜传质系数，$k_{G,i}$ 为气膜传质系数，$c_i^*$ 和 $p_i^*$ 为界面平衡浓度和分压。对于 $k_H$ 大的物种（如 $\text{H}_2\text{O}_2$），液膜阻力可忽略；对于 $k_H$ 小的物种（如 NO），气膜和液膜阻力均需考虑。

### 2.2 液相关键化学反应网络

进入液相后，RONS发生复杂的级联反应。以下为最核心的反应路径：

**（A）过氧化氢的生成路径：**

$$
\begin{aligned}
\text{OH} + \text{OH} &\rightarrow \text{H}_2\text{O}_2 \quad &k = 5.5\times 10^9\,\text{M}^{-1}\text{s}^{-1} \\[4pt]
\text{O}_3 + \text{H}_2\text{O} &\rightarrow \text{H}_2\text{O}_2 + \text{O}_2 \quad &k = 1.0\times 10^{-4}\,\text{s}^{-1} \\[4pt]
\text{HO}_2 + \text{HO}_2 &\rightarrow \text{H}_2\text{O}_2 + \text{O}_2 \quad &k = 8.3\times 10^5\,\text{M}^{-1}\text{s}^{-1}
\end{aligned}
$$

其中OH自由基的二聚反应是 $\text{H}_2\text{O}_2$ 最主要的生成路径，其极高的反应速率常数表明只要液相中有足够的OH浓度，$\text{H}_2\text{O}_2$ 即可快速积累。

**（B）氮物种的转化路径：**

$$
\begin{aligned}
2\text{NO}_2 + \text{H}_2\text{O} &\rightarrow \text{HNO}_2 + \text{HNO}_3 \quad &(\text{歧化反应}) \\[4pt]
\text{NO}_2 + \text{OH} &\rightarrow \text{HNO}_3 \quad &k = 1.0\times 10^{10}\,\text{M}^{-1}\text{s}^{-1} \\[4pt]
\text{NO} + \text{OH} &\rightarrow \text{HNO}_2 \quad &k = 1.0\times 10^{10}\,\text{M}^{-1}\text{s}^{-1}
\end{aligned}
$$

亚硝酸（$\text{HNO}_2$，p$K_a = 3.3$）和硝酸（$\text{HNO}_3$，p$K_a = -1.4$）在溶液中电离为 $\text{NO}_2^-$ 和 $\text{NO}_3^-$，同时释放 $\text{H}^+$，导致PAW呈现酸性（pH通常在 $2\text{--}4$ 范围）。

**（C）过氧亚硝酸盐的形成与分解——RONS耦合的关键节点：**

$$
\begin{aligned}
\text{H}_2\text{O}_2 + \text{HNO}_2 &\rightarrow \text{ONOOH} + \text{H}_2\text{O} \quad &(\text{pH < 5})\\[4pt]
\text{ONOOH} &\rightarrow \text{NO}_3^- + \text{H}^+ \quad &(\text{异构化})\\[4pt]
\text{ONOOH} &\rightarrow \text{NO}_2 + \text{OH} \quad &(\text{均裂，产生活性自由基})
\end{aligned}
$$

过氧亚硝酸（ONOOH，p$K_a = 6.8$）及其共轭碱 ONOO$^-$ 是PAW化学中的"枢纽分子"。它由 $\text{H}_2\text{O}_2$ 和 $\text{HNO}_2$ 的酸催化反应生成，又能在生理pH条件下均裂产生 $\text{NO}_2$ 和 $\text{OH}$ 这两种高活性自由基，是连接 ROS（活性氧）和 RNS（活性氮）两条反应链的关键桥梁。

### 2.3 pH对液相化学的调控

PAW的液相反应化学强烈依赖于pH。pH不仅通过酸碱平衡影响各物种的存在形态：

$$
\begin{aligned}
\text{HNO}_2 &\rightleftharpoons \text{NO}_2^- + \text{H}^+ \quad &\text{p}K_a = 3.3 \\[4pt]
\text{ONOOH} &\rightleftharpoons \text{ONOO}^- + \text{H}^+ \quad &\text{p}K_a = 6.8 \\[4pt]
\text{H}_2\text{O}_2 &\rightleftharpoons \text{HO}_2^- + \text{H}^+ \quad &\text{p}K_a = 11.6
\end{aligned}
$$

还通过影响反应通道的速率决定了活性物种的主导类型。在低pH条件（pH $\sim 2\text{--}3$）下，ONOOH的形成被显著促进，体系中以RNS效应为主；而在近中性条件（pH $\sim 6\text{--}7$）下，ONOO$^-$ 相对稳定，ROS路径贡献更大。这一pH敏感性是PAW在不同应用场景中表现差异的重要原因。

---

## 三、制备技术

### 3.1 等离子体射流法

**结构原理**：工作气体（Ar、He或其与空气的混合气）通过管状电极产生等离子体，形成等离子体射流（Plasma Jet），喷射到液体表面。射流长度可从数毫米到数厘米，功率通常在 $1\text{--}50\,\text{W}$ 范围。

**化学特征**：由于射流周围空气的卷吸作用，实际放电化学以工作气体+空气混合放电为主。使用Ar或He作为工作气体时，电子温度较高（$T_e \sim 2\text{--}5\,\text{eV}$），有利于产生高密度OH自由基（通过电子与 $\text{H}_2\text{O}$ 的碰撞解离）。

**优势与局限**：结构简单、操作灵活，特别适合实验室机理研究；但单支射流处理面积小，规模化困难。

### 3.2 介质阻挡放电法（DBD）

**结构原理**：至少一个电极被介质层（石英、陶瓷、玻璃等）覆盖，在交流或脉冲高压（$5\text{--}50\,\text{kV}$，频率 $50\,\text{Hz}\text{--}100\,\text{kHz}$）驱动下产生均匀或丝状放电。放电间隙 $0.1\text{--}10\,\text{mm}$。

**放电物理特征**：介质层的存在限制了单次放电的电流增长，避免了电弧转变。能量沉积由每周期转移的电荷量决定：

$$
E_{\text{pulse}} = \int_{t_1}^{t_2} V(t) I(t)\,dt
$$

**化学特征**：DBD在空气放电中产生丰富的 $\text{O}_3$ 和 $\text{NO}_x$。由于比表面积大，DBD适合处理大面积液体薄膜或薄液层。**这是目前农业PAW应用中最广泛采用的技术路线。**

### 3.3 等离子体-液体直接接触法

**结构原理**：放电电极直接浸入水中或液面上方极近处，在液相内部或气液界面产生放电。典型形式包括水中脉冲放电、滑弧放电（Gliding Arc）等。

**化学特征**：活性物种浓度高，紫外辐射和冲击波效应可协同作用于液体。但由于电极与液体直接接触，存在严重的电极腐蚀和金属离子溶出问题，可能引入重金属污染。

### 3.4 新兴制备技术

- **远程/离位PAW（Remote/Off-site PAW）**：等离子体先处理气体，再将活化气体通入水中，避免电极污染。
- **微气泡增强PAW**：通过微气泡（直径 $< 100\,\mu\text{m}$）注入增大气液接触面积。比表面积可达到 $10^4\,\text{m}^2/\text{m}^3$ 量级，大幅提升传质效率。

---

## 四、研究现状

### 4.1 农业应用

PAW中同时含有氮源（$\text{NO}_3^-$、$\text{NO}_2^-$）和活性氧（$\text{H}_2\text{O}_2$），对植物兼具营养和信号调控双重作用。

**氮营养效应**：植物通过硝酸还原酶（NR）和亚硝酸还原酶（NiR）将PAW中的无机氮同化为氨基酸：

$$
\text{NO}_3^- \xrightarrow{\text{NR}} \text{NO}_2^- \xrightarrow{\text{NiR}} \text{NH}_4^+ \rightarrow \text{氨基酸}
$$

这使PAW可部分替代传统氮肥。研究表明，PAW处理后的**小麦发芽率提高 $15\text{--}30\%$**，**番茄幼苗生物量增加 $20\text{--}40\%$**，**生菜叶片硝酸盐含量提高的同时抗氧化酶活性也得到增强**。

**ROS信号调控**：低浓度 $\text{H}_2\text{O}_2$（$< 100\,\mu\text{M}$）可作为信号分子激活植物抗氧化防御系统，包括超氧化物歧化酶（SOD）、过氧化氢酶（CAT）和过氧化物酶（POD）的上调表达，增强植物对干旱、盐碱等非生物胁迫的抗逆性。但浓度过高（$> 1\,\text{mM}$）会导致氧化损伤，因此PAW的"施用量窗口"是农业应用的关键参数。

### 4.2 医疗灭菌与生物医学

PAW在医疗领域的核心优势是能够进入**直接CAP无法到达的深层组织和不规则结构**。

**灭菌机理**：PAW中活性物种通过以下多靶点机制破坏微生物：

1. **细胞膜脂质过氧化**：OH和ONOO$^-$ 攻击膜磷脂的不饱和脂肪酸链，增加膜通透性，导致细胞内容物泄漏；
2. **蛋白质氧化损伤**：$\text{H}_2\text{O}_2$ 和RNS氧化半胱氨酸巯基（-SH $\rightarrow$ -SO$_3$H）和甲硫氨酸残基，导致酶失活；
3. **DNA链断裂**：OH自由基引起DNA糖-磷酸骨架的抽氢反应，造成单链和双链断裂。

实验证明PAW对**大肠杆菌、金黄色葡萄球菌、铜绿假单胞菌及耐药菌株（如MRSA）的杀灭率达到 $99.9\%$ 以上**。由于PAW通过多靶点氧化损伤发挥作用，而非单一的受体-配体结合机制，细菌难以产生耐药性——这是相对于传统抗生素的独特优势。

### 4.3 环境治理

PAW用于降解水体中的有机污染物，依赖RONS的氧化还原反应破坏污染物分子结构。

以染料亚甲基蓝（MB）为例，其降解涉及拟一级反应动力学：

$$
\ln\left(\frac{[\text{MB}]_t}{[\text{MB}]_0}\right) = -k_{\text{app}} t
$$

其中表观速率常数 $k_{\text{app}}$ 与液相中 $\text{H}_2\text{O}_2$ 和 $\text{O}_3$ 的稳态浓度相关：

$$
k_{\text{app}} = k_{\text{H}_2\text{O}_2}[\text{H}_2\text{O}_2] + k_{\text{O}_3}[\text{O}_3] + k_{\text{OH}}[\text{OH}]_{\text{ss}}
$$

近年来研究重点已从单一的染料模型扩展到**抗生素（四环素、磺胺类）、全氟和多氟烷基物质（PFAS）、农药残留**等难降解污染物。针对PFAS的降解，$\text{O}_3$ 和 $\text{OH}$ 的协同作用被认为是C-F键断裂的关键。

### 4.4 绿色化工催化

PAW不仅是一种处理介质，也是一种**含活性物种的反应介质**。探索性研究显示：

- **选择性氧化**：利用PAW中的 $\text{H}_2\text{O}_2$ 替代传统化学氧化剂，在温和条件下实现醇→醛/酮的选择性转化；
- **CO$_2$转化**：PAW中活性物种可促进CO$_2$还原为甲酸盐或甲醇，目前效率较低但展示了绿色合成潜力。

---

## 五、当前存在的关键问题

### 5.1 作用机理尚未系统阐明

PAW含有数十种活性组分，且各组分之间存在复杂的生成-消耗耦合关系和pH依赖性。当前研究面临的主要困境是：

- 不同实验条件下观察到的主要杀菌/促生长活性物种不一致——部分研究归因于 $\text{H}_2\text{O}_2$，部分强调 $\text{ONOO}^-$，还有研究表明**酸性pH本身**即可贡献显著的灭菌效果；
- 液相反应网络的"涌现行为"使得难以通过单一物种的浓度变化预测PAW的宏观生物学效应；
- 缺乏从放电参数（电压、频率、气体组成）→气相RONS→液相化学→生物效应的**全链条定量模型**。

### 5.2 标准化缺失

PAW制备涉及多个自由度：放电电压 $V$、频率 $f$、放电功率 $P$、气体种类与流量、初始水质（电导率 $\sigma$、pH、缓冲容量）、处理时间 $t$、液层厚度 $h$ 等。不同实验室在这些参数上的差异导致制备的PAW在 $\text{H}_2\text{O}_2$（$10\,\mu\text{M}\text{--}10\,\text{mM}$）、$\text{NO}_3^-$（$10\,\mu\text{M}\text{--}100\,\text{mM}$）、pH（$2\text{--}7$）等关键指标上跨数个数量级的差异，使实验结果难以横向比较。

### 5.3 储存稳定性

PAW中的活性物种随时间指数衰减。典型的一级衰减动力学：

$$
[\text{species}]_t = [\text{species}]_0 \cdot e^{-k t}
$$

其中 $\text{H}_2\text{O}_2$ 的衰减半衰期 $t_{1/2} = \ln 2/k$ 通常为数小时到数天，$\text{NO}_2^-$ 在酸性条件下的衰减更快（通过歧化反应生成 $\text{NO}_3^-$ 和 NO）。臭氧（$\text{O}_3$）的液相半衰期更短，通常在分钟级。商业化需要至少数周以上的保质期，这要求发展PAW的稳定化策略。

### 5.4 能效经济性

从电能到"有效化学效应"的能量效率是PAW工业化的核心瓶颈。能量效率可定义为单位能耗产生的活性物种量：

$$
\eta = \frac{c_{\text{RONS}} \cdot V}{E_{\text{input}}}
$$

或以G值（每吸收100 eV能量产生的分子数）衡量：

$$
G_{(\text{H}_2\text{O}_2)} = \frac{N_{\text{H}_2\text{O}_2}}{E_{\text{abs}}/100\,\text{eV}}
$$

目前PAW制备中 $\text{H}_2\text{O}_2$ 的典型G值约 $0.1\text{--}1\,\text{分子}/100\,\text{eV}$，远低于理论热力学极限。对于大规模废水处理，PAW的单位能耗仍高于臭氧氧化、Fenton氧化等传统高级氧化工艺。

### 5.5 规模化制造

从实验室烧杯（$50\text{--}500\,\text{mL}$）到工业连续流（$> 1000\,\text{L/h}$），面临以下工程挑战：

- 功率耦合效率随液层厚度增加而下降；
- 活性物种的空间分布不均匀性在放大后加剧；
- 连续流反应器中的停留时间分布需要精确设计和控制。

---

## 六、未来发展趋势

### 6.1 精准PAW调控

从"产生尽可能多的活性物种"转向"**定向产生特定类型的活性物种**"。通过调节工作气体（$\text{O}_2/\text{N}_2$ 比、Ar/O$_2$ 混合比）、放电模式（脉冲 vs. 连续）、功率密度和后处理策略，实现：

- **高 $\text{H}_2\text{O}_2$ 型PAW**：适用于灭菌和氧化降解；
- **高 $\text{NO}_3^-$ 型PAW**：适用于农业氮肥替代；
- **高 $\text{ONOO}^-$ 型PAW**：适用于生物医学（利用其独特的硝化和氧化双重反应性）。

### 6.2 连续流工业化系统

连续流反应器（Continuous Flow Reactor）将成为PAW产业化的核心装备形态。关键设计要素包括：

- 液膜流动模式优化（降膜、螺旋流、射流撞击）；
- 多级串联处理以提高活性物种利用效率；
- 与在线监测（$\text{H}_2\text{O}_2$ 传感器、pH电极、电导率探头）集成，实现反馈控制。

### 6.3 微气泡强化技术

微气泡的引入可将气-液界面面积提升至传统鼓泡的 $10^2\text{--}10^3$ 倍。微气泡内部的气体压力满足Young-Laplace方程：

$$
\Delta P = \frac{2\gamma}{r}
$$

其中 $\gamma$ 为表面张力，$r$ 为气泡半径。微米级气泡（$r \sim 10\,\mu\text{m}$）的内部压力可达 $0.15\,\text{atm}$（相对于大气压），增强了气体向液相的传质驱动力。此外，微气泡的收缩-坍塌过程可能引发局部热点效应，进一步增强化学反应活性。

### 6.4 PAW与催化体系协同

通过引入异相催化剂（如 $\text{TiO}_2$、$\text{Fe}_3\text{O}_4$）或均相Fenton试剂（$\text{Fe}^{2+}/\text{H}_2\text{O}_2$ 体系），可实现PAW中 $\text{H}_2\text{O}_2$ 的催化活化：

$$
\text{Fe}^{2+} + \text{H}_2\text{O}_2 \rightarrow \text{Fe}^{3+} + \text{OH} + \text{OH}^- \quad (\text{Fenton反应})
$$

$$
\text{TiO}_2 + h\nu \rightarrow e^-_{\text{CB}} + h^+_{\text{VB}} \quad (\text{光催化})
$$

催化剂可将PAW中相对温和的 $\text{H}_2\text{O}_2$ 转化为反应性更强的OH自由基，实现"活性升级"。光催化耦合则能利用太阳能驱动PAW体系的持续活化。

### 6.5 数字化与AI优化

PAW制备是典型的多变量、非线性复杂过程，非常适合引入机器学习方法：

- 建立从放电参数（电压、频率、气体流量、液流速率）到PAW关键质量指标（$\text{H}_2\text{O}_2$ 浓度、pH、$\text{NO}_3^-/\text{NO}_2^-$ 比值、杀菌效力）的**代理模型（Surrogate Model）**；
- 使用贝叶斯优化或多目标进化算法在线搜索最优操作窗口；
- 构建**数字孪生（Digital Twin）**，实现PAW制备过程的实时仿真、故障预测和自适应控制。

---

## 七、总结与展望

等离子体活化水技术是常压冷等离子体从"物理-化学现象"走向"工程实用产品"的关键桥梁。其核心科学问题可概括为：**将纳秒至微秒尺度的电子碰撞能量通过级联化学反应"编码"为具有分钟至天级寿命的液相化学信息**。这一跨时间尺度的转换涉及放电物理、界面传质、液相自由基化学和生物响应等多学科交叉。

当前研究已充分证明了PAW在农业增产、灭菌消毒、污染物降解等方面的有效性，正在从"证明有效"（proof-of-concept）阶段向"阐明机理、建立标准、实现工程化"阶段过渡。未来的突破将依赖于：

1. **机理层面**：建立从放电参数→气相化学→液相化学→生物效应的多尺度定量模型；
2. **工程层面**：发展高效连续流反应器和微气泡强化传质技术；
3. **标准层面**：建立以关键活性物种浓度和生物效应为指标的PAW质量评价体系；
4. **智能层面**：AI辅助的过程优化和数字孪生驱动的智能制造。

当这些问题得到系统解决后，PAW有望成为继传统化学消毒剂、液体氮肥和化学氧化剂之后的**新一代绿色液体功能材料**，为可持续农业、公共健康和环境治理提供全新的技术范式。

---

## 参考文献

1. Brandenburg, R. (2024). Plasma-activated water: Mechanisms, applications, and challenges. *Plasma Processes and Polymers*, 21(5), 2400249.

2. Kim, S., Park, J., & Lee, H. (2023). Plasma activated water combined with catalytic systems for organic pollutant degradation. *Chemical Engineering Journal*, 452, 139270.

3. Lukes, P., Locke, B. R., & Brisset, J. L. (2014). Plasma Chemistry and Catalysis in Liquids: Applications for Water Treatment. Springer.

4. Xiong, J., Li, H., & Liu, T. (2024). Off-site plasma-activated water generation: Towards continuous industrial-scale production. *Separation and Purification Technology*, 311, 125952.

5. Zhou, X., Chen, J., & Wang, Y. (2024). Microbubble-enhanced plasma activated water for efficient wastewater treatment. *Journal of Hazardous Materials*, 462, 131893.

6. Bruggeman, P. J., Kushner, M. J., Locke, B. R., et al. (2016). Plasma–liquid interactions: a review and roadmap. *Plasma Sources Science and Technology*, 25(5), 053002.

7. Machala, Z., Tarabová, B., Hensel, K., et al. (2019). Formation of ROS and RNS in Water electro-sprayed through Transient Spark Discharge in Air and their Bactericidal Effects. *Plasma Processes and Polymers*, 16(1), 1800136.

8. Thirumdas, R., Kothakota, A., Annapure, U., et al. (2018). Plasma activated water (PAW): Chemistry, physico-chemical properties, applications in food and agriculture. *Trends in Food Science & Technology*, 77, 21-31.
