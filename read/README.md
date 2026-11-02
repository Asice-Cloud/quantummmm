# PRB 111 研究整理

已有结果的技术推导保留两份主文档，组会建议先阅读独立汇报版，再查阅技术推导。后续控制方向另设讨论文档，区分已证明构造与待研究问题。

| 文档 | 内容 |
|---|---|
| [解析式结果：简要汇报](analytic_continuation_brief.md) | 约 10–15 分钟；完整保真度公式及函数定义、三步推导、图像机制、已有验证和近似范围 |
| [解析延拓保真度工作：完整组会汇报](analytic_continuation_report.md) | 动机、模型、完整公式构造、图像讲解、独立检验与误差控制；含 18 页幻灯片安排、讨论问答和技术附录 |
| [微扰、几何结构与对称性抵消](perturbation.md) | 精确能谱、快模消去、holonomy、二阶与四阶保真度、群平均与九圈校正门；第 10 节讨论新增耦合恢复简并与原控制集内的动态抵消 |
| [解析延拓、具体公式与验证](analytic_continuation.md) | 完整标量函数构造、局部递推与延拓、读出公式、峰位峰高与长时机制、完整 SO(5) 矩阵函数、图像核对、假设和误差预算 |
| [群对称控制：具体协议与研究价值](group_symmetry_control.md) | 超导相位与逻辑旋转的区别、六圈协议及受限最小性、两圈标定停留、高阶抵消障碍、泄漏与公平比较、文献定位和后续研究问题 |

**工作总结：**[2026-11-01：公式、零耦合极限、PRB 105 对照与下一步](summary_2026-11-01.md)。

**控制路线讨论：**[微扰文档第 10 节](perturbation.md#10-后续研究全程简并补偿与有限时间对称抵消) 给出一个依赖额外耦合的精确简并补偿构造，并推导原三段正向路径仅改波形、正时长时的一阶门误差障碍。补偿后的目标 holonomy、控制通道的最小性与新协议的完整误差仍需研究；不属于原保真度公式已经验证的范围。

第 10.9–10.12 节进一步给出**两圈编织加标定停留**的一阶动态校正构造、明确的二阶残余和比例漂移响应。它复用原第二段的解析传播，减少编织次数，但不具有通用群平均序列对独立弱参数的相同鲁棒性，尚未验证完整模型的总误差。

**2026-11-02 检验更新：**[解析延拓文档第 16–17 节](analytic_continuation.md) 已完成全部四阶系数的独立核对、明确的六阶截断余项界，以及相对完整五模模型的可计算物理误差上界。物理界只需脉冲积分，并有适用域证明；它较保守，不等于对整张图保证 0.01 精度。代码与数据同时记录理论上界和已有实际偏差。

**图像入口：**[论文与公式热图](figures/fidelity_formula/paper_vs_formula.pdf)；[用户 solution.png 与公式截面对照](figures/fidelity_formula/solution_vs_scalar_formula.pdf)。论文图用于趋势对照，定量参考是同一 $K=2A$ 约定下的完整模型。

**几何示意：**[平行输运标架图 PDF](figures/parallel_transport/parallel_transport_frame.pdf)（[PNG](figures/parallel_transport/parallel_transport_frame.png)、[读图说明](figures/parallel_transport/README.md)）展示控制参数路径、球面上的随动低能方向及闭环后的标架旋转；已插入简要汇报第 3 节。

**弱耦合作用：**[弱耦合分解图 PDF](figures/parallel_transport/weak_coupling_projection.pdf)（[PNG](figures/parallel_transport/weak_coupling_projection.png)、[读图与分量公式](figures/parallel_transport/README.md#弱耦合怎样作用于剩余自由度)）展示弱耦合连接辅助模式与作用于低能内部的不同分量，并画出三段有效哈密顿量的系数；已插入简要汇报第 3 节。

## 当前目录

```text
read/
  README.md
  perturbation.md
  analytic_continuation.md
  analytic_continuation_brief.md   # 解析式结果的简要汇报
  analytic_continuation_report.md  # 独立组会汇报及技术附录
  group_symmetry_control.md       # 群对称控制方案、研究价值及待验证问题
  summary_2026-11-01.md        # 本次工作总结与验证状态
  code/                       # 公式、完整模型与参数响应模块
  data/                       # 数值数组、参数与验证记录
  figures/fidelity_formula/   # 两张主结果图、通道诊断图与原文面板
  figures/parallel_transport/ # 平行输运标架、弱耦合分解图及说明
  literature/
    2008.md                   # 2020 年 holonomic CNOT 论文汇报
    2008/                     # 文献图与原核对代码
```

[Majorana holonomic CNOT 汇报正文](literature/2008.md) 为约 **10–15 分钟**篇幅的完整报告，包含背景、创新、模型、控制方法与数值结果，并用示意图解释双隧穿机制。论文发表于 2020 年，文件名来自 arXiv:2008.05181；它是理论方案与数值模拟论文。该汇报独立于 PRB 111 的研究推导。[原模型记录](../model.md)、[PRB 111 原论文](../arts/prb111.pdf)、[用户源图](../solution.png) 保留原位。

配套幻灯片为 **13 页主报告 + 2 页备用内容**，包含 Fig. 3(b) 的能谱解释、四模作用到 CNOT 的电路构造，以及附录 A、B 的简述：[Overleaf 上传包](literature/2008/2008_cnot_overleaf.zip) · [TeX 源文件](literature/2008/slides/main.tex) · [编译说明](literature/2008/slides/README.md)。使用 XeLaTeX 编译。

## 复现

在项目根目录执行；依赖和数据说明见解析延拓文档第 13 节。

```bash
MPLCONFIGDIR=/tmp/quantumsss_fidelity_mpl .venv/bin/python read/code/plot_fidelity_formula.py
MPLCONFIGDIR=/tmp/quantumsss_fidelity_mpl .venv/bin/python read/code/compare_solution_cuts.py
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python read/code/check_continuation_fourth_order.py
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python read/code/projection_error_bound.py
```

已有数组会复用；完整模型数组不存在时，`read/code/full_model.py` 可从矩阵函数重算。图片像素不是作者原始数据；JSON 中既有的像素比较保留来源标注。

## 本次整合的范围

原 `r1.md`–`r3.md` 的微扰和当前门控协议下的对称性推导纳入第一篇；`r3.md`–`r8.md` 的收敛级数、标量函数、矩阵函数、延拓和机制纳入第二篇；`r9.md`–`r11.md` 的有效验证、近似与误差分析并入第二篇后半部。原 `read.md` 移至文献区。

已移除对话式重复草稿、论文颜色反读的旧截面图、整页渲染缓存和 Python 缓存；独立的 Heun 归约、额外控制协议、未继续的复时间连接分析等探索分支不再占用主目录。此前 `geometric_gates` 的全程恒定 $t_1$ 协议不与当前门控 $u_1(t)$ 混用；当前协议下已完成的抵消推导保留在微扰文档中。完整矩阵函数的已完成通道分解与参数导数构造保留在解析延拓文档第 10 节。

整理前研究材料已备份到 `/tmp/quantumsss-research-before-consolidation-20261101.tar.gz`，用于本次清理回退；它是临时备份，当前目录中的两份主文档与复现资源才是后续工作入口。
