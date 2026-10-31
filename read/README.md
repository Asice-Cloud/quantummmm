# PRB 111 研究整理

主线现收敛为两份文档，建议按下表顺序阅读。

| 文档 | 内容 |
|---|---|
| [微扰、几何结构与对称性抵消](perturbation.md) | 精确能谱、快模消去、holonomy、二阶与四阶保真度、抵消比例、群平均与九圈二阶校正门 |
| [解析延拓、具体公式与验证](analytic_continuation.md) | 完整标量函数构造、局部递推与延拓、读出公式、峰位峰高与长时机制、完整 SO(5) 矩阵函数、图像核对、假设和误差预算 |

**工作总结：**[2026-11-01：公式、零耦合极限、PRB 105 对照与下一步](summary_2026-11-01.md)。

**图像入口：**[论文与公式热图](figures/fidelity_formula/paper_vs_formula.pdf)；[用户 solution.png 与公式截面对照](figures/fidelity_formula/solution_vs_scalar_formula.pdf)。论文图用于趋势对照，定量参考是同一 $K=2A$ 约定下的完整模型。

## 当前目录

```text
read/
  README.md
  perturbation.md
  analytic_continuation.md
  summary_2026-11-01.md        # 本次工作总结与验证状态
  code/                       # 公式、完整模型与参数响应模块
  data/                       # 数值数组、参数与验证记录
  figures/fidelity_formula/   # 两张主结果图、通道诊断图与原文面板
  literature/
    2008.md                   # 原论文组会解读
    2008/                     # 文献图与原核对代码
```

[2008 年论文组会解读](literature/2008.md) 独立保留，作为文献背景，不与 PRB 111 的研究推导混合。[原模型记录](../model.md)、[原论文](../arts/prb111.pdf)、[用户源图](../solution.png) 保留原位。

## 复现

在项目根目录执行；依赖和数据说明见解析延拓文档第 13 节。

```bash
MPLCONFIGDIR=/tmp/quantumsss_fidelity_mpl .venv/bin/python read/code/plot_fidelity_formula.py
MPLCONFIGDIR=/tmp/quantumsss_fidelity_mpl .venv/bin/python read/code/compare_solution_cuts.py
```

已有数组会复用；完整模型数组不存在时，`read/code/full_model.py` 可从矩阵函数重算。图片像素不是作者原始数据；JSON 中既有的像素比较保留来源标注。

## 本次整合的范围

原 `r1.md`–`r3.md` 的微扰和当前门控协议下的对称性推导纳入第一篇；`r3.md`–`r8.md` 的收敛级数、标量函数、矩阵函数、延拓和机制纳入第二篇；`r9.md`–`r11.md` 的有效验证、近似与误差分析并入第二篇后半部。原 `read.md` 移至文献区。

已移除对话式重复草稿、论文颜色反读的旧截面图、整页渲染缓存和 Python 缓存；独立的 Heun 归约、额外控制协议、未继续的复时间连接分析等探索分支不再占用主目录。此前 `geometric_gates` 的全程恒定 $t_1$ 协议不与当前门控 $u_1(t)$ 混用；当前协议下已完成的抵消推导保留在微扰文档中。完整矩阵函数的已完成通道分解与参数导数构造保留在解析延拓文档第 10 节。

整理前研究材料已备份到 `/tmp/quantumsss-research-before-consolidation-20261101.tar.gz`，用于本次清理回退；它是临时备份，当前目录中的两份主文档与复现资源才是后续工作入口。
