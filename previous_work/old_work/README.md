# 旧工作归档

归档日期：2026-10-11。此前这轮分析的 23 个文档、脚本、图像和数据文件原样移入本目录；文件清单与 SHA-256 校验值见 [archive_manifest.json](archive_manifest.json)。

新的研究入口是 [geometric_state_analysis.md](../geometric_state_analysis.md)。模型与论文来源仍位于 [model.md](../../model.md)、[prb111.pdf](../../arts/prb111.pdf) 和 [pq.txt](../../arts/pq.txt)。旧文档中直接写出的 `model.md`、`arts/...` 路径均按项目根目录理解。

## 阅读索引

| 阶段 | 文档 | 配套文件与适用范围 |
|---|---|---|
| 最初的几何分析 | [213.md](213.md) | 包含特定门控假设及近似分析；不能直接代替论文协议。 |
| 局部保真度展开 | [so5_fidelity_analysis.md](so5_fidelity_analysis.md) | `verify_so5_fidelity.py`、`so5_fidelity_validation.json`；抵消比例只在相应展开条件下成立。 |
| 展开公式的图像 | 同上 | `plot_so5_fidelity_formula.py`、`so5_fidelity_formula_fig1d.*`；图像使用近似公式，其适用范围见原文。 |
| 完整有限参数保真度 | [so5_full_fidelity_analysis.md](so5_full_fidelity_analysis.md) | `so5_full_fidelity.py`、`so5_full_fidelity_fig1d.*`、`so5_full_fidelity_slices.*`，以及完整演化的验证与参数分析 JSON。 |
| 参数响应与局部峰线 | [so5_parameter_pde.md](so5_parameter_pde.md) | `so5_parameter_response.py`、`so5_parameter_response.json`；精确响应方程不等于全空间峰线的解析解。 |
| 标量 PDE 的闭合问题 | [so5_global_scalar_pde.md](so5_global_scalar_pde.md) | `verify_scalar_observability.py`、`so5_scalar_observability.json`；讨论高阶消元和保留态坐标的输运方程。 |

脚本的相互导入仍在本目录内，图像链接也保持有效。可在项目根目录以 `.venv/bin/python previous_work/old_work/脚本名.py` 运行；脚本生成的数据和图像会写入本归档目录，重跑前应注意相应输出会更新。
