# 既有工作归档

归档日期：2026-10-11。本目录集中保存建立固定协议 PDE 之前的研究，包括此前几何分析、验证代码、数据、图像及更早的旧工作，共 72 个已有文件。

当前研究入口是 [fixed_protocol_pde.md](../fixed_protocol_pde.md)。模型与论文来源保留在 [model.md](../model.md)、[prb111.pdf](../arts/prb111.pdf) 和 [pq.txt](../arts/pq.txt)。

| 内容 | 入口 |
|---|---|
| 完整状态量与保真度的几何定义 | [geometric_state_analysis.md](geometric_state_analysis.md) |
| 六段精确演化及接口 | [geometric_full_protocol.md](geometric_full_protocol.md) |
| 相干分量与第一、第三段分析 | [geometric_coherence_dynamics.md](geometric_coherence_dynamics.md) |
| 第二段的精确分块与记忆 | [geometric_second_stage.md](geometric_second_stage.md) |
| 移动标架和慢球面近似 | [geometric_moving_frame.md](geometric_moving_frame.md) |
| 固定输入的换号边界及反馈修正 | [geometric_boundary_validation.md](geometric_boundary_validation.md) |
| 第二周期匹配与反向目标 | [geometric_second_cycle_chain.md](geometric_second_cycle_chain.md) |
| 完整模型的参数区域与峰线 | [geometric_parameter_space.md](geometric_parameter_space.md) |
| 完整周期的角度近似与条带预测 | [geometric_stripe_prediction.md](geometric_stripe_prediction.md) |
| 此前更早的展开、保真度和 PDE 分析 | [old_work/README.md](old_work/README.md) |

几何验证脚本和输出统一保存在 `geometric_validation/`，相互导入及数据相对位置保持不变。从项目根目录可运行，例如：

```sh
.venv/bin/python -B previous_work/geometric_validation/predict_geometric_stripes.py
```

归档移动后已逐文件核对内容一致性，再修正文档链接和复现命令。原位置、归档位置及修改前后的 SHA-256 记录见 [archive_manifest.json](archive_manifest.json)。早期归档内的清单仍保留其历史时点的记录。

当前设定固定每段时长为 100、全部门控及其余模型设置，只研究两个不可控幅度参数 \(E_1,t_1\)。旧文档中其他时间尺度的实验属于历史记录。
