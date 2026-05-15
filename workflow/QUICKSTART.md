# 快速开始指南

## 文件夹内容

```
workflow/
  ├── README.md                      # 完整说明文档
  ├── config.py                      # 全局配置（修改参数在此）
  ├── utils.py                       # 工具函数库
  ├── run_full_workflow.py           # 主程序：运行所有步骤
  ├── step1_eight_vertex.py          # 步骤1：Pauli展开与参数映射
  ├── step2_bloch_rotation.py        # 步骤2：Bloch旋转与子空间分析
  ├── step3_full_chain_bdg.py        # 步骤3：全链BdG与边缘定位
  ├── step4_mzm_abs_criteria.py      # 步骤4：综合MZM/ABS判断
  └── results/                       # 输出结果
      ├── step1_eight_vertex/
      ├── step2_bloch_rotation/
      ├── step3_full_chain/
      └── step4_summary/
```

## 一键运行

```bash
cd /home/asice-cloud/projects/pyyy/quantumsss/workflow
/home/asice-cloud/projects/pyyy/quantumsss/.venv/bin/python run_full_workflow.py
```

或简写（需已激活虚拟环境）：
```bash
python run_full_workflow.py
```

## 单步运行

### 运行指定某一步

```bash
# 步骤1：修改参数并运行 eight-vertex 模型
python step1_eight_vertex.py --u-list "0,0.5,1.57" --delta-list "0,0.05,0.1"

# 步骤2：Bloch 旋转分析
python step2_bloch_rotation.py --u-list "0,1.57" --delta-list "0,0.1"

# 步骤3：全链 BdG 分析，指定链长
python step3_full_chain_bdg.py --u-list "0.78" --delta-list "0" --L-list "40,80,160"

# 步骤4：综合诊断
python step4_mzm_abs_criteria.py --u-list "0,1.57" --delta-list "0,0.1"
```

## 修改参数

编辑 `config.py`：

```python
# 路径参数 u 的取值范围
U_LIST_DEFAULT = [0.0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]

# 耦合参数 δ 的取值范围
DELTA_LIST_DEFAULT = [0.0, 0.015, 0.1]

# 全链长度（用于 E_0(L) 扫描）
L_LIST_DEFAULT = [40, 80, 160, 320]
```

修改后重新运行对应步骤或主程序。

## 输出文件说明

### Step 1: Eight-Vertex 模型
- `eight_vertex_data.npz`：Pauli 系数、d 向量、Kitaev 参数
- 关键指标：
  - `d_z = δ/2`（局部有效子空间的纵向分量）
  - `t, Δ, μ`（全链 BdG 参数，注意 μ 总为 0）

### Step 2: Bloch 旋转
- `bloch_rotation_data.npz`：能级分裂、轨迹几何
- `bloch_rotation.png`：|d| 与能级分裂随 δ 的变化
- 关键观察：δ 增大 → d_z 增大 → 轨迹被抬离 xy 平面

### Step 3: 全链 BdG
- `full_chain_data.npz`：E_0(L)、LDOS、bulk gap 等
- `E0_vs_L.png`：E_0 随 L 的指数衰减（MZM 特征）或振荡（ABS 特征）
- `ldos_*.png`：零能 LDOS 热图（端点峰值 = MZM；内部峰值 = ABS）

### Step 4: MZM/ABS 综合诊断
- `mzm_abs_summary.txt`：每个 (u,δ) 点的诊断结果
- `mzm_abs_summary.png`：汇总表（绿色 = MZM，黄色 = ABS，红色 = 不确定）

## 关键结论

**δ 的作用**：
- 仅影响本地投影子空间的 `d_z` 分量
- **不影响**全链拓扑（因为 μ = 0 恒成立）
- 因此不改变 E_0(L) 的指数衰减，不改变 LDOS 分布

**MZM vs ABS 的区分**：
| 判据 | MZM | ABS |
|------|-----|-----|
| Bulk gap | 满足 \|μ\|<2\|t\| | - |
| E_0(L) 曲线 | 指数趋零 | 有限分裂或振荡 |
| LDOS 集中 | 端点 | 内部特定位置 |
| Bloch d_z | 0（δ=0 时）| 非零（δ≠0 时）|

## 扩展与修改

1. **改变 eight-vertex 模型**：在 `step1_eight_vertex.py` 中修改 `eight_vertex_H4()` 函数。
2. **改变 Kitaev 映射**：在 `utils.py` 中修改 `map_to_kitaev()` 函数。
3. **添加新的诊断指标**：在 `step4_mzm_abs_criteria.py` 中扩展 `diagnose_mzm_abs()` 函数。
4. **自定义可视化**：编辑各 `step*.py` 中的绘图部分。

## 故障排除

**问题**：ModuleNotFoundError: No module named 'numpy'
- 解决：在虚拟环境中运行，或 `pip install numpy scipy matplotlib`

**问题**：输出文件很大
- 解决：在 `config.py` 中减少 `U_LIST_DEFAULT`、`DELTA_LIST_DEFAULT` 或 `L_LIST_DEFAULT` 的元素个数

**问题**：图像没有生成
- 检查：`config.py` 中 `SAVE_FIGURES = True` 且 `results/workflow/` 可写

## 参考文献

- 完整理论推导：见 `../result.md`
- eight-vertex 模型详解：见 `../ver7.md`
- 脚本分层说明：见 `../ver7.md` 第 8 节
