# 完整流程工作区

本文件夹整理了从**路径设计 → eight-vertex 模型 → 全链 BdG → Bloch 旋转 → MZM/ABS 判据**的完整数值流程。

## 目录结构

- `run_full_workflow.py`：总控脚本，一键运行所有步骤。
- `step1_eight_vertex.py`：eight-vertex 模型的 Pauli 展开与映射。
- `step2_bloch_rotation.py`：投影到逻辑子空间，计算 `d` 与 Bloch 旋转。
- `step3_full_chain_bdg.py`：推广到全链 BdG，计算 LDOS、`E_0(L)`、bulk gap。
- `step4_mzm_abs_criteria.py`：综合判断 MZM / ABS 的各项指标。
- `config.py`：统一配置（路径参数、链长、参数扫描范围等）。
- `utils.py`：公共工具函数（Pauli 映射、嵌入、对角化等）。

## 快速开始

```bash
# 运行完整流程（所有参数用默认值）
python run_full_workflow.py

# 运行指定某个步骤
python step1_eight_vertex.py --u-list "0,1.57,3.14" --delta-list "0,0.015,0.1"
python step2_bloch_rotation.py --u 0.0 --delta 0.0
python step3_full_chain_bdg.py --u 0.0 --delta 0.0 --L-list "40,80,160"
python step4_mzm_abs_criteria.py --u 0.0 --delta 0.0
```

## 流程概览

### Step 1: Eight-vertex 模型
- 输入：路径参数 `(u, δ)`
- 输出：Pauli 系数 `h_{αβ}`, 有效两能级的 `d`，映射参数 `(t,Δ,μ)`
- 保存：`results/step1_eight_vertex_*.npz`

### Step 2: Bloch 旋转（投影子空间）
- 输入：有效两能级 `H_eff = d_0 I + d·σ`
- 输出：逻辑子空间的瞬时本征值 `E_eff(u) = ±|d|`，轨迹几何判断
- 保存：`results/step2_bloch_rotation_*.npz` 及图像

### Step 3: 全链 BdG
- 输入：映射参数 `(t,Δ,μ)` 和链长 `L`
- 输出：LDOS 热图、`E_0(L)` 指数衰减、bulk gap、拓扑指示
- 保存：`results/step3_full_chain_*.npz` 及图像

### Step 4: MZM/ABS 综合判断
- 输入：所有前三步的结果
- 输出：综合判断表（子空间 MZM-like/ABS-like vs 全链拓扑/边态判据）
- 保存：`results/step4_mzm_abs_summary.txt` 和对比图

## 输出文件说明

所有结果保存在 `results/` 下：
- `*.npz`：数据文件（用 `np.load()` 读取）
- `*.png`：图像（LDOS 热图、`E0(L)` 曲线、汇总对比表等）
- `*.txt`：文本汇总表

## 主要判据

### 子空间层（投影到 `{|01>, |10>}`）
- **MZM-like**：`δ=0`，轨迹 `d(u)` 在 `xy` 平面绕原点，几何相主导。
- **ABS-like**：`δ≠0`，`d_z` 非零，轨迹被抬离平面，出现能级分裂。

### 全链层（BdG 单粒子谱）
- **MZM**：
  - Bulk gap 满足 `|μ|<2|t|` 和 `gap>0`。
  - 开链 `E_0(L)` 随 `L` 指数趋零。
  - 零能 LDOS 在端点集中。
- **ABS**：
  - `E_0(L)` 有明显有限分裂或随参数振荡。
  - LDOS 在特定内部位置出现峰值。

## 参考文档

- `result.md`：完整理论总结。
- `ver7.md`：eight-vertex 推导与脚本分层说明。

## 扩展建议

若要修改模型或参数：
1. 修改 `config.py` 里的全局参数。
2. 在 `utils.py` 里加新的映射函数。
3. 编写新的 `step*.py` 或修改现有脚本。
4. 再运行 `python run_full_workflow.py`。
