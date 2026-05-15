# Workflow 索引 & 快速导航

## 📋 文件清单

| 文件 | 大小 | 功能 | 快速链接 |
|------|------|------|---------|
| **README.md** | 3.1K | 完整说明：理论背景、流程概览、输出说明 | 开始前必读 |
| **QUICKSTART.md** | 4.3K | 快速开始：常见命令、参数修改、故障排除 | 实际操作参考 |
| **STRUCTURE.md** | 8.0K | 文件夹结构：代码层级、数据流向、关键发现 | 理解架构 |
| **INDEX.md** | - | 本文档：导航与索引 | 当前 |
| | | | |
| **config.py** | 1.5K | 全局配置：参数范围、输出路径、选项 | 修改参数的地方 |
| **utils.py** | 4.9K | 工具库：Pauli、Kitaev 映射、BdG 构建、LDOS | 核心函数集合 |
| | | | |
| **run_full_workflow.py** | 1.9K | 主程序：一键运行 Step 1-4 | `python run_full_workflow.py` |
| **step1_eight_vertex.py** | 3.3K | 模型 → Pauli → Kitaev | 单步运行 |
| **step2_bloch_rotation.py** | 3.8K | 子空间分析与能级分裂 | 单步运行 |
| **step3_full_chain_bdg.py** | 6.7K | BdG 对角化、LDOS、bulk gap | 单步运行 |
| **step4_mzm_abs_criteria.py** | 7.5K | MZM/ABS 综合诊断 | 单步运行 |
| **check_results.py** | 2.0K | 结果检查脚本 | `python check_results.py` |

**总计代码量**: ~46 KB（Python + 文档）

---

## 🚀 常见操作

### 一键运行（完整流程，5×3=15个点，耗时 ~30 sec）
```bash
cd /home/asice-cloud/projects/pyyy/quantumsss/workflow
python run_full_workflow.py
```
**输出**: 
- `results/` 目录下的 4 个子文件夹
- 数据文件 `.npz` + 图像 `.png`
- 诊断汇总 `step4_summary/mzm_abs_summary.txt`

### 查看结果
```bash
python check_results.py          # 摘要
cat results/step4_summary/mzm_abs_summary.txt    # 完整诊断
```

### 修改参数后重新运行
```bash
# 编辑 config.py
nano config.py
# 修改:
# U_LIST_DEFAULT = [...]
# DELTA_LIST_DEFAULT = [...]
# L_LIST_DEFAULT = [...]

# 重新运行
python run_full_workflow.py
```

### 单步调试
```bash
# 例: 仅运行 Step 3 的特定参数
python step3_full_chain_bdg.py --u-list "0,1.57,3.14" --L-list "80,160"
```

---

## 📖 学习路径

### Path A: 快速上手 (~5 min)
1. 读 README.md 的"流程概览"部分
2. 运行 `python run_full_workflow.py`
3. 看 `results/step4_summary/mzm_abs_summary.txt`
4. 查看生成的图像

### Path B: 理解原理 (~20 min)
1. 读 README.md 全文
2. 读 STRUCTURE.md 的"数据流向"部分
3. 浏览 `step1_eight_vertex.py` 的 `eight_vertex_H4()` 函数
4. 读 QUICKSTART.md

### Path C: 修改与扩展 (~1 hour)
1. 完成 Path B
2. 阅读 `utils.py` 中的函数定义
3. 修改 `config.py` 的参数，重新运行
4. 理解 `step3_full_chain_bdg.py` 中 LDOS 的计算
5. 尝试在 `step4_mzm_abs_criteria.py` 中添加新的诊断指标

---

## 🔍 关键概念速查

### δ (delta) 是什么？
- Eight-vertex 模型的耦合参数
- 控制局部有效子空间中 $d_z = δ/2$
- **不影响**全链拓扑（因为全链参数 μ ≡ 0）
- 影响子空间的能级分裂和几何

### 四层架构是什么？
```
Layer 1: Eight-Vertex H_4(u, δ)     ← 四维单粒子 Hamiltonian
Layer 2: Pauli 展开 h_αβ             ← 4×4 Pauli 基展开
Layer 3: 有效 d 向量                 ← 二能级子空间投影
Layer 4: Kitaev 映射 (t,Δ,μ)        ← 全链 BdG 参数
```

### MZM vs ABS 如何判断？
- **MZM**: E_0(L) 指数衰减，LDOS 端点集中，gap > 0
- **ABS**: E_0(L) 有限分裂/振荡，LDOS 内部峰值，任意 gap

### 子空间的 d_z ≠ 0 意味着什么？
- 本地轨迹 $d(u)$ 被抬离 xy 平面
- 出现 ABS-like 的能级分裂
- **但全链拓扑不变**（因为 μ 恒为 0）

---

## 📊 输出文件说明

### Step 1: eight_vertex_data.npz
```python
import numpy as np
data = np.load('results/step1_eight_vertex/eight_vertex_data.npz', allow_pickle=True)
u_list = data['u_list']
delta_list = data['delta_list']
results = data['results'].item()  # Dictionary: (u,δ) → {h, d, d0, t, Delta, mu, ...}
```

### Step 2: bloch_rotation_data.npz
```python
# 同样格式，额外存储:
# E_+, E_-, splitting, |d|
```

### Step 3: full_chain_data.npz
```python
# 包含所有 L 值的数据:
# E_0_data: {L → {E_0, edge_weight, spectrum}}
# ldos: (L_max, NE) 数组
# ldos_energies: 能量轴
```

### Step 4: mzm_abs_summary.txt
```
简单文本格式，每个 (u,δ) 点一行:
u=..., δ=...
  Subspace: MZM-like or ABS-like
  Global:   Majorana Zero Mode or Andreev Bound State
  Evidence: ...
```

---

## ⚙️ 高级用法

### 自定义参数扫描
```bash
# 方法 1: 命令行
python step1_eight_vertex.py --u-list "0,0.5,1,1.5,2,2.5,3.14" --delta-list "0,0.05,0.1,0.15"

# 方法 2: 修改 config.py
# U_LIST_DEFAULT = np.linspace(0, np.pi, 21)
# DELTA_LIST_DEFAULT = np.linspace(0, 0.2, 11)
python run_full_workflow.py
```

### 扩展 BdG 模型
编辑 `utils.py` 的 `map_to_kitaev()` 和 `build_bdg_chain()` 函数

### 添加新的诊断指标
在 `step4_mzm_abs_criteria.py` 的 `diagnose_mzm_abs()` 中添加

---

## 🐛 常见问题

**Q: 为什么 μ 总是 0？**  
A: Eight-vertex 模型中，h_{ZI} + h_{IZ} = δ/2 - δ/2 = 0，所以映射 μ = 4h_{ZZ} - 2(h_{ZI}+h_{IZ}) = 0。

**Q: 为什么改变 δ 不改变 E_0(L) 的行为？**  
A: δ 只影响 d_z（子空间），不影响 (t,Δ,μ)（全链参数）。全链拓扑由 (t,Δ,μ) 决定。

**Q: LDOS 为什么有时端点有峰、有时内部有峰？**  
A: 取决于全链拓扑（由 u 决定的 t,Δ）。δ 无关。

**Q: 能否处理更大的 L 或更多的 (u,δ) 点？**  
A: 可以，但 `full_chain_data.npz` 会很大。建议修改 config.py 中的 NE_LDOS（LDOS 采样点数）以减少文件大小。

---

## 📝 参考与链接

- **完整推导**: 见项目根目录 `result.md` (8 节)
- **eight-vertex 详解**: 见 `ver7.md` (8 节)
- **历史开发笔记**: `ver6.md`, `ver5.md`, ... 

---

## 版本信息

| 部分 | 版本 | 最后更新 | 状态 |
|------|------|--------|------|
| 理论文档 | ver7.md | 完成 | ✓ 验证 |
| 工作流代码 | 1.0 | 刚完成 | ✓ 运行成功 |
| 默认参数 | 5u × 3δ × 4L | - | ✓ 测试通过 |

---

**下一步**: 根据上表选择学习路径，或直接运行 `python run_full_workflow.py`！
