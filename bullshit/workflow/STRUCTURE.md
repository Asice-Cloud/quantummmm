# 完整工作流文件夹结构

## 顶层目录

```
/home/asice-cloud/projects/pyyy/quantumsss/workflow/
│
├── 📄 README.md                     # 完整说明文档（理论背景 + 使用指南）
├── 📄 QUICKSTART.md                 # 快速开始指南（命令示例 + 参数修改）
├── 📄 STRUCTURE.md                  # 本文件 - 文件夹结构说明
│
├── 🐍 config.py                     # 全局配置
│                                    # - 参数范围设定（U_LIST, DELTA_LIST, L_LIST）
│                                    # - 输出目录配置
│                                    # - LDOS 参数、verbose 选项等
│
├── 🐍 utils.py                      # 工具函数库（270 行）
│                                    # - Pauli 矩阵与 Kronecker 积
│                                    # - H_4 到 Pauli 展开：pauli_expand()
│                                    # - Pauli 到 Kitaev 映射：map_to_kitaev()
│                                    # - d 向量提取：extract_d_vector()
│                                    # - BdG 链构建：build_bdg_chain()
│                                    # - 对角化、bulk gap、LDOS 等
│
├── 🐍 step1_eight_vertex.py         # 步骤 1: Eight-Vertex 模型
│                                    # 输入: (u, δ) 参数对
│                                    # 输出: H_4 → h_αβ → d → (t,Δ,μ)
│                                    # 保存: eight_vertex_data.npz
│                                    # 命令: python step1_eight_vertex.py --u-list "..." --delta-list "..."
│
├── 🐍 step2_bloch_rotation.py       # 步骤 2: Bloch 旋转（投影子空间）
│                                    # 输入: step1 的 d 向量
│                                    # 输出: 能级分裂 E_± = d_0 ± |d|
│                                    # 保存: bloch_rotation_data.npz + 图像
│                                    # 图像: |d| vs δ, 能级分裂 vs δ
│
├── 🐍 step3_full_chain_bdg.py       # 步骤 3: 全链 BdG 与边缘定位
│                                    # 输入: (t,Δ,μ) + 链长 L
│                                    # 输出: bulk gap, E_0(L), LDOS
│                                    # 保存: full_chain_data.npz + E0_vs_L.png + LDOS 热图
│                                    # 关键: 指数衰减 vs ABS 分裂判断
│
├── 🐍 step4_mzm_abs_criteria.py     # 步骤 4: 综合 MZM/ABS 诊断
│                                    # 输入: steps 1-3 的所有结果
│                                    # 输出: 每个 (u,δ) 的 MZM/ABS 判断
│                                    # 保存: mzm_abs_summary.txt + 汇总表
│                                    # 逻辑: 子空间性质 + 全链拓扑 + E0 行为 → 综合判断
│
├── 🐍 run_full_workflow.py          # 主程序：一键运行所有步骤
│                                    # 调用: step1() → step2() → step3() → step4()
│                                    # 用法: python run_full_workflow.py
│
├── 🐍 check_results.py              # 结果检查脚本
│                                    # 功能: 列出已生成的文件、大小、摘要
│                                    # 用法: python check_results.py
│
└── 📁 results/                      # 输出文件夹（由脚本自动创建）
    │
    ├── 📁 step1_eight_vertex/
    │   └── eight_vertex_data.npz   # 存储所有 (u,δ) 的 h, d, (t,Δ,μ)
    │
    ├── 📁 step2_bloch_rotation/
    │   ├── bloch_rotation_data.npz  # 存储 E_±, |d|, 分裂等
    │   └── bloch_rotation.png       # 图: |d| vs δ, 分裂 vs δ
    │
    ├── 📁 step3_full_chain/
    │   ├── full_chain_data.npz      # 存储所有 L 的 E_0, LDOS, gap 等（大文件）
    │   ├── E0_vs_L.png              # 图: E_0(L) 网格（每个 (u,δ) 一个子图）
    │   ├── ldos_u0.00_d0.000.png    # 图: LDOS 热图（示例）
    │   └── ...
    │
    └── 📁 step4_summary/
        ├── mzm_abs_summary.txt      # 文本: 诊断结果汇总表
        └── mzm_abs_summary.png      # 图: 彩色诊断汇总表（绿=MZM, 黄=ABS, 红=不确定）
```

## 代码层级关系

```
run_full_workflow.py (主程序，依次调用)
│
├─→ step1_eight_vertex.py
│   ├─ 调用: utils.eight_vertex_H4() → pauli_expand() → map_to_kitaev() → extract_d_vector()
│   └─ 输出: h_αβ, d, (t,Δ,μ)
│
├─→ step2_bloch_rotation.py
│   ├─ 依赖: step1 的 d 向量
│   ├─ 计算: E_± = d_0 ± |d|, 轨迹几何
│   └─ 输出: 能级分裂、Bloch 向量长度
│
├─→ step3_full_chain_bdg.py
│   ├─ 依赖: step1 的 (t,Δ,μ)
│   ├─ 调用: utils.build_bdg_chain() → diagonalize_bdg() → compute_bulk_gap() → compute_ldos()
│   └─ 输出: E_0(L), LDOS, gap, 拓扑判据
│
└─→ step4_mzm_abs_criteria.py
    ├─ 依赖: steps 1-3 的所有结果
    ├─ 逻辑: 子空间类型 + 全链拓扑 + E0 行为 → MZM/ABS 判断
    └─ 输出: 综合诊断表与图像
```

## 数据流向

```
(u, δ) 参数对
    ↓
[Step 1] eight_vertex_H4(u, δ) → Pauli 展开 → d, (t,Δ,μ)
    ├─ 子空间信息: d_z (影响 Bloch 轨迹)
    └─ 全链参数: t, Δ, μ (决定 BdG 拓扑)
    ↓
[Step 2] 有效子空间分析: d → E_± = d_0 ± |d|
    ├─ 判断: δ=0 → MZM-like (平面轨迹) vs δ≠0 → ABS-like (抬升轨迹)
    └─ 注意: 此层为局部投影，不决定全局拓扑
    ↓
[Step 3] 全链 BdG: (t,Δ,μ) + L → E_0(L), LDOS, bulk gap
    ├─ 拓扑判据: |μ| < 2|t| (here μ=0 always)
    ├─ 关键观察: E_0(L) 指数衰减 (MZM) vs 有限分裂 (ABS)
    └─ LDOS 模式: 端点峰值 (MZM) vs 内部峰值 (ABS)
    ↓
[Step 4] 综合诊断: 子空间类型 + 全链拓扑 + E0 行为 → 最终判断
    └─ 输出: Majorana Zero Mode (MZM) / Andreev Bound State (ABS) / 不确定
```

## 关键发现

**δ 的影响范围限制**：
- ✓ 影响: 子空间的 d_z = δ/2 分量 → 本地能级分裂
- ✗ 不影响: 全链 Kitaev 参数 (μ 恒为 0)
- ✗ 不影响: 全链拓扑 → E_0(L) 指数衰减pattern 不变

**MZM vs ABS 区分准则**：

| 层级 | MZM 特征 | ABS 特征 |
|------|---------|---------|
| **子空间 (Step 2)** | δ≈0，|d| 在 xy 平面 | δ≠0，|d| 抬离平面 |
| **全链拓扑 (Step 3)** | gap > 0，μ = 0 | 任意 μ |
| **E_0(L) 行为 (Step 3)** | 指数衰减 ∼ e^{-αL} | 有限分裂或振荡 |
| **LDOS 分布 (Step 3)** | 端点集中 | 内部峰值 |

## 使用建议

### 快速验证（<1 min）
```bash
python run_full_workflow.py          # 用默认参数（5×3=15 个点）
python check_results.py              # 查看摘要
```

### 细致扫描 (custom)
```bash
# 修改 config.py
U_LIST_DEFAULT = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.14]  # 7 个 u 点
DELTA_LIST_DEFAULT = [0, 0.01, 0.02, ..., 0.15]     # 16 个 δ 点
L_LIST_DEFAULT = [40, 80, 160, 320]                  # 4 个链长

# 再运行
python run_full_workflow.py
```

### 调试单个点
```bash
u=1.57; d=0.05; L=160
python step1_eight_vertex.py --u-list "$u" --delta-list "$d"
python step2_bloch_rotation.py --u-list "$u" --delta-list "$d"
python step3_full_chain_bdg.py --u-list "$u" --delta-list "$d" --L-list "$L"
python step4_mzm_abs_criteria.py --u-list "$u" --delta-list "$d" --L-list "$L"
```

## 文件大小预期

- `eight_vertex_data.npz`: 几 KB（参数存储）
- `bloch_rotation_data.npz`: 几 KB
- `full_chain_data.npz`: **大** (~18 MB for 15 points, L=[40,80,160,320])
  - 原因: LDOS 矩阵（L×NE 数组，每点多个 L，15×4=60 个矩阵）
- 图像: 每张 100-500 KB

## 相关文档

见工作区根目录：
- `result.md` - 完整理论总结（8 节）
- `ver7.md` - eight-vertex 推导详解（8 节）
- `ver6.md` 等 - 历史开发笔记

