# 代码 ⇄ 结果映射 (Code → Results)

本文档列出仓库中主要脚本、用途、典型运行命令，以及对应生成的 `results/` 输出文件（图片或 `.npz` 数据）。

主要脚本映射

- `tools/compute_ldos.py`
  - 作用: 从八顶点映射构造有效 Kitaev 链并计算 LDOS(能量 vs 位置)热图与 E=0 空间剖面。
  - 典型命令:

```bash
python tools/compute_ldos.py --u 0 --delta 0.015 --L 160
```
  - 主要输出:
    - `results/ldos/ldos_u{u}_d{delta}_L{L}.npz`（包含 `energies`, `ldos`, `spec`）
    - `results/ldos/ldos_u{u}_d{delta}_L{L}.png`

- `tools/edge_localization.py`
  - 作用: 对最低 BdG 本征态做链长扫描，记录 `E0(L)` 与端点权重，绘制对数尺度的衰减与权重曲线。
  - 典型命令:

```bash
python tools/edge_localization.py --u 0 --delta 0.015 --Ls 40,80,160,320
```
  - 主要输出:
    - `results/edge_localization/edge_loc_u{u:.3f}_d{delta:.3g}.npz`（键: `L`, `E0`, `left_w`, `right_w`）
    - `results/edge_localization/edge_loc_u{u:.3f}_d{delta:.3g}.png`

- `tools/verify_path_braid.py`
  - 作用: 从 ver5.md 中的路径生成 4×4 门，按段积分构造总门，投影到逻辑子空间并计算与理想 braid 门/ YBE 的残差，保存诊断 `.npz` 与汇总图像。
  - 典型命令:

```bash
python tools/verify_path_braid.py --steps 400 --T 1.0 --tc 1.0 --axis x
```
  - 主要输出（默认 `--outdir results/path_braid`）:
    - `results/path_braid/<stem>.npz`（包含 R1,R2,R3,R_total,R_logical, residuals 等）
    - `results/path_braid/<stem>.png`（汇总柱状残差图 + 差异矩阵图像）

- `tools/scan_bloch_grid.py`
  - 作用: 扫描 `d0 × T` 网格，计算最终 Bloch 旋转角并保存热图与 `.npz`。
  - 典型命令:

```bash
python tools/scan_bloch_grid.py
```
  - 主要输出:
    - `results/bloch_scan_d0T.npz`（包含 `d0_list`, `T_list`, `angles`）
    - `results/bloch_scan_heatmap.png`

- `tools/reproduce_figs.py`
  - 作用: 高层驱动脚本，调用仓库内工具复现论文 Fig.1..Fig.5（会调用 `embed_kitaev`, `tetron_path_sim`, `bloch_rotation_sim` 等）。
  - 典型命令:

```bash
python tools/reproduce_figs.py
```
  - 输出目录: 由 `tools/paper_params.py` 中的 `OUTDIR` 控制（默认 `results`）。
  - 常见输出文件示例:
    - `results/reproduce_Fig1_ldos_panel.png`
    - `results/reproduce_tetron_MZM_*.npy` / `results/reproduce_Fig3_abs_eigen_vs_time.png` / `results/reproduce_Fig3_overlap_vs_T.png`
    - `results/reproduce_Fig4_modulated_eigs.png`, `results/reproduce_Fig5_modulation_amplitude.png`

其他常用模块（快速索引）

- `tools/bloch_rotation_sim.py`
  - 导出: `run_for_T(T_step, n_per_step, ...)`，被 `scan_bloch_grid.py` 与 `reproduce_figs.py` 使用。
  - 产生: 中间 Bloch 轨迹图与数据（由调用者保存）。

- `tools/tetron_path_sim.py`
  - 导出: `run_sim(...)` 与 `make_time_grid()` 等，用于 tetron 波函数演化模拟（Fig.2/3 的数据来源）。

- `tools/embed_kitaev.py`
  - 提供: 将门/控制参数映射到 BdG 链的构造与 snapshot LDOS 的工具，`reproduce_figs.py` 中调用 `snapshot_ldos`。

- 校准 / 拟合工具（映射参数 A0/B0/C0）
  - `tools/calibrate_mapping.py`, `tools/fit_mapping_ABC.py`, `tools/optimize_mapping.py`
  - 作用: 用多面板数据联合拟合将 Pauli → Kitaev 的尺度参数确定下来，输出拟合结果与诊断 `.npz/.pkl`。

结果索引（仓库中已存在的部分文件）
- `results/ldos/ldos_u1.571_d0.015_L160.png` 及对应 `.npz`
- `results/edge_localization/edge_loc_summary.npz` 与若干 `edge_loc_u*_d*.png/.npz`
- `results/path_braid/*.png` 与 `*.npz`
- `results/bloch_*`, `results/alpha_scan*.png`, `results/deriveR_report*.pkl/.npz`

提示与下一步建议
- 若要逐图重跑（可复制论文面板），先运行:

```bash
pip install -r requirements.txt
python tools/reproduce_figs.py
```

- 若只想单独重跑某些面板，请按上面脚本示例运行对应脚本。

- 我已把该文档保存在 `docs/code_to_results.md`。


---
（文档是初版；我可以继续把每个脚本的参数说明扩展为完整的 `--help` 摘要，并为每个论文图指定精确命令序列。需要我继续完善吗？）
