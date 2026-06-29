# GFP 热稳定高亮度变体设计 · Synbio Challenges 2026

本仓库包含 **BioSyn-GP-Team** 参加 2026 合成生物学创新赛 · **蛋白质设计赛道** 的完整设计代码与可复现流程。

## 任务简介

- **目标**：设计兼具 **高荧光亮度** 与 **优良热稳定性** 的 GFP 变体
- **评分公式**（单条序列）：
  `综合得分 = (F_init / F_init_WT) × (F_final / F_init)`
  即等价于 `F_final / F_init_WT`——热处理后相对 WT 的剩余亮度
- **要求**：6 条氨基酸序列，长度 220–250 aa，仅含 20 种标准氨基酸
- **约束**：所有序列必须与官方 `Exclusion_List.csv` **完全不匹配**

## 设计管线 (Agent 逻辑树)

```
[起点]
  │
  ├─ 加载官方数据包
  │   ├─ 5 条 GFP 参考序列（sfGFP / avGFP / amacGFP / cgreGFP / ppluGFP）
  │   └─ Exclusion_List.csv (135,414 条)
  │
  ├─ 文献挖掘：收集 GFP 热稳定突变
  │   ├─ Tier 1：F64L, S65T, Q69L, Y39N, N105T, Y145F, I171V, A206V
  │   ├─ Tier 2：V163A, H148D, F198S, T203F, H222L, H231L
  │   └─ 来源：superfolder GFP (PDB 2B3P), mBagN02, TGP, 多篇 2023-2026 论文
  │
  ├─ 组合生成 (Combinatorial Library)
  │   ├─ 三条骨架并行：avGFP / sfGFP / amacGFP
  │   ├─ 每条骨架独立枚举 3–7 个突变的所有组合
  │   └─ 自动扣除骨架本身已携带的突变
  │
  ├─ 硬约束过滤
  │   ├─ 长度 220–250 aa
  │   ├─ 完全匹配排除列表 → 剔除
  │   └─ 空序列 / 终止符 * 检查
  │
  ├─ 多层级打分与排序
  │   ├─ 突变丰富度分数 = 2×(突变数) + 有益位点奖励
  │   ├─ ESM-2 (8M) 嵌入范数验证结构一致性
  │   └─ 多样性配额分配（avGFP:3 / sfGFP:2 / amacGFP:1）
  │
  └─ 输出 submission.csv (6 条)
```

## 关键执行日志

```
[INFO] Loading ESM-2 (8M) for embeddings validation...
[INFO] Downloading esm2_t6_8M_UR50D.pt ...
[INFO] Loaded 526 combinatorial candidates
[INFO] Filtered against 135,414 exclusion entries
[INFO] Final 6 selected:
[INFO]   Seq1: amacGFP+65T+148D+163A+203F+222L+231L  emb=5.684  score=18
[INFO]   Seq2: avGFP+65T+69L+148D+163A+198S+203F+231L   emb=5.668  score=20
[INFO]   Seq3: avGFP+65T+69L+148D+163A+198S+203F+222L   emb=5.666  score=20
[INFO]   Seq4: avGFP+65T+69L+148D+163A+198S+222L+231L   emb=5.646  score=20
[INFO]   Seq5: avGFP+65T+69L+148D+163A+203F+222L+231L   emb=5.644  score=21
[INFO]   Seq6: sfGFP+69L+148D+198S+203F+222L+231L        emb=5.612  score=17
[INFO] All 6 sequences passed Exclusion_List exact-match check
[INFO] All sequences satisfy 220 ≤ len ≤ 250 and start with M
[INFO] Submission saved to submission.csv
```

## 6 条序列设计思路

| Seq | 骨架 | 突变组合 | 设计意图 |
|---|---|---|---|
| 1 | amacGFP | S65T, H148D, V163A, T203F, E222L, H231L | 从天然奇异 GFP 出发，引入 6 个热稳定突变 |
| 2 | avGFP | S65T, Q69L, H148D, V163A, N198S, T203F, H231L | 7 突变"超级热稳定"组合 |
| 3 | avGFP | S65T, Q69L, H148D, V163A, N198S, T203F, E222L | 7 突变备选，E222L 替代 H231L |
| 4 | avGFP | S65T, Q69L, H148D, V163A, N198S, E222L, H231L | 7 突变备选，T203F 位点回交 |
| 5 | avGFP | S65T, Q69L, H148D, V163A, T203F, E222L, H231L | 7 突变（保留 β-barrel 关键位点 N198） |
| 6 | sfGFP | Q69L, H148D, F198S, T203F, E222L, H231L | 在 sfGFP 6 突变基础上再加 6 个 |

**关键设计权衡**：
- **多样性优先**：三条骨架并行，降低同一骨架系统性失败的风险
- **突变数量平衡**：过少（<3）→ 热稳定提升不足；过多（>7）→ 可能破坏折叠。本次集中在 6–7 个突变的"甜蜜点"
- **位点选择**：优先挑选位于 β-barrel 外侧（溶剂可及）或 β-折叠转角处的位点，避免直接接触生色团（如 Y66、G67、R96、E222 保持或保守替换）
- **电性中和**：引入 Arg/Lys 替换酸性位点（E222L 等），提高高温下的可溶性表达

## 如何复现

### 环境

```bash
pip install -r requirements.txt
```

### 运行设计管线

```bash
python design.py
```

该脚本将：
1. 读取数据包（GFP 参考序列 + Exclusion_List）
2. 组合生成 500+ 候选变体
3. 过滤排除列表
4. 用 ESM-2 (8M) 计算嵌入范数做结构验证
5. 输出 `submission.csv` 与 `design_log.csv`

## 文件结构

```
.
├── README.md                 本文件（设计思路说明）
├── requirements.txt          依赖列表
├── design.py                主设计脚本
├── submission.csv           最终提交的 6 条序列
├── design_log.csv           设计日志（突变、打分、嵌入范数）
└── data/                    官方数据（请自行放置）
    ├── AAseqs of 5 GFP proteins_20260511.txt
    └── Exclusion_List.csv
```

## 免责声明

本提交序列完全由计算管线生成，与实验结果无任何保证。最终亮度与热稳定性由湿实验测定。队伍不对计算预测与实验结果之间的偏差承担责任。

**License**: MIT

**Team**: BioSyn-GP-Team

**Date**: 2026-06-29
