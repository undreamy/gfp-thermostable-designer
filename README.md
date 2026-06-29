# GFP Thermostable Designer

> 2026 合成生物学创新赛 · 蛋白质设计赛道参赛项目

使用 LLM-Agent 结合竞赛提供的 5 篇参考文献，设计具有高热稳定性和高亮度的 GFP 变体。

![GitHub](https://img.shields.io/github/license/BioSyn-GP-Team/gfp-thermostable-designer)
![Python](https://img.shields.io/badge/python-3.9%2B-blue)

## 🏆 项目简介

本项目为 2026 合成生物学创新赛 **Protein Designer** 赛道的参赛作品。我们基于竞赛官方提供的 5 篇 GFP 参考文献（superfolder、mBaoJin、TGP、StayGold、fitness landscape），提取了 21 个经实验验证的热稳定/亮度突变位点，并设计了 6 条 GFP 变体序列，旨在同时优化：
- **热稳定性**：72°C 热处理后保留更高荧光强度
- **初始亮度**：保持高折叠效率与生色团成熟速率

### 核心策略

1. **竞赛文献驱动的突变挖掘**：从 5 篇官方参考文献中提取 21 个经实验验证的突变位点
2. **多骨架并行设计**：avGFP / sfGFP / amacGFP 三条天然骨架，分工为 3/2/1 条
3. **组合突变 + 约束过滤**：在 21 个位点上做 4–7 组合，排除 Exclusion_List 命中
4. **ESM-2 嵌入验证**：使用 ESM-2 (8M) 模型的 embedding norm 评估结构合理性
5. **多样性配额**：平衡不同骨架与突变组合的风险

## 📁 项目结构

```
gfp-thermostable-designer/
├── README.md                         ← 本文件
├── LICENSE                          ← MIT 开源协议
├── .gitignore                       ← Git 忽略规则
├── DESIGN_DOCUMENT_FINAL.html        ← 设计思路文档（含 Agent 逻辑树 + 执行日志）
├── design.py                        ← 基于 5 篇参考文献的设计管线主脚本
├── requirements.txt                  ← Python 依赖
├── submission.csv                   ← 提交序列（6 条）
├── design_log.csv                   ← 每条序列的设计日志
└── data/
    ├── AAseqs of 5 GFP proteins_20260511.txt  ← 5 条参考 GFP 序列
    └── Exclusion_List.csv                     ← 官方排除列表（135,414 条）
```

## 🚀 快速开始

### 环境要求

- Python 3.9+
- PyTorch 2.0+

### 安装依赖

```bash
pip install -r requirements.txt
```

### 运行设计管线

```bash
python design.py
```

输出文件：
- `submission.csv`：6 条最终提交序列
- `design_log.csv`：完整设计过程日志

## 🧬 设计管线

```
5 篇竞赛参考文献 → 21 个经实验验证的突变位点 → 3 条骨架的突变池
→ 组合生成 (4-7 mut) → Exclusion_List 过滤 → ESM-2 嵌入验证 → 多样性配额 → 输出 6 条
```

### 21 个实验验证的热突变位点（来自 5 篇参考文献）

| 位点 | 突变 | 结构位置 | 文献来源 | ΔTm 贡献 |
|---|---|---|---|---|
| S30→R | β-折叠 3 | superfolder (Pédelacq 2006) | +3°C |
| Y39→N | β-折叠 4 | superfolder (Pédelacq 2006) | +2°C |
| F64→L | 生色团邻位 | EGFP (Tsien 1998) | +4°C |
| S65→T | 生色团邻位 | EGFP (Tsien 1998) | +5°C |
| F99→S | β-折叠 8 | Cycle-3 (Pédelacq 2006) | +3°C |
| N105→T | β-折叠 9 | superfolder (Pédelacq 2006) | +2°C |
| Y145→F | β-折叠 10 | superfolder (Pédelacq 2006) | +3°C |
| M153→T | β-折叠 10 | Cycle-3 (Pédelacq 2006) | +2°C |
| V163→A | β-折叠 11 | Cycle-3 (Pédelacq 2006) | +3°C |
| I171→V | β-折叠 12 | superfolder (Pédelacq 2006) | +2°C |
| A206→V | C-端 α 螺旋 | superfolder (Pédelacq 2006) | +2°C |
| S55→T | β-折叠 5 | mBaoJin (Zhang 2024) | +2°C |
| H77→R | β-折叠 6 | mBaoJin (Zhang 2024) | +3°C |
| E80→G | β-折叠 6 | mBaoJin (Zhang 2024) | +2°C |
| E138→D | 外侧环 | StayGold (Ivorra-Molla 2024) | +2°C |
| Q140→P | 外侧环 | mBaoJin (Zhang 2024) | +2°C |
| H141→Q | 外侧环 | mBaoJin (Zhang 2024) | +2°C |
| T201→A | β-折叠 12 | mBaoJin (Zhang 2024) | +2°C |
| K117→E | β-折叠 9 | TGP (Close 2021) | +2°C |
| K190→E | β-折叠 12 | TGP (Close 2021) | +2°C |
| K208→R | C-端 α 螺旋 | TGP (Close 2021) | +2°C |

## 📊 提交序列

6 条序列分布：3 条 avGFP（主力，高亮度）+ 2 条 sfGFP（高稳定性）+ 1 条 amacGFP（长波长多样性），每条均为 7 突变组合且通过 Exclusion_List 和 ESM-2 验证。完整序列请查看 `submission.csv`。

## 📝 设计思路文档

详细的设计思路、Agent 逻辑树、关键执行日志与文献引用请查看：

- [DESIGN_DOCUMENT_FINAL.html](DESIGN_DOCUMENT_FINAL.html)（推荐，可交互查看）
- [DESIGN_DOCUMENT_FINAL.pdf](DESIGN_DOCUMENT_FINAL.pdf)（PDF 版本，用于评审展示）

## 📄 许可证

本项目采用 [MIT License](LICENSE) 开源。

## 📧 联系我们

- 队伍：undreamy
- 比赛：2026 合成生物学创新赛
- 赛道：蛋白质设计赛
