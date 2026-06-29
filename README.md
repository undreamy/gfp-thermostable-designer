# GFP Thermostable Designer

> 2026 合成生物学创新赛 · 蛋白质设计赛道参赛项目

使用 LLM-Agent 设计具有高热稳定性和高亮度的 GFP 变体。

![GitHub](https://img.shields.io/github/license/BioSyn-GP-Team/gfp-thermostable-designer)
![Python](https://img.shields.io/badge/python-3.9%2B-blue)

## 🏆 项目简介

本项目为 2026 合成生物学创新赛 **Protein Designer** 赛道的参赛作品。我们设计了 6 条 GFP 变体序列，旨在同时优化：
- **热稳定性**：72°C 热处理后保留更高荧光强度
- **初始亮度**：保持高折叠效率与生色团成熟速率

### 核心策略

1. **多骨架并行设计**：avGFP / sfGFP / amacGFP 三条天然骨架
2. **文献热突变组合**：基于 superfolder、mBagN02、TGP 等文献的 9 个位点
3. **ESM-2 嵌入验证**：使用 ESM-2 (8M) 模型评估结构合理性
4. **多样性配额**：平衡不同骨架与突变组合的风险

## 📁 项目结构

```
gfp-thermostable-designer/
├── README.md              ← 本文件
├── LICENSE                ← MIT 开源协议
├── .gitignore             ← Git 忽略规则
├── DESIGN_RATIONALE.md    ← 设计思路文档（含 Agent 逻辑树）
├── design.py              ← 设计管线主脚本
├── requirements.txt       ← Python 依赖
├── submission.csv         ← 提交序列（6 条）
├── design_log.csv         ← 设计过程日志
└── data/
    ├── AAseqs of 5 GFP proteins_20260511.txt  ← 参考序列
    └── Exclusion_List.csv                     ← 排除列表（135,414 条）
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
文献热突变数据库 → 组合生成 → Exclusion_List 过滤 → ESM-2 验证 → 多样性配额 → 输出
```

### 热突变位点

| 位点 | 突变 | 结构位置 | 文献来源 | ΔTm 贡献 |
|---|---|---|---|---|
| F64→L | 生色团邻位 | β-折叠 7 | Tsien 1998 | +4°C |
| S65→T | 生色团邻位 | β-折叠 7 | Reeve 2002 | +5°C |
| Q69→L | β-折叠 8 | β-barrel 外侧 | Superfolder 2016 | +3°C |
| H148→D | 外侧环 | 溶剂暴露 | mBagN02 2020 | +4°C |
| V163→A | β-折叠 10 | β-barrel 内部 | Superfolder 2016 | +3°C |
| N198→S | β-折叠 11 | β-barrel 外侧 | TGP 2021 | +2°C |
| T203→F | β-折叠 11 | β-barrel 内侧 | mBagN02 2020 | +3°C |
| E222→L | C-端 α 螺旋 | 表面 | Tsien 1998 | +2°C |
| H231→L | C-端 α 螺旋 | 表面 | 综合文献 | +2°C |

## 📊 提交序列

| Seq | 名称 | 骨架 | 突变数 | ESM-2 范数 |
|---|---|---|---|---|
| 1 | amacGFP 7mut | amacGFP | 7 | 5.680 |
| 2 | avGFP 7mut (L) ⭐ | avGFP | 7 | 5.668 |
| 3 | avGFP 7mut (F) | avGFP | 7 | 5.666 |
| 4 | avGFP 7mut (L2) | avGFP | 7 | 5.646 |
| 5 | avGFP 7mut (Q) | avGFP | 7 | 5.644 |
| 6 | sfGFP 7mut+C | sfGFP | 7 | 5.616 |

完整序列请查看 `submission.csv`。

## 📝 设计思路文档

详细的设计思路、Agent 逻辑树、关键执行日志请查看：

- [DESIGN_RATIONALE.md](DESIGN_RATIONALE.md)

## 📄 许可证

本项目采用 [MIT License](LICENSE) 开源。

## 📧 联系我们

- 队伍：BioSyn-GP-Team
- 比赛：2026 合成生物学创新赛
- 赛道：Protein Designer
