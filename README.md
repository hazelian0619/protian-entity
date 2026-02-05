# Protian Entity — Human Protein & RNA Knowledge Graph (Industrial Data Product)

This repository contains curated human Protein and RNA entity datasets and the code, contracts, and QA artifacts required to build, validate, and release them as industrial-grade data products.

Key principles
- Code, contracts, and QA reports are tracked in Git for auditability and reproducibility.
- Large data artifacts (L1 tables) are published via GitHub Releases with a manifest and checksums.

Quick links
- Protein (L1) dataset: `data/processed/protein_master_v6_clean.tsv`
- RNA (L1) dataset: release `rna-l1-v1` (see `pipelines/rna/README.md`)
- Validation tool: `tools/kg_validate_table.py`

Repository layout
- `data/processed/` — final, curated TSV tables (small-to-medium L1 tables are stored here when size permits)
- `pipelines/` — extraction and ETL pipelines (e.g., `pipelines/rna/`)
- `docs/` — design documents, data dictionary, quality gate definitions
- `scripts/`, `tools/` — helper and validation scripts

Data release model
1. Build entity tables using `pipelines/` scripts in a reproducible environment.
2. Produce `manifest.json` (checksums, row counts, git commit, build timestamp).
3. Publish artifacts as a GitHub Release and attach QA reports.

Getting started
1. Clone repository
2. Review `docs/DATA_DICTIONARY.md` and `docs/QUALITY_GATES.md` for schema and validation rules
3. Run validation example:

```bash
python3 tools/kg_validate_table.py --contract pipelines/protein/contracts/protein_master_v6.json \
  --table data/processed/protein_master_v6_clean.tsv --out build/validate/protein_master_v6_report.json
```

Contributing
- See `CONTRIBUTING.md` for development, testing, and release workflow.

License
- MIT License — see `LICENSE` for details.

Contact
- Repository owner: hazelian0619
- Project maintenance: see `CODEOWNERS` or `CONTRIBUTING.md` for maintainers and contact instructions


---

### 🎯 设计原则

- ✅ **一级信息为主**：提取原始数据，不做推断
- ✅ **以蛋白质为主体**：每个 UniProt ID 一行
- ✅ **保留完整原文**：功能描述含证据代码和 PubMed 引用
- ✅ **多源整合**：7 个主要生物数据库

---

### 📋 项目状态

**阶段**：✅ 完成  \
**版本**：v6_clean  \
**更新**：2025-10-27

**数据源**：UniProt | AlphaFold | HGNC | STRING | GO | PDB  \
**时效性**：截止 2025-10-26

---

## 🧬 RNA 实体（L1, v1）

RNA（miRNA + mRNA/transcript）属于 L1 实体表；输出体积较大（单文件 >100MB），所以：

- 数据产物不直接 commit 进仓库（避免触发 GitHub 单文件限制、避免仓库膨胀）
- 统一通过 **GitHub Releases** 发布，并附带可核验的 `manifest.json` 与 QA 报告

- Release: https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1
- 使用说明：`pipelines/rna/README.md`
- 规范：`docs/rna/README.md`
