# 人类知识图谱数据集（Protein + RNA）

这个仓库按“工业级数据产品”的方式组织：

- **代码 / 规范 / 质量报告**进入仓库（可审计、可复现）
- **体积大的数据产物**通过 **GitHub Releases** 发布（可下载、可校验、可回滚）

## 快速入口（给同事看这一段就够）

- **Protein（L1）数据集**：`data/processed/protein_master_v6_clean.tsv`（仓库内可直接下载）
- **RNA（L1, v1）数据集**：Release `rna-l1-v1`（包含 `.tsv.gz` + `manifest.json` + QA 报告）
  - Release: https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1
  - RNA 使用说明：`pipelines/rna/README.md`
  - RNA 规范：`docs/rna/README.md`

---

（原 README 内容已保留为中文）
