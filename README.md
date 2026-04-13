# Protian Entity — Human Protein & RNA Knowledge Graph

Protian Entity is an industrial-grade data product: curated L1 entity tables (Protein, RNA) with reproducible ETL, validation contracts, and QA artifacts.

Quick summary
- Primary entity: human Protein L1 table — `data/processed/protein_master_v6_clean.tsv` (v6 snapshot)
- RNA artifacts are published as release assets (`rna-l1l2-v2`) with `manifest.json` and QA reports
- Validation tool: `tools/kg_validate_table.py`

Badges
- CI: `.github/workflows/data-qa.yml` runs data validation on PRs and pushes validation reports as artifacts

Contents
- `data/processed/` — curated TSV L1 tables
- `pipelines/` — ETL code and contracts
- `docs/` — schema, data dictionary, quality gates, and admin guides
- `tools/` — validation and QA helpers

Quickstart (validate current protein table)
```bash
# run the table validator and write a report
python3 tools/kg_validate_table.py \
  --contract pipelines/protein/contracts/protein_master_v6.json \
  --table data/processed/protein_master_v6_clean.tsv \
  --out build/validate/protein_master_v6_report.json

# open the JSON report
less build/validate/protein_master_v6_report.json
```

CI / Non-interactive usage
- CI runs must be non-interactive. Tools and scripts must support CI mode by either a `--yes` / `--ci` flag or `CI=true` environment variable.
- To run locally in non-interactive mode: `CI=true ./scripts/run_full_build.sh --yes` (if present).
- See `docs/GITHUB_ACTIONS_SETTINGS.md` for repository-level Actions configuration (Admins only).

Data release process
1. Run pipelines to produce L1 artifacts and `manifest.json` (row counts, checksums, git commit, build timestamp).
2. Run all validation contracts and attach reports.
3. Create a GitHub Release and upload artifacts + `manifest.json`.

Governance & contribution
- See `CONTRIBUTING.md` for branching, validation, and release checklists.
- Use PR templates and attach validation report artifacts for any change that modifies tables or contracts.

What stays out of git
- Large raw inputs and big output artifacts should be published as release assets or stored in object storage — do NOT commit >100MB files.

License
- MIT (see `LICENSE`)

Contact
- Repo owner: `@hazelian0619` — open an issue or PR for changes

---

This README provides a concise starting point. For field-level schema and examples, see `docs/DATA_DICTIONARY.md` and `data/processed/README.md`.
