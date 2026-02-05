# Contributing

Thank you for contributing to Protian Entity. This document explains development, validation, and release practices used by the project.

Development workflow
- Branching: create a branch using `chore/` or `feat/` prefixes for non-breaking changes and feature work respectively (e.g., `chore/standardize-docs-ci`).
- Tests & validation: run validation before opening a PR:

```bash
python3 tools/kg_validate_table.py --contract pipelines/protein/contracts/protein_master_v6.json \
  --table data/processed/protein_master_v6_clean.tsv --out build/validate/protein_master_v6_report.json
```

- Commit messages: use clear, imperative messages. Follow conventional commits if possible.

Pull requests
- Open PRs against `main`. Describe the change, test steps, and link to any data releases.
- Include validation reports for any changes to entity tables.

Releases
- Release data artifacts (large L1 tables) via GitHub Releases.
- Attach `manifest.json` with checksums, row counts, git commit SHA, and QA reports.

Contacts
- For maintenance and code ownership see `CODEOWNERS` or raise an issue.
