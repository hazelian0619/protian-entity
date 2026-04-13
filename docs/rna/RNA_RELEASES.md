# RNA Releases

## Current recommended release

- **Tag**: `rna-l1l2-v2`
- **URL**: https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1l2-v2
- **Contents**:
  - 16 compressed RNA tables (`*.tsv.gz`)
  - 15 non-sample validation reports (`*.validation.json`, all PASS)
  - audit artifacts (`mrna_enst_urs_coverage_v2.json`, `rna_trna_features_v2.metrics.json`, etc.)
  - `manifest.json` + `SHA256SUMS.txt`

## Historical release

- **Tag**: `rna-l1-v1`
- **URL**: https://github.com/hazelian0619/protian-entity/releases/tag/rna-l1-v1

## Verification

1. Download `manifest.json` and `SHA256SUMS.txt` from the release page.
2. Verify checksums locally:

```bash
sha256sum -c SHA256SUMS.txt
```

3. Use the corresponding validation JSON reports in the same release for QA traceability.
