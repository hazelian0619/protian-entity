#!/usr/bin/env bash
set -euo pipefail

# Publish molecules v1 to GitHub (code + Releases).
#
# What it does:
# 1) pushes branch `feat/molecules-v1`
# 2) optionally opens a PR to main
# 3) packages release assets (large SQLite -> .zst)
# 4) creates/updates GitHub Releases:
#    - molecules-l1-v1
#    - molecules-psi-l2-v1
#
# Requirements:
# - `gh` CLI authenticated (`gh auth login`)
# - write access to the repo
#
# This script is safe to re-run: it will `upload --clobber` if releases already exist.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

BRANCH="feat/molecules-v1"
BASE_BRANCH="main"

if ! command -v gh >/dev/null 2>&1; then
  echo "ERROR: gh not found. Install GitHub CLI first." >&2
  exit 2
fi

echo "[INFO] gh auth status"
gh auth status >/dev/null

current_branch="$(git branch --show-current)"
if [[ "$current_branch" != "$BRANCH" ]]; then
  echo "[INFO] switching to $BRANCH"
  git checkout "$BRANCH"
fi

echo "[INFO] push branch -> origin/$BRANCH"
git push -u origin "$BRANCH"

if [[ "${SKIP_PR:-}" != "1" ]]; then
  if gh pr view --head "$BRANCH" >/dev/null 2>&1; then
    echo "[INFO] PR already exists for $BRANCH"
  else
    echo "[INFO] creating PR -> $BASE_BRANCH"
    gh pr create --base "$BASE_BRANCH" --head "$BRANCH" \
      --title "Add molecules (L1+PSI L2) pipeline + v1 release metadata" \
      --body "Adds small-molecule L1 (M1+M2) and PSI L2 (M3) industrialized pipeline, QA reports, and release packaging script."
  fi
fi

echo "[INFO] package release assets (dist/molecules)"

# You can override these inputs if your artifacts are in a different location.
# Example:
#   M1_DB=/path/to/molecules_m1.sqlite M3_DB=/path/to/chembl_m3.sqlite bash pipelines/molecules/publish_v1.sh
M1_DB_DEFAULT="data/output/molecules/molecules_m1.sqlite"
M3_DB_DEFAULT="data/output/molecules/chembl_m3.sqlite"
export M1_DB="${M1_DB:-$M1_DB_DEFAULT}"
export M3_DB="${M3_DB:-$M3_DB_DEFAULT}"

bash pipelines/molecules/package_release.sh

target_commit="$(git rev-parse HEAD)"

notes_l1="dist/molecules/RELEASE_NOTES_molecules_l1_v1.md"
cat > "$notes_l1" <<'EOF'
Molecules L1 v1 (M1+M2)

Assets
- `molecules_m1.sqlite.zst`: SQLite DB containing molecule entity (M1) + RDKit descriptors (M2)
- `manifest_molecules_l1_v1.json`: sha256/size manifest
- QC JSON + postmortems for auditability

Source
- Derived from ChEMBL 36 chemreps.
EOF

notes_m3="dist/molecules/RELEASE_NOTES_molecules_psi_l2_v1.md"
cat > "$notes_m3" <<'EOF'
Molecules PSI L2 v1 (M3)

Assets
- `chembl_m3.sqlite.zst`: SQLite DB containing `psi_evidence_v1` and `psi_edges_v1`
- `manifest_molecules_psi_l2_v1.json`: sha256/size manifest
- QC JSON + gates JSON + postmortem for auditability

Source
- Derived from ChEMBL 36 SQLite dump.
EOF

create_or_update_release() {
  local tag="$1"
  local title="$2"
  local notes_file="$3"
  shift 3
  if gh release view "$tag" >/dev/null 2>&1; then
    echo "[INFO] release exists: $tag -> uploading assets (--clobber)"
    gh release upload "$tag" --clobber "$@"
  else
    echo "[INFO] creating release: $tag"
    gh release create "$tag" \
      --title "$title" \
      --notes-file "$notes_file" \
      --target "$target_commit" \
      "$@"
  fi
}

create_or_update_release \
  molecules-l1-v1 \
  "Molecules L1 (M1+M2) v1" \
  "$notes_l1" \
  dist/molecules/molecules_m1.sqlite.zst \
  dist/molecules/manifest_molecules_l1_v1.json \
  pipelines/molecules/reports/molecules_m1_v1.qc.json \
  pipelines/molecules/reports/molecules_m2_v1.qc.json \
  pipelines/molecules/POSTMORTEM_molecules_m1_v1.md \
  pipelines/molecules/POSTMORTEM_molecules_m2_v1.md

create_or_update_release \
  molecules-psi-l2-v1 \
  "Molecules PSI (L2/M3) v1" \
  "$notes_m3" \
  dist/molecules/chembl_m3.sqlite.zst \
  dist/molecules/manifest_molecules_psi_l2_v1.json \
  pipelines/molecules/reports/molecules_m3_smoke_v1.qc.json \
  pipelines/molecules/reports/molecules_m3_full_v1.qc.json \
  pipelines/molecules/reports/molecules_m3_v1.gates.json \
  pipelines/molecules/POSTMORTEM_molecules_m3_v1.md \
  pipelines/molecules/M3_decisions_v1.md

echo "[DONE] Published molecules releases."
echo "- https://github.com/hazelian0619/protian-entity/releases/tag/molecules-l1-v1"
echo "- https://github.com/hazelian0619/protian-entity/releases/tag/molecules-psi-l2-v1"

