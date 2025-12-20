#!/usr/bin/env bash
set -euo pipefail

# Packages molecule artifacts for GitHub Releases.
#
# This script does NOT upload anything; it only creates `dist/molecules/*` assets.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

M1_DB_DEFAULT="data/output/molecules/molecules_m1.sqlite"
M3_DB_DEFAULT="data/output/molecules/chembl_m3.sqlite"
OUTDIR_DEFAULT="dist/molecules"

M1_DB="${M1_DB:-$M1_DB_DEFAULT}"
M3_DB="${M3_DB:-$M3_DB_DEFAULT}"
OUTDIR="${OUTDIR:-$OUTDIR_DEFAULT}"

mkdir -p "$OUTDIR"

echo "[INFO] Backup + compress M1 (L1+M2)"
sqlite3 "$M1_DB" ".backup '$OUTDIR/molecules_m1.sqlite'"
zstd -19 -T0 --rm "$OUTDIR/molecules_m1.sqlite" -o "$OUTDIR/molecules_m1.sqlite.zst"

echo "[INFO] Backup + compress M3 (PSI L2)"
sqlite3 "$M3_DB" ".backup '$OUTDIR/chembl_m3.sqlite'"
zstd -19 -T0 --rm "$OUTDIR/chembl_m3.sqlite" -o "$OUTDIR/chembl_m3.sqlite.zst"

echo "[INFO] Make manifests"
python3 tools/kg_make_manifest.py --data-version molecules-l1-v1 --out "$OUTDIR/manifest_molecules_l1_v1.json" \
  "$OUTDIR/molecules_m1.sqlite.zst" \
  pipelines/molecules/reports/molecules_m1_v1.qc.json \
  pipelines/molecules/reports/molecules_m2_v1.qc.json \
  pipelines/molecules/POSTMORTEM_molecules_m1_v1.md \
  pipelines/molecules/POSTMORTEM_molecules_m2_v1.md

python3 tools/kg_make_manifest.py --data-version molecules-psi-l2-v1 --out "$OUTDIR/manifest_molecules_psi_l2_v1.json" \
  "$OUTDIR/chembl_m3.sqlite.zst" \
  pipelines/molecules/reports/molecules_m3_smoke_v1.qc.json \
  pipelines/molecules/reports/molecules_m3_full_v1.qc.json \
  pipelines/molecules/reports/molecules_m3_v1.gates.json \
  pipelines/molecules/POSTMORTEM_molecules_m3_v1.md \
  pipelines/molecules/M3_decisions_v1.md

echo "[DONE] Release assets -> $OUTDIR"

