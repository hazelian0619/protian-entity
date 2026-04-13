#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

FETCH_DATE = date.today().isoformat()

OUTPUT_COLUMNS = [
    "mapping_id",
    "interaction_type",
    "edge_id",
    "mapping_scope",
    "raw_term",
    "standardized_predicate",
    "ontology_namespace",
    "ontology_id",
    "ontology_label",
    "ontology_uri",
    "mapping_confidence",
    "source",
    "source_version",
    "fetch_date",
]

EDGE_INPUTS = {
    "PPI": Path("data/output/edges/edges_ppi_v1.tsv"),
    "PSI": Path("data/output/edges/drug_target_edges_v1.tsv"),
    "RPI": Path("data/output/edges/rna_protein_edges_v1.tsv"),
}

RPI_FUNCTION_INPUT = Path("data/output/evidence/rpi_function_context_v2.tsv")

GO_RELATION_MAP = {
    "transcription_regulation": (
        "GO:0006355",
        "regulation of DNA-templated transcription",
    ),
    "translation_regulation": (
        "GO:0006417",
        "regulation of translation",
    ),
    "splicing_regulation": (
        "GO:0043484",
        "regulation of RNA splicing",
    ),
    "post_transcriptional_regulation": (
        "GO:0010608",
        "post-transcriptional regulation of gene expression",
    ),
}

SOURCE_URI_MAP = {
    "string": ("SOURCE:STRING", "STRING", "https://identifiers.org/string"),
    "drugbank": ("SOURCE:DRUGBANK", "DrugBank", "https://identifiers.org/drugbank"),
    "starbase": ("SOURCE:STARBASE", "starBase/ENCORI", "https://rnasysu.com/encori"),
    "encori": ("SOURCE:ENCORI", "ENCORI", "https://rnasysu.com/encori"),
    "rnainter": ("SOURCE:RNAINTER", "RNAInter", "http://www.rnainter.org"),
    "npinter": ("SOURCE:NPINTER", "NPInter", "http://bigdata.ibp.ac.cn/npinter4"),
}

URI_RE = re.compile(r"^https?://")


def normalize(v: str) -> str:
    return (v or "").strip()


def lower(v: str) -> str:
    return normalize(v).lower()


def to_float_text(v: float) -> str:
    return f"{float(v):.6f}"


def write_json(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def make_id(prefix: str, payload: str) -> str:
    return prefix + "_" + hashlib.sha1(payload.encode("utf-8")).hexdigest()[:20]


def slugify(s: str) -> str:
    t = re.sub(r"[^a-zA-Z0-9]+", "_", s).strip("_")
    return t.lower() if t else "unknown"


def map_method_to_psi_mi(method_raw: str) -> Tuple[str, str, str, float]:
    m = lower(method_raw)
    # conservative PSI-MI mapping for current method vocabulary
    if "string" in m or "predict" in m:
        return ("MI:0063", "interaction prediction", "http://purl.obolibrary.org/obo/MI_0063", 0.90)
    if "clip" in m or "rip" in m:
        return ("MI:0686", "unspecified method", "http://purl.obolibrary.org/obo/MI_0686", 0.80)
    if any(x in m for x in ["inhib", "agon", "antagon", "induc", "target", "binder", "substrate", "activ"]):
        return ("MI:0915", "physical association", "http://purl.obolibrary.org/obo/MI_0915", 0.70)
    if m:
        return ("MI:0686", "unspecified method", "http://purl.obolibrary.org/obo/MI_0686", 0.65)
    return ("MI:0686", "unspecified method", "http://purl.obolibrary.org/obo/MI_0686", 0.50)


def infer_method_from_edge_source(interaction_type: str, source: str) -> str:
    s = lower(source)
    if interaction_type == "PPI":
        if "string" in s:
            return "string_combined_score"
        return "predicted_interaction"
    if interaction_type == "PSI":
        if "drugbank" in s:
            return "drug_target_curation"
        return "curated_interaction"
    if interaction_type == "RPI":
        if any(x in s for x in ["starbase", "encori"]):
            return "clip_supported_record"
        if "rnainter" in s:
            return "rnainter_record"
        if "npinter" in s:
            return "npinter_record"
        return "rna_protein_record"
    return "interaction_record"


def map_source_to_uri(source: str) -> Tuple[str, str, str, float]:
    s = lower(source)
    for k, val in SOURCE_URI_MAP.items():
        if k in s:
            oid, label, uri = val
            return oid, label, uri, 0.95
    token = slugify(s or "unknown_source")
    return (
        f"SOURCE:{token.upper()}",
        source or "Unknown source",
        f"https://example.org/source/{token}",
        0.60,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build ontology mapping table for interaction evidence/edges")
    p.add_argument(
        "--output",
        type=Path,
        default=Path("data/output/evidence/interaction_ontology_mapping_v2.tsv"),
    )
    p.add_argument(
        "--report",
        type=Path,
        default=Path("pipelines/interaction_ontology_mapping/reports/interaction_ontology_mapping_v2.metrics.json"),
    )
    p.add_argument("--gates-report", type=Path, default=None)
    p.add_argument("--limit-per-type", type=int, default=None)
    p.add_argument("--fetch-date", default=FETCH_DATE)
    p.add_argument("--check-inputs", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    base = {
        "pipeline": "interaction_ontology_mapping_v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "sample_mode": args.limit_per_type is not None,
        "input": {
            "edges": {k: str(v) for k, v in EDGE_INPUTS.items()},
            "rpi_function_context": str(RPI_FUNCTION_INPUT),
        },
    }

    missing = [str(p) for p in EDGE_INPUTS.values() if not p.exists()]

    if args.check_inputs:
        if missing:
            blocked = {
                **base,
                "status": "blocked_missing_inputs",
                "missing_required": missing,
                "message": "缺少必需输入，已按协作约束中断。",
            }
            write_json(args.report, blocked)
            print(f"[BLOCKED] missing inputs: {missing}")
            print(f"[BLOCKED] report -> {args.report}")
            return 2

        ready = {
            **base,
            "status": "inputs_ready",
            "optional_inputs": {
                "rpi_function_context_exists": RPI_FUNCTION_INPUT.exists(),
            },
            "message": "输入检查通过。",
        }
        write_json(args.report, ready)
        print(f"[OK] inputs ready -> {args.report}")
        return 0

    if missing:
        blocked = {
            **base,
            "status": "blocked_missing_inputs",
            "missing_required": missing,
            "message": "缺少必需输入，已按协作约束中断。",
        }
        write_json(args.report, blocked)
        print(f"[BLOCKED] missing inputs: {missing}")
        return 2

    args.output.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    mapped_rows = 0
    uri_valid_rows = 0
    by_scope = Counter()
    by_predicate = Counter()
    by_interaction_type = Counter()

    source_rep_edge: Dict[Tuple[str, str], str] = {}
    selected_edge_ids: Set[str] = set()

    with args.output.open("w", encoding="utf-8", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()

        # method + predicate mapping per edge
        for itype, edge_path in EDGE_INPUTS.items():
            count = 0
            with edge_path.open("r", encoding="utf-8", newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                if reader.fieldnames is None:
                    raise SystemExit(f"[ERROR] missing header: {edge_path}")
                required = {"edge_id", "predicate", "source", "source_version"}
                if not required.issubset(set(reader.fieldnames)):
                    raise SystemExit(f"[ERROR] edge table missing required columns: {edge_path}")

                for row in reader:
                    if args.limit_per_type is not None and count >= args.limit_per_type:
                        break
                    edge_id = normalize(row.get("edge_id", ""))
                    if not edge_id:
                        continue
                    count += 1

                    selected_edge_ids.add(edge_id)
                    by_interaction_type[itype] += 1

                    raw_predicate = normalize(row.get("predicate", "")) or "unknown_predicate"
                    source = normalize(row.get("source", "")) or "Unknown"
                    source_version = normalize(row.get("source_version", "")) or "unknown"

                    # 1) PSI-MI method mapping row
                    raw_method = infer_method_from_edge_source(itype, source)
                    mi_id, mi_label, mi_uri, mi_conf = map_method_to_psi_mi(raw_method)
                    rec_method = {
                        "mapping_id": make_id("iom", f"{itype}|{edge_id}|method|{raw_method}"),
                        "interaction_type": itype,
                        "edge_id": edge_id,
                        "mapping_scope": "psi_mi_method",
                        "raw_term": raw_method,
                        "standardized_predicate": "detected_by",
                        "ontology_namespace": "PSI-MI",
                        "ontology_id": mi_id,
                        "ontology_label": mi_label,
                        "ontology_uri": mi_uri,
                        "mapping_confidence": to_float_text(mi_conf),
                        "source": source,
                        "source_version": source_version,
                        "fetch_date": args.fetch_date,
                    }
                    writer.writerow(rec_method)

                    total_rows += 1
                    by_scope[rec_method["mapping_scope"]] += 1
                    by_predicate[rec_method["standardized_predicate"]] += 1
                    if rec_method["ontology_id"] and rec_method["ontology_uri"]:
                        mapped_rows += 1
                    if URI_RE.match(rec_method["ontology_uri"]):
                        uri_valid_rows += 1

                    # 2) unified predicate standardization row
                    rec_pred = {
                        "mapping_id": make_id("iom", f"{itype}|{edge_id}|predicate|{raw_predicate}"),
                        "interaction_type": itype,
                        "edge_id": edge_id,
                        "mapping_scope": "predicate_standardization",
                        "raw_term": raw_predicate,
                        "standardized_predicate": "standardized_as",
                        "ontology_namespace": "KGPRED",
                        "ontology_id": "KGPRED:participates_in",
                        "ontology_label": "participates_in",
                        "ontology_uri": "https://w3id.org/protian/predicate/participates_in",
                        "mapping_confidence": to_float_text(0.95),
                        "source": source,
                        "source_version": source_version,
                        "fetch_date": args.fetch_date,
                    }
                    writer.writerow(rec_pred)

                    total_rows += 1
                    by_scope[rec_pred["mapping_scope"]] += 1
                    by_predicate[rec_pred["standardized_predicate"]] += 1
                    if rec_pred["ontology_id"] and rec_pred["ontology_uri"]:
                        mapped_rows += 1
                    if URI_RE.match(rec_pred["ontology_uri"]):
                        uri_valid_rows += 1

                    source_key = (itype, source)
                    if source_key not in source_rep_edge:
                        source_rep_edge[source_key] = edge_id

        # 3) source support rows (include supported_by predicate vocabulary)
        for (itype, source), edge_id in sorted(source_rep_edge.items()):
            oid, label, uri, conf = map_source_to_uri(source)
            rec_source = {
                "mapping_id": make_id("iom", f"{itype}|{edge_id}|source|{source}"),
                "interaction_type": itype,
                "edge_id": edge_id,
                "mapping_scope": "source_support",
                "raw_term": source,
                "standardized_predicate": "supported_by",
                "ontology_namespace": "SOURCE",
                "ontology_id": oid,
                "ontology_label": label,
                "ontology_uri": uri,
                "mapping_confidence": to_float_text(conf),
                "source": source,
                "source_version": "derived_from_edge_source",
                "fetch_date": args.fetch_date,
            }
            writer.writerow(rec_source)

            total_rows += 1
            by_scope[rec_source["mapping_scope"]] += 1
            by_predicate[rec_source["standardized_predicate"]] += 1
            if rec_source["ontology_id"] and rec_source["ontology_uri"]:
                mapped_rows += 1
            if URI_RE.match(rec_source["ontology_uri"]):
                uri_valid_rows += 1

        # 4) GO term mapping from RPI function contexts
        go_mapped_rows = 0
        go_input_rows = 0
        if RPI_FUNCTION_INPUT.exists():
            with RPI_FUNCTION_INPUT.open("r", encoding="utf-8", newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                if reader.fieldnames is not None and {"edge_id", "function_relation"}.issubset(set(reader.fieldnames)):
                    for row in reader:
                        edge_id = normalize(row.get("edge_id", ""))
                        if not edge_id:
                            continue
                        if args.limit_per_type is not None and edge_id not in selected_edge_ids:
                            continue

                        go_input_rows += 1
                        raw_rel = normalize(row.get("function_relation", ""))
                        go_id, go_label = GO_RELATION_MAP.get(
                            raw_rel,
                            ("GO:0010608", "post-transcriptional regulation of gene expression"),
                        )
                        go_uri = "http://purl.obolibrary.org/obo/" + go_id.replace(":", "_")
                        rec_go = {
                            "mapping_id": make_id("iom", f"RPI|{edge_id}|go|{raw_rel}|{go_id}"),
                            "interaction_type": "RPI",
                            "edge_id": edge_id,
                            "mapping_scope": "go_term",
                            "raw_term": raw_rel or "post_transcriptional_regulation",
                            "standardized_predicate": "participates_in",
                            "ontology_namespace": "GO",
                            "ontology_id": go_id,
                            "ontology_label": go_label,
                            "ontology_uri": go_uri,
                            "mapping_confidence": to_float_text(0.90),
                            "source": normalize(row.get("source", "")) or "UniProtKB",
                            "source_version": normalize(row.get("source_version", "")) or "unknown",
                            "fetch_date": args.fetch_date,
                        }
                        writer.writerow(rec_go)

                        total_rows += 1
                        go_mapped_rows += 1
                        by_scope[rec_go["mapping_scope"]] += 1
                        by_predicate[rec_go["standardized_predicate"]] += 1
                        if rec_go["ontology_id"] and rec_go["ontology_uri"]:
                            mapped_rows += 1
                        if URI_RE.match(rec_go["ontology_uri"]):
                            uri_valid_rows += 1

    mappable_rate = (mapped_rows / total_rows) if total_rows else 0.0
    uri_valid_rate = (uri_valid_rows / total_rows) if total_rows else 0.0

    checks = {
        "mappable_record_ratio_ge_0_85": mappable_rate >= 0.85,
        "ontology_uri_valid_rate_ge_0_99": uri_valid_rate >= 0.99,
        "three_interaction_types_present": all(by_interaction_type.get(x, 0) > 0 for x in ["PPI", "PSI", "RPI"]),
    }

    metrics = {
        **base,
        "status": "completed",
        "output": str(args.output),
        "row_count": {
            "total_mapping_rows": total_rows,
            "edge_rows_processed_by_type": dict(by_interaction_type),
            "scope_rows": dict(by_scope),
            "go_mapped_rows": go_mapped_rows,
        },
        "coverage": {
            "mappable_record_ratio": mappable_rate,
            "ontology_uri_valid_rate": uri_valid_rate,
        },
        "distribution": {
            "standardized_predicate": dict(by_predicate),
        },
        "gates": {
            "passed": all(checks.values()),
            "checks": checks,
        },
    }

    write_json(args.report, metrics)

    if args.gates_report is not None:
        write_json(
            args.gates_report,
            {
                "pipeline": "interaction_ontology_mapping_v2",
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "status": "PASS" if metrics["gates"]["passed"] else "FAIL",
                "checks": metrics["gates"]["checks"],
                "coverage": metrics["coverage"],
                "row_count": metrics["row_count"],
            },
        )

    print(
        "[OK] "
        f"rows={total_rows} mappable_rate={mappable_rate:.4f} uri_valid_rate={uri_valid_rate:.4f} "
        f"types={dict(by_interaction_type)}"
    )
    print(f"[OK] report -> {args.report}")
    if args.gates_report is not None:
        print(f"[OK] gates -> {args.gates_report}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
