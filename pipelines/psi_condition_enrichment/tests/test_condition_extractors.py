#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipelines.psi_condition_enrichment.scripts.condition_extractors import extract_condition_bundle


class ConditionExtractorTests(unittest.TestCase):
    def test_extract_explicit_ph_and_temperature(self) -> None:
        row = {
            "assay_description": "Assay in 20 mM HEPES buffer (pH 7.4), incubated at 37 C in HEK293 cells",
            "assay_context": "",
            "activity_comment": "",
            "data_validity_comment": "",
            "condition_pH": "",
            "condition_temperature_c": "",
            "condition_system": "",
            "condition_context": "",
        }
        out = extract_condition_bundle(row)
        self.assertEqual(out["condition_pH"], "7.4")
        self.assertEqual(out["condition_temperature_c"], "37")
        self.assertIn("HEPES", out["condition_context_json"])
        self.assertIn("cell_line", out["condition_context"])

    def test_extract_range_and_fahrenheit(self) -> None:
        row = {
            "assay_description": "Binding measured at pH 7-8 and 68 F in Tris buffer",
            "assay_context": "",
            "activity_comment": "",
            "data_validity_comment": "",
            "condition_pH": "",
            "condition_temperature_c": "",
            "condition_system": "",
            "condition_context": "",
        }
        out = extract_condition_bundle(row)
        self.assertEqual(out["condition_pH"], "7-8")
        self.assertEqual(out["condition_temperature_c"], "20")

    def test_infer_from_cell_culture_clues(self) -> None:
        row = {
            "assay_description": "Activity measured in CHO cells in DMEM medium after 16 h incubation",
            "assay_context": "",
            "activity_comment": "",
            "data_validity_comment": "",
            "condition_pH": "",
            "condition_temperature_c": "",
            "condition_system": "",
            "condition_context": "",
        }
        out = extract_condition_bundle(row)
        self.assertEqual(out["condition_pH"], "7.4")
        self.assertEqual(out["condition_temperature_c"], "37")
        self.assertGreaterEqual(float(out["condition_extract_confidence"]), 0.45)

    def test_conflict_detection(self) -> None:
        row = {
            "assay_description": "enzyme assay at pH 7.4",
            "assay_context": "",
            "activity_comment": "assayed again at pH 8.0",
            "data_validity_comment": "",
            "condition_pH": "",
            "condition_temperature_c": "",
            "condition_system": "",
            "condition_context": "",
        }
        out = extract_condition_bundle(row)
        self.assertEqual(out["conflict_flag"], "true")
        self.assertIn("condition_pH", out["conflict_fields"])


if __name__ == "__main__":
    unittest.main()
