import unittest

from pipelines.molecule_pubchem_direct_xref.scripts.build_molecule_xref_pubchem_enhanced_v1 import (
    apply_pubchem_enhancement,
    merge_confidence,
    resolve_pubchem_candidates,
)


class TestPubChemDirectLogic(unittest.TestCase):
    def test_inchikey_unique_has_priority(self):
        resolved = resolve_pubchem_candidates(
            inchikey_cids=["12345"],
            smiles_cids=["99999"],
            has_smiles=True,
        )
        self.assertEqual(resolved["selected_cids"], ["12345"])
        self.assertEqual(resolved["match_strategy"], "pubchem_direct_inchikey_exact_unique")
        self.assertEqual(resolved["confidence"], "high")

    def test_inchikey_multi_can_tiebreak_with_smiles(self):
        resolved = resolve_pubchem_candidates(
            inchikey_cids=["10", "20", "30"],
            smiles_cids=["20"],
            has_smiles=True,
        )
        self.assertEqual(resolved["selected_cids"], ["20"])
        self.assertEqual(resolved["match_strategy"], "pubchem_direct_inchikey_multi_tiebreak_smiles")
        self.assertEqual(resolved["confidence"], "high")
        self.assertEqual(resolved["tie_break"], "resolved_unique_by_smiles_intersection")

    def test_smiles_fallback_multi_is_retained(self):
        resolved = resolve_pubchem_candidates(
            inchikey_cids=[],
            smiles_cids=["200", "100"],
            has_smiles=True,
        )
        self.assertEqual(resolved["selected_cids"], ["100", "200"])
        self.assertEqual(resolved["match_strategy"], "pubchem_direct_smiles_fallback_multi_cid")
        self.assertEqual(resolved["confidence"], "low")

    def test_merge_confidence_keeps_higher(self):
        self.assertEqual(merge_confidence("high", "low"), "high")
        self.assertEqual(merge_confidence("medium", "high"), "high")
        self.assertEqual(merge_confidence("", "medium"), "medium")

    def test_existing_high_confidence_mapping_is_not_overwritten(self):
        original = {
            "inchikey": "AAAAABBBBBCCCC-DDDDDEEEEE-F",
            "pubchem_cid": "777",
            "match_strategy": "seed_v2",
            "confidence": "high",
            "xref_source": "chembl_36.db",
            "source_version": "legacy",
        }
        resolved = {
            "selected_cids": ["123"],
            "match_strategy": "pubchem_direct_inchikey_exact_unique",
            "confidence": "high",
            "tie_break": "none",
            "inchikey_candidates": ["123"],
            "smiles_candidates": [],
        }
        updated, changed, reason = apply_pubchem_enhancement(
            row=original,
            resolved=resolved,
            source_tag="pubchem_pug_rest",
            source_version_tag="PubChem:PUGREST:2026-04-13",
        )
        self.assertFalse(changed)
        self.assertEqual(reason, "locked_existing_high_confidence")
        self.assertEqual(updated["pubchem_cid"], "777")


if __name__ == "__main__":
    unittest.main()

