import unittest

from pipelines.rna_type_features.scripts.build_rna_trna_features_v2 import (
    TRNA_MT_SYMBOL_TO_ANTICODON,
    choose_balanced,
    extract_mt_symbol,
    parse_trna_annotation,
)


class TestBuildRnaTrnaFeaturesV2(unittest.TestCase):
    def test_parse_trna_annotation(self):
        aa, anti, anti_rna, label = parse_trna_annotation(
            "GTRNADB:tRNA-Leu-TAG-3-1:CM000678.2:22195711-22195792",
            "tRNA-Leu-TAG-3-1",
        )
        self.assertEqual(aa, "Leu")
        self.assertEqual(anti, "TAG")
        self.assertEqual(anti_rna, "UAG")
        self.assertTrue(label.startswith("tRNA-Leu-TAG"))

    def test_extract_mt_symbol(self):
        s = extract_mt_symbol("GENECARDS:MT-TP:URS000002176F_9606", "")
        self.assertEqual(s, "MT-TP")
        self.assertEqual(TRNA_MT_SYMBOL_TO_ANTICODON[s], "TGG")

    def test_choose_balanced_resolve_conflict(self):
        anti, source, conf, status, ranked = choose_balanced(
            {
                "TGG": {"score": 3.0, "lines": 3, "dbs": {"GTRNADB"}, "methods": {"current_parser"}, "aa": {}},
                "TGC": {"score": 0.7, "lines": 2, "dbs": {"ENA"}, "methods": {"mt_symbol_map"}, "aa": {}},
            }
        )
        self.assertEqual(anti, "TGG")
        self.assertEqual(status, "weighted_conflict_resolve")
        self.assertEqual(source, "gtrnadb")
        self.assertEqual(conf, "high")
        self.assertGreaterEqual(len(ranked), 2)

    def test_choose_balanced_unresolved(self):
        anti, source, conf, status, ranked = choose_balanced(
            {
                "TGC": {"score": 0.7, "lines": 2, "dbs": {"ENA"}, "methods": {"mt_symbol_map"}, "aa": {}},
                "TGA": {"score": 0.7, "lines": 2, "dbs": {"ENA"}, "methods": {"mt_symbol_map"}, "aa": {}},
            }
        )
        self.assertEqual(anti, "")
        self.assertEqual(source, "none")
        self.assertEqual(conf, "")
        self.assertEqual(status, "conflict_unresolved")
        self.assertEqual(len(ranked), 2)


if __name__ == "__main__":
    unittest.main()
