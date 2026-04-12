import unittest

from pipelines.rna_type_features.scripts.build_rna_type_features import (
    choose_best_urs,
    infer_rrna_class,
    parse_external_locus,
    parse_gtf_attributes,
    parse_trna_annotation,
)


class TestRnaTypeFeaturesHelpers(unittest.TestCase):
    def test_parse_gtf_attributes(self):
        attrs = (
            'gene_id "ENSG00000226803.10"; transcript_id "ENST00000585414.1"; '
            'gene_name "ZNF451-AS1"; transcript_name "ZNF451-AS1-201"; transcript_biotype "lncRNA";'
        )
        parsed = parse_gtf_attributes(attrs)
        self.assertEqual(parsed["gene_id"], "ENSG00000226803.10")
        self.assertEqual(parsed["transcript_id"], "ENST00000585414.1")
        self.assertEqual(parsed["transcript_biotype"], "lncRNA")

    def test_parse_external_locus_ena(self):
        loc = parse_external_locus("CM034956.1:56983721..56986952:ncRNA")
        self.assertIsNotNone(loc)
        self.assertEqual(loc["seq_accession"], "CM034956.1")
        self.assertEqual(loc["start"], "56983721")
        self.assertEqual(loc["end"], "56986952")
        self.assertEqual(loc["strand"], ".")
        self.assertEqual(loc["feature"], "ncRNA")

    def test_parse_external_locus_gtrnadb(self):
        loc = parse_external_locus("GTRNADB:tRNA-Thr-TGT-2-1:CM000663.1:222638347-222638419")
        self.assertIsNotNone(loc)
        self.assertEqual(loc["seq_accession"], "CM000663.1")
        self.assertEqual(loc["start"], "222638347")
        self.assertEqual(loc["end"], "222638419")
        self.assertEqual(loc["feature"], "tRNA")

    def test_parse_trna_annotation(self):
        aa, anti, anti_rna, label = parse_trna_annotation(
            "GTRNADB:tRNA-Arg-CCT-3-1:CM000678.2:3152900-3152972",
            "tRNA-Arg-CCT-3-1",
        )
        self.assertEqual(aa, "Arg")
        self.assertEqual(anti, "CCT")
        self.assertEqual(anti_rna, "CCU")
        self.assertEqual(label, "tRNA-Arg-CCT-3-1")

    def test_infer_rrna_class(self):
        self.assertEqual(infer_rrna_class("RNA5-8SN5", "NR_003285", "rRNA"), "5.8S_rRNA")
        self.assertEqual(infer_rrna_class("RNA5S1", "NR_023363", "rRNA"), "5S_rRNA")
        self.assertEqual(infer_rrna_class("", "18S ribosomal RNA", "rRNA"), "18S_rRNA")
        self.assertEqual(infer_rrna_class("", "unclassified", "rRNA"), "rRNA")

    def test_choose_best_urs(self):
        choice = choose_best_urs({"URS1_9606": 2, "URS2_9606": 5, "URS3_9606": 1})
        self.assertEqual(choice, "URS2_9606")


if __name__ == "__main__":
    unittest.main()
