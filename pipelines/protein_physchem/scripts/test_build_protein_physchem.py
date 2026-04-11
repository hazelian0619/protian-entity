import unittest

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from pipelines.protein_physchem.scripts.build_protein_physchem import (
    calculate_numeric_parse_rates,
    compute_physchem,
    sanitize_sequence,
)


class TestProteinPhyschemHelpers(unittest.TestCase):
    def test_sanitize_sequence(self):
        self.assertEqual(sanitize_sequence(" acdEFg "), "ACDEFG")

    def test_compute_physchem_matches_biopython(self):
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        expected = ProteinAnalysis(sequence)

        result = compute_physchem(sequence)

        self.assertEqual(result["sequence_len"], len(sequence))
        self.assertAlmostEqual(result["mass_recalc"], expected.molecular_weight(), places=6)
        self.assertAlmostEqual(result["isoelectric_point"], expected.isoelectric_point(), places=6)
        self.assertAlmostEqual(result["gravy"], expected.gravy(), places=6)
        self.assertAlmostEqual(result["aromaticity"], expected.aromaticity(), places=6)
        self.assertAlmostEqual(result["instability_index"], expected.instability_index(), places=6)

    def test_calculate_numeric_parse_rates(self):
        df = pd.DataFrame(
            {
                "sequence_len": ["3", "4"],
                "mass_recalc": ["123.1", "456.2"],
                "isoelectric_point": ["5.2", "7.4"],
                "gravy": ["-0.5", "0.1"],
                "aromaticity": ["0.1", "0.2"],
                "instability_index": ["20.1", "33.2"],
            }
        )

        rates = calculate_numeric_parse_rates(df)
        self.assertTrue(all(rate == 1.0 for rate in rates.values()))


if __name__ == "__main__":
    unittest.main()
