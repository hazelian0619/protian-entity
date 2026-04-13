import unittest

from pipelines.molecule_physchem_descriptors.scripts.build_molecule_physchem_descriptors_v1 import (
    choose_best_record,
    merge_descriptors,
    parse_pubchem_property_entry,
)


class TestMoleculePhyschemDescriptorLogic(unittest.TestCase):
    def test_choose_best_record_prefers_more_non_empty_fields(self):
        records = [
            {"molecular_weight": "", "logp": "", "tpsa": "12.2", "hbd": "", "hba": "", "rotatable_bonds": ""},
            {"molecular_weight": "180.16", "logp": "1.2", "tpsa": "63.6", "hbd": "1", "hba": "4", "rotatable_bonds": "3"},
        ]
        best = choose_best_record(records)
        self.assertEqual(best["molecular_weight"], "180.16")
        self.assertEqual(best["rotatable_bonds"], "3")

    def test_merge_descriptors_chembl_has_priority_with_pubchem_fill(self):
        chembl = {
            "molecular_weight": "300.12",
            "logp": "",
            "tpsa": "45.1",
            "hbd": "1",
            "hba": "4",
            "rotatable_bonds": "2",
        }
        pubchem = {
            "molecular_weight": "301.00",
            "logp": "2.5",
            "tpsa": "46.0",
            "hbd": "2",
            "hba": "5",
            "rotatable_bonds": "3",
        }
        merged = merge_descriptors(chembl, pubchem)
        self.assertEqual(merged["molecular_weight"], "300.12")
        self.assertEqual(merged["logp"], "2.5")
        self.assertEqual(merged["tpsa"], "45.1")

    def test_parse_pubchem_property_entry(self):
        entry = {
            "CID": 2244,
            "MolecularWeight": "180.16",
            "XLogP": 1.2,
            "TPSA": 63.6,
            "HBondDonorCount": 1,
            "HBondAcceptorCount": 4,
            "RotatableBondCount": 3,
        }
        out = parse_pubchem_property_entry(entry)
        self.assertEqual(out["molecular_weight"], "180.16")
        self.assertEqual(out["logp"], "1.2")
        self.assertEqual(out["hbd"], "1")
        self.assertEqual(out["hba"], "4")
        self.assertEqual(out["rotatable_bonds"], "3")


if __name__ == "__main__":
    unittest.main()
