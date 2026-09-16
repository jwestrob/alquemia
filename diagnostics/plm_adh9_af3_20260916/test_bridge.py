#!/usr/bin/env python3
"""Synthetic adapter tests; no protein preparation or quantum calculation."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import gemmi

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("af3_bridge_under_test", HERE / "prepare_af3.py")
bridge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bridge)


def summary():
    return {"chain_pair_iptm": [[.8, .93, .95], [.93, None, .51], [.95, .51, .8]], "iptm": .9}


def synthetic_structure(path, residue_name="GLU"):
    structure = gemmi.Structure()
    model = gemmi.Model(1)
    for chain_name, name, atom_name, element in (("A", residue_name, "CA", "C"),
            ("B", "LA", "LA", "La"), ("C", "PQQ", "O5", "O")):
        chain = gemmi.Chain(chain_name)
        residue = gemmi.Residue()
        residue.name, residue.seqid = name, gemmi.SeqId(1, " ")
        atom = gemmi.Atom()
        atom.name, atom.element, atom.pos = atom_name, gemmi.Element(element), gemmi.Position(1, 2, 3)
        residue.add_atom(atom)
        chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    structure.make_mmcif_document().write_file(str(path))


class BridgeTests(unittest.TestCase):
    def test_native_null_diagonal_preserved(self):
        value = summary()
        before = copy.deepcopy(value)
        self.assertEqual(bridge.validate_confidence(value, ["A", "B", "C"]), .93)
        self.assertEqual(value, before)

    def test_missing_protein_metal_confidence_is_not_imputed(self):
        value = summary(); value["chain_pair_iptm"][0][1] = None
        with self.assertRaises(ValueError):
            bridge.validate_confidence(value, ["A", "B", "C"])

    def test_nan_is_rejected(self):
        value = summary(); value["chain_pair_iptm"][0][1] = float("nan")
        with self.assertRaises(ValueError):
            bridge.validate_confidence(value, ["A", "B", "C"])

    def test_confidence_chain_order_mismatch_is_rejected(self):
        with self.assertRaises(ValueError):
            bridge.validate_confidence(summary(), ["A", "C", "B"])

    def test_role_annotation_removal_does_not_alter_selector(self):
        value = {"anchor_glutamate": {"chain": "A", "resnum": 1, "icode": "", "resname": "GLU", "note": "fixture"}}
        self.assertEqual(bridge.clean_roles(value), {"anchor_glutamate": {"chain": "A", "resnum": 1, "icode": "", "resname": "GLU"}})

    def fixture(self, root):
        source = root / "synthetic_model.cif"
        synthetic_structure(source)
        confidence = root / "synthetic_summary_confidences.json"
        confidence.write_text(json.dumps(summary(), indent=4) + "\n")
        roles = {"anchor_glutamate": {"chain": "A", "resnum": 1, "icode": "", "resname": "GLU"}}
        case = {"case_id": "synthetic_only", "target_id": "synthetic_only", "rank": 0,
            "source_cif": bridge.record(source), "summary_confidences": bridge.record(confidence),
            "protein_chain": "A", "metal": {"chain": "B", "resname": "LA", "resnum": 1, "icode": "", "atom": "LA"},
            "pqq": {"chain": "C", "resname": "PQQ", "resnum": 1, "icode": ""}, "roles": roles,
            "confidence_chain_order": ["A", "B", "C"]}
        prior = {"source_cif": bridge.record(source), "roles": copy.deepcopy(roles)}
        selection = root / "synthetic_selection.json"
        selection.write_text("{}\n")
        return case, prior, bridge.record(selection)

    def test_adapter_changes_names_only(self):
        with tempfile.TemporaryDirectory(prefix="plm_af3_adapter_test_") as tmp:
            root = Path(tmp); case, prior, selection = self.fixture(root)
            target, receipt = bridge.adapt_case(case, prior, root, selection)
            self.assertEqual(target["source_cif"]["sha256"], case["source_cif"]["sha256"])
            self.assertEqual(target["summary_confidence"]["sha256"], case["summary_confidences"]["sha256"])
            self.assertEqual(bridge.read(target["summary_confidence"]["path"]), summary())
            self.assertFalse(receipt["confidence_values_changed"])
            self.assertFalse(list(root.rglob("*.inp")))

    def test_changed_protein_sequence_rejected(self):
        with tempfile.TemporaryDirectory(prefix="plm_af3_adapter_test_") as tmp:
            root = Path(tmp); case, prior, selection = self.fixture(root)
            altered = root / "altered.cif"; synthetic_structure(altered, "ASP")
            case["source_cif"] = bridge.record(altered)
            with self.assertRaisesRegex(ValueError, "protein sequence/numbering"):
                bridge.adapt_case(case, prior, root, selection)

    def test_changed_role_rejected(self):
        with tempfile.TemporaryDirectory(prefix="plm_af3_adapter_test_") as tmp:
            root = Path(tmp); case, prior, selection = self.fixture(root)
            case["roles"]["anchor_glutamate"]["resnum"] = 2
            with self.assertRaisesRegex(ValueError, "role selectors differ"):
                bridge.adapt_case(case, prior, root, selection)

    def test_implementation_pins_pass(self):
        bridge.verify_pins(HERE / "implementation_pins.json")


if __name__ == "__main__":
    unittest.main(verbosity=2)
