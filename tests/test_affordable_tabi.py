"""Real archived parser/preparation checks; no manufactured solver success.

These tests execute no electrostatic solver. Native integration is separately
reported when the approved scientific tasks have actually run.
"""
import copy
import importlib.util
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from affordable_common import InvalidArtifact, read_json, record, write_new
import affordable_tabi as tabi

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "workspaces/affordable_challenger_20260915"
BIN = ARCHIVE / "software/apbs-3.4.1/APBS-3.4.1.Linux/bin"
EXAMPLES = BIN.parent / "share/apbs/examples/bem"
NUMERICS = ROOT / "diagnostics/global_electrostatic_20260916/NUMERICS.md"
STATES = ARCHIVE / "solver_completion"


class RealTABIParser(unittest.TestCase):
    def test_loaded_identity_survives_later_source_file_edit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "actual_adapter_copy.py"
            shutil.copyfile(Path(tabi.__file__), path)
            initial = record(path)
            specification = importlib.util.spec_from_file_location("tabi_identity_integrity_check", path)
            module = importlib.util.module_from_spec(specification)
            specification.loader.exec_module(module)
            # Technical source-copy edit, not a solver run or scientific result.
            with path.open("a") as stream:
                stream.write("\n# Explicit later-file-edit fixture for identity integrity.\n")
            self.assertEqual(module.LOADED_IMPLEMENTATION, initial)
            self.assertNotEqual(record(path)["sha256"], module.LOADED_IMPLEMENTATION["sha256"])

    def test_archived_real_bem_outputs(self):
        for name in ("451c_order1.out", "451c_order5.out"):
            with self.subTest(name=name):
                result = tabi.parse_output(EXAMPLES / name)
                self.assertFalse(result["coulomb_added_to_descriptor"])
                self.assertLessEqual(result["gmres_residual"], 1e-4)
                self.assertAlmostEqual(result["reaction_field_kJ_mol"] / 4.184,
                                       result["reaction_field_kcal_mol"], places=12)
                self.assertNotEqual(result["reaction_field_kJ_mol"], result["audit_only_free_energy_kJ_mol"])

    def test_exact_archived_components(self):
        result = tabi.parse_output(EXAMPLES / "451c_order5.out")
        self.assertEqual(result["reaction_field_kJ_mol"], -4.678186467206e3)
        self.assertEqual(result["coulomb_kJ_mol"], -7.777800511228e4)
        self.assertAlmostEqual(result["reaction_field_kJ_mol"] + result["coulomb_kJ_mol"],
                               result["audit_only_free_energy_kJ_mol"], places=5)

    def test_identified_corrupted_real_outputs_rejected(self):
        original = (EXAMPLES / "451c_order5.out").read_text()
        corruptions = {
            "missing_convergence": original.replace("GMRES completed.", "CORRUPTED convergence marker"),
            "duplicate_output": original + original,
            "nonfinite": original.replace("Solvation energy = -4678.186467", "Solvation energy = nan"),
            "accounting": original.replace("Free energy = -82456.191579", "Free energy = -92456.191579"),
            "nonconvergence": original.replace("9.333179e-05 residual.", "9.333179e-02 residual."),
        }
        with tempfile.TemporaryDirectory() as directory:
            for name, text in corruptions.items():
                path = Path(directory) / (name + ".out")
                path.write_text(text)
                with self.subTest(name=name), self.assertRaises(InvalidArtifact):
                    tabi.parse_output(path)

    def test_empty_csv_header_is_explicit_failure(self):
        original = EXAMPLES / "451c_order5.out"
        with tempfile.TemporaryDirectory() as directory:
            # Explicit malformed-input copy of a real archived output. It is
            # truncated to an empty header, never fabricated as successful CSV.
            broken_header = Path(directory) / "corrupted_headers.csv"
            shutil.copyfile(original, broken_header)
            broken_header.write_bytes(b"")
            with self.assertRaisesRegex(InvalidArtifact, "empty/truncated"):
                tabi.parse_output(original, original, broken_header)


class RealTABIStatePreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.workspace = Path(cls.directory.name)
        cls.source_modes = {name: (BIN / name).stat().st_mode for name in tabi.BINARY_HASHES}
        tabi.install_software(BIN, cls.workspace / "software")
        cls.software = cls.workspace / "software/software_manifest.json"

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def state_path(self, partition="qm33", metal="Ca"):
        if partition == "qm33" and metal == "La":
            return ARCHIVE / "solver_recovery/1h4i_qm33_La/primary/state.json"
        return STATES / f"1h4i_{partition}_{metal}/primary/state.json"

    def prepare(self, name, state=None, **kwargs):
        return tabi.prepare(state or self.state_path(), self.workspace / name,
                            self.software, numerics=NUMERICS, **kwargs)

    def test_software_copy_does_not_modify_archive(self):
        software = read_json(self.software)
        for name, expected in tabi.BINARY_HASHES.items():
            self.assertEqual(record(BIN / name)["sha256"], expected)
            self.assertEqual((BIN / name).stat().st_mode, self.source_modes[name])
            self.assertTrue(os.access(software["executables"][name]["path"], os.X_OK))

    def test_real_partition_and_metal_mesh_inputs_identical(self):
        manifests = [self.prepare(f"mesh_{partition}_{metal}", self.state_path(partition, metal))
                     for partition in ("qm33", "qm36") for metal in ("La", "Ca")]
        self.assertEqual(len({m["physical_xyzr"]["sha256"] for m in manifests}), 1)
        self.assertEqual(len({m["physical_boundary_hash"] for m in manifests}), 1)
        self.assertEqual(len({m["charges"]["sha256"] for m in manifests}), 4)
        self.assertTrue(all(m["physical_atom_count"] == 9141 for m in manifests))

    def test_component_charge_changes_cannot_change_mesh(self):
        manifests = [self.prepare(f"component_{name}", component=name) for name in ("total", "core", "environment", "zero")]
        self.assertEqual(len({m["physical_xyzr"]["sha256"] for m in manifests}), 1)
        self.assertEqual(len({m["charges"]["sha256"] for m in manifests}), 4)
        self.assertEqual(manifests[-1]["total_charge_e"], 0.)

    def test_rigid_transforms_move_real_source_and_charge_atoms_together(self):
        state = read_json(self.state_path())
        source = next(a for a in state["physical_atoms"] if a["id"] == "metal")
        charge = next(a for a in state["core_atoms"] if a["id"] == "metal")
        for name in ("translated", "rotated"):
            definition = tabi.transform_definition(name)
            a = tabi.transformed([source], definition)[0]
            b = tabi.transformed([charge], definition)[0]
            self.assertEqual(a["xyz_A"], b["xyz_A"])
            self.assertNotEqual(a["xyz_A"], source["xyz_A"])
            self.prepare(f"transform_{name}", transform=name)

    def test_caps_remain_charge_sites_but_not_production_mesh(self):
        state = read_json(self.state_path())
        caps = [a for a in state["core_atoms"] if a.get("kind") == "cap"]
        self.assertTrue(caps)
        self.assertFalse({a["id"] for a in caps} & {a["id"] for a in state["physical_atoms"]})
        broken = copy.deepcopy(state)
        broken["physical_atoms"].append(caps[0])
        with self.assertRaisesRegex(InvalidArtifact, "synthetic cap"):
            tabi.validate_state(broken)

    def test_isolated_reduction_role_is_explicit(self):
        state = read_json(self.state_path())
        state["physical_atoms"] = copy.deepcopy(state["core_atoms"])
        state["environment_atoms"] = []
        state["expected_environment_charge_e"] = 0.
        with self.assertRaisesRegex(InvalidArtifact, "synthetic cap"):
            tabi.validate_state(state)
        state["cavity_role"] = "isolated_reduction"
        tabi.validate_state(state)
        path = self.workspace / "isolated_state.json"
        write_new(path, state)
        self.prepare("isolated", state=path)

    def test_invalid_physics_and_charge_ownership_rejected(self):
        state = read_json(self.state_path())
        bad = copy.deepcopy(state); bad["settings"]["solute_dielectric"] = 4.
        with self.assertRaisesRegex(InvalidArtifact, "physical setting"):
            tabi.validate_state(bad)
        bad = copy.deepcopy(state); bad["environment_atoms"].append(bad["core_atoms"][0])
        with self.assertRaisesRegex(InvalidArtifact, "force-field charge"):
            tabi.validate_state(bad)
        bad = copy.deepcopy(state); bad["core_atoms"][0]["charge_e"] += .2
        with self.assertRaisesRegex(InvalidArtifact, "charge closure"):
            tabi.validate_state(bad)

    def test_settings_change_invalidates_cache(self):
        a = self.prepare("cache_primary")
        b = self.prepare("cache_refined", level="refined")
        c = self.prepare("cache_tree", level="tree")
        self.assertEqual(len({m["cache_key"] for m in (a, b, c)}), 3)
        self.assertEqual(a["physical_xyzr"]["sha256"], c["physical_xyzr"]["sha256"])
        self.assertNotEqual(a["input"]["sha256"], c["input"]["sha256"])
        tabi._verified_manifest(self.workspace / "cache_primary/tabi_manifest.json")

    def test_refuses_overwrite_and_missing_numerics(self):
        self.prepare("exclusive")
        with self.assertRaisesRegex(InvalidArtifact, "existing preparation"):
            self.prepare("exclusive")
        with self.assertRaisesRegex(InvalidArtifact, "NUMERICS"):
            tabi.prepare(self.state_path(), self.workspace / "missing_numerics", self.software)

    def test_execution_requires_allocation_without_running_solver(self):
        with patch.dict(os.environ, {}, clear=True), self.assertRaisesRegex(InvalidArtifact, "SLURM"):
            tabi.execute("unopened_manifest.json")


if __name__ == "__main__":
    unittest.main()
