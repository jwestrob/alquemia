"""Parser/contract tests on real pinned parameters and archived molecular outputs.

Malformed cases are explicitly corrupted copies. No fabricated successful ORCA
output or scientific energy is used. These tests never execute a chemistry tool.
"""
import copy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("native_gfn2_parameters", ROOT / "scripts/native_gfn2_parameters.py")
g = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(g)
PARAM = g.DATA / "corrected_la.xtb.json"
POLICY = g.DATA / "la_shell_three_arrays_v1.json"
ORIGINAL_POLICY = g.DATA / "original_diagnostic_v1.json"
ARCHIVE = Path("/home/jwestrob/jwestrob/obsidian-vault/projects/lanthanide-binding/orca_parameter_audit_2026-10-12/la_followup/native_reciprocal_controls/results")


class ParameterTests(unittest.TestCase):
    def test_real_la_nine_changes_and_shell_identity(self):
        audit = g.audit_parameters(PARAM, POLICY)
        self.assertEqual(len(audit["differences"]), 9)
        self.assertEqual(audit["elements"]["La"]["shells"], ["5d", "6s", "6p"])
        for element in ("La", "Ca", "H", "C", "N", "O", "S"):
            self.assertTrue(audit["elements"][element]["parameter_consistent"])
            self.assertFalse(audit["elements"][element]["runtime_qualified"])

    def test_no_silent_qualification_of_other_metals(self):
        audit = g.audit_parameters(PARAM, POLICY)
        for element in ("Dy", "Cu", "Mn", "Fe", "Co", "Ni"):
            self.assertFalse(audit["elements"][element]["parameter_consistent"])
        self.assertFalse(audit["elements"]["Cu"]["fields"]["refocc"]["passed"])

    def test_original_detects_la_mismatch(self):
        audit = g.audit_parameters(g.ORIGINAL, ORIGINAL_POLICY)
        self.assertEqual(audit["differences"], [])
        self.assertFalse(audit["elements"]["La"]["parameter_consistent"])
        self.assertTrue(audit["elements"]["La"]["fields"]["refocc"]["passed"])

    def corrupt(self, mutation, policy_mutation=None):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / "explicitly_corrupted_real_parameters.json"
        data = g.load_json(PARAM)
        mutation(data)
        path.write_text(json.dumps(data))
        policy = g.load_json(POLICY)
        policy["candidate_sha256"] = g.sha256(path)
        if policy_mutation:
            policy_mutation(policy)
        return path, policy

    def test_unapproved_global_change_rejected(self):
        path, policy = self.corrupt(lambda d: d["dispersion"]["d4"].update(a1=0.521))
        with self.assertRaisesRegex(g.GuardError, "unapproved"):
            g.audit_parameters(path, policy)

    def test_tiny_unrelated_change_is_not_roundoff_exception(self):
        path, policy = self.corrupt(lambda d: d["element"]["Ca"]["kcn"].__setitem__(0, 0.057093 + 1e-14))
        with self.assertRaisesRegex(g.GuardError, "unapproved"):
            g.audit_parameters(path, policy)

    def test_wrong_la_permutation_rejected(self):
        path, policy = self.corrupt(lambda d: d["element"]["La"]["kcn"].reverse())
        with self.assertRaisesRegex(g.GuardError, "does not match author"):
            g.audit_parameters(path, policy)

    def test_global_shell_reorder_rejected(self):
        path, policy = self.corrupt(lambda d: d["element"]["Ca"]["shells"].reverse())
        with self.assertRaisesRegex(g.GuardError, "unapproved"):
            g.audit_parameters(path, policy)

    def test_missing_full_json_section_rejected(self):
        path, policy = self.corrupt(lambda d: d.pop("dispersion"))
        with self.assertRaisesRegex(g.GuardError, "keys changed"):
            g.audit_parameters(path, policy)

    def test_unsupported_correction_field_rejected(self):
        policy = g.load_json(POLICY)
        policy["corrections"]["La"].append("shells")
        with self.assertRaisesRegex(g.GuardError, "unsupported correction"):
            g.audit_parameters(PARAM, policy)

    def test_wrong_full_file_digest_rejected(self):
        policy = g.load_json(POLICY)
        policy["candidate_sha256"] = g.ORIGINAL_SHA
        with self.assertRaisesRegex(g.GuardError, "full-file pin"):
            g.audit_parameters(PARAM, policy)

    def test_duplicate_and_nonfinite_json_rejected(self):
        for text in ('{"a":1,"a":2}', '{"a":NaN}'):
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "malformed_input.json"
                path.write_text(text)
                with self.assertRaises(g.GuardError):
                    g.load_json(path)


@unittest.skipUnless((ARCHIVE / "La_water_corrected/receipt.json").exists(),
                     "real archived ORCA controls unavailable; no fabricated replacement")
class TaskContractTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.task = Path(self.tmp.name)
        self.source = ARCHIVE / "La_water_corrected"
        self.receipt = g.load_json(self.source / "receipt.json")
        for name in ("core.xyz", "endpoint.inp"):
            shutil.copyfile(self.source / name, self.task / name)
        self.executable = Path(self.receipt["binary"]["path"])
        if not self.executable.exists():
            self.skipTest("actual recorded ORCA executable unavailable")
        # Test binding records actual files/environment; it does not execute ORCA.
        self.binding = {
            "schema_version": 1, "protocol_id": "parameter-contract-real-fixture-test",
            "parameter_filename": "control.xtb.json", "parameter_sha256": g.sha256(PARAM),
            "input": g.pin(self.task / "endpoint.inp"), "xyz": g.pin(self.task / "core.xyz"),
            "executable": g.pin(self.executable), "runtime_files": [],
            "source_inputs": [g.pin(self.source / "core.xyz")],
            "environment": {"OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS")},
            "model": {"method": g.METHOD, "medium": "ALPB(Water)"},
            "state_convention": g.STATE_CONVENTION,
            "state": {"charge": -5, "multiplicity": 1},
            "elements": g.xyz_elements(self.task / "core.xyz"),
        }

    def stage(self, **kwargs):
        return g.stage(self.binding, PARAM, POLICY, mode="diagnostic", **kwargs)

    def change_input(self, change):
        path = self.task / "endpoint.inp"
        path.write_text(change(path.read_text()))
        self.binding["input"] = g.pin(path)

    def test_real_staging_preflight_and_archived_postflight(self):
        staged = self.stage()
        self.assertFalse(staged["runtime_qualified"])
        manifest = staged["manifest"]["path"]
        self.assertTrue(g.preflight(manifest)["parameter_preflight_passed"])
        result = g.postflight(manifest, self.source / "endpoint.out", 0)
        self.assertTrue(result["parameter_postflight_passed"])
        self.assertIsNone(result["scientific_result_accepted"])
        self.assertEqual((self.task / "control.xtb.json").read_bytes(), PARAM.read_bytes())

    def test_research_without_runtime_approval_rejected(self):
        with self.assertRaisesRegex(g.GuardError, "runtime approval required"):
            g.stage(self.binding, PARAM, POLICY)

    def test_original_cannot_enter_research(self):
        with self.assertRaisesRegex(g.GuardError, "diagnostic mode"):
            g.stage(self.binding, g.ORIGINAL, ORIGINAL_POLICY, variant="original")

    def test_do_not_overwrite_existing_stage(self):
        self.stage()
        with self.assertRaisesRegex(g.GuardError, "overwrite"):
            self.stage()

    def test_parameter_changes_after_staging_rejected(self):
        staged = self.stage()
        with (self.task / "control.xtb.json").open("a") as stream:
            stream.write("\n")
        with self.assertRaisesRegex(g.GuardError, "digest changed"):
            g.preflight(staged["manifest"]["path"])

    def test_geometry_changes_after_staging_rejected(self):
        staged = self.stage()
        with (self.task / "core.xyz").open("a") as stream:
            stream.write("\n")
        with self.assertRaisesRegex(g.GuardError, "digest changed"):
            g.preflight(staged["manifest"]["path"])

    def test_environment_changes_rejected(self):
        staged = self.stage()
        with patch.dict(os.environ, {"OMP_NUM_THREADS": "corrupted-value"}):
            with self.assertRaisesRegex(g.GuardError, "environment changed"):
                g.preflight(staged["manifest"]["path"])

    def test_echo_only_import_is_not_runtime_evidence(self):
        staged = self.stage()
        corrupted = self.task / "explicitly_corrupted_real_stdout.txt"
        text = (self.source / "endpoint.out").read_text()
        corrupted.write_text("\n".join(line for line in text.splitlines()
                                      if not line.strip().startswith("Reading XTB parameter file")))
        with self.assertRaisesRegex(g.GuardError, "runtime import evidence"):
            g.postflight(staged["manifest"]["path"], corrupted, 0)

    def test_execution_failure_rejected(self):
        staged = self.stage()
        with self.assertRaisesRegex(g.GuardError, "exit successfully"):
            g.postflight(staged["manifest"]["path"], self.source / "endpoint.out", 1)

    def test_missing_import_directive_rejected(self):
        self.change_input(lambda s: s.replace(" ReadXTBParam true", " ReadXTBParam false"))
        with self.assertRaisesRegex(g.GuardError, "ReadXTBParam"):
            self.stage()

    def test_parameter_write_enabled_rejected(self):
        self.change_input(lambda s: s.replace(" WriteXTBParam false", " WriteXTBParam true"))
        with self.assertRaisesRegex(g.GuardError, "WriteXTBParam"):
            self.stage()

    def test_nonlocal_parameter_path_rejected(self):
        self.binding["parameter_filename"] = "../control.xtb.json"
        with self.assertRaisesRegex(g.GuardError, "task-local"):
            self.stage()

    def test_state_mismatch_rejected(self):
        self.binding["state"]["charge"] = -4
        with self.assertRaisesRegex(g.GuardError, "charge/multiplicity"):
            self.stage()

    def test_unannounced_element_rejected(self):
        self.binding["elements"].remove("La")
        with self.assertRaisesRegex(g.GuardError, "actual XYZ"):
            self.stage()

    def test_spin_polarized_method_not_silent_substitution(self):
        self.change_input(lambda s: s.replace("Native-GFN2-xTB", "Native-spGFN2-xTB"))
        with self.assertRaisesRegex(g.GuardError, "only Native-GFN2"):
            self.stage()

    def test_duplicate_directive_rejected(self):
        self.change_input(lambda s: s.replace(" ReadXTBParam true", " ReadXTBParam true\n ReadXTBParam true"))
        with self.assertRaisesRegex(g.GuardError, "ReadXTBParam"):
            self.stage()

    def test_cache_changes_for_scientific_bindings(self):
        # Direct digest property uses actual binding with explicit changed metadata.
        original = g.digest(self.binding)
        for key, value in (("protocol_id", "different-protocol"),
                           ("state", {"charge": -4, "multiplicity": 2}),
                           ("model", {"method": g.METHOD, "medium": "vacuum"}),
                           ("parameter_sha256", g.ORIGINAL_SHA)):
            changed = copy.deepcopy(self.binding)
            changed[key] = value
            self.assertNotEqual(original, g.digest(changed))

    def test_runtime_approval_reusable_across_new_geometry_bindings(self):
        # A new task binding is not a new runtime; neither object is scored here.
        changed = copy.deepcopy(self.binding)
        changed["protocol_id"] = "separate-research-task"
        changed["source_inputs"] = [g.pin(self.source / "receipt.json")]
        self.assertNotEqual(g.digest(self.binding), g.digest(changed))
        self.assertEqual(g.runtime_identity(self.binding, g.sha256(PARAM)),
                         g.runtime_identity(changed, g.sha256(PARAM)))

    def test_actual_original_import_diagnostic(self):
        source = ARCHIVE / "La_water_original"
        self.binding["parameter_sha256"] = g.ORIGINAL_SHA
        staged = g.stage(self.binding, g.ORIGINAL, ORIGINAL_POLICY,
                         mode="diagnostic", variant="original")
        result = g.postflight(staged["manifest"]["path"], source / "endpoint.out", 0)
        self.assertFalse(result["runtime_qualified"])

    def test_actual_native_default_export_diagnostic(self):
        source = ARCHIVE / "La_water_default"
        shutil.copyfile(source / "endpoint.inp", self.task / "endpoint.inp")
        self.binding.update(input=g.pin(self.task / "endpoint.inp"), parameter_filename="endpoint.xtb.json",
                            parameter_sha256=g.ORIGINAL_SHA)
        staged = g.stage(self.binding, g.ORIGINAL, ORIGINAL_POLICY,
                         mode="diagnostic", variant="native-default")
        g.preflight(staged["manifest"]["path"])
        # Copy the ACTUAL archived export, not a fabricated runtime artifact.
        shutil.copyfile(source / "endpoint.xtb.json", self.task / "endpoint.xtb.json")
        result = g.postflight(staged["manifest"]["path"], source / "endpoint.out", 0)
        self.assertFalse(result["runtime_qualified"])

    def test_old_control_receipts_do_not_qualify_unpinned_runtime(self):
        helper = self.executable.parent / "orca_leanscf"
        self.assertTrue(helper.exists(), "actual native helper is required by this fixture")
        self.binding["runtime_files"] = [dict(g.pin(helper), role="native_helper")]
        receipts = [g.pin(ARCHIVE / f"La_{medium}_corrected/receipt.json") for medium in ("vacuum", "water")]
        registry = {"schema_version": 1, "status": "qualified_for_research",
                    "runtime_identity": g.runtime_identity(self.binding, g.sha256(PARAM)),
                    "qualified_elements": self.binding["elements"],
                    "qualified_media": ["vacuum", "ALPB(Water)"], "qualified_multiplicities": [1], "receipts": receipts}
        path = self.task / "test_registry_from_real_archived_receipts.json"
        path.write_text(json.dumps(registry))
        # The old numerical controls are real, but do not pin today's libraries.
        with self.assertRaisesRegex(g.GuardError, "current runtime binding"):
            g.validate_approval(g.pin(path), self.binding, g.audit_parameters(PARAM, POLICY))

    def current_binding(self):
        runtime = ROOT / "workspaces/orca_parameter_audit_20261012/runtime_qualification_v1/RUNTIME_BINDING.json"
        registry = ROOT / "diagnostics/orca_parameter_audit_20261012/parameter_guard/runtime_registry_ca_la_corrected_20261012.json"
        if not runtime.exists() or not registry.exists():
            self.skipTest("actual current-runtime control receipts not available")
        approval = g.load_json(registry)
        binding = g.build_binding(input_path=self.task / "endpoint.inp", xyz=self.task / "core.xyz",
                                  parameters=PARAM, runtime_binding=runtime,
                                  observed_runtime=approval["observed_runtime"]["path"],
                                  source_inputs=[self.source / "core.xyz"], protocol_id="real-control-contract-test",
                                  charge=-5, multiplicity=1, medium="ALPB(Water)", parameter_filename="control.xtb.json")
        return binding, registry

    def test_current_real_runtime_research_preflight_postflight(self):
        binding, registry = self.current_binding()
        with patch.dict(os.environ, {k: v for k, v in binding["environment"].items() if v is not None}):
            for key, value in binding["environment"].items():
                if value is None:
                    os.environ.pop(key, None)
            staged = g.stage(binding, PARAM, POLICY, approval=g.pin(registry))
            self.assertTrue(staged["runtime_qualified"])
            g.preflight(staged["manifest"]["path"])
            actual_stdout = ROOT / "workspaces/orca_parameter_audit_20261012/runtime_qualification_v1/results/La_water_corrected/endpoint.out"
            result = g.postflight(staged["manifest"]["path"], actual_stdout, 0)
            self.assertTrue(result["parameter_postflight_passed"])
            self.assertIsNone(result["scientific_result_accepted"])

    def test_runtime_approval_rejects_changed_libraries(self):
        binding, registry = self.current_binding()
        binding["runtime_files"][0]["sha256"] = g.ORIGINAL_SHA
        with self.assertRaisesRegex(g.GuardError, "runtime approval identity differs"):
            g.validate_approval(g.pin(registry), binding, g.audit_parameters(PARAM, POLICY))

    def test_runtime_approval_rejects_unqualified_element(self):
        binding, registry = self.current_binding()
        binding["elements"].append("Cu")
        with self.assertRaisesRegex(g.GuardError, "element lacks runtime"):
            g.validate_approval(g.pin(registry), binding, g.audit_parameters(PARAM, POLICY))

    def test_runtime_approval_rejects_unqualified_multiplicity(self):
        binding, registry = self.current_binding()
        binding["state"]["multiplicity"] = 3
        with self.assertRaisesRegex(g.GuardError, "multiplicity bookkeeping"):
            g.validate_approval(g.pin(registry), binding, g.audit_parameters(PARAM, POLICY))

    def test_research_approval_rejects_corrupted_control_acceptance(self):
        binding, registry = self.current_binding()
        approval = g.load_json(registry)
        receipt = g.load_json(approval["receipts"][0]["path"])
        receipt["output"]["accepted"] = False
        corrupted = self.task / "explicitly_corrupted_real_receipt.json"
        corrupted.write_text(json.dumps(receipt))
        approval["receipts"][0] = g.pin(corrupted)
        bad_approval = self.task / "approval_pointing_to_corrupted_real_receipt.json"
        bad_approval.write_text(json.dumps(approval))
        with self.assertRaisesRegex(g.GuardError, "state/population checks failed"):
            g.validate_approval(g.pin(bad_approval), binding, g.audit_parameters(PARAM, POLICY))

    def test_actual_ca_corrected_import_uses_same_research_registry(self):
        source = ROOT / "workspaces/orca_parameter_audit_20261012/runtime_qualification_v1/results/Ca_water_corrected_file"
        if not (source / "receipt.json").exists():
            self.skipTest("actual corrected-full-file Ca control unavailable")
        for name in ("core.xyz", "endpoint.inp"):
            shutil.copyfile(source / name, self.task / name)
        binding, registry = self.current_binding()
        binding["state"]["charge"] = -6
        binding["source_inputs"] = [g.pin(source / "core.xyz")]
        with patch.dict(os.environ, {k: v for k, v in binding["environment"].items() if v is not None}):
            for key, value in binding["environment"].items():
                if value is None:
                    os.environ.pop(key, None)
            staged = g.stage(binding, PARAM, POLICY, approval=g.pin(registry))
            self.assertTrue(staged["runtime_qualified"])
            g.preflight(staged["manifest"]["path"])
            result = g.postflight(staged["manifest"]["path"], source / "endpoint.out", 0)
            self.assertTrue(result["parameter_postflight_passed"])


if __name__ == "__main__":
    unittest.main()
