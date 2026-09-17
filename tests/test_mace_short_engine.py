"""Actual short-readout execution and separation from charged total energies."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
from mace_hybrid import accepted_attempt
from mace_short_engine import validate,collect_short
MANIFEST=ROOT/'workspaces/mace_short_engine_20260916/pilot_v1/manifest.json'
COLLECTION=MANIFEST.with_name('collection_job_1200736.json')


@unittest.skipUnless(COLLECTION.exists(),'requires actual executed short engine validation')
class ShortEngineTests(unittest.TestCase):
    def test_actual_energy_gradient_gates_and_no_invented_density(self):
        self.assertEqual(validate(MANIFEST)['tasks'],34)
        c=read_json(COLLECTION)
        self.assertTrue(c['numerical_checks_pass']);self.assertEqual(len(c['checks']),50)
        for r in c['rows'].values():
            self.assertIsNone(r['density_coefficients']);self.assertIsNone(r['charge_check'])
            self.assertEqual(r['energy_component'],'interaction_energy')
            self.assertEqual(r['force_definition'],'negative_Cartesian_gradient_of_interaction_energy')
            self.assertTrue(np.isfinite(np.load(verify(r['forces']))).all())

    def test_short_receipt_cannot_satisfy_total_energy_task(self):
        m=read_json(MANIFEST);t=m['tasks'][0]
        attempt=MANIFEST.parent/'execution'/t['task_id']/'attempt_001'
        self.assertIsNotNone(accepted_attempt(attempt,t,MANIFEST))
        corrupted=dict(t);corrupted.pop('energy_component')
        self.assertIsNone(accepted_attempt(attempt,corrupted,MANIFEST))

    def test_missing_output_is_not_a_passing_gate(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'manifest.json';p.write_bytes(MANIFEST.read_bytes())
            c=collect_short(p);self.assertEqual(c['status'],'incomplete')
            self.assertFalse(c['numerical_checks_pass']);self.assertIsNone(c['relaxation_correction_kcal_mol'])


if __name__=='__main__':unittest.main()
