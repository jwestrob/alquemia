"""Parser/algebra checks on pinned, actually computed scout; no molecular calls."""
import unittest, json, hashlib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
WORK=ROOT/'workspaces/orca_parameter_audit_20261012'

def read(path):return json.loads(path.read_text())

class RealScout(unittest.TestCase):
    def test_historical_absolute_gate_preserved(self):
        scout=read(WORK/'context_gfn2_v1/SCOUT.json')
        self.assertEqual(scout['declared'],4);self.assertEqual(scout['accepted'],0)
        for row in scout['rows']:
            self.assertFalse(row['accepted']);self.assertIsNone(row['energy_hartree'])
            self.assertGreater(abs(row['parity']['energy_difference_hartree']),1e-6)
            self.assertLess(abs(row['reference_energy_difference_hartree']),1e-6)

    def test_component_is_direct_same_geometry_algebra(self):
        coll=read(WORK/'context_gfn2_component_v2/SCOUT_COMPONENT_COLLECTION.json')
        self.assertEqual(coll['executed_inputs'],4);self.assertEqual(coll['accepted_components'],2)
        for row in coll['rows']:
            if not row['accepted']:continue
            real={}
            for medium,pin in row['inputs'].items():
                path=Path(pin['path']);self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),pin['sha256'])
                real[medium]=read(path)
            self.assertEqual(real['water']['task']['xyz'],real['vacuum']['task']['xyz'])
            primary=real['water']['backends']['xtb']['energy_hartree']-real['vacuum']['backends']['xtb']['energy_hartree']
            self.assertEqual(primary,row['primary_solvation_hartree'])
            self.assertLess(abs(row['component_error_hartree']),1e-6)
            self.assertFalse(row['original_absolute_energy_gate_pass'])

    def test_missing_is_unavailable(self):
        rows=read(WORK/'context_gfn2_component_v2/SCOUT_COMPONENT_COLLECTION.json')['rows']
        missing=[r for r in rows if r.get('failure_reason')=='not_executed']
        self.assertEqual(len(missing),38)
        for row in missing:
            self.assertIsNone(row['primary_solvation_hartree']);self.assertFalse(row['accepted'])

    def test_s9_causal_component_preserves_populations(self):
        rows=read(WORK/'context_gfn2_diagnosis_v1/RESULT.json')['rows']
        self.assertEqual(len(rows),2)
        for row in rows:
            self.assertLess(abs(row['energy_difference_Eh']),1e-6)
            self.assertLess(row['max_population_difference_e'],1e-5)
            self.assertGreater(abs(row['removed_ATM_energy_Eh']['xtb']-row['removed_ATM_energy_Eh']['tblite']),1e-6)

if __name__=='__main__':unittest.main()
