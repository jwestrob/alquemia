"""Real QM/forcefield boundary invariants and executed solvent comparisons."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,write_new
from mace_omol_solvent import validate,ligand_ledger,TOL
W=ROOT/'workspaces/mace_omol_20260917';M=W/'full_boundary_GB_v1/manifest.json'

@unittest.skipUnless(M.exists(),'real solvent preparation not available')
class FullBoundary(unittest.TestCase):
    def test_no_forcefield_charge_on_QM_and_local_boundary_closure(self):
        m=read_json(M)
        for pin in m['states'].values():
            s=read_json(verify(pin));ids=[a['id'] for a in s['physical_atoms']];lookup={k:i for i,k in enumerate(ids)}
            support=set(s['projection_support_ids']);env=np.array(s['environment_charges_e'])
            self.assertTrue(all(env[lookup[k]]==0. for k in support))
            self.assertLess(abs(env.sum()-s['environment_charge']),TOL['boundary_charge_e'])
            for ledger in s['boundary_ledger']:
                for a,b in ledger['recipient_bonds']:
                    self.assertIn(a,support);self.assertNotIn(b,support)
                    self.assertIn(b,ledger['recipient_ids']);self.assertEqual(b.rsplit('/',1)[0],ledger['residue'])
                self.assertAlmostEqual(sum(env[lookup[k]] for k in ids if k.rsplit('/',1)[0]==ledger['residue']),ledger['exterior_target_e'],places=8)
            for metal,e in s['endpoints'].items():
                np.testing.assert_array_equal(np.array(e['QM_charges_e'])+env,e['full_charges_e'])
                self.assertLess(abs(sum(e['full_charges_e'])-e['charge']),TOL['charge_e'])

    def test_same_physical_GGR_boundary_and_exact_solver_inventory(self):
        m=read_json(M);a,b=(read_json(verify(m['states'][k])) for k in ('GGR_extended','GGR_connected'))
        self.assertEqual(a['physical_atoms'],b['physical_atoms'])
        self.assertNotEqual(a['projection_support_ids'],b['projection_support_ids'])
        self.assertEqual(validate(M)['tasks'],48)
        self.assertEqual(sum(t['variant']=='identity' for t in m['tasks']),8)
        self.assertEqual(sum(t['platform']=='Reference' for t in m['tasks']),2)

    def test_corrupted_actual_charged_group_and_cavity_are_rejected(self):
        c=read_json(W/'matched_H_prepared_v1/preparation.json');case=read_json(verify(c['cases']['GGR_extended']))
        sources={a['physical_id'] for a in case['mapping'] if a['kind']=='source'};sources.remove('A/134//OD1')
        with self.assertRaisesRegex(InvalidArtifact,'incomplete QM charged'):ligand_ledger(case,[],sources)
        m=copy.deepcopy(read_json(M));m['model']['radii_A']['La']=2.
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'corrupted_actual_manifest.json';write_new(f,m)
            with self.assertRaisesRegex(InvalidArtifact,'scientific settings'):validate(f)

    def test_actual_scientific_results_have_independent_numerical_and_accuracy_gates(self):
        p=W/'full_boundary_GB_report_v1/result.json'
        if not p.exists():self.skipTest('real solvent execution/report not complete')
        r=read_json(p);self.assertEqual(r['status'],'complete');self.assertEqual(len(r['rows']),48)
        self.assertEqual(r['numerical_gate_pass'],all(c['pass'] for c in r['checks']))
        self.assertEqual(r['partition_gate_pass'],abs(r['partition_shift_kcal_scale'])<=2)
        self.assertEqual(r['ordering_gate_pass'],r['numerical_gate_pass'] and all(c['difference_kcal_scale']>.02 for c in r['contrasts']))
        self.assertEqual(len(r['contrasts']),4)
        for c in r['cases'].values():
            self.assertAlmostEqual(c['GB_R_kcal_mol'],c['GB_QM_R_kcal_mol']+c['GB_cross_R_kcal_mol'],places=7)
        self.assertIsNone(r['aqueous_affinity_score']);self.assertIsNone(r['combined_gradient']);self.assertFalse(r['baseline_changed'])

if __name__=='__main__':unittest.main()
