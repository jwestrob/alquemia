"""Real derived geometries/receipts; no fabricated predictions or fitted labels."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify,xyz
from mace_hybrid import accepted_attempt,write_xyz,EV_TO_KCAL
from mace_omol_coordination import validate,collect
W=ROOT/'workspaces/mace_omol_20260917'
M=W/'coordination_qualification_v1/manifest.json'


@unittest.skipUnless(M.exists(),'requires real prepared coordination qualification')
class CoordinationTests(unittest.TestCase):
    def test_only_metal_moves_and_charge_spin_inventory_are_preserved(self):
        m=read_json(M);self.assertEqual(len(m['tasks']),10)
        for task in m['tasks']:
            source=xyz(verify(task['source_xyz']));actual=xyz(verify(task['xyz']))
            self.assertEqual([r[0] for r in source],[r[0] for r in actual])
            pos=np.array([r[1:] for r in actual])
            self.assertGreater(np.linalg.norm(pos[1:]-pos[0],axis=1).min(),6.)
            self.assertEqual(task['charge'],0 if task['metal']=='La' else -1)
            self.assertEqual(task['spin_multiplicity'],1)
            if task['variant']!='rotate':self.assertEqual(source[1:],actual[1:])

    def test_bound_receipt_cannot_satisfy_detached_task(self):
        m=read_json(M);task=m['tasks'][0]
        path=W/'benchmark_v1/execution'/task['source_task_id']/'attempt_001'
        self.assertTrue(path.exists())
        self.assertIsNone(accepted_attempt(path,task,M))

    def test_corrupted_real_ligand_displacement_is_rejected(self):
        m=copy.deepcopy(read_json(M));task=m['tasks'][0]
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);rows=xyz(verify(task['xyz']));row=list(rows[1]);row[1]+=.01;rows[1]=tuple(row)
            p=d/'corrupted_real_geometry.xyz';write_xyz(p,rows);task['xyz']=record(p)
            p=d/'corrupted_manifest.json';p.write_text(json.dumps(m))
            with self.assertRaises(InvalidArtifact):validate(p)

    def test_executed_numerical_qualification(self):
        paths=sorted(M.parent.glob('collection_job_*.json'))
        if not paths:self.skipTest('actual detached-reference scientific execution unavailable')
        actual=collect(M);saved=read_json(paths[-1])
        self.assertEqual(actual,saved)
        self.assertTrue(actual['numerical_gate_pass'])
        self.assertEqual(len(actual['checks']),19)
        self.assertTrue(all(r['input_state_check']['metal_neighbor_edge_count']==0 for r in actual['rows'].values()))

    def test_actual_four_endpoint_algebra_has_no_inherited_aquo_offset(self):
        paths=sorted((W/'coordination_report_v1').glob('result.json'))
        if not paths:self.skipTest('actual matched coordination benchmark report unavailable')
        r=read_json(paths[0]);rows=r['canonical_rows']+list(r['nonPQQ_rows'].values())
        self.assertEqual(len(rows),32)
        for row in rows:
            if row['status']!='computed':continue
            c,l=(row['endpoints'][m] for m in ('Ca','La'))
            expected=((c['bound']['energy_eV']-c['detached']['energy_eV'])-
                      (l['bound']['energy_eV']-l['detached']['energy_eV']))*EV_TO_KCAL
            self.assertEqual(row['R_kcal_mol'],expected)
        self.assertIsNone(r['binding_free_energy_kcal_mol']);self.assertFalse(r['baseline_changed'])


if __name__=='__main__':unittest.main()
