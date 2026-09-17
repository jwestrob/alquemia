"""Pinned real preparations/receipts; no fabricated scientific outputs."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,cache_key,read_json,verify,xyz
from mace_omol import accepted_state
from mace_omol_intact import preparations,expected_tasks,geometry,collect,validate
W=ROOT/'workspaces/mace_omol_20260917'
P=ROOT/'workspaces/mace_global_benchmark_20260916/prepared_v1/preparation_manifest.json'


@unittest.skipUnless(P.exists(),'requires pinned real intact physical preparations')
class IntactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.physical=preparations(P)
        cls.tasks=expected_tasks('intact_qualification',[],cls.physical)

    def test_real_metal_mapping_and_only_metal_moves(self):
        self.assertEqual(len(self.tasks),14)
        for t in self.tasks:
            source=xyz(verify(t['source_xyz']));actual=geometry(t);i=t['metal_index']
            self.assertGreater(i,0);self.assertEqual(source[i][0],t['metal'])
            self.assertEqual(len(source),1932)
            self.assertEqual(t['charge'],-4 if t['metal']=='La' else -5)
            self.assertEqual(t['spin_multiplicity'],1)
            if t['variant']!='rotate':
                self.assertEqual(source[:i]+source[i+1:],actual[:i]+actual[i+1:])
            if t['position']=='detached':
                pos=np.array([r[1:] for r in actual])
                self.assertGreater(np.linalg.norm(np.delete(pos,i,axis=0)-pos[i],axis=1).min(),6.)

    def test_paired_coordinates_and_native_inventory(self):
        for name,p in self.physical.items():
            c=p['data'];self.assertFalse(any('cap' in a['kind'] for a in c['physical_atoms']))
            la,ca=[xyz(verify(c['endpoints'][m]['xyz'])) for m in ('La','Ca')]
            self.assertEqual([r[1:] for r in la],[r[1:] for r in ca])
            self.assertEqual(c['endpoints']['La']['charge']-c['endpoints']['Ca']['charge'],1)
            self.assertEqual([i for i,(a,b) in enumerate(zip(la,ca)) if a[0]!=b[0]],[p['metal_index']])

    def test_energy_only_does_not_claim_or_satisfy_force_task(self):
        paths=sorted((W/'intact_core_v1').glob('collection_job_*.json'))
        if not paths:self.skipTest('actual energy-only inference not yet available')
        c=read_json(paths[-1]);m=read_json(verify(c['manifest']))
        self.assertTrue(c['numerical_gate_pass'])
        for t in m['tasks']:
            r=c['rows'][t['task_id']];self.assertTrue(accepted_state(r,t))
            self.assertIsNone(r['forces']);self.assertFalse(r['gradient_computation_enabled'])
            force_task=copy.deepcopy(t);force_task.pop('energy_only')
            self.assertFalse(accepted_state(r,force_task))
            corrupted=copy.deepcopy(r);corrupted['input_state_check']['selected_metal_index']=1
            self.assertFalse(accepted_state(corrupted,t))

    def test_actual_full_numerical_qualification(self):
        paths=sorted((W/'intact_qualification_v1').glob('collection_job_*.json'))
        if not paths:self.skipTest('actual intact qualification not yet available')
        c=read_json(paths[-1]);self.assertEqual(c,collect(verify(c['manifest'])))
        if c['status']!='complete':
            self.assertFalse(c['numerical_gate_pass'])
            self.skipTest('native intact execution incomplete after real A5000 OOM; H200 recovery pending')
        self.assertTrue(c['numerical_gate_pass']);self.assertEqual(len(c['rows']),14)
        for key,r in c['rows'].items():
            if '_detached_' in key:self.assertEqual(r['input_state_check']['metal_neighbor_edge_count'],0)

    def test_real_OOM_has_no_invented_energy_or_pass(self):
        p=W/'intact_qualification_v1/collection_job_1200808.json'
        if not p.exists():self.skipTest('requires actual failed native A5000 attempt')
        c=read_json(p);self.assertEqual(c['status'],'incomplete');self.assertFalse(c['numerical_gate_pass'])
        self.assertTrue(all(r['energy_eV'] is None for r in c['rows'].values()))
        self.assertEqual(len(c['attempts']),1);attempt=c['attempts'][0];self.assertFalse(attempt['accepted'])
        receipt=read_json(verify(attempt['receipt']));failure=read_json(verify(receipt['result']))
        self.assertEqual(failure['status'],'failed');self.assertIn('out of memory',failure['reason'])

    def test_unrecorded_execution_change_is_rejected_even_with_new_cache_key(self):
        p=W/'intact_core_v1/manifest.json'
        if not p.exists():self.skipTest('requires real prepared core bridge')
        m=copy.deepcopy(read_json(p));t=m['tasks'][0];t['execution_device']='cpu'
        payload={k:v for k,v in t.items() if k!='cache_key'}
        t['cache_key']=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']})
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'corrupted_real_manifest.json';p.write_text(json.dumps(m))
            with self.assertRaises(InvalidArtifact):validate(p)


if __name__=='__main__':unittest.main()
