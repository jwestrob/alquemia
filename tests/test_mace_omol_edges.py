"""Exact adapter tests on pinned real native inputs and executed receipts."""
import copy
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
from mace_omol import accepted_state
from mace_omol_edge_run import tasks,collect
W=ROOT/'workspaces/mace_omol_20260917'
N=W/'intact_core_v1/manifest.json'


@unittest.skipUnless(N.exists(),'requires pinned real native energy-only bridge')
class EdgeTests(unittest.TestCase):
    def test_exact_native_scientific_inputs_at_both_chunk_sizes(self):
        native=read_json(N);expected=tasks(native,'edge_core');self.assertEqual(len(expected),8)
        sources={t['task_id']:t for t in native['tasks']}
        for t in expected:
            s=sources[t['reference_task_id']]
            for key in ('xyz','source_xyz','charge','state','spin_multiplicity','metal_index','energy_component'):
                self.assertEqual(t[key],s[key])
            self.assertIn(t['edge_adapter']['chunk_size'],(1024,2048))

    def test_native_receipt_does_not_satisfy_adapter_task(self):
        n=read_json(N);c=read_json(W/'intact_core_v1/collection_job_1200807.json')
        for t in tasks(n,'edge_core'):
            r=c['rows'][t['reference_task_id']];self.assertFalse(accepted_state(r,t))

    def test_CPU_tasks_keep_geometry_and_reject_GPU_receipts(self):
        n=read_json(N);c=read_json(W/'intact_core_v1/collection_job_1200807.json');ts=tasks(n,'cpu_core')
        self.assertEqual(len(ts),4)
        source={t['task_id']:t for t in n['tasks']}
        for t in ts:
            self.assertEqual(t['execution_device'],'cpu');self.assertNotIn('edge_adapter',t)
            s=source[t['reference_task_id']]
            for key in ('xyz','charge','spin_multiplicity','state','metal_index'):self.assertEqual(t[key],s[key])
            self.assertFalse(accepted_state(c['rows'][t['reference_task_id']],t))

    def test_executed_core_equivalence_and_complete_edge_inventory(self):
        files=sorted((W/'edge_core_v1').glob('collection_job_*.json'))
        if not files:self.skipTest('actual edge-batched scientific execution not yet available')
        c=read_json(files[-1]);m=read_json(verify(c['manifest']))
        self.assertEqual(c,collect(verify(c['manifest'])));self.assertTrue(c['numerical_gate_pass'])
        self.assertEqual(len(c['checks']),20)
        for t in m['tasks']:
            r=c['rows'][t['task_id']];self.assertTrue(accepted_state(r,t))
            native=copy.deepcopy(t);native.pop('edge_adapter');self.assertFalse(accepted_state(r,native))
            changed=copy.deepcopy(r);changed['execution_adapter']['layers'][0]['processed_edges']-=1
            self.assertFalse(accepted_state(changed,t))


if __name__=='__main__':unittest.main()
