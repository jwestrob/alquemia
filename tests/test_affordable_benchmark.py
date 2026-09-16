"""Benchmark checks on real pinned preparations and archived aquo outputs."""
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_benchmark import audit_reference,collect
from affordable_common import read_json,verify,energy,contrast,classify_raw,digest
from affordable_workflow import dry_run

class BenchmarkRegression(unittest.TestCase):
    def test_archived_reference_algebra(self):
        r=audit_reference(ROOT)
        self.assertAlmostEqual(r['delta_E_aquo_hartree'],-646.0775458314704,places=10)
        self.assertEqual(len(r['artifacts']),4)

    def test_prepared_tasks_and_preserved_bytes(self):
        mp=ROOT/'workspaces/baseline_benchmark_20260915/run_v1/manifest.json'
        if not mp.exists(): self.skipTest('real prepared benchmark not available')
        m=read_json(mp)
        self.assertEqual(dry_run(mp)['tasks'],26)
        self.assertEqual(len(m['cases']),13)
        for t in m['tasks']:
            self.assertEqual(verify(t['input']).read_bytes(),verify(t['source_input']).read_bytes())
            self.assertEqual(verify(t['xyz']).read_bytes(),verify(t['source_xyz']).read_bytes())
        for c in m['cases']:
            self.assertEqual(c['paired_invariants']['status'],'pass')
        ae=[c['case'] for c in m['cases'] if c['biological_group']=='aequorin' and c['lane']=='repaired']
        self.assertEqual(ae,['aequorin_1sl8_EF1','aequorin_1sl8_EF3','aequorin_1sl8_EF4'])

    def test_actual_energies_do_not_transfer_threshold(self):
        release=read_json(ROOT/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json')
        e={m:energy(verify(r['output'])) for m,r in release['scores'][0]['artifacts'].items()}
        s=contrast(e['Ca'],e['La'])
        self.assertIsNone(s['S_kcal_mol'])
        self.assertEqual(classify_raw(s['R_kcal_mol'],release,'generic_peptide_amide_vertical_native_r2scan3c_v3'),'uncalibrated_protocol')

    def test_real_receipt_rejects_partial_or_changed_output(self):
        from run_orca_task_manifest import load_manifest_tasks,_completed_attempt_is_valid
        mp=ROOT/'workspaces/baseline_benchmark_20260915/run_v1/manifest.json'
        if not mp.exists(): self.skipTest('real benchmark unavailable')
        m,tasks=load_manifest_tasks(mp)
        done=[t for t in tasks if Path(str(t['output'])+'.execution.json').exists()]
        if not done: self.skipTest('no completed scientific receipt yet')
        t=done[0];op=t['output'];rp=Path(str(op)+'.execution.json')
        args=dict(manifest_sha256=digest(mp),task=t,runner_identity=m['execution_policy']['task_runner'],runtime_renderer_identity=m['execution_policy']['runtime_renderer'])
        self.assertTrue(_completed_attempt_is_valid(rp,op,**args))
        with tempfile.TemporaryDirectory() as d:
            bad=Path(d)/'explicitly_corrupted_real_output.out'
            bad.write_text(op.read_text().replace('ORCA TERMINATED NORMALLY','REMOVED BY MALFORMED-INPUT TEST'))
            self.assertFalse(_completed_attempt_is_valid(rp,bad,**args))
            self.assertFalse(_completed_attempt_is_valid(Path(d)/'absent_receipt.json',op,**args))

if __name__=='__main__':unittest.main()
