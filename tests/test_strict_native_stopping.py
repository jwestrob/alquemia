"""Real-source recipe/state/two-start checks for the strict native diagnostic."""
import sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import strict_native_stopping as s
from affordable_common import read_json,verify,HA_TO_KCAL
from native_pool_continuation import recipe as old_recipe
D=ROOT/'workspaces/strict_native_stopping_20260923/run_v1';M=D/'manifest.json'


class StrictNativeTests(unittest.TestCase):
    def test_exact_twenty_real_starts_and_recipe_change(self):
        v=s.validate(M);self.assertEqual(v['tasks'],20);m=read_json(M)
        self.assertEqual(m['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':20})
        for t in m['tasks']:
            body=verify(t['input']).read_text()
            self.assertEqual(body.replace(' TolE 1e-10\n',''),old_recipe(t['charge'],1,t['medium']))
            self.assertEqual(verify(t['xyz']).read_bytes(),verify(t['source']['xyz']).read_bytes())
            self.assertFalse(t['gradient_requested'])

    def test_each_state_has_exact_cold_and_pass2_seed(self):
        m=read_json(M)
        for cid,q,medium in s.CELLS:
            for z in ('Ca','La'):
                rows={t['seed_kind']:t for t in m['tasks'] if (t['case_id'],t['candidate'],t['metal'],t['medium'])==(cid,q,z,medium)}
                self.assertEqual(set(rows),{'cold','pass2'})
                for field,k in [('immutable_seed','xtbw'),('immutable_gbw','gbw')]:
                    self.assertEqual(verify(rows['cold'][field]).read_bytes(),verify(rows['cold']['source'][k]).read_bytes())
                    after='preserved_after' if k=='xtbw' else 'preserved_gbw_after'
                    self.assertEqual(verify(rows['pass2'][field]).read_bytes(),verify(rows['pass2']['old_pass2']['seed_after'][after]).read_bytes())

    def test_actual_outputs_require_effective_tolerance_and_preserve_failures(self):
        path=D/'COLLECTION.json'
        if not path.exists():self.skipTest('strict-native molecular calls unrun')
        c=read_json(path);self.assertEqual(len(c['rows']),20);self.assertEqual(len(c['cells']),10)
        for r in c['rows']:
            if r['status']=='confirmed_restart':
                self.assertEqual(r['observed_TolE_hartree'],1e-10)
                self.assertIn('INITIAL GUESS: XTBRESTART',verify(r['actual']['output']).read_text())
                self.assertTrue(r['details']['native_mixer_observed'])
            else:self.assertIsNone(r['energy_hartree'])
        for r in c['cells']:
            if r['status']=='complete':
                delta=(r['from_pass2_hartree']-r['from_cold_hartree'])*HA_TO_KCAL
                self.assertAlmostEqual(delta,r['seed_difference_kcal_mol'],places=10)
                self.assertEqual(r['seed_agreement_pass'],abs(delta)<=.1)
            else:self.assertIsNone(r['seed_difference_kcal_mol'])
        self.assertIsNone(c['classifier_or_reference'])


if __name__=='__main__':unittest.main()
