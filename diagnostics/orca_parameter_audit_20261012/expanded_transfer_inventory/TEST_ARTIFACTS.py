"""Read-only real-artifact tests; no molecular calls."""
import importlib.util,json,unittest
from pathlib import Path
S=Path(__file__).with_name('REEVALUATE.py')
spec=importlib.util.spec_from_file_location('expanded_recheck',S)
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
class ArchiveTests(unittest.TestCase):
    def test_released_pool_replay_and_missing_candidate(self):
        cfg=json.loads(S.with_name('CONFIG.json').read_text());count=0
        for source in cfg['collections']:
            for case in r.f.read(source)['cases']:
                matrix={z:{q:case['matrix'][z][q].get('components') for q in r.ORDER} for z in ('Ca','La')}
                got=r.choose(matrix)
                if got['status']!='available':continue
                count+=1
                for policy in ('mathematical','operational'):
                    self.assertAlmostEqual(got[policy]['R_model_kcal_mol'],case['pool'][policy]['composite_R_model_kcal_mol'],places=10)
                original=got['origin']['R_model_kcal_mol']
                # Deliberately unavailable copy of an actual archived cell.
                matrix['La']['adaptive_Ca']=None
                missing=r.choose(matrix)
                self.assertIsNone(missing['operational']);self.assertIsNone(missing['mathematical'])
                self.assertEqual(missing['origin']['R_model_kcal_mol'],original)
        self.assertEqual(count,134)
    def test_frozen_population_and_cost_scope(self):
        m=r.f.read('workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3/manifest.json')
        self.assertEqual(len(m['cases']),135);self.assertEqual(len(m['tasks']),1620)
        self.assertEqual(sum(t['reuse'] is not None for t in m['tasks']),24)
        self.assertEqual(sum(t['stage']=='canonical' and t['reuse'] is None for t in m['tasks']),276)
        self.assertEqual({c['matrix']['Ca']['origin']['charge'] for c in m['cases']},{-4,-3,-2})
        self.assertEqual(sum(t['status']!='prepared' for t in m['triples']),6)
        self.assertEqual(sum(c['role']=='canonical_calibration' and c['expected_class']=='Ca' for c in m['cases']),14)
if __name__=='__main__':unittest.main(verbosity=2)
