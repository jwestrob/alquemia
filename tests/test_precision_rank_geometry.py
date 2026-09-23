"""Real scalar source copies and paired-rank finite scope; no invented energies."""
from pathlib import Path
import sys,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
import precision_rank_geometry as experiment

class RankGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=ROOT/'workspaces/precision_rank_geometry_20260923/run_v1';cls.ready=read_json(cls.root/'READY.json')

    def test_actual_finite_shards_and_contained_paths(self):
        tasks=[];cpu=0
        for pin in self.ready['manifests']:
            path=verify(pin);m=read_json(path);v=experiment.validate(path)
            self.assertEqual(v['new_calls'],2);self.assertEqual(v['new_physical_geometries'],0)
            cpu+=m['execution_resources']['mpi_ranks']*m['execution_resources']['concurrent_tasks']
            for t in m['tasks']:
                tasks.append((t['case_id'],t['geometry_version'],m['rank_count']))
                for k in ('input','xyz'):self.assertTrue(verify(t[k]).is_relative_to(path.parent))
                self.assertTrue(Path(t['output_path']).is_relative_to(path.parent))
        self.assertEqual((len(tasks),len(set(tasks)),cpu),(4,4,18))
        self.assertEqual({(g,r) for c,g,r in tasks},{('old',1),('new',8)})

    def test_exact_actual_source_bytes_and_fresh_recipe(self):
        for pin in self.ready['manifests']:
            for t in read_json(verify(pin))['tasks']:
                for k in ('input','xyz'):self.assertEqual(verify(t[k]).read_bytes(),verify(t['actual_source_task'][k]).read_bytes())
                body=verify(t['input']).read_text()
                self.assertIn('NoAutostart',body);self.assertIn('MaxIter 500',body)
                self.assertNotIn('EnGrad',body);self.assertEqual((t['charge'],t['multiplicity']),(-1,1))

    def test_prior_inventory_has_no_exact_missing_cell_reuse(self):
        inv=read_json(verify(self.ready['inventory']))
        self.assertEqual(sum(s['task_count'] for s in inv['searched_scopes']),616)
        self.assertFalse(any(s['exact_geometry_state_matches'] for s in inv['searched_scopes']))
        self.assertEqual(len(inv['missing_cells']),4)
        for target in inv['targets']:
            self.assertEqual(target['versions']['old']['actual_ranks'],8)
            self.assertEqual(target['versions']['new']['actual_ranks'],1)
            self.assertEqual(target['versions']['old']['original_recipe'],target['versions']['new']['original_recipe'])

    @unittest.skipUnless((ROOT/'workspaces/precision_rank_geometry_20260923/run_v1/COLLECTION.json').exists(),
                         'four actual molecular evaluations have not completed')
    def test_actual_complete_matrix_arithmetic(self):
        from affordable_common import HA_TO_KCAL
        r=read_json(self.root/'COLLECTION.json');self.assertEqual((r['fresh_denominator'],r['fresh_complete']),(4,4))
        for t in r['targets']:
            self.assertEqual(t['status'],'complete');m=t['matrix'];f=t['effects']
            for g in ('old','new'):
                d=(m[g]['1']['energy_hartree']-m[g]['8']['energy_hartree'])*HA_TO_KCAL
                self.assertEqual(f['rank1_minus_rank8_kcal_mol'][g],d)
                self.assertEqual(f['rank_gate_pass'][g],abs(d)<=.1)
            for rank in ('1','8'):
                d=(m['new'][rank]['energy_hartree']-m['old'][rank]['energy_hartree'])*HA_TO_KCAL
                self.assertEqual(f['new_minus_old_geometry_kcal_mol'][rank],d)
            self.assertAlmostEqual(f['interaction_kcal_mol'],f['new_minus_old_geometry_kcal_mol']['1']-f['new_minus_old_geometry_kcal_mol']['8'],places=10)

if __name__=='__main__':unittest.main()
