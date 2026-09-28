"""Actual prepared-geometry invariants; no molecular evaluations."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import read_json,verify,xyz
ROOT=Path(__file__).resolve().parents[1]
class CompactHydrogenTransferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=read_json(ROOT/'workspaces/metal_environment_response_20260926/compact_CH_repair_v1/manifest.json')
    def test_real_all_source_fixed_and_paired_coordinates(self):
        self.assertEqual({x['source_id'] for x in self.m['sources']},{'Hans8DQ2','Hans8FNR','Mex8FNS'})
        for r in self.m['sources']:
            self.assertTrue(all(r['checks'].values()));new={k:xyz(verify(v)) for k,v in r['xyz'].items()};np.testing.assert_array_equal([a[1:] for a in new['La']],[a[1:] for a in new['Dy']])
            for metal in ('La','Dy'):
                oldlines=verify(r['original_xyz'][metal]).read_text().splitlines();lines=verify(r['xyz'][metal]).read_text().splitlines()
                for i in range(len(new[metal])):
                    if i not in r['repaired_H_indices']:self.assertEqual(lines[i+2],oldlines[i+2])
    def test_real_graph_and_admission_and_count(self):
        expected={'Hans8DQ2':14,'Hans8FNR':14,'Mex8FNS':10}
        for r in self.m['sources']:
            rev=read_json(verify(r['reviewed_admission']));self.assertTrue(rev['admitted']);ex=read_json(verify(rev['parent_export']));atoms={a['id']:a for a in read_json(verify(ex['atoms']))};bonds={frozenset(b) for b in read_json(verify(ex['bonds']))};self.assertEqual(len(r['repaired_H_indices']),expected[r['source_id']])
            for a in r['transferred_atoms']:
                self.assertEqual(atoms[a['source_id']]['element'],'H');self.assertEqual(atoms[a['parent_id']]['element'],'C');self.assertIn(frozenset((a['source_id'],a['parent_id'])),bonds);self.assertIn(a['source_id'],rev['free_H_ids'])
if __name__=='__main__':unittest.main()
