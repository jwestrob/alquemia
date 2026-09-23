"""Actual-source rank diagnostics; no successful scientific fixture synthesis."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import native_gfn2_rank_scaling as a
DESIGN=ROOT/'workspaces/native_gfn2_ranks_20260923/run_v1/design.json'

class NativeRanks(unittest.TestCase):
    def test_exact_24_sources_and_rank_only_changes(self):
        d=a.read_json(DESIGN);self.assertEqual(d['new_calls'],24);observed=[]
        for pin in d['manifests']:
            mp=a.verify(pin);m=a.read_json(mp);self.assertEqual(a.validate(mp)['tasks'],8)
            observed.append(m['ranks'])
            for t in m['tasks']:
                self.assertEqual(t['candidate'],'adaptive_Ca');self.assertEqual(t['multiplicity'],1)
                self.assertIn('NoAutostart',a.verify(t['input']).read_text())
        self.assertEqual(observed,[1,4,8])

    def test_actual_paired_geometries_and_archived_complete_receipts(self):
        m=a.read_json(a.verify(a.read_json(DESIGN)['manifests'][0]))
        for cid in a.CASES:
            ca,la=[a.xyz(a.verify(next(t['xyz'] for t in m['tasks'] if t['case_id']==cid and t['metal']==z and t['medium']=='vacuum'))) for z in ('Ca','La')]
            self.assertEqual(ca[1:],la[1:]);self.assertEqual(ca[0][1:],la[0][1:])
        for t in m['tasks']:
            rec=a.read_json(a.verify(t['archived']['receipt']))
            self.assertTrue(rec['normal_termination']);self.assertTrue(rec['scf_converged']);self.assertEqual(rec['parallelism']['nprocs'],8)

    def test_corrupted_actual_rank_scope_rejected(self):
        m=copy.deepcopy(a.read_json(a.verify(a.read_json(DESIGN)['manifests'][0])));m['execution_resources']['mpi_ranks']=2
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'corrupted_actual_manifest.json';a.write_new(p,m)
            with self.assertRaises(a.InvalidArtifact):a.validate(p)

    def test_actual_final_comparison_when_executed(self):
        p=DESIGN.parent/'COMPARISON.json'
        if not p.exists():self.skipTest('native rank experiment has not completed')
        r=a.read_json(p);self.assertEqual([x['ranks'] for x in r['rank_rows']],[1,4,8])
        for row in r['rank_rows']:
            self.assertEqual(len(row['cells']),8);self.assertEqual(len(row['pairs']),2)
            for cell in row['cells']:
                if cell['status']=='available':self.assertEqual(cell['passes_component_gate'],abs(cell['delta_vs_fresh8_kcal_mol'])<=.1)
            for pair in row['pairs']:
                c=pair['components']
                if c:self.assertAlmostEqual(c['native_R_model_kcal_mol']+c['ALPB_R_kcal_mol']-c['vacuum_R_kcal_mol'],c['composite_R_model_kcal_mol'],places=9)
        self.assertFalse(r['production_changed'])

if __name__=='__main__':unittest.main()
