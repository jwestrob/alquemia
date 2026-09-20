"""Read actual DFT receipts and immutable MACE triples; no molecular calls."""
import json
from pathlib import Path
import statistics
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import accommodation_fold_dft_triples as triples


class ExactTripleComparator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base=ROOT/'workspaces/accommodation_fold_DFT_20260920'
        cls.dft=cls.base/'triples_v1/dft_partial_collection.json'
        cls.mace=ROOT/'workspaces/accommodation_goal_20260920/folds_v1/triples_v1/result.json'
        cls.plan=ROOT/'diagnostics/accommodation_fold_DFT_20260920/THREE_FOLD_DFT_PLAN.md'
        cls.result=cls.base/'triples_v1/comparison_with_units/result.json'
        if not all(p.exists() for p in (cls.dft,cls.mace,cls.result)):
            raise unittest.SkipTest('real partial DFT and completed MACE triple artifacts unavailable')
        cls.actual=json.loads(cls.result.read_text())

    def test_exact100_members_and_unchanged_mace(self):
        original=json.loads(self.mace.read_text())
        old={(r['root_case_id'],tuple(sorted(r['members'])),r['representation']):r for r in original['rows']}
        self.assertEqual(len(self.actual['rows']),100)
        self.assertEqual(len(self.actual['proteins']),25)
        for r in self.actual['rows']:
            for rep in ('core','context'):
                prior=old[r['root_case_id'],tuple(sorted(r['members'])),rep]
                for method in ('native','composite'):
                    self.assertEqual(r['methods'][f'MACE_{method}_{rep}'],prior['methods'][method])

    def test_real_completed_energy_algebra_and_strict_median(self):
        d=json.loads(self.dft.read_text());index={r['case_id']:r for r in d['rows']};available=0
        self.assertGreater(len(self.actual['verified_DFT_endpoint_receipts']),0)
        for row in self.actual['rows']:
            e=row['methods']['DFT_baseline'];values=[]
            for cid in row['members']:
                r=index[cid]
                if r['R_kcal_mol'] is None:continue
                value=(r['endpoints']['Ca']['energy_hartree']-r['endpoints']['La']['energy_hartree'])*627.509474
                self.assertEqual(value,r['R_kcal_mol']);values.append(value)
            if len(values)==3:
                available+=1;self.assertEqual(e['median_R_kcal_mol'],statistics.median(values))
            else:
                self.assertIsNone(e['median_R_kcal_mol']);self.assertEqual(e['outcome'],'unavailable')
                self.assertEqual(len(e['missing_members']),3-len(values))
        self.assertGreater(available,0)

    def test_actual_prelaunch_snapshot_stays_entirely_unavailable(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=triples.compare(self.base/'run_v1/prelaunch_collection.json',self.mace,self.plan,Path(tmp)/'output')
            self.assertEqual(r['counts']['DFT_baseline']['triples'],{'correct':0,'wrong':0,'inconclusive':0,'unavailable':100})
            self.assertEqual(r['counts']['MACE_native_core'],self.actual['counts']['MACE_native_core'])

    def test_explicitly_corrupted_real_band_record_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad=json.loads(self.dft.read_text());bad['bands']['Ca_max']+=1.
            path=Path(tmp)/'corrupted_real_collection.json';path.write_text(json.dumps(bad))
            with self.assertRaisesRegex(triples.InvalidArtifact,'DFT bands differ'):
                triples.compare(path,self.mace,self.plan,Path(tmp)/'should_not_exist')

    def test_explicitly_corrupted_real_triple_membership_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad=json.loads(self.mace.read_text());bad['rows'].pop()
            path=Path(tmp)/'corrupted_real_MACE_triples.json';path.write_text(json.dumps(bad))
            with self.assertRaisesRegex(triples.InvalidArtifact,'triple denominator changed'):
                triples.compare(self.dft,path,self.plan,Path(tmp)/'should_not_exist')


if __name__=='__main__':unittest.main()
