"""Real frozen225 selectors; no generated scientific values."""
from pathlib import Path
import sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import read_json,verify,write_new,InvalidArtifact
import union_adaptive_transfer as transfer
ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'workspaces/union_adaptive_20260923/transfer225_v1'

class Selection(unittest.TestCase):
    def test_four_disjoint_actual_shards(self):
        top=read_json(RUN/'SELECTION.json');ids=[]
        for i in range(4):
            p=RUN/f'INPUTS_shard_{i}.json';r=transfer.validate_selection(p,verify(top['calibration']))
            self.assertEqual(r['new_sources'],51);ids.extend(c['case_id'] for c in read_json(p)['cases'])
        self.assertEqual((len(ids),len(set(ids))),(204,204))
        self.assertFalse(set(ids)&set(top['reuse_case_ids']))
        self.assertEqual(len(top['all_cases']),225)
        self.assertEqual(sum(c['union_origin_status']!='complete' for c in top['all_cases']),17)

    def test_no_canonical_folds_or_label_replacement(self):
        top=read_json(RUN/'SELECTION.json');source=read_json(verify(top['sources']))
        original={c['case_id']:c for c in source['cases']}
        for row in top['all_cases']:
            actual=original[row['case_id']]
            self.assertTrue(actual['primary_evaluation_pool'])
            self.assertFalse(actual['canonical_coordinate_match'])
            self.assertEqual(row['known_class'],actual['expected_class'])
        self.assertEqual(set(top['reuse_case_ids']),{r['case_id'] for r in top['all_cases'] if ('a0a3f2yly8' in r['case_id'] and any(r['case_id'].endswith('conditioned_Ca__seed-1_sample-'+i) for i in ('1','3'))) or ('a0acd6b9f2' in r['case_id'] and r['case_id'].endswith('sample-4'))})

    def test_corrupted_real_shard_rejected(self):
        p=RUN/'INPUTS_shard_0.json';data=read_json(p);top=read_json(verify(data['selection']))
        data['cases'].pop()
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'corrupted_missing_source.json';write_new(out,data)
            with self.assertRaises(InvalidArtifact):transfer.validate_selection(out,verify(top['calibration']))

    def test_real_comparator_join_keeps_all_populations(self):
        from union_adaptive_transfer_compare import add_ledger
        base=read_json(ROOT/'workspaces/nikasha_recovery_20260922/proposal_comparison.json')
        native=read_json(ROOT/'workspaces/adaptive_minimal_pool_20260923/COMPARISON_recovered_v1.json')
        standalone=read_json(ROOT/'workspaces/standalone_xtb_transfer_20260923/run_v1/COMPARISON_v1.json')
        for kind,count in (('rows',225),('pools',75),('triples',100)):
            joined=add_ledger(base[kind],native[kind],kind,{'minimal_recovered':'native_minimal_recovered'})
            joined=add_ledger(joined,standalone[kind],kind,{'standalone_static':'standalone_static'})
            self.assertEqual(len(joined),count)
            for before,after in zip(base[kind],joined):
                self.assertEqual(before['methods']['DFT'],after['methods']['DFT'])
                self.assertEqual(before['methods']['context_composite'],after['methods']['context_composite'])
        corrupted=[dict(r) for r in standalone['rows']];corrupted[0]['expected_class']='corrupted_label'
        with self.assertRaises(InvalidArtifact):add_ledger(base['rows'],corrupted,'rows',{'standalone_static':'standalone_static'})

if __name__=='__main__':unittest.main()
