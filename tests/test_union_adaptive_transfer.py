"""Real frozen225 selectors; no generated scientific values."""
from pathlib import Path
import statistics,sys,tempfile,unittest
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


class FinalTransfer(unittest.TestCase):
    """Actual terminal results only; absent scientific outputs are explicit skips."""
    @classmethod
    def setUpClass(cls):
        path=RUN/'COMPARISON_v1.json'
        if not path.exists():raise unittest.SkipTest('actual final225 comparison has not completed')
        cls.result=read_json(path)

    def test_actual_full_replay_and_denominators(self):
        from union_adaptive_transfer_compare import compare
        r=self.result
        self.assertEqual([len(r[k]) for k in ('rows','pools','triples')],[225,75,100])
        self.assertEqual(len({x['case_id'] for x in r['rows']}),225)
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'actual_comparison_replay.json'
            compare(selection=verify(r['selection']),collections=[verify(p) for p in r['collections']],
                **{k:verify(p) for k,p in r['comparators'].items()},output=out)
            replay=read_json(out)
        for k in ('rows','pools','triples','bands','counts','matched','reference'):
            self.assertEqual(replay[k],r[k])

    def test_actual_raw_components_both_selection_policies(self):
        from affordable_common import HA_TO_KCAL
        from mace_hybrid import EV_TO_KCAL
        r=self.result;index={x['case_id']:x for x in r['rows']};seen=set()
        collections=[read_json(verify(p)) for p in r['collections']]
        selection=read_json(verify(r['selection']))
        collections.append(read_json(verify(selection['pilot_collection'])))
        for collection in collections:
            for case in collection['cases']:
                cid=case['case_id']
                if cid not in index:continue
                self.assertNotIn(cid,seen);seen.add(cid)
                row=index[cid]
                if case['pool']['status']!='available':
                    for name in ('union_adaptive','union_adaptive_mathematical'):
                        self.assertIsNone(row['methods'][name]['R'])
                    continue
                candidates=[x['id'] for x in case['candidates']];chosen={}
                for z in ('Ca','La'):
                    raw={q:case['matrix'][z][q]['components'] for q in candidates};origin=raw['origin']
                    work={q:(v['MACE_eV']-origin['MACE_eV'])*EV_TO_KCAL+
                        ((v['GFN2_ALPB_hartree']-origin['GFN2_ALPB_hartree'])-
                         (v['GFN2_vacuum_hartree']-origin['GFN2_vacuum_hartree']))*HA_TO_KCAL for q,v in raw.items()}
                    minimum=min(candidates,key=work.get)
                    chosen[z]={'mathematical':minimum,'operational':minimum if work[minimum]<-.1 else 'origin'}
                    for mode in ('operational','mathematical'):
                        self.assertEqual(chosen[z][mode],row['selected_candidates'][mode][z])
                for mode,name in (('operational','union_adaptive'),('mathematical','union_adaptive_mathematical')):
                    ca=case['matrix']['Ca'][chosen['Ca'][mode]]['components'];la=case['matrix']['La'][chosen['La'][mode]]['components']
                    direct=(ca['MACE_eV']-la['MACE_eV'])*EV_TO_KCAL+((ca['GFN2_ALPB_hartree']-ca['GFN2_vacuum_hartree'])-
                        (la['GFN2_ALPB_hartree']-la['GFN2_vacuum_hartree']))*HA_TO_KCAL
                    self.assertAlmostEqual(row['methods'][name]['R'],direct,places=7)
        self.assertEqual(len(seen),208)

    def test_actual_strict_missing_and_independent_aggregation(self):
        r=self.result;index={x['case_id']:x for x in r['rows']}
        for group in r['pools']+r['triples']:
            for name in ('union_adaptive','union_adaptive_mathematical'):
                actual=group['methods'][name]
                missing=sorted(x for x in group['members'] if index[x]['methods'][name]['R'] is None)
                self.assertEqual(sorted(actual['missing_members']),missing)
                if missing:
                    self.assertIsNone(actual['R']);self.assertEqual(actual['outcome'],'unavailable');continue
                values=[index[x]['methods'][name]['R'] for x in group['members']]
                if group.get('pool')=='balanced':
                    value=sum(statistics.median(index[x]['methods'][name]['R'] for x in group['members']
                        if index[x]['source_conditioning_metal']==z) for z in ('Ca','La'))/2
                else:value=statistics.median(values)
                self.assertEqual(actual['R'],value)

    def test_actual_failures_and_changed_decision_ids_remain_visible(self):
        r=self.result;rows=r['rows'];name='union_adaptive'
        self.assertEqual(r['counts']['all225'][name],{'denominator':225,'correct':205,'wrong':0,'inconclusive':1,'unavailable':19})
        self.assertEqual([x['case_id'] for x in rows if x['methods'][name]['outcome']=='inconclusive'],
            ['c5axv8-pqq-la_model__conditioned_La__seed-1_sample-3'])
        failed={x['case_id'] for x in rows if any(o and o['optimizer'] and not o['optimizer']['success'] for o in x['optimizers'].values())}
        self.assertEqual(failed,{'q92wy9-pqq-la_model__conditioned_Ca__seed-1_sample-2',
            'q60ar6-pqq-la_model__conditioned_La__seed-1_sample-0'})
        selection=read_json(verify(r['selection']))
        excluded={x['case_id'] for x in selection['all_cases'] if x['union_origin_status']!='complete'}
        self.assertEqual(len(excluded),17)
        self.assertEqual({x['case_id'] for x in rows if x['methods'][name]['R'] is None},excluded|failed)
        self.assertEqual({x['case_id'] for x in rows if x['methods']['context_composite']['outcome']=='wrong'
            and x['methods'][name]['outcome']=='correct'},
            {'a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1',
             'a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4'})

if __name__=='__main__':unittest.main()
