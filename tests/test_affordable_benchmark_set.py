"""Integrity checks on the actual evidence/structure/task release; no synthetic science."""
from pathlib import Path
import sys,json,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,paired,InvalidArtifact,record
from affordable_benchmark_set import assemble
from run_orca_task_manifest import load_manifest_tasks

class BenchmarkSetIntegrity(unittest.TestCase):
    def setUp(self):
        self.base=ROOT/'workspaces/benchmark_set_20260915'
        if not (self.base/'structure_inventory_v4.json').exists():self.skipTest('pinned real structure inventory unavailable')

    def test_all_metals_preserved_in_source_inventory(self):
        rows={r['pdb_id']:r for r in read_json(self.base/'structure_inventory_v4.json')['records']}
        self.assertEqual([s['metal']['element'] for s in rows['6MI5']['metal_sites_first_model']],['Y']*3)
        self.assertEqual(sum(s['metal']['element']=='La' for s in rows['9VY7']['metal_sites_first_model']),5)
        self.assertIn('Mn',{s['metal']['element'] for s in rows['3CNA']['metal_sites_first_model']})
        self.assertEqual([s['metal']['element'] for s in rows['9B1U']['metal_sites_first_model']],['Na','Na'])

    def test_prepared_physical_pairs_and_output_locations(self):
        m,tasks=load_manifest_tasks(self.base/'ready_tasks_v4/manifest.json')
        self.assertEqual(len(tasks),10)
        for t in tasks:self.assertEqual(t['output'].parent,t['input'].parent)
        for c in m['cases']:
            source=read_json(verify(c['source_manifest']))
            endpoints=source['outputs']
            self.assertEqual(paired(verify(endpoints['La']['xyz']),verify(endpoints['Ca']['xyz']),endpoints['La']['charge'],endpoints['Ca']['charge'])['status'],'pass')
            self.assertLessEqual(c['source_heavy_max_displacement_A'],.002)
        hans=[c for c in m['cases'] if c['target_id']=='hans_lanm_wt']
        self.assertEqual(len(hans),3)
        self.assertEqual({c['biological_group'] for c in hans},{'lanmodulin_family'})
        self.assertTrue(all(not c['explicit_water_inventory'] for c in hans))
        alpha=[c for c in m['cases'] if c['target_id']=='ALACTA_BOVINE_STRONG_SITE']
        self.assertEqual(len(alpha),2)
        self.assertNotEqual(len(alpha[0]['explicit_water_inventory']),len(alpha[1]['explicit_water_inventory']))

    def test_changed_real_capture_rejected(self):
        path=self.base/'evidence/lanthanide/evidence_rows.json'
        row=read_json(path)[0]
        with tempfile.TemporaryDirectory() as td:
            base=Path(td);lane=base/'evidence/lane';lane.mkdir(parents=True)
            actual=path;copy=base/'explicitly_corrupted_real_capture.json';copy.write_bytes(actual.read_bytes()+b'\nCORRUPTED FOR INTEGRITY TEST\n')
            row['primary_sources']=[{'capture':dict(record(actual),path=str(copy))}]
            (lane/'evidence_rows.json').write_text(json.dumps([row]))
            with self.assertRaises(InvalidArtifact):assemble(base/'evidence',self.base/'structure_inventory_v4.json',base/'release')

    def test_release_preserves_exposure_and_unresolved_labels(self):
        manifest=self.base/'release_v1/benchmark_manifest.json'
        if not manifest.exists():self.skipTest('assembled real evidence ledger unavailable')
        rows={r['target_id']:r for r in read_json(manifest)['rows']}
        self.assertEqual(len(rows['GGR_1GLG']['existing_scores']),2)
        self.assertEqual(len(rows['AEQUORIN_1SL8_VECTOR']['existing_scores']),6)
        for target in ('CALBINDIN_N','CALBINDIN_C','H19_FAM1_XoxF5'):
            self.assertEqual(rows[target]['direction'],'unresolved')
        self.assertEqual(rows['PQQT_WT']['family_holdout_group'],rows['PQQT_K142D']['family_holdout_group'])
        self.assertEqual(rows['hans_lanm_wt']['biological_group'],rows['mex_lanm_wt']['biological_group'])
        self.assertFalse(any(row['new_calculation_performed'] for row in rows.values()))

    def test_archived_protonation_replay_preserves_actual_coordinates(self):
        for source in ('1f6s','6ip9'):
            old=self.base/f'prepared/alacta_{source}_v1/strong_site/amide_v3'
            new=self.base/f'prepared/alacta_{source}_v2/strong_site/amide_v3'
            for original in old.glob('*.xyz'):
                self.assertEqual(original.read_bytes(),(new/original.name).read_bytes())

    def test_unrun_real_manifest_collects_missing_scores_as_grouped_nulls(self):
        from affordable_benchmark import collect
        manifest=self.base/'ready_tasks_v4/manifest.json'
        if any(Path(t['output_path']).exists() for t in read_json(manifest)['tasks']):
            self.skipTest('the unrun fixture has since been executed')
        result=collect(manifest)
        self.assertEqual(result['status'],'incomplete')
        self.assertEqual(result['completed_endpoints'],0)
        groups={g['observation_group']:g for g in result['ordered_vectors']}
        self.assertEqual(groups['hans_lanm_wt']['S_vector_kcal_mol'],[None,None,None])
        self.assertEqual(groups['ALACTA_BOVINE_STRONG_SITE']['S_vector_kcal_mol'],[None,None])

if __name__=='__main__':unittest.main()
