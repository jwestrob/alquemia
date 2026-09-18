"""Real source-backed structural transfer; no fabricated scientific outputs."""
from pathlib import Path
import copy,json,sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,xyz,InvalidArtifact
from mace_omol_hybrid_transfer import normalized_core,prepared_states,validate,validate_quantum
BASE=ROOT/'workspaces/mace_omol_hybrid_transfer_20260918'
PREP=BASE/'prepared_v2/preparation.json'
OLD=ROOT/'workspaces/mace_omol_20260917/matched_H_prepared_v1'

@unittest.skipUnless((OLD/'preparation.json').exists(),'requires real archived normalized GGR cores')
class TransferTests(unittest.TestCase):
    def test_uniform_mapping_reproduces_both_archived_1GLG_cores(self):
        for representation in ('extended','connected'):
            c=read_json(OLD/('GGR_'+representation)/'mapping.json')
            original=read_json(verify(c['source_mapping']))
            repair=read_json(verify(original['source_preparation']))
            ends,mapping,checks=normalized_core(repair,read_json(verify(c['normalized_global_preparation'])))
            self.assertTrue(checks['heavy_coordinates_unchanged'])
            for metal,e in ends.items():
                archive=xyz(verify(c['endpoints'][metal]['xyz']))
                self.assertEqual([a[0] for a in e['atoms']],[a[0] for a in archive])
                np.testing.assert_allclose([a[1:] for a in e['atoms']],[a[1:] for a in archive],atol=5e-10,rtol=0)
                self.assertEqual(e['charge'],c['endpoints'][metal]['charge'])

    @unittest.skipUnless(PREP.exists(),'requires real new transfer preparations')
    def test_transfer_counts_pairs_and_source_graphs(self):
        p,cases=prepared_states(PREP);self.assertEqual(len(cases),4)
        for name,c in cases.items():
            repair=read_json(verify(c['source_preparation']));w=read_json(verify(c['whole_preparation']))
            replay,mapping,checks=normalized_core(repair,w)
            self.assertEqual(mapping,c['mapping']);self.assertEqual(checks,c['normalization_replay'])
            expected=111 if name.endswith('connected') else 58
            for metal,e in replay.items():
                rows=xyz(verify(c['endpoints'][metal]['xyz']));self.assertEqual(len(rows),expected)
                np.testing.assert_allclose([a[1:] for a in rows],[a[1:] for a in e['atoms']],atol=5e-10,rtol=0)
            self.assertEqual(replay['Ca']['atoms'][1:],replay['La']['atoms'][1:])

    @unittest.skipUnless((BASE/'initial_v1/initial.json').exists(),'requires actual finite manifests')
    def test_finite_tasks_and_corrupted_real_charge_rejected(self):
        initial=read_json(BASE/'initial_v1/initial.json');total=0
        for ref in initial['MACE_manifests']:total+=validate(verify(ref))['tasks']
        self.assertEqual(total,156);self.assertEqual(validate_quantum(verify(initial['quantum_manifest']))['status'],'dry_run_pass')
        m=read_json(verify(initial['MACE_manifests'][0]));m['tasks'][0]['charge']+=2
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted_actual_charge.json';p.write_text(json.dumps(m))
            with self.assertRaises(InvalidArtifact):validate(p)

    @unittest.skipUnless((BASE/'minimum_v1/preparation.json').exists(),'requires actual transfer grids and fixed native points')
    def test_actual_predictions_and_native_source_coordinates(self):
        from mace_omol_hybrid_transfer_minimum import inputs
        ap=BASE/'assessment_v1/result.json';a,p,cases,selected=inputs(ap)
        top=read_json(BASE/'minimum_v1/preparation.json')
        self.assertEqual(len(selected),8);self.assertEqual(top['excluded'],{})
        self.assertEqual(validate(BASE/'minimum_v1/mace/manifest.json')['tasks'],16)
        self.assertEqual(validate_quantum(BASE/'minimum_v1/quantum/manifest.json')['status'],'dry_run_pass')
        for key,state in top['states'].items():
            name,metal=state['case_id'],state['metal'];c=cases[name];u=np.array(selected[key]['prediction']['displacement_A'])
            self.assertLessEqual(np.linalg.norm(u),.20+1e-10)
            positions=[]
            for kind,e in [('core',c['endpoints'][metal]),('full',p['fulls'][c['global_id']][metal])]:
                old=xyz(verify(e['xyz']));new=xyz(verify(state['endpoints'][kind]['xyz']));i=e['metal_index']
                self.assertEqual(old[:i]+old[i+1:],new[:i]+new[i+1:])
                np.testing.assert_allclose(np.array(new[i][1:])-old[i][1:],u,atol=5e-10,rtol=0)
                positions.append(new[i][1:])
            np.testing.assert_allclose(*positions,atol=5e-10,rtol=0)
        corrupted=copy.deepcopy(a);corrupted['rows']['GGR_2FW0_connected']['Ca']['prediction']['displacement_A'][0]+=.01
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'corrupted_actual_prediction.json';path.write_text(json.dumps(corrupted))
            with self.assertRaises(InvalidArtifact):inputs(path)

    @unittest.skipUnless((BASE/'report_v1/result.json').exists(),'requires all actual native transfer outputs')
    def test_actual_transfer_report_algebra_and_complete_denominator(self):
        from mace_omol_hybrid_transfer_minimum import report
        archived=read_json(BASE/'report_v1/result.json')
        with tempfile.TemporaryDirectory() as td:
            output=Path(td)/'replay'
            report(BASE/'minimum_v1/preparation.json',BASE/'minimum_v1/quantum/manifest.json',BASE/'minimum_v1/mace/manifest.json',output)
            actual=read_json(output/'result.json')
        for key in ('rows','scores','partition','contrasts'):
            self.assertEqual(actual[key],archived[key])
        self.assertEqual(len(actual['contrasts']),12)
        self.assertEqual(sum(c['reused_comparison'] for c in actual['contrasts']),4)
        self.assertEqual(actual['biological_groups'],2)
        reference=read_json(verify(actual['reference_report']))
        for name,s in actual['scores'].items():
            ca,la=[actual['rows'][name+'_'+metal] for metal in ('Ca','La')]
            self.assertEqual(s['actual_delta_R_kcal_scale'],ca['actual_energy_change_kcal_mol']-la['actual_energy_change_kcal_mol'])
            self.assertIsNone(s['calibrated_class'])
        for c in actual['contrasts']:
            y=(reference['scores']['GGR_'+c['GGR'].rsplit('_',1)[1]] if c['reused_comparison'] else actual['scores'][c['GGR']])
            x=reference['scores'][c['alpha']]
            self.assertEqual(c['actual_margin_kcal'],x['actual_R_kcal_scale']-y['actual_R_kcal_scale'])
            self.assertIsNone(c['qualified_margin_kcal'])
        self.assertFalse(actual['all_twelve_qualified_directions_pass'])
        self.assertIsNone(actual['aqueous_score'])
        self.assertIsNone(actual['entropy_correction'])

if __name__=='__main__':unittest.main()
