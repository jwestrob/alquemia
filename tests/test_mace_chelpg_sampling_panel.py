"""Real archived NoIter outputs and pinned scientific input guards."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,write_new,verify,xyz
from mace_chelpg_sampling_panel import parse_sampling,SETTINGS,validate,source_inventory,population_row
W=ROOT/'workspaces/mace_omol_20260917'

class SamplingPanel(unittest.TestCase):
    def test_actual_default_charge_parser_and_wrong_setting(self):
        m=read_json(W/'chelpg_sampling_qualification_v1/manifest.json')
        pop=read_json(W/'chelpg_sampling_population_recollection_v2/result.json')
        for t in m['tasks']:
            text=Path(t['output_path']).read_text().split('CHELPG CHARGES GENERATION',1)[1].split('CHELPG charges calculated...',1)[0]+'CHELPG charges calculated...'
            elements=[a[0] for a in xyz(verify(t['xyz']))]
            q=parse_sampling(text,elements,t['charge'],SETTINGS['original'])
            np.testing.assert_array_equal(q,pop['rows'][t['task_id']]['charge_e'])
            with self.assertRaisesRegex(InvalidArtifact,'sampling control'):parse_sampling(text,elements,t['charge'],SETTINGS['finer'])

    def test_real_panel_and_exact_qualification_reuse(self):
        p=W/'chelpg_sampling_panel_v1/manifest.json';r=validate(p)
        self.assertEqual(r['tasks'],38)
        m=read_json(p);self.assertEqual({r['task_id'] for r in m['reused']},{'GGR_2FVY_Ca_original','GGR_2FVY_La_original'})
        self.assertEqual(m['planned_SCF_optimizations'],0)

    def test_corrupted_real_input_cache_rejected(self):
        m=read_json(W/'chelpg_sampling_panel_v1/manifest.json');m['tasks'][0]['sampling']['grid_A']=.4
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_panel.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'changed sampling'):validate(p)

    def test_old_and_new_source_schemas_keep_same_physical_projection(self):
        c=read_json(W/'chelpg_sampling_panel_config_v1/config.json');sources,_,_=source_inventory(c)
        self.assertEqual(len(sources),10)
        for source in sources.values():
            proj=read_json(verify(source['observation_task']['projection']))
            q=np.array(source['density_row']['charge_e']);expected=source['density_row']['projected_charge_e']
            np.testing.assert_allclose(np.array(proj['weights'])@q,expected,rtol=0,atol=1e-12)
            self.assertTrue(source['exterior_indices'])

    def test_actual_native_variant_outputs(self):
        mp=W/'chelpg_sampling_panel_v1/manifest.json';m=read_json(mp)
        tasks=[t for t in m['tasks'] if t['case_id']=='GGR_extended' and t['metal']=='Ca']
        if not all(Path(t['output_path']+'.execution.json').exists() for t in tasks):self.skipTest('native sampling variants not yet available')
        for t in tasks:
            r=population_row(t,mp)
            self.assertEqual(r['status'],'native_NoIter_population_completed')
            self.assertTrue(r['preliminary_identity_pass'])
            self.assertFalse(r['SCF_optimization_performed'])
            self.assertLessEqual(abs(sum(r['charge_e'])-t['charge']),5e-5)

    def test_actual_panel_identity_and_sensitivity_algebra(self):
        p=W/'chelpg_sampling_report_v1/result.json'
        if not p.exists():self.skipTest('actual 40-result comparison not yet available')
        r=read_json(p);self.assertTrue(r['complete']);self.assertTrue(r['all_density_identities_pass'])
        self.assertEqual(len(r['rows']),40);self.assertEqual(len(r['pair_quality']),20)
        for case,settings in r['paired_source_self_R_kcal'].items():
            for setting,value in settings.items():
                ca=r['rows'][case+'_Ca_'+setting]['source_self_kcal'];la=r['rows'][case+'_La_'+setting]['source_self_kcal']
                self.assertEqual(ca-la,value)
        for item in r['sensitivity_screens']:
            if item['primary']:self.assertEqual(item['pass_'],abs(item['change_kcal'])<=.5)
        self.assertIsNone(r['new_score']);self.assertIsNone(r['biological_classification'])

if __name__=='__main__':unittest.main()
