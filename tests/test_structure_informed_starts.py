"""Real pinned structures and explicit corruption only; no scientific mock energies."""
import copy
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import structure_informed_starts as starts
from affordable_common import InvalidArtifact, read_json, verify, write_new
from mace_site_kinematics import Kinematics

PREP=ROOT/'workspaces/structure_informed_starts_20260922/prepared_v3/PREPARATION.json'

class ActualStructures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=read_json(PREP)
        cls.data,cls.sources,cls.templates=starts.sources(verify(cls.p['inputs']))

    def test_exact_declared_denominators_and_no_execution(self):
        p=self.p
        self.assertEqual((p['case_denominator'],p['start_denominator'],p['available_starts']),(8,16,4))
        self.assertEqual((p['maximum_optimizer_starts'],p['maximum_new_GFN2_singlepoints']),(8,32))
        self.assertEqual(p['new_molecular_calls'],0)
        self.assertEqual([r['case_id'] for r in p['rows']],[r['case_id'] for r in self.data['cases']])
        verify(p['implementation'])

    def test_actual_transplants_replay_with_same_paired_graph(self):
        for source,row in zip(self.sources,self.p['rows']):
            ca,la=source['tasks']
            self.assertEqual(ca['active_indices'],la['active_indices'])
            self.assertEqual(read_json(verify(ca['mapping'])),read_json(verify(la['mapping'])))
            for old in row['starts']:
                got=starts.transplant(source,self.templates[old['template_case_id']],old['template_class'])
                # JSON converts descriptive quartet tuples to lists.
                import json
                self.assertEqual(json.loads(json.dumps(got)),old)
                if old['status']=='start_available':
                    self.assertLessEqual(old['dihedral_replay_max_error_radian'],1e-8)
                    self.assertTrue(old['geometry_checks']['mapping_checks']['pass'])
                    self.assertLessEqual(old['maximum_heavy_displacement_A'],.8+1e-7)

    def test_true_source_quartet_can_extend_template_local_map(self):
        src=self.templates['1H4I'];kin=Kinematics(read_json(verify(src['tasks'][0]['mapping']))['context'])
        with self.assertRaises(InvalidArtifact):starts.indices(src,kin,'catalytic_aspartate',('N','CA','CB','CG'))
        coords,provenance=starts.template_quartet(src,kin,'catalytic_aspartate',('N','CA','CB','CG'))
        self.assertEqual(coords.shape,(4,3));verify(provenance['source'])
        self.assertTrue(np.isfinite(starts.dihedral(coords)))

    def test_absent_asp_is_not_replaced_or_charge_changed(self):
        case=next(r for r in self.p['rows'] if r['case_id']=='4MAE')
        ca=next(s for s in case['starts'] if s['template_class']=='Ca')
        self.assertEqual(ca['status'],'start_available')
        self.assertEqual(len(ca['mapped']),2);self.assertEqual(len(ca['unmapped']),2)
        self.assertTrue(all(x['role']=='extra_acidic_ligand_homolog' and x['delta_radian']==0 for x in ca['unmapped']))
        self.assertEqual(case['tasks'][1]['charge']-case['tasks'][0]['charge'],1)

    def test_difficult_sources_remain_exact_outside_domain(self):
        for row in self.p['rows'][4:]:
            self.assertTrue(all(s['status']=='unsupported' for s in row['starts']))
            self.assertTrue(all('domain' in s['reason'] for s in row['starts']))
            self.assertTrue(all(s['maximum_heavy_displacement_A']>.8 for s in row['starts']))
        q9=self.p['rows'][2]['starts'][0]
        self.assertIn('heavy displacement',q9['reason']) # Full-source atom fix preserves genuine rejection.

    def test_actual_carboxylate_equivalence_has_state_gate(self):
        src=self.templates['1H4I'];template=self.templates['4MAE']
        kin=Kinematics(read_json(verify(src['tasks'][0]['mapping']))['context'])
        tk=Kinematics(read_json(verify(template['tasks'][0]['mapping']))['context'])
        self.assertEqual(starts.terminal_period(src,template,kin,tk,'anchor_glutamate','GLU','chi3'),np.pi)
        corrupt=copy.deepcopy(template);corrupt['fragments']['anchor_glutamate']['formal_charge']=0
        self.assertEqual(starts.terminal_period(src,corrupt,kin,tk,'anchor_glutamate','GLU','chi3'),2*np.pi)
        self.assertEqual(starts.terminal_period(src,template,kin,tk,'anchor_asparagine','ASN','chi2'),2*np.pi)

    def test_corrupted_real_role_is_explicit_failure(self):
        corrupt=copy.deepcopy(self.sources[0]);corrupt['roles']['anchor_glutamate']['resnum']=-12345
        r=starts.transplant(corrupt,self.templates['4MAE'],'La')
        self.assertEqual(r['status'],'unsupported');self.assertIn('canonical role',r['reason'])

    def test_eight_actual_native_searches_return_existing_basins(self):
        root=ROOT/'workspaces/structure_informed_starts_20260922/searches_v1'
        manifest=read_json(root/'manifest.json');self.assertEqual(len(manifest['tasks']),8)
        for t in manifest['tasks']:
            new=read_json(root/'proposals'/t['task_id']/'result.json')
            old=read_json(Path(t['original_manifest']['path']).parent/'proposals'/t['original_task_id']/'result.json')
            self.assertEqual(new['status'],'proposal_available')
            self.assertLessEqual(new['proposal']['MACE_eV']-new['start']['MACE_eV'],1e-7)
            self.assertLess(abs(new['proposal']['MACE_eV']-old['proposal']['MACE_eV']),3e-10)
            self.assertLess(max(abs(a-b) for a,b in zip(new['proposal']['full_q'],old['proposal']['full_q'])),3e-5)

    def test_exact_repeats_preserve_state_recipe_and_no_restart(self):
        path=ROOT/'workspaces/structure_informed_starts_20260922/numerical_v2/manifest.json'
        self.assertEqual(starts.numerical_validate(path)['tasks'],8)
        m=read_json(path)
        for t in m['tasks']:
            self.assertEqual(verify(t['xyz']).read_bytes(),verify(t['original_task']['xyz']).read_bytes())
            text=verify(t['input']).read_text()
            self.assertEqual(text,verify(t['original_task']['input']).read_text())
            self.assertIn('UseXTBMixer true',text);self.assertIn('NoAutostart',text)
            self.assertNotIn('Convergence Tight',text);self.assertIsNone(t['restart'])

    def test_real_repeats_diagnose_coordinate_sensitivity_without_reclassification(self):
        result=read_json(ROOT/'workspaces/structure_informed_starts_20260922/numerical_v2/collection_1210115.json')
        self.assertEqual(result['complete'],8);self.assertIsNone(result['full_Ca_La_score'])
        for geometry in ('old_adaptive','new_template'):
            for medium in ('vacuum','alpb'):
                rows=[r for r in result['rows'] if r['geometry']==geometry and r['medium']==medium]
                self.assertEqual(rows[0]['energy_hartree'],rows[1]['energy_hartree'])
                for row in rows:
                    details=starts.scf_details(verify(row['actual']['output']).read_text())
                    self.assertEqual(details,row['details']);self.assertEqual(details['electrons'],438)
                    self.assertAlmostEqual(sum(details['printed_orbital_occupations']),438,places=3)
        delta=(result['La_transfer_kcal_mol']['new_template']['primary_repeat_1']-
               result['La_transfer_kcal_mol']['old_adaptive']['primary_repeat_1'])
        self.assertLess(delta,-4.8)

    def test_absent_convergence_diagnostic_is_not_a_pass(self):
        r=read_json(ROOT/'workspaces/structure_informed_starts_20260922/numerical_v2/collection_1210115.json')['rows'][0]
        text=verify(r['actual']['output']).read_text()
        corrupted='\n'.join(line for line in text.splitlines() if 'Last MAX-Density change' not in line)
        with self.assertRaises(InvalidArtifact):starts.scf_details(corrupted)

    def test_real_comparison_keeps_unadmitted_sources_unavailable(self):
        r=read_json(ROOT/'workspaces/structure_informed_starts_20260922/COMPARISON_v1.json')
        self.assertEqual((r['denominator'],r['available']),(8,4));self.assertIsNone(r['new_calibration'])
        for row in r['rows'][4:]:
            self.assertEqual(row['branch_status'],'unavailable')
            self.assertIsNone(row['variants']['operational']['alternative_start_R'])
        q9=next(row for row in r['rows'] if row['case_id']=='q9z4j7-pqq-la_model')['variants']['operational']
        self.assertEqual(q9['alternative_start_adaptive_band_transfer'],'inconclusive')
        self.assertLess(abs(q9['delta_R_model_kcal_mol']),1e-5) # No padding a real boundary result.

if __name__=='__main__':unittest.main()
