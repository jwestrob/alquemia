"""Real source probes, unchanged solver protocol and actual response components."""
from pathlib import Path
import sys,unittest,copy,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,write_new
from mace_responsive_charges import probe_inventory,MODEL
from mace_responsive_solvent import experiment_specs,MODEL as GB_MODEL
from mace_omol_solvent import MODEL as OLD_GB_MODEL
W=ROOT/'workspaces/mace_omol_20260917'

class ResponsiveComponents(unittest.TestCase):
    def test_real_combined_probes_preserve_each_source_order(self):
        c=read_json(W/'normalized_charge_v1/manifest.json');s=read_json(W/'explicit_field_short_v1/manifest.json')
        for old in c['tasks']:
            field=next(t for t in s['tasks'] if t['task_id']==old['task_id']);pts,groups=probe_inventory(old,field)
            np.testing.assert_array_equal(pts[groups['exterior_indices']],np.loadtxt(verify(old['points']),skiprows=1))
            np.testing.assert_array_equal(pts[groups['environment_indices']],np.loadtxt(verify(field['points']),skiprows=1))
            self.assertEqual(groups['environment_weights'],field['weights'])
        self.assertEqual(MODEL['scheme'],'native_ORCA6.1.1_CHELPG')
        self.assertFalse(MODEL['fit_charge_renormalization'])

    def test_exact_solver_inventory_and_unchanged_physical_settings(self):
        specs=experiment_specs();self.assertEqual(len(specs),44)
        self.assertTrue(all(s[2]!='environment' for s in specs))
        self.assertEqual(sum(s[3]=='identity' for s in specs),8)
        self.assertEqual(sum(s[3]=='reference' for s in specs),2)
        allowed={'vacuum_parent','charge_representation','energy_definition'}
        for k,v in OLD_GB_MODEL.items():
            if k not in allowed:self.assertEqual(GB_MODEL[k],v)

    def test_real_embedded_charge_manifest_and_missing_receipts(self):
        from mace_responsive_charges import validate,report
        m=W/'responsive_charges_v3/manifest.json'
        if not m.exists():self.skipTest('actual charge preparation unavailable')
        self.assertEqual(validate(m)['tasks'],8)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'manifest.json';write_new(p,read_json(m))
            r=report(p,Path(tmp)/'report')
            self.assertEqual(r['status'],'incomplete');self.assertFalse(r['coupling_gate_pass'])
            self.assertFalse(r['projection_gate_pass'])

    def test_actual_worker_roundoff_and_corrupted_cap_anchor(self):
        from mace_responsive_charges import same_projection
        from mace_omol_charges import projection
        from affordable_common import xyz
        t=read_json(W/'normalized_charge_v1/manifest.json')['tasks'][0]
        p=read_json(verify(t['projection']));self.assertTrue(same_projection(p,p))
        import json
        for line in (W/'projection_replay_debug_1200985.log').read_text().splitlines():
            row=json.loads(line)
            for cap in row['changed_caps']:
                a={'weights':p['weights'],'caps':[cap['actual']]};b={'weights':p['weights'],'caps':[cap['saved']]}
                # Only the actual recorded cap comparison is tested; no scientific geometry/energy is invented.
                self.assertTrue(same_projection(a,b))
                b['caps'][0]['retained']='corrupted_real_source_anchor';self.assertFalse(same_projection(a,b))

    def test_corrupted_real_density_reference_rejected(self):
        from mace_responsive_charges import validate
        m=W/'responsive_charges_v3/manifest.json'
        if not m.exists():self.skipTest('actual charge preparation unavailable')
        data=read_json(m);data['tasks'][0]['files']['gbw']=data['tasks'][1]['files']['gbw']
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'manifest.json';write_new(p,data)
            with self.assertRaisesRegex(InvalidArtifact,'wavefunction/density changed'):validate(p)

    @unittest.skipUnless((W/'responsive_charge_report_v2/result.json').exists(),'actual embedded utility outputs unavailable')
    def test_actual_charge_representation(self):
        r=read_json(W/'responsive_charge_report_v2/result.json');self.assertEqual(r['status'],'complete')
        self.assertEqual(len(r['coupling_checks']),5)
        for row in r['rows'].values():
            self.assertTrue(row['projection_conservation_pass']);self.assertAlmostEqual(sum(row['charge_e']),sum(row['projected_charge_e']),places=9)
            for k in ('exact','fitted','projected'):self.assertTrue(np.isfinite(row['direct_coupling_kcal_mol'][k]))

    def test_prepared_responsive_GB_changes_only_QM_distribution(self):
        from mace_responsive_solvent import validate
        path=W/'responsive_GB_v1/manifest.json'
        if not path.exists():self.skipTest('actual responsive GB preparation unavailable')
        m=read_json(path);old=read_json(W/'full_boundary_GB_v1/manifest.json')
        self.assertEqual(validate(path)['tasks'],44)
        for name,pin in m['states'].items():
            a=read_json(verify(pin));b=read_json(verify(old['states'][name]))
            self.assertEqual({k:v for k,v in a.items() if k!='endpoints'},{k:v for k,v in b.items() if k!='endpoints'})
            for metal,e in a['endpoints'].items():
                excluded={'QM_charges_e','full_charges_e','actual_charge_receipt'}
                self.assertEqual({k:v for k,v in e.items() if k not in excluded},{k:v for k,v in b['endpoints'][metal].items() if k not in excluded})
                np.testing.assert_array_equal(np.array(e['QM_charges_e'])+a['environment_charges_e'],e['full_charges_e'])
                self.assertNotEqual(e['actual_charge_receipt'],b['endpoints'][metal]['actual_charge_receipt'])

    @unittest.skipUnless((W/'responsive_GB_result_v1.json').exists(),'actual responsive GB outputs unavailable')
    def test_actual_hybrid_component_accounting(self):
        r=read_json(W/'responsive_GB_result_v1.json');self.assertEqual(r['status'],'complete');self.assertTrue(r['numerical_gate_pass'])
        self.assertEqual(len(r['contrasts']),4)
        for c in r['cases'].values():
            self.assertAlmostEqual(c['hybrid_R_kcal_scale']-c['frozen_hybrid_R_kcal_scale'],c['response_R_kcal_mol']+c['GB_update_R_kcal_mol'],places=7)
        for k in ('reference','calibrated_class','combined_gradient','relaxation_correction'):self.assertIsNone(r[k])

if __name__=='__main__':unittest.main()
