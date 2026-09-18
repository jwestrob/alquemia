"""Saved real embedded states; no substitute scientific backend or output."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,energy,read_json,verify
from density_embedding import parse_potential
from mace_density_gk_boundary import validate as validate_boundary,comparable
from mace_density_gk_hybrid import attach_environment_reuse,load_trial_terms,pipeline_preflight,TRIAL_PROTOCOL

W=ROOT/'workspaces/mace_omol_20260917/trial_density_gk_source_audit_v1.json'


class ActualTrialSource(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not W.exists():raise unittest.SkipTest('actual trial-density source audit absent')
        cls.r=read_json(W)

    def test_old_generating_field_removed_before_current_environment(self):
        for row in self.r['rows'].values():
            task=row['source_task'];qt=row['quantum_task'];groups=read_json(verify(task['probe_groups']))
            pc=np.loadtxt(verify(qt['pointcharges']),skiprows=1)
            points=np.loadtxt(verify(task['points']),skiprows=1)
            phi=parse_potential(verify(row['existing_center_potential']),points)
            ids=groups['environment_indices']
            np.testing.assert_allclose(points[ids]*BOHR_TO_A,pc[:,1:],atol=1e-11,rtol=0)
            direct=np.dot(pc[:,0],phi[ids]);embedded=energy(verify(task['source_output']))
            self.assertAlmostEqual(embedded-direct,row['intrinsic_core_hartree'],places=10)
            self.assertNotEqual(embedded,row['intrinsic_core_hartree'])
            self.assertGreaterEqual(row['intrinsic_polarization_cost_kcal'],-.05)

    def test_actual_wavefunctions_and_exact_nuclei_are_available_without_new_calls(self):
        self.assertEqual(len(self.r['rows']),8)
        for row in self.r['rows'].values():
            self.assertTrue(row['physical_state_exact']);self.assertTrue(row['core_coordinates_exact'])
            self.assertTrue(row['wavefunctions_available'])
            for key,pin in row['source_task']['files'].items():
                self.assertTrue(verify(pin).is_file())
                self.assertEqual(pin['sha256'],row['source_task']['source_wavefunctions'][key]['sha256'])
        self.assertIsNone(self.r['numerical_score'])
        self.assertTrue(all(self.r[k]==0 for k in ('new_DFT_calls','new_charge_fits','new_potential_utilities','new_native_solves')))

    def test_current_trial_boundary_keeps_all_environment_references_identical(self):
        folder=W.parent/'trial_density_gk_boundary_v1'
        if not (folder/'manifest.json').exists():self.skipTest('actual trial native preparation unavailable')
        m=validate_boundary(folder/'manifest.json');new=read_json(folder/'result.json')
        self.assertTrue(new['gates_pass']);self.assertEqual(new['actual_native_initializations'],16)
        old=read_json(W.parent/'density_gk_boundary_v1/result.json')
        for t in m['tasks']:
            if t['mode']!='zero_source':continue
            a=read_json(verify(new['tasks'][t['task_id']]['result']))
            b=read_json(verify(old['tasks'][t['task_id']]['result']))
            self.assertEqual(comparable(a),comparable(b))
        terms=load_trial_terms(W,m,W.parent/'explicit_field_short_result_v1.json')
        for tid,term in terms.items():
            self.assertEqual(term['intrinsic_core_hartree'],self.r['rows'][tid]['intrinsic_core_hartree'])
            self.assertNotEqual(term['intrinsic_core_hartree'],term['embedded_energy_hartree'])

    def test_reuse_only_for_unchanged_real_environment_inputs(self):
        # Prepared metadata copied only to exercise compatibility checks; this
        # is not a scientific output or an executed trial-density manifest.
        m=copy.deepcopy(read_json(verify(self.r['parent_manifest'])))
        m['protocol']=TRIAL_PROTOCOL;attach_environment_reuse(m,W)
        reused=[t for t in m['static_tasks']+m['response_tasks'] if 'reuse_component' in t]
        self.assertEqual(len(reused),37);self.assertTrue(all(t['state']=='environment' for t in reused))
        self.assertEqual((m['new_energy_calls'],m['new_field_queries'],m['new_response_solves']),(34,32,40))
        bad=copy.deepcopy(m);t=next(t for t in bad['static_tasks'] if t['state']=='environment')
        t['overrides']['sha256']='explicitly-corrupted-real-charge-file-pin'
        with self.assertRaisesRegex(InvalidArtifact,'scientific inputs changed'):
            attach_environment_reuse(bad,W)

    def test_frozen_dependent_pipeline_has_declared_inventory(self):
        p=W.parent/'trial_density_gk_pipeline_v1/config.json'
        if not p.exists():self.skipTest('actual dependent configuration unavailable')
        c=pipeline_preflight(p);self.assertEqual(len(c['task_ids']),111)
        self.assertEqual(c['reused_components'],37)
        with tempfile.TemporaryDirectory() as d:
            damaged=copy.deepcopy(c);damaged['new_calls']['energy']+=1
            q=Path(d)/'explicitly_corrupted_config.json';q.write_text(json.dumps(damaged))
            with self.assertRaisesRegex(InvalidArtifact,'inventory changed'):
                pipeline_preflight(q)


if __name__=='__main__':unittest.main()
