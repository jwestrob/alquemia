"""Pinned real observation inputs; explicitly corrupted configurations only."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,BOHR_TO_A
from mace_density_boundary_panel import validate as boundary
from mace_density_observations import validate,collect

BASE=ROOT/'workspaces/mace_omol_20260917'
OBS=BASE/'trial_gk_expansion_observations_v2/manifest.json'
ENV=BASE/'trial_gk_expansion_environment_boundary_v1/manifest.json'


class DensityObservationTests(unittest.TestCase):
    def test_real_zero_source_reference_does_not_invent_QM_charges(self):
        m=boundary(ENV);self.assertEqual(len(m['tasks']),4)
        for t in m['tasks']:
            self.assertIsNone(t['QM_actual_charge']);self.assertIsNone(t['source_charge_execution'])
            for i in t['frozen_indices']:self.assertEqual(t['charge_e'][i-1],0.)
        r=read_json(ENV.parent/'result.json');self.assertTrue(r['gates_pass'])

    def test_all_real_generating_field_charges_are_observed_at_actual_coordinates(self):
        m=validate(OBS);self.assertEqual(len(m['tasks']),4)
        for t in m['tasks']:
            p=read_json(verify(t['probes']));pc=np.loadtxt(verify(t['source_task']['pointcharges']),skiprows=1)
            observed=np.c_[p['old_generating_charges_e'],np.array(p['positions_bohr'])*BOHR_TO_A]
            observed=observed[observed[:,0]!=0]
            order=lambda a:a[np.lexsort(a[:,::-1].T)]
            np.testing.assert_allclose(order(pc),order(observed),atol=2e-12,rtol=0)
            self.assertEqual(p['exterior_start'],49*len(p['physical_ids']))

    def test_changed_real_charge_model_invalidates_query(self):
        m=read_json(OBS);m['model']['grid_A']=0.4
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_charge_model.json';p.write_text(json.dumps(m))
            with self.assertRaisesRegex(InvalidArtifact,'model/tolerances'):validate(p)

    def test_unexecuted_real_queries_remain_missing(self):
        m=read_json(OBS)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'copied_real_manifest_without_execution.json';p.write_text(json.dumps(m))
            r=collect(p,Path(d)/'report');self.assertEqual(r['status'],'incomplete')
            self.assertIsNone(r['full_hybrid_score']);self.assertFalse(r['numerical_pass'])
            self.assertTrue(all(x['status']=='unavailable' for x in r['rows'].values()))


if __name__=='__main__':unittest.main()
