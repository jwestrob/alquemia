"""Field algebra on real charges/coordinates; actual native integration when present."""
from pathlib import Path
import json
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,HA_TO_KCAL
from mace_qm_field import validate,point_field,offset_points,derivative,diagonal_response,accepted
from density_embedding import parse_potential

W=ROOT/'workspaces/mace_omol_20260917/qm_electric_field_v1'


class RealFieldInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('real field preparation unavailable')
        cls.manifest=validate(W/'manifest.json')

    def test_charge_field_sign_units_and_derivative_on_real_geometry(self):
        # Analytical consistency only; point-charge potentials are not substitutes
        # for native quantum outputs or evidence of field quality.
        task=self.manifest['tasks'][0];data=read_json(verify(task['probes']))
        probes=np.array(data['positions_bohr']);src=np.array(data['QM_positions_bohr']);q=np.array(data['charge_e'])
        idx=np.argsort(data['nearest_QM_or_cap_A'])[:12];sample=probes[idx]
        pp=offset_points(sample)
        potentials=np.sum(q[None,:]/np.linalg.norm(pp[:,None,:]-src[None,:,:],axis=2),axis=1)
        finite=derivative(potentials,len(sample));analytic=point_field(q,src,sample)
        np.testing.assert_allclose(finite[1],analytic,atol=1e-7,rtol=1e-5)
        alpha=np.array(data['alpha_bohr3'])[idx]
        expected=-.5*sum(a*np.dot(e,e) for a,e in zip(alpha,analytic))*HA_TO_KCAL
        self.assertAlmostEqual(diagonal_response(analytic,alpha),expected,places=9)

    def test_actual_pair_states_and_all_environment_probes(self):
        for task in self.manifest['tasks']:
            data=read_json(verify(task['probes']));state=read_json(verify(task['state']))
            expected=[a['id'] for a in state['physical_atoms'] if a['id'] not in state['projection_support_ids']]
            self.assertEqual(data['physical_ids'],expected)
            self.assertEqual(len(np.loadtxt(verify(task['points']),skiprows=1)),12*len(expected))

    def test_changed_field_step_invalidates_manifest(self):
        damaged=json.loads(json.dumps(self.manifest));damaged['steps_bohr'][0]=.002
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'damaged_real_manifest.json';p.write_text(json.dumps(damaged))
            with self.assertRaisesRegex(InvalidArtifact,'undeclared field configuration'):validate(p)

    def test_native_fields_only_from_actual_receipts(self):
        reports=list(W.glob('report_job_*/result.json'))
        if not reports:raise unittest.SkipTest('native utility outputs/report not yet available')
        report=read_json(reports[-1]);self.assertTrue(report['complete'])
        for task in self.manifest['tasks']:
            r=report['rows'][task['task_id']];receipt=read_json(verify(r['execution_receipt']))
            self.assertIsNotNone(accepted(Path(r['execution_receipt']['path']).parent,task,W/'manifest.json'))
            data=read_json(verify(task['probes']));p=np.loadtxt(verify(task['points']),skiprows=1)
            v=parse_potential(verify(receipt['potential']),p);fields=derivative(v,len(data['physical_ids']))
            with np.load(verify(r['arrays'])) as saved:
                np.testing.assert_array_equal(saved['exact'],fields[1])
        self.assertIsNone(report['environment_correction'])
        self.assertIsNone(report['numerical_score'])


if __name__=='__main__':unittest.main()
