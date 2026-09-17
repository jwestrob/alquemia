"""Real source mappings and executed fitted moments; no mock quantum outputs."""
from pathlib import Path
import json
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_distributed_source import validate,evaluate,accepted
from mace_qm_field import point_field
from density_embedding import parse_potential

W=ROOT/'workspaces/mace_omol_20260917/distributed_source_v1'


class SourceInputs(unittest.TestCase):
    workspace=W
    @classmethod
    def setUpClass(cls):
        if not (cls.workspace/'manifest.json').exists():raise unittest.SkipTest('real distributed-source preparation unavailable')
        cls.manifest=validate(cls.workspace/'manifest.json')

    def test_exact_reduction_to_real_projected_monopoles(self):
        for t in self.manifest['tasks']:
            data=read_json(verify(t['probes']));q=data['projected_charge_e'];source=data['projected_positions_bohr']
            probes=np.array(data['positions_bohr'])[::37]
            phi,e=evaluate(source,probes,q,np.zeros((len(q),3)))
            np.testing.assert_allclose(e,point_field(q,source,probes),atol=1e-14,rtol=1e-12)
            direct=np.sum(np.array(q)[None,:]/np.linalg.norm(probes[:,None,:]-np.array(source)[None,:,:],axis=2),axis=1)
            np.testing.assert_allclose(phi,direct,atol=1e-14,rtol=1e-12)

    def test_real_source_identity_and_no_new_atoms(self):
        for t in self.manifest['tasks']:
            state=read_json(verify(t['state']))
            self.assertEqual(set(t['source_ids']),set(state['projection_support_ids']))
            self.assertEqual(len(t['source_ids']),len(set(t['source_ids'])))

    def test_changed_regularization_rejected(self):
        damaged=json.loads(json.dumps(self.manifest));damaged['model']['regularization_fraction']=1e-3
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'damaged_real_manifest.json';p.write_text(json.dumps(damaged))
            with self.assertRaisesRegex(InvalidArtifact,'configuration differs'):validate(p)


class SourceIntegration(unittest.TestCase):
    workspace=W
    @classmethod
    def setUpClass(cls):
        reports=list(cls.workspace.glob('report_job_*/result.json'))
        if not reports:raise unittest.SkipTest('actual fits/native validation not yet available')
        cls.result=read_json(reports[-1]);cls.manifest=validate(cls.workspace/'manifest.json')
        if not cls.result['complete']:raise unittest.SkipTest('actual scientific attempts incomplete')

    def test_native_receipts_and_no_scoring_substitution(self):
        for t in self.manifest['tasks']:
            row=self.result['rows'][t['task_id']];rp=verify(row['execution_receipt'])
            self.assertIsNotNone(accepted(rp.parent,t,self.workspace/'manifest.json'))
            fit=read_json(verify(row['fit']));self.assertLessEqual(abs(fit['charge_residual_e']),1e-9)
        self.assertIsNone(self.result['numerical_score'])
        self.assertIsNone(self.result['environment_correction'])

    def test_actual_fitted_dipole_field_sign_and_tensor_order(self):
        for t in self.manifest['tasks']:
            data=read_json(verify(t['probes']));fit=read_json(verify(self.result['rows'][t['task_id']]['fit']))
            probes=np.array(data['positions_bohr'])[np.argsort(data['nearest_QM_or_cap_A'])[:12]]
            q,mu=fit['charge_e'],fit['dipole_e_bohr'];source=data['projected_positions_bohr']
            _,analytic=evaluate(source,probes,q,mu);finite=np.empty_like(analytic);h=1e-4
            for axis in range(3):
                shift=np.eye(3)[axis]*h
                plus,_=evaluate(source,probes+shift,q,mu);minus,_=evaluate(source,probes-shift,q,mu)
                finite[:,axis]=-(plus-minus)/(2*h)
            np.testing.assert_allclose(analytic,finite,atol=1e-8,rtol=1e-6)


class SourceSamplingInputs(SourceInputs):
    workspace=W.with_name('distributed_source_v2_recovery_v1')

    def test_additional_training_uses_actual_native_potentials(self):
        for t in self.manifest['tasks']:
            data=read_json(verify(t['probes']));n=len(data['positions_bohr'])
            points=np.loadtxt(verify(t['additional_training_points']),skiprows=1)
            actual=parse_potential(verify(t['additional_training_potential']),points)
            with np.load(verify(t['training'])) as training:
                self.assertEqual(training['potential'].shape,(2*n,))
                np.testing.assert_array_equal(training['potential'][n:],actual[:n])
                np.testing.assert_array_equal(training['positions_bohr'][n:],points[:n])
            validation=np.loadtxt(verify(t['points']),skiprows=1)[:n]
            self.assertFalse(np.allclose(validation,points[:n],atol=1e-8,rtol=0))


class SourceSamplingIntegration(SourceIntegration):
    workspace=W.with_name('distributed_source_v2_recovery_v1')


if __name__=='__main__':unittest.main()
