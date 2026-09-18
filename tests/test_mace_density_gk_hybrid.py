"""Pinned source/algebra tests, with scientific integration only from real runs."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,read_json,verify
from mace_density_gk_hybrid import derive,key,validate,accepted,parse,load_terms
from mace_native_field_input import derive as prior_derive

BASE=ROOT/'workspaces/mace_omol_20260917'
W=BASE/'density_gk_hybrid_v1'
S=BASE/'density_gk_hybrid_software_v1'


class NativeSource(unittest.TestCase):
    def test_unrounded_residual_diagnostic_preserves_solver_equations(self):
        if not (S/'receipt.json').exists():self.skipTest('real native build unavailable')
        r=read_json(S/'receipt.json');text=verify(r['original_source']).read_text()
        original,derived=derive(text);old_original,old_derived=prior_derive(text)
        self.assertEqual(original,old_original)
        restored=derived.replace('subroutine alquemia_induce0c (source_fields,eps_final)',
                                 'subroutine alquemia_induce0c (source_fields)',1)
        restored=restored.replace('real*8 source_fields(3,n,4),eps_final','real*8 source_fields(3,n,4)',1)
        restored=restored.replace('      eps_final = eps\n','',1)
        self.assertEqual(restored,old_derived)
        self.assertEqual(derived,verify(r['derived_routine']).read_text())


class RealPreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('real hybrid preparation unavailable')
        cls.m=validate(W/'manifest.json')

    def test_reused_quantum_and_short_terms_are_actual_compatible_receipts(self):
        boundary=read_json(verify(self.m['sources']['boundary_manifest']))
        terms=load_terms(verify(self.m['sources']['short_report']),boundary)
        self.assertEqual(terms,read_json(verify(self.m['terms'])))
        self.assertEqual(len(terms),8)
        # Independent archived Hartree contrast, unit conversion exactly once.
        report=read_json(verify(self.m['sources']['short_report']))
        for case in boundary['cases']:
            ca,la=terms[case+'_Ca'],terms[case+'_La']
            old=report['cases'][case]['endpoints']
            actual=(ca['DFT_vacuum_hartree']-la['DFT_vacuum_hartree'])*HA_TO_KCAL
            expected=old['Ca']['DFT_vacuum_kcal_mol']-old['La']['DFT_vacuum_kcal_mol']
            self.assertAlmostEqual(actual,expected,places=7)

    def test_cache_changes_for_actual_source_cavity_rotation_and_solver(self):
        t=self.m['response_tasks'][0];original=key(t,self.m)
        changes=[('poleps',1e-8),('metal_GK_diameter','3.65'),('state','La')]
        for field,value in changes:
            changed=copy.deepcopy(t);changed[field]=value
            self.assertNotEqual(key(changed,self.m),original)
        changed=copy.deepcopy(t);changed['density_arrays']['sha256']='corrupted-real-pin'
        self.assertNotEqual(key(changed,self.m),original)
        changed=copy.deepcopy(t);changed['rotation'][0][0]+=.01
        self.assertNotEqual(key(changed,self.m),original)

    def test_dependency_failure_is_not_cached_scientific_success(self):
        t=self.m['response_tasks'][0]
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'attempt_0001';p.mkdir()
            # Explicit runner failure metadata, not fabricated scientific output.
            (p/'result.json').write_text(json.dumps(dict(cache_key=t['cache_key'],status='failed',receipt=None)))
            self.assertIsNone(accepted(t,d))
            damaged=read_json(p/'result.json');damaged['cache_key']='corrupted-real-key'
            (p/'result.json').write_text(json.dumps(damaged))
            with self.assertRaisesRegex(InvalidArtifact,'incompatible cached'):
                accepted(t,d)

    def test_real_supplied_fields_use_exact_density_and_native_reaction_difference(self):
        tasks=[t for t in self.m['response_tasks'] if t['variant']=='primary' and t['poleps']==1e-9]
        if any(not accepted(t,W/'tasks'/t['task_id']) for t in tasks):self.skipTest('real full-model response calls incomplete')
        for t in tasks:
            result,_=accepted(t,W/'tasks'/t['task_id'])
            deps=[read_json(verify(p)) for p in result['dependencies']]
            q,e=deps;n=q['parameters']['inventory'][0]
            Q=np.array([q['fields'][str(i)] for i in range(1,n+1)])
            E=np.array([e['fields'][str(i)] for i in range(1,n+1)])
            actual=np.loadtxt(verify(result['supplied_fields']),skiprows=1)[:,1:]
            if t['state']=='environment':np.testing.assert_array_equal(actual,E);continue
            meta=read_json(verify(t['boundary']));ids={p:i for i,p in enumerate(meta['physical_ids'])}
            probes=read_json(verify(t['density_task']['probes']));ix=[ids[p] for p in probes['physical_ids']]
            with np.load(verify(t['density_task']['field_arrays'])) as arr:field=arr['exact']/BOHR_TO_A**2
            np.testing.assert_allclose(actual[ix,1:4]-E[ix,1:4],field,rtol=1e-12,atol=1e-14)
            for vacuum,solvent in ((1,7),(4,10)):
                expected=E[ix,solvent:solvent+3]+field+(Q[ix,solvent:solvent+3]-Q[ix,vacuum:vacuum+3])-(E[ix,solvent:solvent+3]-E[ix,vacuum:vacuum+3])
                np.testing.assert_allclose(actual[ix,solvent:solvent+3],expected,rtol=1e-12,atol=1e-14)

    def test_actual_native_parsing_and_no_false_completion(self):
        reports=sorted(W.glob('collection_job_*.json'))
        if not reports:self.skipTest('actual hybrid collection unavailable')
        report=read_json(reports[-1]);self.assertEqual(report['complete'],all(v['status']=='computed_native_hybrid_component' for v in report['tasks'].values()))
        self.assertIsNone(report['reference']);self.assertIsNone(report['calibrated_class'])
        self.assertIsNone(report['combined_gradient']);self.assertFalse(report['baseline_changed'])
        for mode in ('static','identity','solve'):
            t=next(t for t in self.m['static_tasks']+self.m['response_tasks'] if t['mode']==mode)
            item=accepted(t,W/'tasks'/t['task_id'])
            if not item:continue
            result,_=item;text=verify(result['receipt']['log']).read_text()
            parsed=parse(text,t);self.assertEqual(parsed['energies_kcal_mol'],result['energies_kcal_mol'])
            with self.assertRaisesRegex(InvalidArtifact,'incomplete/nonconverged'):
                parse(text.replace('ALQUEMIA_DENSITY_GK_COMPLETE','CORRUPTED_COMPLETION'),t)


if __name__=='__main__':unittest.main()
