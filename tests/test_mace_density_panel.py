"""Frozen source-backed panel checks; all energies come from real executions."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import HA_TO_KCAL,InvalidArtifact,read_json,record,verify
from mace_density_panel import validate,inputs
from mace_density_gk_hybrid import collect,key,accepted,parse
BASE=ROOT/'workspaces/mace_omol_20260917'
W=BASE/'trial_gk_expansion_native_v1'

class Panel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('real prepared panel unavailable')
        cls.m=validate(W/'manifest.json')

    def test_actual_preparation_receipts_and_intrinsic_subtraction(self):
        m=self.m;c=read_json(verify(m['sources']['config']))
        data=inputs(**{k:verify(v) for k,v in c['inputs'].items()})
        terms=data[-1];self.assertEqual(terms,read_json(verify(m['terms'])))
        self.assertEqual(len(terms),4)
        for t in terms.values():
            self.assertEqual(t['intrinsic_core_kcal'],t['embedded_energy_hartree']*HA_TO_KCAL-t['old_field_coupling_kcal'])
            verify(t['source_receipt']);verify(t['source_output'])

    def test_same_physical_boundary_and_exact_paired_coordinates(self):
        for case in self.m['case_ids']:
            for variant in ('primary','rigid','radius_minus','radius_plus'):
                tasks=[t for t in self.m['static_tasks'] if t['case_id']==case and t['variant']==variant]
                self.assertEqual(len(tasks),3)
                xyz=[np.loadtxt(verify(t['xyz']),skiprows=1,usecols=(2,3,4)) for t in tasks]
                for v in xyz[1:]:np.testing.assert_array_equal(v,xyz[0])
                meta=[read_json(verify(t['boundary'])) for t in tasks]
                self.assertTrue(all(v==meta[0] for v in meta))

    def test_cache_invalidates_source_radius_and_convergence(self):
        m=self.m;t=m['response_tasks'][0]
        for field,value in [('poleps',1e-8),('metal_GK_diameter','3.8'),('state','La')]:
            damaged=copy.deepcopy(t);damaged[field]=value
            self.assertNotEqual(key(damaged,m),t['cache_key'])
        changed=copy.deepcopy(m);changed['sources']['density']['sha256']='corrupted-real-pin'
        self.assertNotEqual(key(t,changed),t['cache_key'])

    def test_missing_outputs_are_unavailable(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'real_manifest_without_results.json';p.write_text(json.dumps(self.m))
            r=collect(p,Path(d)/'collection.json')
            self.assertFalse(r['complete']);self.assertFalse(r['numerical_pass'])
            self.assertFalse(r['ordering_pass']);self.assertIsNone(r['calibrated_class'])
            self.assertTrue(all(v['status']=='unavailable' for v in r['tasks'].values()))

    def test_corrupted_actual_inventory_rejected(self):
        m=copy.deepcopy(self.m);m['static_tasks'].pop()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_manifest.json';p.write_text(json.dumps(m))
            with self.assertRaisesRegex(InvalidArtifact,'inventory'):validate(p)

    def test_legacy_trial_collector_reproduces_frozen_result(self):
        old=BASE/'trial_density_gk_hybrid_v1/collection_job_1201137.json';r=read_json(old)
        with tempfile.TemporaryDirectory() as d:
            actual=collect(verify(r['manifest']),Path(d)/'collection.json')
            # Four rigid dipole contractions differ by <=3.6e-15 between hosts.
            # Require exact scores/decisions and <=1e-12 only for those contractions.
            expected=copy.deepcopy(r)
            for label,v in actual['variants'].items():
                for case,row in v['cases'].items():
                    for metal,e in row['endpoints'].items():
                        a=e['direct_components_kcal']['dipole']
                        b=expected['variants'][label]['cases'][case]['endpoints'][metal]['direct_components_kcal']['dipole']
                        self.assertLessEqual(abs(a-b),1e-12)
                        e['direct_components_kcal']['dipole']=b
            self.assertEqual(actual,expected)

    def test_real_panel_execution_and_direct_algebra(self):
        reports=sorted(W.glob('collection_job_*.json'))
        if not reports:self.skipTest('new native panel has not actually run')
        r=read_json(reports[-1]);self.assertTrue(r['complete']);self.assertEqual(len(r['checks']),18)
        self.assertEqual((r['actual_energy_calls'],r['actual_field_queries'],r['actual_response_solves']),(27,24,30))
        self.assertIsNone(r['partition_pass']);self.assertIsNone(r['reference'])
        for c in r['variants']['primary']['cases'].values():
            self.assertAlmostEqual(c['R_kcal'],c['endpoints']['Ca']['total_kcal']-c['endpoints']['La']['total_kcal'],places=7)
        for mode in ('static','solve','identity'):
            t=next(t for t in self.m['static_tasks']+self.m['response_tasks'] if t['mode']==mode)
            item=accepted(t,W/'tasks'/t['task_id']);self.assertIsNotNone(item)
            out,_=item;text=verify(out['receipt']['log']).read_text()
            self.assertEqual(parse(text,t)['energies_kcal_mol'],out['energies_kcal_mol'])
            with self.assertRaisesRegex(InvalidArtifact,'incomplete/nonconverged'):
                parse(text.replace('ALQUEMIA_DENSITY_GK_COMPLETE','CORRUPTED_COMPLETION'),t)

if __name__=='__main__':unittest.main()
