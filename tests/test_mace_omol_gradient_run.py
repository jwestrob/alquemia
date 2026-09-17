"""Real gradient-pilot geometry, finite inventory and qualification guards."""
from pathlib import Path
import copy
import csv
import json
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,write_new,xyz,record
from mace_file_checks import cached_file_checks
from mace_omol_gradient_run import prepare,validate,geometry,collect,check_parent,STAGES


class RealGradientPilot(unittest.TestCase):
    def setUp(self):
        self.work=ROOT/'workspaces/mace_omol_20260917'
        self.manifest=self.work/'masked_gradient_core_v4/manifest.json'
        if not self.manifest.exists():self.skipTest('requires real pinned1H4I gradient manifest')

    @cached_file_checks
    def test_actual_core_inventory_and_only_selected_metal_displaced(self):
        self.assertEqual(validate(self.manifest)['tasks'],10)
        m=read_json(self.manifest)
        for t in m['tasks']:
            source=xyz(t['source_xyz']['path']);new=geometry(t)
            self.assertEqual(len(new),73)
            self.assertEqual([a[0] for a in source],[a[0] for a in new])
            if t['gradient_variant'] in ('positive','negative'):
                moved=[i for i,(a,b) in enumerate(zip(source,new)) if a!=b]
                self.assertEqual(moved,[t['metal_index']])
                self.assertAlmostEqual(np.linalg.norm(np.array(source[moved[0]][1:])-new[moved[0]][1:]),.01,places=12)

    def test_full_stage_requires_actual_passing_core_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(InvalidArtifact,'passing actual core-gradient qualification'):
                prepare(self.work/'charge_ablation_development_v2/collection_job_1200828.json',
                        ROOT/'diagnostics/mace_omol_20260917/MASKED_GRADIENT_PLAN.md',
                        Path(tmp)/'unexecuted',STAGES[1])

    @cached_file_checks
    def test_corrupted_real_manifest_cannot_change_displacement_direction(self):
        m=copy.deepcopy(read_json(self.manifest));m['tasks'][0]['direction'][0]+=.1
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'corrupted_real_manifest.json';write_new(path,m)
            with self.assertRaisesRegex(InvalidArtifact,'fixed inventory'):validate(path)

    def test_actual_complete_core_checks_are_serializable_and_pass(self):
        if not (self.work/'masked_gradient_core_report_v5/result.json').exists():
            self.skipTest('requires actual ten-forward gradient qualification')
        r=collect(self.manifest)
        json.dumps(r,allow_nan=False)
        self.assertEqual(r['status'],'complete')
        self.assertTrue(r['numerical_gate_pass'])
        self.assertEqual(len(r['checks']),23)

    def test_actual_gradient_units_and_negative_gradient_convention(self):
        from mace_omol_gradient_worker import accepted
        p=self.work/'masked_gradient_core_v4/execution/Ca_checkpointed_center/attempt_001/result.json'
        if not p.exists():self.skipTest('requires actual checkpointed gradient')
        r=read_json(p);m=read_json(self.manifest)
        t=next(t for t in m['tasks'] if t['task_id']==r['task_id'])
        self.assertTrue(accepted(r,t))
        np.testing.assert_array_equal(np.load(r['gradient']['path']),-np.load(r['forces']['path']))
        damaged=copy.deepcopy(r);damaged['gradient_unit']='Hartree_per_Bohr'
        self.assertFalse(accepted(damaged,t))

    def test_roundoff_tolerance_cannot_change_actual_qualification_decisions(self):
        full=self.work/'masked_gradient_full_v3/manifest.json'
        if not full.exists():self.skipTest('requires declared full-stage manifest and actual core report')
        m=read_json(full);r=read_json(m['core_report']['path'])
        with tempfile.TemporaryDirectory() as tmp:
            # Explicitly perturbed copies of real diagnostic fields, never model outputs.
            rounded=copy.deepcopy(r);rounded['projections']['R']['gradient_model_kcal_per_A']+=1e-13
            path=Path(tmp)/'roundoff_copy_of_actual_report.json';write_new(path,rounded)
            changed=copy.deepcopy(m);changed['core_report']=record(path);check_parent(changed)
            damaged=copy.deepcopy(r);damaged['checks'][0]['pass']=False
            path=Path(tmp)/'corrupted_actual_report_decision.json';write_new(path,damaged)
            changed['core_report']=record(path)
            with self.assertRaisesRegex(InvalidArtifact,'criterion or pass/fail'):check_parent(changed)

    def test_actual_whole_protein_source_map_and_paired_gradient(self):
        path=self.work/'masked_gradient_full_report_v1/result.json'
        if not path.exists():self.skipTest('requires the actual six-call whole-GGR gradient result')
        r=read_json(path);self.assertTrue(r['numerical_gate_pass'])
        self.assertEqual(len(r['checks']),11)
        with Path(r['source_mapped_gradients']['path']).open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        atoms=read_json(r['preparation']['path'])['physical_atoms']
        self.assertEqual([v['source_id'] for v in rows],[v['id'] for v in atoms])
        gradient=np.load(r['paired_gradient']['path'])
        self.assertEqual(gradient.shape,(4698,3))
        for i,row in enumerate(rows):
            actual=[float(row[f'R_d{axis}_model_kcal_per_A']) for axis in ('x','y','z')]
            np.testing.assert_array_equal(actual,gradient[i])
            difference=[float(row[f'Ca_d{axis}_model_kcal_per_A'])-float(row[f'La_d{axis}_model_kcal_per_A'])
                        for axis in ('x','y','z')]
            np.testing.assert_allclose(actual,difference,rtol=0,atol=1e-12)


if __name__=='__main__':unittest.main()
