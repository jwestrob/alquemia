"""Actual refined outputs and immutable source fixtures; no invented SCF success."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify
from mace_omol_vacuum import embedded_input,parse_endpoint,terminal_scf_residuals,STRICT_SCF_LIMITS
import metal_environment_scf_check as sc
BASE=ROOT/'workspaces/metal_environment_response_20260926/grid_check_v1'
SOURCE=BASE/'FINAL_COLLECTION.json'
PLAN=ROOT/'diagnostics/metal_environment_response_20260926/SCF_CHECK_PLAN.md'


class StrictSCFTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not SOURCE.exists():raise unittest.SkipTest('Real refined collection unavailable')
        cls.c=read_json(SOURCE);cls.m=read_json(verify(cls.c['manifest']))

    def test_all_four_actual_density_residuals_exceed_thresholds(self):
        for metal,env in sc.inventory():
            row=self.c['rows'][f'{metal}_{env}_refined'];r=terminal_scf_residuals(verify(row['output']).read_text())
            for key in ('MAX-Density change','RMS-Density change'):
                self.assertGreater(abs(r[key]['value']),STRICT_SCF_LIMITS[key])

    def test_default_refined_parser_unchanged_and_strict_rejects_loose(self):
        task=next(t for t in self.m['tasks'] if t['task_id']=='Ca_A_refined');row=self.c['rows'][task['task_id']]
        result=parse_endpoint(task,verify(row['output']),verify(row['engrad']),permanent_field=True,grid_profile=sc.GRID)
        self.assertEqual(result['energy_hartree'],row['energy_hartree'])
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'strict.inp';path.write_text(embedded_input(task['charge'],grid_profile=sc.GRID,scf_profile=sc.SCF))
            changed=copy.deepcopy(task);changed['input']=record(path)
            with self.assertRaisesRegex(InvalidArtifact,'actual forced all-criteria'):
                parse_endpoint(changed,verify(row['output']),verify(row['engrad']),permanent_field=True,grid_profile=sc.GRID,scf_profile=sc.SCF)
        with self.assertRaisesRegex(InvalidArtifact,'unsupported'):
            embedded_input(-3,grid_profile=sc.GRID,scf_profile='arbitrary_tolerance')

    def test_four_cell_real_prepare_and_inventory_failure(self):
        if not PLAN.exists():self.skipTest('Root agreement unavailable')
        with tempfile.TemporaryDirectory() as directory:
            result=sc.prepare(SOURCE,PLAN,Path(directory)/'prepared',workers=4,mpi_ranks=86)
            path=verify(result['manifest']);m=read_json(path)
            self.assertEqual(m['execution_resources'],{'mpi_ranks':86,'concurrent_tasks':4})
            for t in m['tasks']:
                src=next(x for x in self.m['tasks'] if x['task_id']==f"{t['metal']}_{t['environment']}_refined")
                for k in ('xyz','pointcharges'):self.assertEqual(t[k]['sha256'],src[k]['sha256'])
            sc.validate(path);m['tasks'].pop();path.write_text(json.dumps(m))
            with self.assertRaisesRegex(InvalidArtifact,'four-cell inventory'):sc.validate(path)


if __name__=='__main__':unittest.main()
