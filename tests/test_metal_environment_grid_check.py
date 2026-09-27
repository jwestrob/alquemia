"""Real-artifact numerical-profile/geometry tests, not refined molecular results."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,record
from mace_omol_vacuum import embedded_input,parse_endpoint
import metal_environment_grid_check as gc
BASE=ROOT/'workspaces/metal_environment_response_20260926'
INPUTS=BASE/'preparation/scout_v3/INPUTS.json'
SOURCE=BASE/'reference_scout_v1/FINAL_COLLECTION_1219207.json'
FORCE=BASE/'force_checks_v1/COLLECTION_1219319.json'


class GridCheckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not all(p.exists() for p in (INPUTS,SOURCE,FORCE)):raise unittest.SkipTest('Actual completed fixtures unavailable')
        cls.config=read_json(INPUTS);cls.source=read_json(SOURCE)

    def test_eight_cells_and_unchanged_ab_geometry(self):
        self.assertEqual(len(gc.inventory()),8);self.assertEqual(len(set(gc.inventory())),8)
        for metal in ('Ca','La'):
            for env in ('A','B'):
                x,p=gc.rendered_geometry(self.config,metal,env)
                self.assertEqual(x,verify(self.config['endpoints'][metal]['xyz']).read_text())
                self.assertEqual(p,verify(self.config['environments'][env]['pointcharges']).read_text())

    def test_translation_and_rigid_preserve_real_distances(self):
        _,_,x=gc.xyz_data(verify(self.config['endpoints']['Ca']['xyz']))
        _,q,y=gc.pc_data(verify(self.config['environments']['A']['pointcharges']))
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory)
            for env in ('translation','rigid'):
                a,b=gc.rendered_geometry(self.config,'Ca',env);(d/'x.xyz').write_text(a);(d/'e.pc').write_text(b)
                _,_,xx=gc.xyz_data(d/'x.xyz');_,qq,yy=gc.pc_data(d/'e.pc')
                np.testing.assert_array_equal(q,qq)
                np.testing.assert_allclose(np.linalg.norm(xx[:,None,:]-yy[None,:15,:],axis=2),np.linalg.norm(x[:,None,:]-y[None,:15,:],axis=2),atol=6e-14,rtol=0)
                if env=='translation':np.testing.assert_allclose(xx-x,np.broadcast_to((.173,.117,.231),x.shape),atol=1e-14,rtol=0)

    def test_default_parser_unchanged_refined_rejects_old_grid(self):
        manifest=read_json(verify(self.source['manifest']));task=next(t for t in manifest['tasks'] if t['task_id']=='Ca_A');row=self.source['rows']['Ca_A']
        result=parse_endpoint(task,verify(row['output']),verify(row['engrad']),permanent_field=True)
        self.assertEqual(result['energy_hartree'],row['energy_hartree'])
        with self.assertRaisesRegex(InvalidArtifact,'unsupported'):
            embedded_input(-3,grid_profile='arbitrary_override')
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'refined.inp';path.write_text(embedded_input(task['charge'],grid_profile=gc.PROFILE))
            changed=copy.deepcopy(task);changed['input']=record(path)
            # Real old-grid output must never qualify as refined. No fake output.
            with self.assertRaisesRegex(InvalidArtifact,'actual refined grid/header'):
                parse_endpoint(changed,verify(row['output']),verify(row['engrad']),permanent_field=True,grid_profile=gc.PROFILE)

    def test_prepare_dryrun_and_missing_cell_rejection(self):
        agreement=ROOT/'diagnostics/metal_environment_response_20260926/FORCE_CHECK_PLAN.md'
        with tempfile.TemporaryDirectory() as directory:
            result=gc.prepare(verify(self.source['manifest']),SOURCE,FORCE,agreement,Path(directory)/'prepared',workers=8,mpi_ranks=43)
            path=verify(result['manifest']);m=read_json(path)
            self.assertEqual(m['execution_resources'],{'mpi_ranks':43,'concurrent_tasks':8})
            gc.validate(path)
            m['tasks'].pop();path.write_text(json.dumps(m))
            with self.assertRaisesRegex(InvalidArtifact,'eight-cell'):gc.validate(path)


if __name__=='__main__':unittest.main()
