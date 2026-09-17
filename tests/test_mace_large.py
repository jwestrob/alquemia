"""The same independent kernel tests, with the actual large checkpoint parameters."""
from pathlib import Path
import unittest
import numpy as np
import test_mace_analytic as reference_tests
from affordable_common import read_json,verify,xyz

ROOT=Path(__file__).resolve().parents[1]
MODEL=ROOT/'workspaces/mace_large_20260916/software_v1/MACE-POLAR-1-L.model'

@unittest.skipUnless(MODEL.exists(),'requires pinned released large checkpoint')
class LargeAnalyticMultipoleTests(reference_tests.AnalyticMultipoleTests):
    @classmethod
    def setUpClass(cls):
        import torch
        from mace.calculators import mace_polar
        from mace_realspace_compat import configure
        from mace_analytic import configure as analytic
        cls.torch=torch;torch.set_num_threads(1);torch.set_default_dtype(torch.float64)
        cls.c=read_json(reference_tests.COLLECTION);cls.m=read_json(verify(cls.c['manifest']))
        cls.calc=mace_polar(model=str(MODEL),device='cpu',default_dtype='float64')
        configure(cls.calc)
        cls.old_state={k:v.clone() for k,v in cls.calc.models[0].state_dict().items()}
        cls.fixtures=[]
        for t in cls.m['tasks']:
            if t['kind']=='core':
                cls.fixtures.append((t['task_id'],np.array([a[1:] for a in xyz(verify(t['xyz']))]),
                    np.load(verify(cls.c['rows'][t['task_id']]['density_coefficients']))))
        analytic(cls.calc,tile=17,edge_tile=128,node_tile=17)

if __name__=='__main__':unittest.main()
