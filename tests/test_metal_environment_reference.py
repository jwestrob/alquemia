"""Archived real native artifacts and corrupted copies; no fabricated energies."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import verify, read_json, InvalidArtifact
from metal_environment_reference import read_pcgrad
from mace_omol_vacuum import parse_endpoint


class ReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.archive=ROOT/'workspaces/mace_omol_20260917/responsive_quantum_v2'
        m=read_json(cls.archive/'manifest.json')
        cls.task=next(t for t in m['tasks'] if t['task_id']=='GGR_extended_Ca')
        cls.pc=cls.archive/'GGR_extended_Ca/endpoint.runtime.pcgrad'
        cls.count=int(verify(cls.task['pointcharges']).read_text().splitlines()[0])

    def test_native_analytic_gradient_and_external_inventory(self):
        t=self.task
        result=parse_endpoint(t,t['output_path'],t['engrad_path'],permanent_field=True)
        g=read_pcgrad(self.pc,self.count)
        self.assertEqual(g.shape,(self.count,3))
        self.assertEqual(result['quantity'],'gradient_not_force')
        self.assertTrue(all(result['native_components'].values()))

    def test_corrupted_real_pcgrad_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'corrupted.pcgrad'
            path.write_text('\n'.join(self.pc.read_text().splitlines()[:-1])+'\n')
            with self.assertRaises(InvalidArtifact): read_pcgrad(path,self.count)

    def test_pointcharge_inventory_mismatch(self):
        with self.assertRaises(InvalidArtifact): read_pcgrad(self.pc,self.count+1)


if __name__=='__main__': unittest.main()
