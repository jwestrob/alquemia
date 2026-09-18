"""Real-source invariants for the numerical resolution experiment."""
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,cache_key,read_json,verify,write_new
from mace_ddx_resolution import preflight
W=ROOT/'workspaces/mace_omol_20260917'


class Resolution(unittest.TestCase):
    def test_actual_source_and_physical_method_remain_identical(self):
        m=preflight(W/'ddx_resolution_v1/manifest.json');parent=read_json(verify(m['parent']))
        source=next(g['source'] for g in parent['groups'] if g['group_id']=='GGR_2FW0_primary_primary')
        self.assertEqual(m['model'],parent['model']);self.assertEqual(m['tolerances'],parent['tolerances'])
        self.assertEqual(sum(len(g['new_states']) for g in m['groups']),8)
        for g in m['groups']:
            self.assertEqual(g['source'],source);self.assertEqual(g['states'],['Ca','La'])
        self.assertEqual(m['grids']['integration24']['n_lebedev'],m['grids']['basis30']['n_lebedev'])
        self.assertEqual(m['grids']['basis30']['lmax'],m['grids']['integration30']['lmax'])

    def test_corrupted_real_cavity_parameter_rejected_despite_new_hash(self):
        m=read_json(W/'ddx_resolution_v1/manifest.json');m['model']['eta']=.2
        m['cache_key']=cache_key({k:v for k,v in m.items() if k!='cache_key'})
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_cavity.json';write_new(p,m)
            with self.assertRaises(InvalidArtifact):preflight(p)


if __name__=='__main__':unittest.main()
