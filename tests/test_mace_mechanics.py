"""Real source-mapped coupled geometries and quantum manifest invariants."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,xyz
from affordable_response import source_key
from affordable_workflow import dry_run
from mace_mechanics import GRID,physical_id
MANIFEST=ROOT/'workspaces/mace_mechanics_20260916/prepared_v2/manifest.json'


@unittest.skipUnless(MANIFEST.exists(),'requires real coupled preparation')
class MechanicsPreparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=read_json(MANIFEST)

    def test_exact_quantum_inventory_and_analytic_centers(self):
        q=read_json(verify(self.m['DFT_tasks']));self.assertEqual(dry_run(verify(self.m['DFT_tasks']))['tasks'],36)
        self.assertEqual(sum(t['task_type']=='analytic_gradient' for t in q['tasks']),4)
        self.assertEqual(len(self.m['reused_GGR']),20)
        for t in q['tasks']:
            text=verify(t['input']).read_text()
            self.assertIn('r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF',text)
            self.assertNotIn('NumGrad',text)
            self.assertEqual('EnGrad' in text,t['task_type']=='analytic_gradient')

    def test_source_atoms_match_full_motion_and_caps_stay_fixed(self):
        for ref in self.m['cases'].values():
            c=read_json(verify(ref));physical=read_json(verify(c['physical_atoms']))
            indices={a['source_key']:i for i,a in enumerate(physical)}
            full=self.m['physical_cases'][c['global_id']]
            center=xyz(verify(c['grids']['center']['endpoints']['La']['xyz']))
            self.assertEqual(set(c['grids']),set(GRID))
            for point,grid in c['grids'].items():
                core=xyz(verify(grid['endpoints']['La']['xyz']))
                whole=xyz(verify(full['grids'][point]['endpoints']['La']['xyz']))
                np.testing.assert_allclose(core[0][1:],whole[indices['metal']][1:],atol=1e-9,rtol=0)
                for atom in c['source_graph']['source_to_qm']:
                    i=atom['qm_index']
                    if atom['kind']=='source':
                        j=indices[source_key(atom['source'])]
                        self.assertEqual(core[i][0],whole[j][0]);np.testing.assert_allclose(core[i][1:],whole[j][1:],atol=1e-9,rtol=0)
                    else:np.testing.assert_allclose(core[i][1:],center[i][1:],atol=1e-9,rtol=0)
                ca=xyz(verify(grid['endpoints']['Ca']['xyz']))
                self.assertEqual(core[0][1:],ca[0][1:]);self.assertEqual(core[1:],ca[1:])
                self.assertLessEqual(grid['validation']['maximum_heavy_displacement_A'],.05)

    def test_physical_selection_and_common_GGR_scaffold(self):
        cases={k:read_json(verify(v)) for k,v in self.m['cases'].items()}
        for name,c in cases.items():
            self.assertEqual(c['selected_coordinate']['oxygen']['resnum'],140 if name.startswith('GGR') else 84)
            self.assertEqual(c['selected_coordinate']['amide_nitrogen']['resnum'],141 if name.startswith('GGR') else 85)
            jac=np.load(verify(c['core_jacobians']));self.assertEqual(jac.shape[0],2)
            caps=[r['qm_index'] for r in c['source_graph']['source_to_qm'] if r['kind']!='source']
            np.testing.assert_array_equal(jac[:,caps],0)
        self.assertEqual(cases['GGR_extended']['physical_atoms'],cases['GGR_connected']['physical_atoms'])
        self.assertEqual(cases['GGR_extended']['specs'],cases['GGR_connected']['specs'])


if __name__=='__main__':unittest.main()
