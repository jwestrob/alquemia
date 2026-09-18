"""Matched physical centers and source-mapped 3D grids from real hybrid artifacts."""
from pathlib import Path
import copy,json,sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,xyz,InvalidArtifact
from mace_omol_hybrid_response import sources,validate
PREP=ROOT/'workspaces/mace_omol_hybrid_response_20260918/prepared_v1/preparation.json'

@unittest.skipUnless(PREP.exists(),'requires actual matched-H response preparation')
class HybridResponseTests(unittest.TestCase):
    def test_exact_real_centers_and_vacuum_gradients_reused(self):
        top=read_json(PREP);p,r,gr,gm,states,fulls,reused=sources(verify(top['source_preparation']),verify(top['static_report']),verify(top['ggr_gradients']))
        self.assertEqual(len(reused),2);self.assertEqual(len(states),4);self.assertEqual(len(fulls),3)
        self.assertEqual(states['GGR_extended']['Ca']['full'],states['GGR_connected']['Ca']['full'])
        for name,row in states.items():
            for metal,s in row.items():
                self.assertEqual(s['DFT_source']['energy_scope'],'isolated_vacuum_endpoint')
                self.assertEqual(s['DFT_gradient_kcal_mol_A'],s['DFT_source']['gradient_kcal_mol_per_A'][0])
                self.assertEqual(xyz(verify(s['core']['xyz']))[0][0],metal)
        self.assertEqual(set(reused),{'GGR_1GLG_Ca','GGR_1GLG_La'})

    def test_finite_grid_counts_and_fixed_scaffold(self):
        top=read_json(PREP);total=0;gradient=0
        for pin in top['jobs']:
            mp=verify(pin);m=read_json(mp);self.assertEqual(validate(mp)['status'],'pass')
            total+=len(m['tasks'])
            for t in m['tasks']:
                old=xyz(verify(t['source_xyz']));new=xyz(verify(t['xyz']));i=t['metal_index']
                self.assertEqual(old[:i]+old[i+1:],new[:i]+new[i+1:])
                np.testing.assert_allclose(np.array(new[i][1:])-old[i][1:],t['displacement_A'],atol=5e-10,rtol=0)
                gradient+=not t['energy_only']
        self.assertEqual(total,228);self.assertEqual(gradient,12)

    def test_corrupted_real_state_and_hidden_grid_extension_rejected(self):
        top=read_json(PREP);m=read_json(verify(top['jobs'][0]))
        for mutation in ('charge','displacement'):
            damaged=copy.deepcopy(m)
            if mutation=='charge':damaged['tasks'][0]['charge']+=1
            else:damaged['tasks'][0]['displacement_A'][0]+=.01
            with tempfile.TemporaryDirectory() as td:
                path=Path(td)/'corrupted_real_manifest.json';path.write_text(json.dumps(damaged))
                with self.assertRaises(InvalidArtifact):validate(path)

if __name__=='__main__':unittest.main()
