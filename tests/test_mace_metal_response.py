"""Pinned real-structure/gradient checks; no invented scientific values."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,xyz,HA_TO_KCAL,BOHR_TO_A,InvalidArtifact
from affordable_response import read_engrad
from mace_metal_response import grid,shift,validate
from mace_hybrid import EV_TO_KCAL
PREP=ROOT/'workspaces/mace_metal_response_20260918/prepared_v2/preparation.json'

@unittest.skipUnless(PREP.exists(),'requires the real prepared 3D response fixtures')
class MetalResponseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.p=read_json(PREP)

    def test_exact_pairs_only_physical_metal_moves(self):
        self.assertEqual(len(grid()),37)
        for name,c in self.p['cases'].items():
            base=xyz(verify(c['states']['center']['La']['xyz']))
            for point,ends in c['states'].items():
                la,ca=(xyz(verify(ends[m]['xyz'])) for m in ('La','Ca'))
                self.assertEqual(la[0][1:],ca[0][1:]);self.assertEqual(la[1:],ca[1:])
                self.assertEqual(la[1:],base[1:])
                np.testing.assert_allclose(np.array(la[0][1:])-base[0][1:],grid()[point],atol=1e-9,rtol=0)
                self.assertEqual(ends['La']['charge']-ends['Ca']['charge'],1)
                full=xyz(verify(self.p['physical'][c['global_id']][point]['La']['xyz']))
                _,i=shift(full,'La',[0,0,0]);np.testing.assert_allclose(full[i][1:],la[0][1:],atol=1e-9,rtol=0)
        self.assertEqual(len(self.p['physical']),3)

    def test_archived_DFT_and_short_force_sign_units(self):
        for name,ends in self.p['centers'].items():
            for metal,c in ends.items():
                g=read_engrad(verify(c['DFT_artifacts']['engrad']))
                np.testing.assert_allclose(c['DFT_gradient_kcal_mol_A'],np.array(g['gradient_Ha_per_bohr'][0])*HA_TO_KCAL/BOHR_TO_A,atol=1e-10)
                gj=(-np.load(verify(c['full_short']['forces']))[c['full_metal_index']]+np.load(verify(c['core_short']['forces']))[0])*EV_TO_KCAL
                np.testing.assert_allclose(c['J_gradient_kcal_mol_A'],gj,atol=1e-10)
                self.assertEqual(c['DFT_energy_Ha'],g['energy_Ha'])
                if name.startswith('ALPHA'):self.assertLess(c['old_axis_projection_fraction'],.32)

    def test_all_finite_manifests_dry_run(self):
        jobs=read_json(PREP.with_name('jobs.json'))['jobs'];self.assertEqual(len(jobs),7)
        for ref in jobs:self.assertEqual(validate(verify(ref))['tasks'],72)

    def test_corrupted_real_metal_identity_rejected(self):
        c=self.p['cases']['GGR_extended']['states']['center']['La'];rows=xyz(verify(c['xyz']))
        corrupted=list(rows);corrupted[0]=('C',*rows[0][1:])
        with self.assertRaises(InvalidArtifact):shift(corrupted,'La',[.01,0,0])

if __name__=='__main__':unittest.main()
