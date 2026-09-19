"""Software checks on pinned real hydration fixtures; no scientific backend."""
import copy
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import hydration_square as h


class RealHydrationFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        config=ROOT/'diagnostics/hydration_square_20260918/SOURCE_CONFIG.json'
        if not config.exists():raise unittest.SkipTest('real archived hydration fixtures unavailable')
        cls.config=h.read_json(config)
        if not all(Path(i['preparation']['path']).exists() for i in cls.config['parents']):
            raise unittest.SkipTest('real archived preparation unavailable')
        cls.parents=[h.checked_parent(i) for i in cls.config['parents']]
        cls.reference=h.xyz(h.verify(cls.config['water_reference_geometry']))

    def test_water_identity_and_mapping(self):
        self.assertEqual([len(p[2]) for p in self.parents],[2,3])
        self.assertEqual([[w['source']['resnum'] for w in p[2]] for p in self.parents],[[211,212],[310,322,326]])
        for parent,atoms,groups in self.parents:
            for water in groups:
                left=[a for i,a in enumerate(atoms) if i not in water['indices']]
                self.assertEqual(len(left),len(atoms)-3)
                self.assertEqual(h.Counter(a[0] for a in atoms)-h.Counter(a[0] for a in left),h.Counter({'H':2,'O':1}))
                self.assertEqual((parent['outputs']['La']['all_electron_count_for_parity']-10)%2,0)

    def test_internal_geometry_preserves_heavy_atoms_plane_bisector(self):
        ref=self.reference;o=np.array(next(a[1:] for a in ref if a[0]=='O'))
        refv=np.array([a[1:] for a in ref if a[0]=='H'])-o
        target=np.linalg.norm(refv,axis=1)
        targetcos=refv[0]@refv[1]/np.prod(target)
        for parent,atoms,groups in self.parents:
            repaired=h.reference_geometry(atoms,groups,ref)
            altered={i for w in groups for i in w['hydrogen_indices']}
            for i,a in enumerate(atoms):
                if i not in altered:self.assertEqual(a,repaired[i])
            for w in groups:
                origin=np.array(atoms[w['oxygen_index']][1:])
                old=np.array([atoms[i][1:] for i in w['hydrogen_indices']])-origin
                new=np.array([repaired[i][1:] for i in w['hydrogen_indices']])-origin
                norms=np.linalg.norm(new,axis=1)
                np.testing.assert_allclose(norms,target,rtol=0,atol=1e-13)
                self.assertAlmostEqual(new[0]@new[1]/np.prod(norms),targetcos,places=13)
                u=old/np.linalg.norm(old,axis=1)[:,None]
                np.testing.assert_allclose(np.cross(old[0],old[1])/np.linalg.norm(np.cross(old[0],old[1])),np.cross(new[0],new[1])/np.linalg.norm(np.cross(new[0],new[1])),atol=1e-13)
                np.testing.assert_allclose(u.sum(axis=0)/np.linalg.norm(u.sum(axis=0)),new.sum(axis=0)/np.linalg.norm(new.sum(axis=0)),atol=1e-13)
            np.testing.assert_allclose(np.array([a[1:] for a in h.reference_geometry(repaired,groups,ref)]),np.array([a[1:] for a in repaired]),atol=1e-13,rtol=0)

    def test_corrupted_real_water_fails(self):
        parent,atoms,groups=self.parents[0]
        corrupted=copy.deepcopy(parent)
        corrupted['explicit_water_inventory'].append(copy.deepcopy(corrupted['explicit_water_inventory'][0]))
        with self.assertRaises(h.InvalidArtifact):h.water_groups(corrupted,atoms)
        corrupted=copy.deepcopy(parent)
        corrupted['atom_graph']['source_to_qm'][-1]['qm_index']=1
        with self.assertRaises(h.InvalidArtifact):h.water_groups(corrupted,atoms)

    def test_real_endpoint_components_and_receipts(self):
        for item in self.config['parents']:
            for metal,archived in item['endpoints'].items():
                e=h.endpoint(archived['output'],archived['receipt'])
                self.assertEqual(e['energy_hartree'],archived['energy_hartree'])
                self.assertEqual(e['orca_version'],'6.1.1')
                self.assertTrue(all(v is not None for v in e['components_hartree'].values()))
                c=e['components_hartree']
                self.assertAlmostEqual(c['SCF']+c['dispersion']+c['gCP'],e['energy_hartree'],places=10)


if __name__=='__main__':unittest.main()
