"""Physical invariants of the H operation on the actual pinned 1H4I state."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from mace_hydrogen import repair, load_bonds, validate
from mace_hybrid import rotation, collect

AUDIT=ROOT/'workspaces/mace_response_trace_20260916/preparation_audit_v1/result.json'
MANIFEST=ROOT/'workspaces/mace_hydrogen_20260916/pilot_v2/medium/manifest.json'


@unittest.skipUnless(AUDIT.exists(),'requires real bonded protein audit')
class HydrogenGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        a=read_json(AUDIT)
        cls.physical=read_json(verify(a['source_state']))['physical_atoms']
        cls.bonds=load_bonds(a)

    def test_physical_identity_fixed_atoms_and_bond_lengths(self):
        original=copy.deepcopy(self.physical)
        atoms,moves=repair(self.physical,self.bonds)
        self.assertEqual(self.physical,original)
        changed={m['physical_index'] for m in moves}
        self.assertEqual(len(changed),4467)
        for i,(old,new) in enumerate(zip(original,atoms)):
            self.assertEqual({k:v for k,v in old.items() if k!='xyz_A'},
                             {k:v for k,v in new.items() if k!='xyz_A'})
            if i not in changed:self.assertEqual(old,new)
        for m in moves:
            x=np.array(atoms[m['heavy_index']]['xyz_A'])
            new=np.array(atoms[m['physical_index']]['xyz_A'])-x
            old=np.array(m['before_A'])-x
            self.assertAlmostEqual(np.linalg.norm(new),m['target_bond_length_A'],places=12)
            self.assertLess(np.linalg.norm(np.cross(new,old)),1e-12)
            self.assertGreater(np.dot(new,old),0)

    def test_rigid_transform_covariance(self):
        matrix=rotation(); shift=np.array([10.,-7.,3.])
        transformed=copy.deepcopy(self.physical)
        for a in transformed:a['xyz_A']=(np.array(a['xyz_A'])@matrix.T+shift).tolist()
        base,_=repair(self.physical,self.bonds)
        moved,_=repair(transformed,self.bonds)
        expected=np.array([a['xyz_A'] for a in base])@matrix.T+shift
        np.testing.assert_allclose([a['xyz_A'] for a in moved],expected,atol=1e-12,rtol=0)

    def test_corrupted_real_bond_inventory_rejected(self):
        index=next(i for i,b in enumerate(self.bonds) if 'H' in b['elements'].split('-'))
        duplicate=self.bonds+[self.bonds[index]]
        with self.assertRaisesRegex(InvalidArtifact,'multiple'):repair(self.physical,duplicate)
        missing=self.bonds[:index]+self.bonds[index+1:]
        with self.assertRaisesRegex(InvalidArtifact,'every protein hydrogen'):repair(self.physical,missing)

    @unittest.skipUnless(MANIFEST.exists(),'requires prepared real H pilot')
    def test_manifest_and_missing_results(self):
        m=read_json(MANIFEST);validate(m)
        bad=copy.deepcopy(m);bad['model']['external_field']=[1.,0.,0.]
        with self.assertRaisesRegex(InvalidArtifact,'electronic model'):validate(bad)
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'manifest.json';p.write_bytes(MANIFEST.read_bytes())
            c=collect(p)
            self.assertEqual(c['status'],'incomplete')
            self.assertIsNone(c['direct']);self.assertIsNone(c['hybrid']);self.assertIsNone(c['S_kcal_mol'])


if __name__=='__main__':unittest.main()
