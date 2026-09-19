"""Geometry/algebra on real prepared alpha structures; no fake energies."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import chemical_states_proton as proton
from affordable_common import read_json,verify,xyz

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'workspaces/chemical_states_20260919/proton_paths_v1/manifest.json'


class ProtonPathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest=read_json(MANIFEST)
        cls.prep=read_json(verify(cls.manifest['preparation']))

    def test_complete_predeclared_coverage(self):
        self.assertEqual(len(self.prep['paths']),6)
        self.assertEqual(len(self.manifest['tasks']),24)
        self.assertEqual(self.manifest['unsupported_endpoints'],0)
        self.assertEqual(proton.validate(MANIFEST)['tasks'],24)
        for key in self.prep['paths']:
            self.assertEqual(sorted(t['fraction'] for t in self.manifest['tasks'] if t['path_id']==key),[.25,.5,.75,1.])

    def test_balanced_atoms_charge_and_frozen_scaffold(self):
        for t in self.manifest['tasks']:
            path=self.prep['paths'][t['path_id']];c=self.prep['centers'][path['center_id']]
            before=xyz(verify(c['xyz']));after=xyz(verify(t['xyz']));h=path['contact']['hydrogen_index']
            self.assertEqual([a[0] for a in before],[a[0] for a in after])
            self.assertEqual(t['charge'],c['charge'])
            self.assertEqual([a for i,a in enumerate(before) if i!=h],[a for i,a in enumerate(after) if i!=h])
        self.assertEqual(self.prep['reaction']['delta_H_atoms'],0)
        self.assertEqual(self.prep['reaction']['proton_reservoir_coefficient'],0)

    def test_real_geometry_endpoint_and_analytic_mapping(self):
        for path in self.prep['paths'].values():
            c=self.prep['centers'][path['center_id']];atoms=xyz(verify(c['xyz']));contact=path['contact']
            h=contact['hydrogen_index'];a=contact['acceptor_index']
            np.testing.assert_allclose(proton.path_position(atoms,contact,0)[0],atoms[h][1:],atol=1e-13,rtol=0)
            final=proton.path_position(atoms,contact,1)[0]
            self.assertAlmostEqual(np.linalg.norm(final-np.array(atoms[a][1:])),.98,places=12)
            for t in [.25,.5,.75]:
                eps=1e-5
                numeric=(proton.path_position(atoms,contact,t+eps)[0]-proton.path_position(atoms,contact,t-eps)[0])/(2*eps)
                np.testing.assert_allclose(proton.path_position(atoms,contact,t)[1],numeric,atol=1e-8,rtol=0)

    def test_actual_source_energy_and_gradient_receipts(self):
        _,centers=proton.source_centers(verify(self.prep['source_collection']))
        self.assertEqual(centers,self.prep['centers'])
        self.assertEqual({c['case'] for c in centers.values()},{'1F6S','6IP9'})

    def test_no_water_oxygen_crossing_on_real_curved_paths(self):
        for path in self.prep['paths'].values():
            c=self.prep['centers'][path['center_id']];atoms=xyz(verify(c['xyz']));contact=path['contact']
            od=np.array(atoms[contact['donor_oxygen_index']][1:])
            for t in np.linspace(0,1,101):
                h,_=proton.path_position(atoms,contact,float(t))
                self.assertGreater(np.linalg.norm(h-od),.95)

    def test_actual_partial_collection_preserves_unavailable_states(self):
        c=read_json(ROOT/'workspaces/chemical_states_20260919/proton_paths_v1/collection_partial_v3.json')
        self.assertEqual(len(c['rows']),24)
        complete=[r for r in c['rows'] if r['status']=='complete']
        self.assertEqual(len(complete),4)
        for row in c['rows']:
            if row['status']!='complete':
                self.assertIsNone(row['energy_change_kcal_mol'])
                self.assertIsNone(row['path_gradient_kcal_mol'])
                continue
            t=next(t for t in self.manifest['tasks'] if t['task_id']==row['task_id'])
            source=self.prep['centers'][self.prep['paths'][row['path_id']]['center_id']]
            actual=proton.endpoint(row['result']['output'],row['result']['receipt'],t['xyz'],t['input'])
            self.assertEqual(actual,row['result'])
            self.assertAlmostEqual(row['energy_change_kcal_mol'],(actual['energy_hartree']-source['result']['energy_hartree'])*proton.HA_TO_KCAL,places=12)

    def test_actual_balanced_metal_difference_sign_and_single_conversion(self):
        c=read_json(ROOT/'workspaces/chemical_states_20260919/proton_paths_v1/collection_partial8_v3.json')
        end={r['metal']:r for r in c['rows'] if r['case']=='1F6S' and r['fraction']==1.}
        source={metal:self.prep['centers'][self.prep['paths'][r['path_id']]['center_id']]['result']['energy_hartree'] for metal,r in end.items()}
        direct=((end['Ca']['result']['energy_hartree']-end['La']['result']['energy_hartree'])-(source['Ca']-source['La']))*proton.HA_TO_KCAL
        split=end['Ca']['energy_change_kcal_mol']-end['La']['energy_change_kcal_mol']
        self.assertAlmostEqual(direct,split,places=7)
        self.assertAlmostEqual(split,3.638224171429037,places=10)
        self.assertGreater(split,0.)

    def test_all_actual_native_receipts_and_gradients(self):
        c=read_json(ROOT/'workspaces/chemical_states_20260919/proton_paths_v1/collection_1202444_recovered_v4.json')
        self.assertEqual(c['status'],'complete');self.assertEqual(len(c['rows']),24)
        for r in c['rows']:
            t=next(t for t in self.manifest['tasks'] if t['task_id']==r['task_id'])
            self.assertEqual(r['status'],'complete')
            self.assertEqual(proton.endpoint(r['result']['output'],r['result']['receipt'],t['xyz'],t['input']),r['result'])
            gradient=proton.checked_gradient(r['gradient'],r['result'],xyz(verify(t['xyz'])))
            self.assertTrue(np.isfinite(gradient).all())


if __name__=='__main__':unittest.main()
