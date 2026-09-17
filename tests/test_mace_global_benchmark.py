"""Physical preparation invariants on five pinned real structures, no dummy science."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_global_prepare import CASES,terminal_position
from mace_global_benchmark import physical_preparations,validate,collect_global,compare_report
from mace_hybrid import rotation
PREP=ROOT/'workspaces/mace_global_benchmark_20260916/prepared_v1/preparation_manifest.json'
MANIFEST=ROOT/'workspaces/mace_global_benchmark_20260916/mace_v1/medium/manifest.json'


@unittest.skipUnless(PREP.exists(),'real five-protein preparations required')
class GlobalPreparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=read_json(PREP)
        cls.cases={r['case_id']:read_json(verify(r['preparation'])) for r in cls.p['rows']}

    def test_real_physical_pair_mapping_and_formal_charges(self):
        physical_preparations(PREP)
        self.assertEqual(tuple(self.cases),CASES)
        for c in self.cases.values():
            a,b=(c['endpoints'][m] for m in ('La','Ca'))
            self.assertEqual(a['charge']-b['charge'],1)
            self.assertEqual(a['state']['all_electron_count']-b['state']['all_electron_count'],36)
            for x,y in zip(xyz(verify(a['xyz'])),xyz(verify(b['xyz']))):self.assertEqual(x[1:],y[1:])

    def test_existing_heavy_atoms_fixed_and_H_radial_targets(self):
        for c in self.cases.values():
            byid={a['id']:a for a in c['physical_atoms']}
            details=c['preparation_details']
            for old in details['original_protein_atoms']:
                if old['element']!='H':self.assertEqual(byid[old['id']]['xyz_A'],old['xyz_A'])
            for move in details['protein_H_moves']:
                parent=np.array(byid[move['heavy_source_id']]['xyz_A'])
                a=np.array(move['before_A'])-parent;b=np.array(byid[move['source_id']]['xyz_A'])-parent
                self.assertAlmostEqual(np.linalg.norm(b),move['target_bond_length_A'],places=12)
                self.assertLess(np.linalg.norm(np.cross(a,b)),1e-12)
                self.assertGreater(np.dot(a,b),0)

    def test_only_expected_terminal_completions_and_covariance(self):
        mat=rotation();shift=np.array([10.,-7.,3.])
        for name,c in self.cases.items():
            additions=c['preparation_details']['terminal_additions']
            self.assertEqual(len(additions),1 if name.startswith('PQQ') else 0)
            atoms={a['id']:a for a in c['physical_atoms']}
            for addition in additions:
                carbon=np.array(atoms[addition['bond_to']]['xyz_A'])
                stem=addition['bond_to'].rsplit('/',1)[0]
                ca=np.array(atoms[stem+'/CA']['xyz_A']);o=np.array(atoms[stem+'/O']['xyz_A'])
                expected=terminal_position(carbon,ca,o)
                np.testing.assert_array_equal(expected,addition['xyz_A'])
                moved=terminal_position(carbon@mat.T+shift,ca@mat.T+shift,o@mat.T+shift)
                np.testing.assert_allclose(moved,expected@mat.T+shift,atol=1e-12,rtol=0)
                self.assertAlmostEqual(np.linalg.norm(expected-carbon),1.25,places=12)
                self.assertGreater(addition['metal_distance_A'],40)

    def test_frozen_PQQ_and_water_oxygen_inventory(self):
        for c in self.cases.values():
            core=xyz(verify(c['source_audit_row']['endpoints']['La']['xyz']))
            atoms=c['physical_atoms'];byid={a['id']:a for a in atoms}
            pqq=[a for a in atoms if a['kind']=='frozen_PQQ']
            if pqq:self.assertEqual([(a['element'],*a['xyz_A']) for a in pqq],core[1:28])
            waters=[a for a in atoms if a['kind']=='retained_site_water']
            self.assertEqual(len(waters),3*len(c['explicit_waters']))
            for a in waters:
                if a['element']=='O':self.assertEqual(tuple(a['xyz_A']),core[a['source_qm_index']][1:])
            for move in c['water_H_moves']:
                o=np.array(byid[move['id'].rsplit('/',1)[0]+'/O']['xyz_A'])
                before=np.array(move['before_A'])-o;after=np.array(move['after_A'])-o
                self.assertAlmostEqual(np.linalg.norm(after),.9572,places=12)
                self.assertLess(np.linalg.norm(np.cross(before,after)),1e-12)

    @unittest.skipUnless(MANIFEST.exists(),'requires prepared MACE execution manifest')
    def test_manifest_and_corrupted_case_rejection(self):
        self.assertEqual(validate(MANIFEST)['tasks'],14)
        m=read_json(MANIFEST)
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'manifest.json';p.write_bytes(MANIFEST.read_bytes())
            c=collect_global(p)
            self.assertEqual(c['status'],'incomplete');self.assertFalse(c['development_ordering_screen_pass'])
            self.assertTrue(all(s['R_kcal_mol'] is None for s in c['scores'].values()))
            altered=copy.deepcopy(m);altered['tasks'][0]['charge']+=1;p.write_text(json.dumps(altered))
            with self.assertRaises(InvalidArtifact):validate(p)

    def test_actual_completed_failure_cannot_become_success_in_reporting(self):
        medium=ROOT/'workspaces/mace_global_benchmark_20260916/gb_v1/medium/collection_job_1200711.json'
        large=ROOT/'workspaces/mace_global_benchmark_20260916/gb_v1/large/collection_job_1200712.json'
        if not medium.exists() or not large.exists():self.skipTest('requires actually completed global GB panel')
        with tempfile.TemporaryDirectory() as folder:
            r=compare_report(medium,large,Path(folder)/'report')
            d=read_json(verify(r['comparison']))
            self.assertFalse(d['primary_candidate_ordering_screen_pass'])
            self.assertEqual(d['numerical_checks'],{'medium':True,'large':True})
            self.assertAlmostEqual(d['contrasts']['medium'][0]['difference_kcal_mol'],-17.628167761780787,places=8)
            self.assertTrue(all(c['expected_order'] is False for cs in d['contrasts'].values() for c in cs))
            self.assertIsNone(d['calibrated_class']);self.assertFalse(d['goal_completion_claimed'])
            with self.assertRaises(InvalidArtifact):compare_report(large,medium,Path(folder)/'swapped')


if __name__=='__main__':unittest.main()
