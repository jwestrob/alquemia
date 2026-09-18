"""Actual multisite AMOEBA/native parameter receipts, no fabricated energies."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_amoeba_capability import topology
from mace_tinker_framework_solver import check_parameters,parse_parameters
BASE=ROOT/'workspaces/mace_omol_20260917'

class Multisite(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=read_json(BASE/'multisite_amoeba_capability_v3/result.json')
        cls.m=read_json(verify(cls.r['manifest']));cls.names=tuple(cls.m['cases'])

    def test_five_real_native_frameworks_retain_all_atoms_and_ions(self):
        self.assertTrue(self.r['complete']);self.assertEqual(len(self.r['rows']),5)
        for row in self.r['rows']:
            r=row['audit'];p=read_json(verify(r['source_preparation']));mapping=read_json(verify(r['mapping']))
            self.assertEqual(mapping['physical_atoms'],p['physical_atoms'])
            self.assertEqual(len(mapping['system_atom_ids'])+1,len(p['physical_atoms']))
            self.assertEqual(mapping['background_metals'],p['background_metals'])
            self.assertEqual(r['framework_formal_charge_e'],p['protein_charge_e']+2*len(p['background_metals']))
            self.assertIsNone(r['gk_settings']);self.assertIsNone(self.r['numerical_score'])
            n=r['native_GK'];receipt=read_json(verify(n['receipt']));params=parse_parameters(verify(receipt['log']).read_text())
            self.assertEqual(params,read_json(verify(n['parameters'])))
            self.assertTrue(check_parameters(params,n['task'],mapping))
            self.assertLess(n['maximum_charge_polarizability_difference'],1e-10)
            for ca in n['background_calcium']:
                a=ca['parameters'];self.assertEqual(a['charge_e'],2)
                self.assertAlmostEqual(a['polarizability_A3'],.55,places=12)
                self.assertAlmostEqual(a['radius_A'],1.82485,places=12)
                self.assertTrue(a['response_allowed'])

    def test_source_acetyl_and_water_inventory_are_preserved(self):
        for row in self.r['rows']:
            r=row['audit'];p=read_json(verify(r['source_preparation']));m=read_json(verify(r['mapping']))
            self.assertEqual(m['explicit_waters'],p['explicit_waters'])
            self.assertEqual(len(m['retained_water_bonds']),2*len(p['explicit_waters']))
            if p['case_id'].startswith('PARV'):
                self.assertIn(['A/0//C','A/1//N'],r['source_covalent_connections'])
                self.assertTrue(any(x['template']=='ACE' for x in m['residue_templates']))

    def test_unapproved_background_element_fails_explicitly(self):
        p=copy.deepcopy(read_json(verify(self.r['rows'][0]['audit']['source_preparation'])))
        a=next(a for a in p['physical_atoms'] if a['kind']=='background_metal');a['element']='La'
        for b in p['background_metals']:
            if b['id']==a['id']:b['element']='La'
        with self.assertRaisesRegex(InvalidArtifact,'species differ|background calcium'):
            topology(p,self.names,allow_background_calcium=True)

    def test_old_single_metal_default_still_rejects_multisite(self):
        p=read_json(verify(self.r['rows'][0]['audit']['source_preparation']))
        with self.assertRaisesRegex(InvalidArtifact,'unsupported cofactor'):
            topology(p,self.names)

if __name__=='__main__':unittest.main()
