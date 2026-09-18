"""Actual multisite physical/core/source bridge and deliberate corruptions."""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_multisite_density import audit_alias
from mace_density_inputs import normalize
BASE=ROOT/'workspaces/mace_omol_20260917'

class MultisiteInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=read_json(BASE/'multisite_density_prepared_v2/preparation.json')
        cls.cases={k:read_json(verify(v)) for k,v in cls.p['cases'].items()}

    def test_actual_paired_geometry_and_complete_core_water_membership(self):
        self.assertEqual(len(self.cases),5)
        for case in self.cases.values():
            p=read_json(verify(case['normalized_global_preparation']));audit_alias(p)
            ca,la=(xyz(verify(case['endpoints'][m]['xyz'])) for m in ('Ca','La'))
            self.assertEqual(ca[1:],la[1:]);self.assertEqual(ca[0][1:],la[0][1:])
            c=read_json(verify(case['source_core']));support=set(read_json(verify(case['projection']))['physical_ids'])
            for w in c['explicit_water_inventory']:
                ids={a['id'] for a in p['physical_atoms'] if a['id'].startswith(f"{w['chain']}/{w['resnum']}/{w['insertion_code']}/")}
                self.assertEqual(len(ids),3);self.assertTrue(ids<=support)
            self.assertEqual([a for a in ca if a[0]!='H'],[a for a in xyz(verify(c['outputs']['Ca']['xyz'])) if a[0]!='H'])

    def test_fixed_field_charge_closure_and_no_charge_on_QM_support(self):
        expected={'AEQ_1SL8_EF1':-4,'AEQ_1SL8_EF3':-3,'AEQ_1SL8_EF4':-3,'PARV_4CPV_CD':-1,'PARV_4CPV_EF':-1}
        for name,c in self.cases.items():
            s=read_json(verify(c['state']));support=set(s['projection_support_ids']);p=read_json(verify(c['normalized_global_preparation']))
            charge=dict(zip([a['id'] for a in s['physical_atoms']],s['environment_charges_e']))
            self.assertTrue(all(charge[i]==0 for i in support));self.assertAlmostEqual(sum(charge.values()),expected[name],places=10)
            for a in p['background_metals']:self.assertEqual(charge[a['id']],2)
            for w in p['explicit_waters']:
                ids=[a['id'] for a in p['physical_atoms'] if a['id'].rsplit('/',1)[0]==w['oxygen_id'].rsplit('/',1)[0]]
                self.assertAlmostEqual(sum(charge[i] for i in ids),0,places=12)
            self.assertIsNone(s['QM_charges']);self.assertIsNone(s['environment_correction'])

    def test_corrupted_real_alias_cannot_drop_a_background_ion(self):
        c=next(iter(self.cases.values()));p=copy.deepcopy(read_json(verify(c['normalized_global_preparation'])))
        a=next(a for a in p['physical_atoms'] if a['kind']=='background_metal');p['physical_atoms'].remove(a)
        with self.assertRaisesRegex(InvalidArtifact,'changed a source state'):audit_alias(p)

    def test_dry_single_metal_interface_remains_explicit(self):
        c=next(iter(self.cases.values()))
        with self.assertRaisesRegex(InvalidArtifact,'dry single-metal'):
            normalize(verify(c['normalized_global_preparation']),verify(c['source_core']))

if __name__=='__main__':unittest.main()
