#!/usr/bin/env python3
"""Saved-only controls for balanced sign, missing data and artifact rejection."""
import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('collector',Path(__file__).with_name('COLLECT.py'))
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)


class Controls(unittest.TestCase):
    def setUp(self):
        self.g={};self.data={m:{} for m in c.METHODS}
        for metal in ('La','Dy'):
            for prefix,symbols in [('wet',[metal,'O','H','H','N']),('dry',[metal,'N']),('aq',[metal,'O','H','H'])]:
                gid=prefix+metal;self.g[gid]={'symbols':symbols,'charge':3}
                self.data['MACE'][gid]={'energy_kcal_mol':0.,'result':{'fixture':gid}}
        for gid,e in {'wetDy':40.,'wetLa':30.,'dryDy':9.,'dryLa':3.,'aqDy':80.,'aqLa':100.}.items():
            self.data['MACE'][gid]['energy_kcal_mol']=e

    def test_balanced_sign_and_common_reference_cancellation(self):
        wet=c.cycle(self.data,'MACE',['wetDy','aqLa','wetLa','aqDy'],self.g)
        dry=c.cycle(self.data,'MACE',['dryDy','aqLa','dryLa','aqDy'],self.g)
        direct=c.cycle(self.data,'MACE',['wetDy','dryLa','wetLa','dryDy'],self.g)
        self.assertEqual(wet['value_kcal_mol'],30.)
        self.assertEqual(dry['value_kcal_mol'],26.)
        self.assertEqual(wet['value_kcal_mol']-dry['value_kcal_mol'],direct['value_kcal_mol'])
        self.assertEqual(direct['value_kcal_mol'],4.)
        self.data['MACE']['aqLa']['energy_kcal_mol']+=1e5
        self.assertEqual(c.cycle(self.data,'MACE',['wetDy','aqLa','wetLa','aqDy'],self.g)['value_kcal_mol']-
                         c.cycle(self.data,'MACE',['dryDy','aqLa','dryLa','aqDy'],self.g)['value_kcal_mol'],4.)

    def test_missing_is_not_zero(self):
        del self.data['MACE']['dryLa']
        r=c.cycle(self.data,'MACE',['wetDy','dryLa','wetLa','dryDy'],self.g)
        self.assertIsNone(r['value_kcal_mol']);self.assertEqual(r['status'],'unavailable')
        self.assertEqual(r['missing_terms'],['dryLa'])

    def test_unbalanced_charge_rejected(self):
        self.g['dryLa']['charge']=2
        with self.assertRaisesRegex(ValueError,'Unbalanced'):c.cycle(self.data,'MACE',['wetDy','dryLa','wetLa','dryDy'],self.g)

    def test_unbalanced_atoms_rejected(self):
        self.g['dryLa']['symbols'].append('H')
        with self.assertRaisesRegex(ValueError,'Unbalanced'):c.cycle(self.data,'MACE',['wetDy','dryLa','wetLa','dryDy'],self.g)

    def test_force_sign_and_all_ligand_atoms(self):
        self.assertEqual(c.radial([[100,100,100],[1,2,3],[4,5,6]],{'ligand_indices':[1,2],'radial_direction':[0,1,0]}),-7.)

    def test_changed_artifact_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x';p.write_text('original');record=c.pin(p);p.write_text('changed')
            with self.assertRaisesRegex(ValueError,'Changed artifact'):c.verify(record)

    def test_nonfinite_endpoint_rejected(self):
        g={'symbols':['La'],'xyz':{},'charge':3,'physical_multiplicity':1}
        with self.assertRaisesRegex(ValueError,'Nonfinite'):c.add(self.data,'MACE','bad',float('nan'),[[0,0,0]],g,{})

    def test_incomplete_force_rejected(self):
        g={'symbols':['La','O'],'xyz':{},'charge':3,'physical_multiplicity':1}
        with self.assertRaisesRegex(ValueError,'Force dimensions'):c.add(self.data,'MACE','bad',0.,[[0,0,0]],g,{})


if __name__=='__main__':unittest.main()
