"""Frozen real Hans preparation and explicitly corrupted copies; no fake energies."""
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify
from metal_environment_lanm_reference import check_config
from metal_environment_electronic_state import describe

CONFIG=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2/INPUTS.json'
CONFIG_SHA='c5886704ee541752ec2adfe7ada0b5f8b2c0df6014d5ed8676b74447d1141c3f'

class ActualLanMPreparation(unittest.TestCase):
    def setUp(self):
        self.assertEqual(hashlib.sha256(CONFIG.read_bytes()).hexdigest(),CONFIG_SHA)
        self.config=read_json(CONFIG)

    def test_actual_preparation_and_states(self):
        c=check_config(self.config)
        self.assertEqual(c['source_id'],'Hans8DQ2')
        for label in ('A','B'):
            for metal,explicit in [('La',806),('Dy',833)]:
                e=c['configurations'][label]['endpoints'][metal]
                state=describe(verify(e['xyz']),metal,e['charge'],e['multiplicity'])
                self.assertEqual(state['explicit_electrons'],explicit)
        self.assertEqual(c['configurations']['A']['pointcharges'],c['configurations']['B']['pointcharges'])

    def test_singlet_dy_rejected(self):
        self.config['configurations']['B']['endpoints']['Dy']['multiplicity']=1
        with self.assertRaises(InvalidArtifact):check_config(self.config)

    def test_charge_parity_rejected(self):
        self.config['configurations']['A']['endpoints']['Dy']['charge']=0
        with self.assertRaises(InvalidArtifact):check_config(self.config)

    def test_false_state_metadata_rejected(self):
        self.config['configurations']['A']['endpoints']['Dy']['all_electron_count']=852
        with self.assertRaisesRegex(InvalidArtifact,'metadata'):check_config(self.config)

    def test_corrupt_paired_coordinates_rejected(self):
        endpoint=self.config['configurations']['B']['endpoints']['Dy']
        lines=verify(endpoint['xyz']).read_text().splitlines()
        fields=lines[3].split();fields[1]=str(float(fields[1])+.01);lines[3]=' '.join(fields)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'explicitly_corrupt_Dy_B.xyz';p.write_text('\n'.join(lines)+'\n')
            endpoint['xyz']=record(p)
            with self.assertRaisesRegex(InvalidArtifact,'paired coordinates'):check_config(self.config)

    def test_corrupt_environment_charge_inventory_rejected(self):
        conf=self.config['configurations']['B'];lines=verify(conf['pointcharges']).read_text().splitlines()
        fields=lines[1].split();fields[0]=str(float(fields[0])+.01);lines[1]=' '.join(fields)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'explicitly_corrupt_environment.pc';p.write_text('\n'.join(lines)+'\n')
            conf['pointcharges']=record(p)
            with self.assertRaisesRegex(InvalidArtifact,'charge inventory'):check_config(self.config)

    def test_changed_file_without_new_pin_rejected(self):
        conf=self.config['configurations']['B'];original=conf['pointcharges']
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'explicitly_corrupt_environment.pc';p.write_text(verify(original).read_text()+'\n')
            conf['pointcharges']={**original,'path':str(p)}
            with self.assertRaisesRegex(InvalidArtifact,'changed artifact'):check_config(self.config)

if __name__=='__main__':unittest.main()
