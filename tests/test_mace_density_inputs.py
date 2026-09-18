"""Real source replay and explicitly corrupted copies; no scientific energies."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_density_inputs import normalize,environment,pointcharges
from mace_omol_charges import projection
from mace_omol_solvent import forcefield_background

BASE=ROOT/'workspaces/mace_omol_20260917'
OLD=BASE/'matched_H_prepared_v1/GGR_extended/mapping.json'


class RealDensityInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old=read_json(OLD)
        cls.state=read_json(BASE/'full_boundary_GB_v1/states/GGR_extended/state.json')
        cls.physical=verify(cls.old['normalized_global_preparation'])
        cls.core=verify(cls.state['ligand_ledger_source'])

    def test_reproduces_executed_1glg_coordinates_and_generating_field(self):
        c=normalize(self.physical,self.core);c['physical_atoms']=self.old['physical_atoms']
        for metal in ('Ca','La'):
            self.assertEqual(c['endpoints'][metal]['atoms'],xyz(verify(self.old['endpoints'][metal]['xyz'])))
        p=projection(c,c['endpoints']['Ca']['atoms'])
        s=environment(c,p,forcefield_background(read_json(self.physical)))
        for key in ('physical_atoms','projection_support_ids','environment_charges_e','environment_charge','boundary_ledger','ligand_ledger'):
            self.assertEqual(s[key],self.state[key],key)
        oldq=read_json(BASE/'responsive_quantum_v2/manifest.json')
        task=next(t for t in oldq['tasks'] if t['task_id']=='GGR_extended_Ca')
        self.assertEqual(pointcharges(s),verify(task['pointcharges']).read_text())

    def test_real_alternative_structures_keep_paired_geometry_and_source_heavy_atoms(self):
        for tag in ('2fw0','2fvy'):
            pp=BASE/f'ggr_source_bridge_{tag}_v1/preparation.json'
            cp=ROOT/f'workspaces/ggr_mechanism_20260915/stage_b_prepared_v1/{tag}/alpha_caps/preparation_manifest.json'
            c=normalize(pp,cp);a,b=(c['endpoints'][m]['atoms'] for m in ('Ca','La'))
            self.assertEqual(a[1:],b[1:]);self.assertEqual(a[0][1:],b[0][1:])
            old=xyz(verify(read_json(cp)['outputs']['Ca']['xyz']))
            self.assertEqual([x for x in a if x[0]!='H'],[x for x in old if x[0]!='H'])
            self.assertTrue(c['geometry_checks']['chemical_state_unchanged'])

    def test_explicitly_corrupted_duplicate_mapping_fails(self):
        c=read_json(self.core);c['atom_graph']['source_to_qm'][1]['qm_index']=c['atom_graph']['source_to_qm'][0]['qm_index']
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_mapping.json';p.write_text(json.dumps(c))
            with self.assertRaisesRegex(InvalidArtifact,'duplicate or invalid'):normalize(self.physical,p)

    def test_explicitly_corrupted_unsupported_cofactor_fails(self):
        p=read_json(self.physical);p['physical_atoms'][0]['kind']='unsupported_cofactor'
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'corrupted_real_inventory.json';path.write_text(json.dumps(p))
            with self.assertRaisesRegex(InvalidArtifact,'unsupported nonstandard'):normalize(path,self.core)

    def test_explicitly_corrupted_missing_charged_donor_fails(self):
        c=normalize(self.physical,self.core);c['physical_atoms']=self.old['physical_atoms']
        p=projection(c,c['endpoints']['Ca']['atoms'])
        c['mapping']=[a for a in c['mapping'] if a.get('physical_id')!='A/134//OD1']
        with self.assertRaisesRegex(InvalidArtifact,'incomplete charged source'):
            environment(c,p,forcefield_background(read_json(self.physical)))


if __name__=='__main__':unittest.main()
