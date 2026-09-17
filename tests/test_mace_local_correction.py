"""Actual core/source and archived DFT invariants; no generated scientific data."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,energy,HA_TO_KCAL,read_json,verify,xyz
from mace_local_correction import source_mapping,validate,collect_local
MANIFEST=ROOT/'workspaces/mace_local_correction_20260916/mace_v2/manifest.json'


@unittest.skipUnless(MANIFEST.exists(),'requires prepared real local-core panel')
class LocalCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=read_json(MANIFEST)

    def test_all_exact_archives_and_source_H_only_changes(self):
        self.assertEqual(validate(MANIFEST)['tasks'],20)
        for case,ref in self.m['cases'].items():
            data=read_json(verify(ref));prep=read_json(verify(data['global_preparation']))
            old=xyz(verify(data['source_endpoints']['La']['xyz']))
            new,mapping,moves=source_mapping(prep,old)
            full={a['id']:a for a in prep['physical_atoms']}
            changed={m['core_index'] for m in moves}
            for i,(a,b,m) in enumerate(zip(old,new,mapping)):
                if i not in changed:self.assertEqual(a,b)
                else:
                    self.assertEqual(a[0],'H');self.assertEqual(tuple(full[m['physical_id']]['xyz_A']),b[1:])
            caps={m['core_index'] for m in mapping if m['kind']=='synthetic_cap'}
            self.assertTrue(caps.isdisjoint(changed))
            for t in self.m['tasks']:
                if t['case_id']==case and t['local_state']=='archived':
                    self.assertEqual(t['xyz'],data['source_endpoints'][t['metal']]['xyz'])

    def test_real_archived_DFT_units_and_new_geometry_unavailability(self):
        for ref in self.m['cases'].values():
            d=read_json(verify(ref));end=d['archived_DFT']
            la,ca=(energy(verify(end[m]['output'])) for m in ('La','Ca'))
            self.assertEqual((ca-la)*HA_TO_KCAL,d['archived_DFT_R_kcal_mol'])
            self.assertIsNone(d['global_H_DFT']);self.assertEqual(d['global_H_DFT_status'],'not_computed')
        d=read_json(verify(self.m['cases']['PQQ_4MAE']))
        self.assertAlmostEqual(d['archived_DFT_R_kcal_mol'],-405381.6912148747,places=8)

    def test_corrupted_real_missing_source_rejected(self):
        d=read_json(verify(self.m['cases']['GGR_1GLG']));prep=read_json(verify(d['global_preparation']))
        core=xyz(verify(d['source_endpoints']['La']['xyz']))
        removed=d['mapping'][1]['physical_id']
        prep['physical_atoms']=[a for a in prep['physical_atoms'] if a['id']!=removed]
        with self.assertRaisesRegex(InvalidArtifact,'missing'):source_mapping(prep,core)

    def test_missing_results_and_changed_model_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'manifest.json';p.write_bytes(MANIFEST.read_bytes())
            c=collect_local(p);self.assertEqual(c['status'],'incomplete');self.assertIsNone(c['hybrid'])
            for s in c['scores'].values():
                self.assertIsNone(s['environment_contribution_kcal_mol']);self.assertIsNone(s['geometry_effect_kcal_mol'])
            bad=copy.deepcopy(self.m);bad['model']['solvent']='water';p.write_text(json.dumps(bad))
            with self.assertRaisesRegex(InvalidArtifact,'electronic model'):validate(p)


if __name__=='__main__':unittest.main()
