"""Actual mapped multisite runner inputs and reusable native parameter outputs."""
import copy
from pathlib import Path
import sys
import tempfile
import json
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_density_frameworks import validate as frameworks
from mace_density_short import validate as short
from mace_hybrid import check_atoms,dry_run
BASE=ROOT/'workspaces/mace_omol_20260917'

class MultisiteRunners(unittest.TestCase):
    def test_real_native_framework_import_keeps_every_source_atom(self):
        m=frameworks(BASE/'multisite_density_frameworks_v1/manifest.json')
        self.assertEqual(m['reused_native_parameter_reads'],5);self.assertEqual(m['new_parameter_reads'],0)
        for t in m['tasks']:
            mapping=read_json(verify(t['mapping']));old=read_json(verify(mapping['source_mapping']))
            alias=mapping['selected_metal_identity_alias'];source=copy.deepcopy(old['physical_atoms'])
            selected=next(a for a in source if a['id']==alias['source_id']);selected['id']='metal'
            self.assertEqual(source,mapping['physical_atoms'])
            self.assertEqual(old['system_atom_ids'],mapping['system_atom_ids'])
            receipt=read_json(verify(t['preflight_receipt']));self.assertEqual(receipt['returncode'],0);verify(receipt['log'])

    def test_all_twenty_short_inputs_have_explicit_background_ions(self):
        path=BASE/'multisite_density_short_v2/manifest.json';self.assertEqual(dry_run(path)['tasks'],20)
        m=read_json(path)
        for t in m['tasks']:
            c=read_json(verify(t['source_mapping']));p=read_json(verify(c['normalized_global_preparation']))
            expected=[i for i,a in enumerate(p['physical_atoms']) if a['kind']=='background_metal'] if t['kind']=='full' else []
            self.assertEqual(t.get('background_calcium_indices',[]),expected)
            check_atoms(xyz(verify(t['xyz'])),t['charge'],background_calcium_indices=expected)

    def test_missing_background_mapping_rejected_before_execution(self):
        m=read_json(BASE/'multisite_density_short_v2/manifest.json');t=next(t for t in m['tasks'] if t['kind']=='full')
        del t['background_calcium_indices'] # Explicitly corrupted real manifest.
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_manifest.json';p.write_text(json.dumps(m))
            with self.assertRaisesRegex(InvalidArtifact,'background calcium'):short(p)

if __name__=='__main__':unittest.main()
