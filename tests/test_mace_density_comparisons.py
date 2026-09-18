"""Frozen real-reference comparison scope, before new full model scores."""
import copy
from pathlib import Path
import sys
import tempfile
import json
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record
from mace_density_comparisons import validate_spec
BASE=ROOT/'workspaces/mace_omol_20260917'

class Comparisons(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=BASE/'multisite_density_comparisons_v1/spec.json'
        cls.cases=list(read_json(BASE/'multisite_density_prepared_v2/preparation.json')['cases'])
        cls.software=read_json(BASE/'trial_density_gk_hybrid_v1/manifest.json')['software']

    def test_all_reference_structures_and_ordered_sites_retained(self):
        s,refs=validate_spec(record(self.path),self.cases,self.software)
        self.assertEqual(len(s['comparisons']),6)
        self.assertEqual({r['left_new_case'] for r in s['comparisons']},{'PARV_4CPV_CD','PARV_4CPV_EF'})
        self.assertEqual({r['right_case'] for r in s['comparisons']},{'GGR_extended','GGR_2FW0','GGR_2FVY'})
        self.assertEqual(s['ordered_vectors']['aequorin'],['AEQ_1SL8_EF1','AEQ_1SL8_EF3','AEQ_1SL8_EF4'])
        self.assertFalse(refs['GGR_transfer']['ordering_pass']) # The inconvenient real reference remains included.

    def test_omitted_real_aequorin_site_rejected(self):
        s=read_json(self.path);s['ordered_vectors']['aequorin'].pop()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_spec.json';p.write_text(json.dumps(s))
            with self.assertRaisesRegex(InvalidArtifact,'omit or duplicate'):validate_spec(record(p),self.cases,self.software)

    def test_changed_real_reference_backend_rejected(self):
        altered=copy.deepcopy(self.software);altered['sha256']='corrupted-real-backend'
        with self.assertRaisesRegex(InvalidArtifact,'backend differs'):validate_spec(record(self.path),self.cases,altered)

if __name__=='__main__':unittest.main()
