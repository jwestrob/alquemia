"""Real canonical/crystal inputs and an explicitly corrupted selector copy."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'workspaces/mace_pqq_utility_20260918/interface_source_v4/implementation'
sys.path.insert(0,str(SOURCE))
from mace_omol_prepared import audit_preparation
from affordable_common import InvalidArtifact,read_json,verify

PANEL=ROOT/'workspaces/mace_omol_20260917/intact_panel_prepared_v2/preparation_manifest.json'


class RealSourceSelectorTests(unittest.TestCase):
    def test_selected_metal_on_separate_chain_is_valid(self):
        r=read_json(PANEL)['rows'][0];p=verify(r['preparation'])
        self.assertEqual(read_json(p)['source_audit_row']['selected_site']['chain'],'B')
        self.assertEqual(audit_preparation(p)['status'],'pass')

    def test_existing_crystal_chain_stays_valid(self):
        r=next(r for r in read_json(PANEL)['rows'] if r['case_id']=='1H4I')
        self.assertEqual(audit_preparation(verify(r['preparation']))['status'],'pass')

    def test_legacy_insertion_code_field_stays_valid(self):
        p=ROOT/'workspaces/mace_global_benchmark_20260916/prepared_v1/GGR_1GLG/preparation.json'
        self.assertEqual(audit_preparation(p)['status'],'pass')

    def test_corrupted_real_selected_residue_is_rejected(self):
        r=read_json(PANEL)['rows'][0];p=copy.deepcopy(read_json(verify(r['preparation'])))
        p['source_audit_row']['selected_site']['resnum']=1000001
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'corrupted_real_selector.json';path.write_text(json.dumps(p))
            with self.assertRaisesRegex(InvalidArtifact,'explicitly selected source metal'):
                audit_preparation(path)


if __name__=='__main__':unittest.main()
