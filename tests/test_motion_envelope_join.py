"""Actual strict origins and declared same-source baselines; zero molecular calls."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,write_new
import motion_envelope_pool as pool
import motion_envelope_compare as compare
BASE=ROOT/'workspaces/motion_envelope_20260923'
class EnvelopeJoin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.origins=BASE/'pilot34_scalar_origins_v1/COLLECTION.json'
        cls.source=read_json(BASE/'pilot34_searches_v1/manifest.json')
    def test_actual_all68_join_with_missing_composite_preserved_in_source(self):
        joined=pool.joined_origins(self.origins,self.source)
        self.assertEqual(len(joined),34)
        self.assertTrue(all(c['status']=='complete'for row in joined.values()for c in row.values()))
        self.assertTrue(all(t['q0']['low']is None for t in self.source['tasks']))
        for row in joined.values():
            for c in row.values():self.assertEqual(set(c['components']),{'MACE_eV','GFN2_vacuum_hartree','GFN2_ALPB_hartree'})
    def test_corrupted_actual_scalar_collection_rejected(self):
        bad=read_json(self.origins);bad['rows'][0]['energy_hartree']+=1
        with tempfile.TemporaryDirectory()as td:
            path=Path(td)/'corrupted.json';write_new(path,bad)
            with self.assertRaisesRegex(InvalidArtifact,'differs from actual output'):pool.joined_origins(path,self.source)
    def test_all34_actual_original_core_identities(self):
        b=compare.strict_baselines(ROOT/'workspaces/strict_native_pool_20260923/run_v1/COMPARISON.json',ROOT/'workspaces/strict_native_comparator_20260923/run_v1/COLLECTION.json')
        native=read_json(BASE/'pilot34_origins_v1/NATIVE_COLLECTION_v1.json');inventory=read_json(verify(native['inventory']));count=0
        for r,s in zip(native['rows'],inventory['pilot']):
            src=r['source'];cid=src['source'].get('root_case_id',s['case_id'])if s['role']=='canonical_calibration'else s['case_id']
            self.assertEqual(compare.core_match(src,b[cid]),{'Ca':0.,'La':0.});count+=1
        self.assertEqual(count,34)
if __name__=='__main__':unittest.main()
