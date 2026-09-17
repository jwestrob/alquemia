"""Pinned physical preparations and actual outputs, without mock energies."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,write_new
from mace_omol_vacuum_hybrid import collect,validate,PROTOCOL

W=ROOT/'workspaces/mace_omol_20260917'
M=W/'vacuum_hybrid_v1/manifest.json'


@unittest.skipUnless(M.exists(),'real vacuum-hybrid preparation unavailable')
class VacuumHybrid(unittest.TestCase):
    def test_exact_six_original_H_endpoints(self):
        self.assertEqual(validate(M)['tasks'],6)
        m=read_json(M)
        self.assertEqual(m['protocol_id'],PROTOCOL)
        self.assertEqual({t['case_id'] for t in m['tasks']},{'GGR_1GLG','ALPHA_1F6S','ALPHA_6IP9'})
        self.assertTrue(all(t['energy_only'] and t['capture_native_readout'] for t in m['tasks']))

    def test_corrupted_real_charge_changes_cache_and_is_rejected(self):
        m=copy.deepcopy(read_json(M));m['tasks'][0]['charge']+=2
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'corrupted_actual_manifest.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'task or cache differs'):
                validate(p)

    def test_actual_results_preserve_separate_partition_and_ordering_gates(self):
        p=W/'vacuum_hybrid_report_v1/result.json'
        if not p.exists():self.skipTest('actual scientific inference not yet executed/reported')
        r=read_json(p);actual=collect(M)
        for k,v in actual.items():self.assertEqual(r[k],v)
        self.assertTrue(r['numerical_gate_pass'])
        self.assertEqual(len(r['rows']),6)
        self.assertEqual(len(r['contrasts']),4)
        self.assertEqual(r['ordering_gate_pass'],all(x['difference_kcal_scale']>.02 for x in r['contrasts']))
        self.assertIsNone(r['calibrated_class'])
        self.assertIsNone(r['aqueous_score'])


if __name__=='__main__':unittest.main()
