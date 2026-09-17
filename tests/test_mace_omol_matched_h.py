"""Real source H normalization and actual matched quantum/ML outputs."""
import copy
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_omol_matched_h import (COUNTS,normalized_physical,validate_prepared,
                                validate_quantum,validate_mace,collect_quantum,collect_mace)
W=ROOT/'workspaces/mace_omol_20260917'
P=W/'matched_H_prepared_v1/preparation.json'
Q=W/'matched_H_quantum_v1/manifest.json'
M=W/'matched_H_mace_v1/manifest.json'


@unittest.skipUnless(P.exists(),'real normalized H preparation unavailable')
class MatchedH(unittest.TestCase):
    def test_only_source_H_changes_and_all_caps_heavy_atoms_survive(self):
        self.assertEqual(validate_prepared(P)['cases'],4)
        p=read_json(P)
        for name,pin in p['cases'].items():
            c=read_json(verify(pin))
            for metal,e in c['endpoints'].items():
                old=xyz(verify(e['original_xyz']));new=xyz(verify(e['xyz']))
                caps={a['qm_index'] for a in c['mapping'] if a['kind']=='sigma_link_H'}
                moved={a['qm_index'] for a in e['H_moves']}
                self.assertEqual(len(moved),COUNTS[name])
                self.assertFalse(caps&moved)
                self.assertEqual([a[0] for a in old],[a[0] for a in new])
                for i,(a,b) in enumerate(zip(old,new)):
                    if a[0]!='H' or i in caps:self.assertEqual(a,b)

    def test_actual_whole_reuse_and_duplicate_selection(self):
        p=read_json(P)
        for case,ends in p['whole_reuses'].items():
            for metal,ref in ends.items():
                self.assertEqual(ref['selected']['task_id'],min(ref['all_matching_task_ids']))
                self.assertLessEqual(ref['duplicate_spread_model_kcal'],.01)
                self.assertEqual(ref['selected']['result']['status'],'computed')
        self.assertEqual(sum(len(v) for v in p['whole_reuses'].values()),6)

    def test_corrupted_real_heavy_coordinate_is_rejected(self):
        p=read_json(P);c=read_json(verify(p['cases']['GGR_extended']))
        prep=copy.deepcopy(read_json(verify(c['normalized_global_preparation'])))
        original=read_json(verify(c['physical_atoms']))
        prep['physical_atoms'][0]['xyz_A'][0]+=.001
        with self.assertRaisesRegex(InvalidArtifact,'heavy/metal coordinate changed'):
            normalized_physical(prep,original)

    def test_quantum_and_MACE_share_all_eight_exact_coordinates(self):
        if not Q.exists() or not M.exists():self.skipTest('manifest preparation not complete')
        self.assertEqual(validate_quantum(Q)['tasks'],8)
        self.assertEqual(validate_mace(M)['tasks'],8)
        q=read_json(Q);m=read_json(M);qm={t['task_id']:t for t in q['tasks']}
        for t in m['tasks']:
            self.assertEqual(t['xyz']['sha256'],qm[t['task_id']]['xyz']['sha256'])
            self.assertEqual(t['charge'],qm[t['task_id']]['charge'])

    def test_actual_scientific_outputs_and_independent_gates(self):
        path=W/'matched_H_report_v1/result.json'
        if not path.exists():self.skipTest('actual quantum/ML outputs are not yet executed/reported')
        r=read_json(path)
        self.assertEqual(r['quantum'],collect_quantum(Q));self.assertEqual(r['MACE'],collect_mace(M))
        self.assertEqual(r['status'],'complete');self.assertTrue(r['numerical_gate_pass'])
        self.assertEqual(len(r['contrasts']),4)
        self.assertEqual(r['ordering_gate_pass'],all(x['difference_kcal_scale']>.02 for x in r['contrasts']))
        self.assertEqual(r['partition_gate_pass'],abs(r['partition_shift_kcal_scale'])<=2.)
        self.assertIsNone(r['aqueous_score']);self.assertFalse(r['baseline_changed'])


if __name__=='__main__':unittest.main()
