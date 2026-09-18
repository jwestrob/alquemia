"""Collective substitutions of actual pinned multisite structures."""
from pathlib import Path
import copy,json,sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_collective import joint_state,validate
from mace_hybrid import check_atoms
P=ROOT/'workspaces/mace_collective_20260918/prepared_v1/manifest.json'

@unittest.skipUnless(P.exists(),'real collective preparation not available')
class CollectiveTests(unittest.TestCase):
    def test_actual_full_occupancy_and_unchanged_single_site_guard(self):
        m=read_json(P)
        for t in m['tasks']:
            new=xyz(verify(t['xyz']));old=xyz(verify(t['source_xyz']));ids=t['joint_metal_indices']
            self.assertEqual(joint_state(new,t['charge'],ids,t['metal']),t['state'])
            self.assertEqual(len(ids),m['groups'][t['case_id']]['metal_count'])
            if t['variant']=='La_primary':
                self.assertEqual([a[1:] for a in new],[a[1:] for a in old])
                self.assertEqual(t['charge']-m['groups'][t['case_id']]['task']['charge'],len(ids))
                for i,(a,b) in enumerate(zip(new,old)):
                    if i not in ids:self.assertEqual(a,b)
                with self.assertRaises(InvalidArtifact):check_atoms(new,t['charge'])
                broken=list(new);i=ids[0];broken[i]=('Ca',*broken[i][1:])
                with self.assertRaises(InvalidArtifact):joint_state(broken,t['charge'],ids,'La')

    def test_real_manifest_and_corrupted_collective_inventory(self):
        self.assertEqual(validate(P)['tasks'],8);m=read_json(P);broken=copy.deepcopy(m)
        broken['tasks'][0]['joint_metal_indices'].pop()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_manifest.json';p.write_text(json.dumps(broken))
            with self.assertRaises(InvalidArtifact):validate(p)

    @unittest.skipUnless(P.parent.parent.joinpath('report_v1/result.json').exists(),'actual collective outputs not yet reported')
    def test_actual_report_replay(self):
        from mace_collective import report
        old=read_json(P.parent.parent/'report_v1/result.json')
        with tempfile.TemporaryDirectory() as d:
            report(P,Path(d)/'replay');new=read_json(Path(d)/'replay/result.json')
            for key in ('scores','checks','domain_comparisons','parvalbumin_supporting_comparisons','numerical_gate_pass','domain_pass_count'):
                self.assertEqual(old[key],new[key])

if __name__=='__main__':unittest.main()
