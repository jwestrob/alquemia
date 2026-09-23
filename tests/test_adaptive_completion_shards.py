"""Actual410-task manifest partition tests, without optimizer/model execution."""
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import adaptive_completion_folds as folds
from affordable_common import InvalidArtifact,read_json


class AdaptiveCompletionShards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=ROOT/'workspaces/adaptive_completion_20260922'
        cls.old=read_json(cls.root/'primary225_v1/manifest.json')
        cls.new=read_json(cls.root/'primary225_v2_sharded/manifest.json')

    def test_exact_scientific_task_replay(self):
        for key in ('tasks','cases','declared_case_ids','settings','protocol_id','source_manifest',
                    'base_collection','diagnostic','model','software','orca','endpoint_denominator'):
            self.assertEqual(self.old[key],self.new[key],key)
        self.assertEqual(self.new['execution_partition'],{'rule':'task_index_modulo','shards':4})
        for module in ('adaptive_completion.py','adaptive_angular_proposals.py','adaptive_force_diagnostic.py',
                       'mace_site_kinematics.py','accommodation_proposals.py','accommodation_nonlinear.py',
                       'accommodation_torsion_profiles.py','affordable_common.py','mace_hybrid.py','mace_omol.py'):
            self.assertEqual(self.old['implementation'][module]['sha256'],
                             self.new['implementation'][module]['sha256'],module)

    def test_shards_cover_exact410_once_in_frozen_order(self):
        selections=[folds.shard_tasks(self.new,i) for i in range(4)]
        self.assertEqual(list(map(len,selections)),[103,103,102,102])
        self.assertEqual(sorted(i for selected in selections for i,_ in selected),list(range(410)))
        for shard,selected in enumerate(selections):
            for index,task in selected:
                self.assertEqual(index%4,shard)
                self.assertEqual(task,self.new['tasks'][index])
        self.assertEqual(len({t['task_id'] for selected in selections for _,t in selected}),410)
        self.assertEqual(self.new['endpoint_denominator'],450)

    def test_outside_shard_range_rejected_and_old_serial_remains_explicit(self):
        for value in (-1,4):
            with self.assertRaises(InvalidArtifact):folds.shard_tasks(self.new,value)
        self.assertEqual(len(folds.shard_tasks(self.old,0)),410)
        with self.assertRaises(InvalidArtifact):folds.shard_tasks(self.old,1)


if __name__=='__main__':unittest.main()
