"""Real pinned GGR mappings and executed readouts; corrupted copies test rejection."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify
from mace_omol import validate
from mace_omol_readout import physical_mapping, checked_arrays

W=ROOT/'workspaces/mace_omol_20260917'
M=W/'readout_v1/manifest.json'


@unittest.skipUnless(M.exists(),'requires actual prepared four-endpoint GGR readout manifest')
class ReadoutTests(unittest.TestCase):
    def test_four_unchanged_endpoints_preserve_source_atoms_and_caps(self):
        m=read_json(M);parent=read_json(verify(read_json(verify(m['source_collection']))['manifest']))
        self.assertEqual(len(m['tasks']),4)
        for task in m['tasks']:
            original=next(t for t in parent['tasks'] if t['task_id']==task['task_id'])
            self.assertEqual({k:v for k,v in task.items() if k not in ('cache_key','capture_native_readout')},
                             {k:v for k,v in original.items() if k!='cache_key'})
            self.assertNotEqual(task['cache_key'],original['cache_key'])
        maps,shared=physical_mapping(m['tasks'])
        self.assertEqual([len(maps[k]) for k in ('GGR_extended','GGR_connected')],[58,111])
        self.assertEqual([sum(r['kind']=='sigma_link_H' for r in maps[k]) for k in maps],[9,7])
        self.assertEqual(len(shared),52)

    def test_duplicate_source_index_in_corrupted_real_preparation_is_rejected(self):
        tasks=copy.deepcopy(read_json(M)['tasks'])
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);mapping=read_json(verify(tasks[0]['preparation']))
            prep=read_json(verify(mapping['source_preparation']))
            prep['atom_graph']['source_to_qm'][0]['qm_index']=2
            p=d/'corrupted_preparation.json';p.write_text(json.dumps(prep))
            mapping['source_preparation']=record(p)
            p=d/'corrupted_mapping.json';p.write_text(json.dumps(mapping))
            for task in tasks[:2]:task['preparation']=record(p)
            with self.assertRaises(InvalidArtifact):physical_mapping(tasks)

    def test_changed_scientific_task_cannot_reuse_original_cache(self):
        m=read_json(M)
        with tempfile.TemporaryDirectory() as d:
            m['tasks'][0]['charge']+=1
            p=Path(d)/'corrupted_charge.json';p.write_text(json.dumps(m))
            with self.assertRaises(InvalidArtifact):validate(p)

    def test_actual_readouts_sum_to_native_energy(self):
        m=read_json(M)
        for task in m['tasks']:
            paths=sorted((M.parent/'execution'/task['task_id']).glob('attempt_*/result.json'))
            valid=[read_json(p) for p in paths if read_json(p)['status']=='computed']
            if not valid:self.skipTest('actual scientific readout execution unavailable')
            arrays=checked_arrays(valid[-1],task)
            self.assertTrue(all(np.isfinite(a).all() for a in arrays.values()))
            damaged=copy.deepcopy(valid[-1]);damaged['native_readout']['arrays']['sha256']='0'*64
            with self.assertRaises(InvalidArtifact):checked_arrays(damaged,task)


if __name__=='__main__':unittest.main()
