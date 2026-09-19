"""Real archived GGR structures, unchanged expansion policy; no dummy energies."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,xyz,record
from second_shell_context import POLICY
from environment_replica_context import validate,endpoint


class StaticReplicas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=ROOT/'workspaces/environment_replicas_20260919/prepared_v1/manifest.json'
        cls.manifest=read_json(cls.path)

    def test_unchanged_real_source_expansion(self):
        result=validate(self.path)
        self.assertEqual(result['tasks'],4)
        self.assertEqual(self.manifest['policy'],POLICY)
        self.assertEqual([r['expanded_atoms'] for r in result['inventory']],[115,115])
        self.assertEqual([r['added_formal_charge'] for r in result['inventory']],[-1,-1])

    def test_preserved_original_cores_and_actual_energy_receipts(self):
        for s in self.manifest['states']:
            prep=read_json(verify(s['preparation']))
            for metal,old in s['source']['endpoints'].items():
                t=next(t for t in self.manifest['tasks'] if t['case']==s['case'] and t['metal']==metal)
                original=xyz(verify(old['xyz']));expanded=xyz(verify(t['xyz']))
                for i,j in prep['core_to_context'].items():self.assertEqual(original[int(i)],expanded[j])
                actual=endpoint(old['output'],old['receipt'],old['xyz'],old['input'])
                self.assertEqual(old['energy_hartree'],actual['energy_hartree'])

    def test_paired_geometry_and_explicit_unavailable_reference(self):
        self.assertIsNone(self.manifest['reference']);self.assertIsNone(self.manifest['calibrated_decision'])
        for case in ('2FW0','2FVY'):
            tasks={t['metal']:t for t in self.manifest['tasks'] if t['case']==case}
            ca,la=[xyz(verify(tasks[m]['xyz'])) for m in ('Ca','La')]
            self.assertEqual(ca[1:],la[1:]);self.assertEqual(ca[0][1:],la[0][1:])
            self.assertEqual(tasks['La']['charge']-tasks['Ca']['charge'],1)


if __name__=='__main__':unittest.main()
