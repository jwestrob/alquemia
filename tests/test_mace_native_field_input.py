"""Pinned native-source transformation and real response replay evidence."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_native_field_input import derive,validate,accepted,parse

W=ROOT/'workspaces/mace_omol_20260917/native_field_input_v1'
S=ROOT/'workspaces/mace_omol_20260917/native_field_input_software_v2'


class NativeTransformation(unittest.TestCase):
    def test_only_declared_native_source_changes(self):
        if not (S/'receipt.json').exists():self.skipTest('actual native adapter build unavailable')
        r=read_json(S/'receipt.json');original,derived=derive(verify(r['original_full_source']).read_text())
        self.assertEqual(original,verify(r['original_routine']).read_text())
        self.assertEqual(derived,verify(r['derived_routine']).read_text())
        restored=derived.replace('      subroutine alquemia_induce0c (source_fields)\n','      subroutine induce0c\n',1)
        restored=restored.replace('      real*8 source_fields(3,n,4)\n','',1)
        restored=restored.replace('      field(:,:) = source_fields(:,:,1)\n      fieldp(:,:) = source_fields(:,:,2)\n      fields(:,:) = source_fields(:,:,3)\n      fieldps(:,:) = source_fields(:,:,4)',
            '      call dfield0d (field,fieldp,fields,fieldps)',1)
        self.assertEqual(restored,original)


class RealReplay(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('native input preparation unavailable')
        cls.m=validate(W/'manifest.json')

    def test_supplied_fields_are_exact_actual_native_arrays(self):
        for t in self.m['tasks']:
            source=read_json(verify(t['field_result']))
            lines=verify(t['supplied_fields']).read_text().splitlines();self.assertEqual(int(lines[0]),len(lines)-1)
            for i,line in enumerate(lines[1:],1):
                f=line.split();self.assertEqual(int(f[0]),i)
                self.assertEqual(list(map(float,f[1:])),source['fields_e_A2'][str(i)])

    def test_native_responses_and_corrupted_completion(self):
        items=[(t,accepted(t,W/'tasks'/t['task_id'])) for t in self.m['tasks']]
        if any(not r for _,r in items):self.skipTest('six actual native solves incomplete')
        for t,(r,_) in items:
            text=verify(r['receipt']['log']).read_text();p=parse(text,t)
            self.assertIsNone(p['numerical_score'])
            self.assertEqual(p['solver_summaries'],r['solver_summaries'])
            self.assertTrue(all(np.all(np.array(p['response_eA'][i][1:])==0) for i in t['frozen_indices']))
            with self.assertRaisesRegex(InvalidArtifact,'incomplete/nonconverged'):
                parse(text.replace('ALQUEMIA_RESPONSE_COMPLETE','DAMAGED_COMPLETION'),t)

    def test_actual_report_retains_failures_or_validated_replay(self):
        reports=sorted(W.glob('collection_job_*.json'))
        if not reports:self.skipTest('actual response report unavailable')
        r=read_json(reports[-1]);self.assertIsNone(r['numerical_score']);self.assertFalse(r['full_model_qualified'])
        self.assertEqual(r['gates_pass'],r['complete'] and len(r['checks'])==9 and all(c['pass_'] for c in r['checks']))


if __name__=='__main__':unittest.main()
