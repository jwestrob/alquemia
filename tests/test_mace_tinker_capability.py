"""Executed real Tinker parameter output; no energy or force substitutes."""
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_tinker_capability import inspect_output,validate

W=ROOT/'workspaces/mace_omol_20260917/tinker_capability_v1'


class TinkerCapability(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'result.json').exists():
            raise unittest.SkipTest('actual native parameter-only outputs unavailable')
        cls.manifest=validate(W/'manifest.json')
        cls.results=read_json(W/'result.json')
        cls.mapping=read_json(verify(cls.manifest['tasks'][0]['mapping']))
        cls.text=verify(cls.results['cases'][0]['native_log']).read_text()

    def test_real_roundtrip_and_observed_native_freeze_bug(self):
        for task, result in zip(self.manifest['tasks'], self.results['cases']):
            mapping=read_json(verify(task['mapping']))
            parsed=inspect_output(verify(result['native_log']).read_text(),mapping)
            self.assertTrue(parsed['coordinates_exact'])
            self.assertTrue(parsed['connectivity_exact'])
            self.assertEqual(parsed['differences']['max_charge_error_e'],0.)
            self.assertEqual(parsed['differences']['axis_mismatches'],0)
            self.assertEqual(parsed['native_keyword_reenabled_sources'],len(mapping['frozen_indices']))
            self.assertTrue(parsed['adapter_mask_pass'])
            self.assertFalse(parsed['full_model_qualified'])

    def test_truncated_actual_native_output_is_not_completion(self):
        text='\n'.join(line for line in self.text.splitlines() if 'ALQUEMIA_PARAMETER_ONLY_COMPLETE' not in line)
        with self.assertRaisesRegex(InvalidArtifact,'completion missing'):
            inspect_output(text,self.mapping)

    def test_corrupted_actual_freeze_flag_rejected(self):
        lines=self.text.splitlines()
        for i,line in enumerate(lines):
            fields=line.split()
            if fields and fields[0]=='FROZEN_AFTER':
                fields[2]='T';lines[i]=' '.join(fields);break
        with self.assertRaisesRegex(InvalidArtifact,'freeze adapter'):
            inspect_output('\n'.join(lines),self.mapping)

    def test_corrupted_real_damping_is_not_a_source_freeze(self):
        lines=self.text.splitlines()
        for i,line in enumerate(lines):
            fields=line.split()
            if fields and fields[0]=='FROZEN_AFTER':
                fields[-1]='0.0';lines[i]=' '.join(fields);break
        with self.assertRaisesRegex(InvalidArtifact,'freeze adapter'):
            inspect_output('\n'.join(lines),self.mapping)


if __name__=='__main__':unittest.main()
