"""Real pinned core/ORCA fixtures; no fabricated scientific results."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify
from mace_omol_vacuum import collect, parse_endpoint, scientific_input, validate

WORK=ROOT/'workspaces/mace_omol_20260917'
MANIFEST=WORK/'matched_vacuum_v1/manifest.json'


@unittest.skipUnless(MANIFEST.exists(),'real prepared vacuum fixtures unavailable')
class MatchedVacuum(unittest.TestCase):
    def test_only_solvent_removed_with_identical_paired_sources(self):
        self.assertEqual(validate(MANIFEST)['tasks'],8)
        m=read_json(MANIFEST)
        for t in m['tasks']:
            old=verify(t['CPCM_gradient']['artifacts']['input']).read_text()
            self.assertEqual(verify(t['input']).read_text(),old.replace('CPCM(Water) ',''))
            self.assertEqual(t['xyz']['sha256'],t['source_xyz']['sha256'])

    def test_real_CPCM_output_cannot_be_certified_as_vacuum(self):
        t=read_json(MANIFEST)['tasks'][0];a=t['CPCM_gradient']['artifacts']
        with self.assertRaisesRegex(InvalidArtifact,'contains solvent'):
            parse_endpoint(t,verify(a['output']),verify(a['engrad']))

    def test_corrupted_real_recipe_is_rejected(self):
        t=dict(read_json(MANIFEST)['tasks'][0]);a=t['CPCM_gradient']['artifacts']
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'corrupted_real_input.inp'
            p.write_text(verify(t['input']).read_text().replace('EnGrad','NumGrad'))
            t['input']=record(p)
            with self.assertRaisesRegex(InvalidArtifact,'exact vacuum EnGrad'):
                parse_endpoint(t,verify(a['output']),verify(a['engrad']))

    def test_actual_vacuum_scope_and_component_closure(self):
        p=WORK/'matched_vacuum_report_v1/result.json'
        if not p.exists():self.skipTest('new scientific outputs have not been executed/collected')
        r=read_json(p);actual=collect(MANIFEST)
        self.assertEqual(actual['status'],'complete')
        self.assertEqual(actual['rows'],r['rows'])
        self.assertTrue(all(v['pass'] for v in r['checks']))
        self.assertEqual(len(r['rows']),8)
        self.assertTrue(all(v['energy_scope']=='isolated_vacuum_endpoint' for v in r['rows'].values()))
        self.assertIsNone(r['aqueous_score'])
        self.assertFalse(r['conditional_whole_inference_eligible'])


if __name__=='__main__':unittest.main()
