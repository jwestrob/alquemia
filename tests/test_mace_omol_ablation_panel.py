"""Conditional extension guard tests using actual saved scientific artifacts."""
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact
from mace_file_checks import cached_file_checks
from mace_omol_ablation_panel import inputs

WORK=ROOT/'workspaces/mace_omol_20260917'
PREP=WORK/'intact_panel_prepared_v2/preparation_manifest.json'


@unittest.skipUnless(PREP.exists(),'requires actual strict whole-chain preparation')
class AblationPanelTests(unittest.TestCase):
    @cached_file_checks
    def test_native_success_does_not_authorize_modified_model_extension(self):
        # A real successful native qualification cannot qualify a changed model.
        with self.assertRaisesRegex(InvalidArtifact,'own completed development protocol'):
            inputs(PREP,WORK/'product_intact_v1/collection_job_1200818.json')

    def test_legacy_preparation_cannot_silently_replace_strict_version(self):
        with self.assertRaisesRegex(InvalidArtifact,'strict v2 preparation'):
            inputs(WORK/'intact_panel_prepared_v1/preparation_manifest.json',
                   WORK/'product_intact_v1/collection_job_1200818.json')

    @cached_file_checks
    def test_partial_descriptor_collection_cannot_unlock_extension(self):
        path=WORK/'charge_ablation_development_v2/progress_first.json'
        if not path.exists():self.skipTest('requires actual saved incomplete collection')
        with self.assertRaisesRegex(InvalidArtifact,'all actual ablation development gates'):
            inputs(PREP,path)


if __name__=='__main__':unittest.main()
