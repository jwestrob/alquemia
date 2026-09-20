"""Real prepared PQQ fixtures and independent archived donor measurements."""
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from accommodation_inventory import geometry
from affordable_common import InvalidArtifact

REFERENCE = ROOT / 'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/prepared/01_a0a3f2yly8-pqq-la_model/pmdh_fc_a0a3f2yly8-pqq-la_model_carve_manifest.json'
PLM = ROOT / 'workspaces/plm_xoxf_all_fixed_core_20260916/retry_1200785/prepared_batch/prepared/PQQSEQ_026611884af742b7c6cc_AF3_sample1_carve_manifest.json'


class GeometryTests(unittest.TestCase):
    def test_both_real_naming_conventions_match_archived_direct_donors(self):
        for path in (REFERENCE, PLM):
            m = json.loads(path.read_text()); result = geometry(path)
            key = 'prepared_typed_direct_donors_ordered_by_distance' if path == PLM else 'source_typed_direct_donors_ordered_by_distance'
            donors = m['coordination'][key]
            name = 'N6' if path == PLM else 'N2'
            archived = next(d['distance_A'] for d in donors if d['resname'] == 'PQQ' and d['atom'] == name)
            self.assertLess(abs(result['features']['pqq_cofactor.N6'] - archived), .00051)
            self.assertEqual(result['charge_La'] - result['charge_Ca'], 1)

    def test_corrupted_copy_of_real_mapping_rejected(self):
        m = json.loads(PLM.read_text())
        m['qm_fragments'][0]['atom_records'][0]['xyz_A'][0] += .1
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / 'explicitly_corrupted_real_carve.json'; p.write_text(json.dumps(m))
            with self.assertRaisesRegex(InvalidArtifact, 'mapping'):
                geometry(p)

    def test_extreme_contact_from_actual_prepared_coordinates(self):
        p = PLM.parent / 'PQQSEQ_83440678cbbd658047c9_AF3_sample0_carve_manifest.json'
        r = geometry(p)
        self.assertAlmostEqual(r['features']['extra_acidic_ligand_homolog.nearest'], 1.766154296770245)
        self.assertEqual(r['formula'], {'C': 27, 'H': 31, 'N': 6, 'O': 15})


if __name__ == '__main__': unittest.main()
