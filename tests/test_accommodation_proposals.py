"""Real archived inputs/energies only; no fabricated scientific successes."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import accommodation_proposals as proposal
from affordable_common import read_json, verify, xyz, HA_TO_KCAL
from mace_hybrid import EV_TO_KCAL


class ProposalPreflight(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = ROOT / 'workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json'
        if not cls.manifest.exists():
            raise unittest.SkipTest('real prepared proposal manifest unavailable')
        cls.m = read_json(cls.manifest)

    def test_actual_sixty_origin_sources_and_common_modes(self):
        check = proposal.validate(self.manifest)
        self.assertEqual(check['starts'], 60)
        self.assertEqual(check['q0_available'], 60)
        self.assertEqual(sum(len(t['active_indices']) == 2 for t in self.m['tasks']), 28)
        for t in self.m['tasks']:
            for medium, low in t['q0']['low'].items():
                self.assertEqual(xyz(verify(low['task']['xyz'])), xyz(verify(t['xyz'])))
                self.assertEqual(low['task']['medium'], medium)
            self.assertLessEqual(t['executed_q0_difference_A'], 1e-12)

    def test_existing_score_algebra_exactly_reproduces_released_reference(self):
        reference = read_json(verify(self.m['frozen_comparison']))
        for c in self.m['cases']:
            if c['case_id'] in proposal.PLM:
                continue
            ca, la = [next(t for t in self.m['tasks'] if t['case_id'] == c['case_id'] and t['metal'] == z) for z in ('Ca', 'La')]
            score = proposal.score(ca['q0']['components'], la['q0']['components'])
            old = next(r for r in reference['rows'] if r['case_id'] == c['case_id'] and r['representation'] == 'context')
            self.assertEqual(score['composite_R_model_kcal_mol'], old['composite_R_model_kcal_mol'])

    def test_selection_and_signs_on_all_real_torsion_points(self):
        data = read_json(verify(self.m['primary_result'])); choices = set()
        keys = ('MACE_eV', 'GFN2_ALPB_hartree', 'GFN2_vacuum_hartree')
        for p in data['points']:
            if p['status'] != 'complete':
                continue
            o = next(r for r in data['points'] if r['case_id'] == p['case_id'] and r['metal'] == p['metal'] and r['point'] == 'origin')
            original = {k: o[k] for k in keys}; moved = {k: p[k] for k in keys}
            result = proposal.select_endpoint(original, moved)
            expected = (p['MACE_eV'] - o['MACE_eV']) * EV_TO_KCAL + ((p['GFN2_ALPB_hartree'] - o['GFN2_ALPB_hartree']) - (p['GFN2_vacuum_hartree'] - o['GFN2_vacuum_hartree'])) * HA_TO_KCAL
            self.assertEqual(result['work']['composite_kcal_mol'], expected)
            self.assertEqual(result['selected'], 'proposal' if expected < -.1 else 'origin')
            choices.add(result['selected'])
        self.assertEqual(choices, {'origin', 'proposal'})

    def test_missing_member_never_falls_back(self):
        origin = self.m['tasks'][0]['q0']['components']
        result = proposal.select_endpoint(origin, None)
        self.assertEqual(result['status'], 'unavailable')
        self.assertIsNone(result['selected'])

    def test_explicitly_corrupted_real_source_state_is_rejected(self):
        t = self.m['tasks'][0]; low = t['q0']['low']['vacuum']
        corrupted = dict(low['task']); corrupted['charge'] += 2
        with self.assertRaisesRegex(proposal.InvalidArtifact, 'electronic state differs'):
            proposal.check_low_source(verify(low['receipt']['manifest']), corrupted,
                                      {'xyz': t['xyz'], 'charge': t['charge'], 'multiplicity': 1}, self.m['orca'])

    def test_actual_prelaunch_collection_retains_missing_results(self):
        if any(self.manifest.parent.glob('proposals/*/result.json')):
            self.skipTest('the real proposal run has started; preserve original prelaunch test receipt')
        with tempfile.TemporaryDirectory() as tmp:
            result = proposal.collect(self.manifest, Path(tmp) / 'collection.json')
            self.assertEqual(result['case_denominator'], 30)
            self.assertEqual(result['available_cases'], 0)
            self.assertEqual(result['available_endpoints'], 0)
            self.assertEqual(result['MACE_calls_started'], 0)


if __name__ == '__main__': unittest.main()
