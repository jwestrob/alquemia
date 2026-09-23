"""Checks of actual saved 250-source union preparations; no scientific calls."""
from collections import defaultdict
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from affordable_common import read_json, verify, xyz
from consistent_context import fragment_id, same_geometry


class ConsistentContextFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prep = read_json(ROOT/'workspaces/consistent_context_20260922/prepared_v1/PREPARATION.json')
        cls.old = read_json(verify(cls.prep['source_preparation']))

    def test_full_denominator_and_same_state_without_coverage_loss(self):
        self.assertEqual((self.prep['denominator'], self.prep['supported']), (250, 233))
        self.assertEqual((self.prep['pilot_denominator'], self.prep['pilot_supported']), (61, 50))
        old_by_id = {r['case_id']: r for r in self.old['cases']}; groups = defaultdict(list)
        for row in self.prep['cases']:
            self.assertEqual(row['status'] == 'prepared', old_by_id[row['case_id']]['status'] == 'prepared')
            groups[row['source']['root_case_id']].append(row)
        self.assertEqual(len(groups), 25)
        for rows in groups.values():
            self.assertEqual(len(rows), 10)
            for z in ('Ca', 'La'): self.assertEqual(sum(r['source']['source_conditioning_metal'] == z for r in rows), 5)
            admitted = [r for r in rows if r['status'] == 'prepared']
            self.assertEqual(len({r['state_key'] for r in admitted}), 1)
            self.assertEqual(len({r['protein_key'] for r in admitted}), 1)
            self.assertEqual(sum(r['source']['canonical_coordinate_match'] for r in admitted), 1)

    def test_union_is_every_prior_supported_fragment_without_score_selection(self):
        for group in self.prep['groups']:
            actual = read_json(verify(group['union'])); wanted = set(); contributors = defaultdict(set)
            for row in self.old['cases']:
                if row['source']['root_case_id'] != group['root_case_id'] or row['status'] != 'prepared': continue
                audit = read_json(verify(row['representations']['context']['preparation']))
                for fragment in audit['added_fragments']:
                    key = fragment_id(fragment); wanted.add(key); contributors[key].add(row['case_id'])
            self.assertEqual({fragment_id(f) for f in actual['fragments']}, wanted)
            self.assertEqual({fragment_id(r['fragment']): set(r['source_cases']) for r in actual['contributors']}, dict(contributors))

    def test_all_existing_context_replays_and_original_core_coordinates(self):
        for row in self.prep['cases']:
            if row['status'] != 'prepared': continue
            self.assertLessEqual(row['original_selection_replay_max_A'], 1e-12)
            rep = row['representations']['context']; audit = read_json(verify(rep['preparation']))
            self.assertFalse(audit['new_protons_or_waters'])
            self.assertEqual(audit['paired_differing_nonmetal_indices'], [])
            endpoints = {z: xyz(verify(rep['endpoints'][z]['xyz'])) for z in ('Ca', 'La')}
            self.assertEqual(endpoints['Ca'][1:], endpoints['La'][1:])
            for z in ('Ca', 'La'):
                old = xyz(verify(row['original_core']['endpoints'][z]['xyz']))
                for oi, ni in audit['core_to_context'].items(): self.assertEqual(old[int(oi)], endpoints[z][ni])
            sig = read_json(verify(row['state_signature']))
            self.assertEqual(sig, read_json(verify(row['canonical_signature'])))
            self.assertEqual(sig['endpoint_charges']['La']-sig['endpoint_charges']['Ca'], 1)

    def test_three_singleton_crystals_are_identity_controls(self):
        m = read_json(ROOT/'workspaces/consistent_context_20260922/crystal_controls_v1/CRYSTALS.json')
        self.assertEqual({r['case_id'] for r in m['cases']}, {'1H4I', '4MAE', '1KB0'})
        for row in m['cases']:
            self.assertEqual(row['status'], 'prepared')
            self.assertEqual(row['atom_count'], row['original_atom_count'])
            self.assertEqual(row['cap_count'], row['original_caps'])
            self.assertEqual(row['new_La_charge'], row['original_La_charge'])
            for z in ('Ca', 'La'):
                self.assertTrue(same_geometry(xyz(verify(row['representations']['context']['endpoints'][z]['xyz'])),
                                              xyz(verify(row['old_context']['endpoints'][z]['xyz']))))

    def test_actual_reuse_comparison_rejects_displacement_and_atom_swap(self):
        # Explicitly corrupted in-memory copies of an actual crystal input.
        m = read_json(ROOT/'workspaces/consistent_context_20260922/crystal_controls_v1/CRYSTALS.json')
        x = xyz(verify(m['cases'][0]['representations']['context']['endpoints']['Ca']['xyz']))
        moved = list(x); moved[1] = (x[1][0], x[1][1]+.01, *x[1][2:])
        self.assertFalse(same_geometry(x, moved))
        changed = list(x); changed[0] = ('La', *x[0][1:])
        self.assertFalse(same_geometry(x, changed))


if __name__ == '__main__': unittest.main()
